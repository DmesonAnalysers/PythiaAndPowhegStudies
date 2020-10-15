'''
Python script for the evaluation of the difference between the rapidity distribution
of the beauty quark and the corresponding beauty hadron
'''

import sys
import argparse
import numpy as np
import uproot
from root_numpy import fill_hist
from ROOT import TH2F, TH1F, TFile, TCanvas, TLegend, TGaxis # pylint: disable=import-error,no-name-in-module
from ROOT import gStyle, kRainBow, kRed, kAzure, kBlack, kFullCircle # pylint: disable=import-error,no-name-in-module

gStyle.SetOptStat(0)
gStyle.SetPadLeftMargin(0.12)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadRightMargin(0.12)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleSize(0.05, 'xy')
gStyle.SetTitleOffset(0.9, 'y')
gStyle.SetLabelSize(0.045, 'xy')
gStyle.SetPalette(kRainBow)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
TGaxis.SetMaxDigits(3)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inFileName', help='input root file name')
parser.add_argument('--yMin', type=float, default=-0.5, help='minimum rapidity')
parser.add_argument('--yMax', type=float, default=0.5, help='maximum rapidity')
parser.add_argument('--charm', action='store_true', default=False,
                    help='activate calculation for charm quarks')
parser.add_argument('--beauty', action='store_true', default=False,
                    help='activate calculation for beauty quarks')
args = parser.parse_args()

if args.charm and args.beauty:
    print('ERROR: only calculation for one flavour can be performed at a time! Exit')
    sys.exit()
if not args.charm and not args.beauty:
    print('ERROR: calculation for charm or beauty should be enabled! Exit')
    sys.exit()

if args.charm:
    quarkName = 'c'
    quark = 4
elif args.beauty:
    quarkName = 'b'
    quark = 5

kineTree = uproot.open(args.inFileName)['treeEvents']
kineDf = kineTree.pandas.df()
dfQuarks = kineDf.query(f'pdg == {quark}')
dfAntiQuarks = kineDf.query(f'pdg == -{quark}')
# do not consider higher mass states
if args.charm:
    dfHadrons = kineDf.query('pdg in [411, 421, 431, 4122]')
    dfAntiHadrons = kineDf.query('pdg in [-411, -421, -431, -4122]')
elif args.beauty:
    dfHadrons = kineDf.query('pdg in [511, 521, 531, 5122]')
    dfAntiHadrons = kineDf.query('pdg in [-511, -521, -531, -5122]')

dfPart = dfQuarks.merge(dfAntiHadrons, left_on='event_number', right_on='event_number', suffixes=('Q', 'H'))
dfAntiPart = dfAntiQuarks.merge(dfHadrons, left_on='event_number', right_on='event_number', suffixes=('Q', 'H'))

hRapidityQuarkVsHadron = TH2F('hRapidityQuarkVsHadron', f';#it{{y}}({quarkName});#it{{y}}(H_{{{quarkName}}})',
                              400, -8., 8., 400, -8., 8.)
hRapidityAntiQuarkVsHadron = TH2F('hRapidityAntiQuarkVsHadron', f';#it{{y}}({quarkName});#it{{y}}(H_{{quarkName}})',
                                  400, -8., 8., 400, -8., 8.)
fill_hist(hRapidityQuarkVsHadron, np.array([dfPart['yQ'].values, dfPart['yH'].values]).transpose())
fill_hist(hRapidityAntiQuarkVsHadron, np.array([dfAntiPart['yQ'].values, dfAntiPart['yH'].values]).transpose())
hRapidityQuarkVsHadron.Add(hRapidityAntiQuarkVsHadron)
hRapidityQuark = hRapidityQuarkVsHadron.ProjectionX()
hRapidityQuark.SetLineColor(kAzure+4)
hRapidityQuark.SetLineWidth(2)
hRapidityQuark.GetXaxis().SetTitle('#it{y}')
hRapidityQuark.GetYaxis().SetTitle('d#sigma/d#it{y} (a.u)')
hRapidityHadron = hRapidityQuarkVsHadron.ProjectionY()
hRapidityHadron.SetLineColor(kRed+1)
hRapidityHadron.SetLineWidth(2)
hRapidityHadron.GetXaxis().SetTitle('#it{y}')
hRapidityHadron.GetYaxis().SetTitle('d#sigma/d#it{y} (a.u)')

cCorr = TCanvas('cCorr', '', 500, 500)
cCorr.SetLogz()
hRapidityQuarkVsHadron.Draw('colz')

leg = TLegend(0.6, 0.75, 0.8, 0.9)
leg.SetTextSize(0.045)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(hRapidityQuark, f'{quarkName}-quark', 'l')
leg.AddEntry(hRapidityHadron, f'{quarkName}-hadron', 'l')

cY = TCanvas('cY', '', 500, 500)
hRapidityQuark.Draw()
hRapidityHadron.Draw('same')
leg.Draw()

yBinMin = hRapidityHadron.GetXaxis().FindBin(args.yMin)
yBinMax = hRapidityHadron.GetXaxis().FindBin(args.yMax)
integralQuark = hRapidityQuark.Integral(yBinMin, yBinMax)
integralHadron = hRapidityHadron.Integral(yBinMin, yBinMax)
integralQuarkUnc = np.sqrt(integralQuark)
integralHadronUnc = np.sqrt(integralHadron)
if integralHadron > 0:
    ratio = integralQuark / integralHadron
    ratioUnc = np.sqrt((integralQuarkUnc/integralQuark)**2 + (integralHadronUnc/integralHadron)**2) * ratio
else:
    ratio = 0
    ratioUnc = 0

hCorrFactor = TH1F('hCorrFactorRapidityPythia',
                   f';;d#sigma_{{{quarkName}}}/d#it{{y}}_{{{quarkName}}}|_{{|#it{{y}}_{{{quarkName}}}| < 0.5}} '
                   f'/ d#sigma_{{H_{{{quarkName}}}}}/d#it{{y}}_{{H_{{{quarkName}}}}}|_{{|#it{{y}}_{{H_{{{quarkName}}}}}}}| < 0.5}}', 1, 0.5, 1.5)
hCorrFactor.SetBinContent(1, ratio)
hCorrFactor.SetBinError(1, ratioUnc)
hCorrFactor.SetLineColor(kBlack)
hCorrFactor.SetMarkerColor(kBlack)
hCorrFactor.SetMarkerStyle(kFullCircle)
hCorrFactor.SetLineWidth(2)
hCorrFactor.GetYaxis().SetRangeUser(0, 1.15)
hCorrFactor.GetYaxis().SetDecimals()

outFile = TFile(args.inFileName.replace('.root', '_corrFactor_rapidity.root'), 'recreate')
cCorr.Write()
cY.Write()
hRapidityQuark.Write()
hRapidityHadron.Write()
hRapidityQuarkVsHadron.Write()
hCorrFactor.Write()
outFile.Close()

input('Press enter to exit')
