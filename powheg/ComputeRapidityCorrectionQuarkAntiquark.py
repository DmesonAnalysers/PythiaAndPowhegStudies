'''
Python script for the evaluation of the difference between the rapidity distribution of single quark and quark-antiquark pairs
'''

import sys
import argparse
import numpy as np
import pandas as pd
from root_numpy import fill_hist
from ROOT import TH2F, TH1F, TCanvas, TLegend, TGaxis, TFile, gStyle, kRainBow, kRed, kAzure, kGreen, kBlack, kFullCircle # pylint: disable=import-error,no-name-in-module
from ReadLHE import ReadLHEFile

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
parser.add_argument('inFileName', help='input lhe file name')
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

# load powheg simulation output
dfEvents = ReadLHEFile(args.inFileName)

# compute rapidity
pd.eval("y = 1./2 * log((dfEvents.E+dfEvents.pz)/(dfEvents.E-dfEvents.pz))", target=dfEvents, inplace=True)

dfQ = dfEvents.query(f'pdg == {quark}')
dfAQ = dfEvents.query(f'pdg == -{quark}')

dfQAQ = dfQ.merge(dfAQ, left_on='event', right_on='event', suffixes=('Q', 'AQ'))

# compute rapidity of quark-antiquark pair
pd.eval("pzQAQ = dfQAQ.pzQ + dfQAQ.pzAQ", target=dfQAQ, inplace=True)
pd.eval("EQAQ = dfQAQ.EQ + dfQAQ.EAQ", target=dfQAQ, inplace=True)
pd.eval("yQAQ = 1./2 * log((dfQAQ.EQAQ+dfQAQ.pzQAQ)/(dfQAQ.EQAQ-dfQAQ.pzQAQ))", target=dfQAQ, inplace=True)

hRapidityQuarkVsAntiQuark = TH2F('hRapidityQuarkVsAntiQuark', f';#it{{y}}({quarkName});#it{{y}}(#bar{{{quarkName}}})',
                                 400, -8., 8., 400, -8., 8.)
fill_hist(hRapidityQuarkVsAntiQuark, np.array([dfQAQ['yQ'].values, dfQAQ['yAQ'].values]).transpose())

hRapidityQuark = hRapidityQuarkVsAntiQuark.ProjectionY('hRapidityQuark')
hRapidityQuark.SetLineColor(kRed+1)
hRapidityQuark.SetLineWidth(2)
hRapidityQuark.GetXaxis().SetTitle('#it{y}')
hRapidityQuark.GetYaxis().SetTitle('d#sigma/d#it{y} (a.u)')
hRapidityAntiQuark = hRapidityQuarkVsAntiQuark.ProjectionX('hRapidityAntiQuark')
hRapidityAntiQuark.SetLineColor(kAzure+4)
hRapidityAntiQuark.SetLineWidth(2)
hRapidityAntiQuark.GetXaxis().SetTitle('#it{y}')
hRapidityAntiQuark.GetYaxis().SetTitle('d#sigma/d#it{y} (a.u)')

hRapidityAverageQuarkAntiQuark = hRapidityAntiQuark.Clone('hRapidityAverageQuarkAntiQuark')
hRapidityAverageQuarkAntiQuark.Add(hRapidityAntiQuark)
hRapidityAverageQuarkAntiQuark.Scale(1./2)

hRapidityPair = TH1F('hRapidityPair', f'd#sigma/d#it{{y}} (a.u);#it{{y}}({quarkName}#bar{{{quarkName}}})', 400, -8., 8.)
fill_hist(hRapidityPair, dfQAQ['yQAQ'].values)
hRapidityPair.SetLineColor(kGreen+2)
hRapidityPair.SetLineWidth(2)

cCorr = TCanvas('cCorr', '', 500, 500)
cCorr.SetLogz()
hRapidityQuarkVsAntiQuark.Draw('colz')

leg = TLegend(0.65, 0.75, 0.85, 0.9)
leg.SetTextSize(0.045)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(hRapidityQuark, f'{quarkName}', 'l')
leg.AddEntry(hRapidityAntiQuark, f'#bar{{{quarkName}}}', 'l')
leg.AddEntry(hRapidityPair, f'{quarkName}#bar{{{quarkName}}}', 'l')

cY = TCanvas('cY', '', 500, 500)
hRapidityQuark.GetYaxis().SetRangeUser(0., hRapidityQuark.GetMaximum()*1.1)
hRapidityQuark.Draw()
hRapidityAntiQuark.Draw('same')
hRapidityPair.Draw('same')
leg.Draw()

yBinMin = hRapidityAverageQuarkAntiQuark.GetXaxis().FindBin(args.yMin)
yBinMax = hRapidityAverageQuarkAntiQuark.GetXaxis().FindBin(args.yMax)

integralQuark = hRapidityAverageQuarkAntiQuark.Integral(yBinMin, yBinMax)
integralQuarkUnc = np.sqrt(integralQuark)
integralPair = hRapidityPair.Integral(yBinMin, yBinMax)
integralPairUnc = np.sqrt(integralPair)
ratio = integralPair/integralQuark
ratioUnc = np.sqrt((integralPairUnc/integralPair)**2 + (integralQuarkUnc/integralQuark)**2) * ratio

hCorrFactor = TH1F('hCorrFactorRapidityPOWHEG',
                   f';;d#sigma_{{{quarkName}#bar{{{quarkName}}}}}/d#it{{y}}_{{{quarkName}#bar{{{quarkName}}}}}|_{{|#it{{y}}_{{{quarkName}#bar{{{quarkName}}}}}| < 0.5}} '
                   f'/ d#sigma_{{{quarkName}}}/d#it{{y}}_{{{quarkName}}}|_{{|#it{{y}}_{{{quarkName}}}| < 0.5}}', 1, 0.5, 1.5)
hCorrFactor.SetBinContent(1, ratio)
hCorrFactor.SetBinError(1, ratioUnc)
hCorrFactor.SetLineColor(kBlack)
hCorrFactor.SetMarkerColor(kBlack)
hCorrFactor.SetMarkerStyle(kFullCircle)
hCorrFactor.SetLineWidth(2)
hCorrFactor.GetYaxis().SetRangeUser(0, 1.15)
hCorrFactor.GetYaxis().SetDecimals()

outFile = TFile(args.inFileName.replace('.lhe', '_corrFactor_rapidity_POWHEG.root'), 'recreate')
cCorr.Write()
cY.Write()
hRapidityQuark.Write()
hRapidityAntiQuark.Write()
hRapidityPair.Write()
hCorrFactor.Write()
outFile.Close()

input('Press enter to exit')
