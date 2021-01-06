'''
Script for the computation of the Hc/D0 ratios from Pythia8 output trees
'''

import sys
import argparse
import pandas as pd
import numpy as np
import uproot
from ROOT import TH1F, TCanvas, TGaxis, gStyle, kBlack, kFullCircle # pylint: disable=import-error,no-name-in-module

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inFileName', help='input root file name')
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

gStyle.SetOptStat(0)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadRightMargin(0.035)
gStyle.SetPadTopMargin(0.035)
gStyle.SetTitleSize(0.05, 'xy')
gStyle.SetTitleOffset(1.4, 'y')
gStyle.SetLabelSize(0.045, 'y')
gStyle.SetLabelSize(0.055, 'x')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
TGaxis.SetMaxDigits(3)

if args.charm:
    labels = {'D0': 'D^{0}', 'Dplus': 'D^{+}', 'Dstar': 'D^{*+}', 'Ds': 'D_{s}^{+}', 'Lc': '#Lambda_{c}^{+}',
            'Sigmac0': '#Sigma_{c}^{0}', 'Sigmacplus': '#Sigma_{c}^{+}', 'Sigmacplusplus': '#Sigma_{c}^{++}',
            'Xic0': '#Xi_{c}^{0}', 'Xicplus': '#Xi_{c}^{+}', 'Omegac': '#Omega_{c}^{0}', 'Jpsi': 'J/#psi'}

    pdgCodes = {'D0': 421, 'Dplus': 411, 'Dstar': 413, 'Ds': 431,
                'Lc': 4122, 'Sigmac0': 4112, 'Sigmacplus': 4212, 'Sigmacplusplus': 4222,
                'Xic0': 4132, 'Xicplus': 4232, 'Omegac': 4332, 'Jpsi': 443}
    refPart = 'D0'
    nonStrMes = ['D0', 'Dplus']
    strangeMes = 'Ds'
else:
    labels = {'B0': 'B^{0}', 'Bplus': 'B^{+}', 'Bs': 'B_{s}^{0}', 'Lb': '#Lambda_{b}^{0}',
            'Sigmab0': '#Sigma_{b}^{0}', 'Sigmabplus': '#Sigma_{b}^{+}', 'Sigmabminus': '#Sigma_{b}^{-}',
            'Xib0': '#Xi_{b}^{0}', 'Xibminus': '#Xi_{b}^{-}', 'Omegab': '#Omega_{b}^{-}', 'Upsilon1S': '#Upsilon(1S)'}

    pdgCodes = {'Bplus': 511, 'B0': 521, 'Bs': 531,
                'Lb': 5122, 'Sigmabminus': 5112, 'Sigmab0': 5212, 'Sigmabplus': 5222,
                'Xibminus': 5132, 'Xib0': 5232, 'Omegab': 5332, 'Upsilon1S': 553}
    refPart = 'Bplus'
    nonStrMes = ['B0', 'Bplus']
    strangeMes = 'Bs'

hRatioToRefPart = TH1F(f'hRatioTo{refPart}', f';;Ratio to {labels[refPart]}', len(labels)-1, 0.5, len(labels)-0.5)
hRatioToRefPart.SetLineColor(kBlack)
hRatioToRefPart.SetMarkerColor(kBlack)
hRatioToRefPart.SetLineWidth(2)
hRatioToRefPart.SetMarkerStyle(kFullCircle)

kineTree = uproot.open(args.inFileName)['treeEvents']
kineDf, counts, unc, ratio, uncRatio = ({} for _ in range(5))
kineDf['all'] = kineTree.pandas.df().query('abs(y) < 0.5') # select only prompt and midrapidity
for iSpecie, specie in enumerate(pdgCodes):
    if 'Jpsi' not in specie:
        kineDf[specie] = kineDf['all'].query(f'origin == 4 and abs(pdg) == {pdgCodes[specie]}')
    else:
        kineDf[specie] = kineDf['all'].query(f'abs(pdg) == {pdgCodes[specie]}')
    counts[specie] = len(kineDf[specie])
    unc[specie] = np.sqrt(counts[specie])
    if refPart not in specie:
        ratio[f'{specie}Over{refPart}'] = counts[specie] / counts[refPart]
        uncRatio[f'{specie}Over{refPart}'] = \
            np.sqrt(1./counts[specie] + 1./counts[refPart]) * counts[specie] / counts[refPart]
        hRatioToRefPart.GetXaxis().SetBinLabel(iSpecie, labels[specie])
        hRatioToRefPart.SetBinContent(iSpecie, ratio[f'{specie}Over{refPart}'])
        hRatioToRefPart.SetBinError(iSpecie, uncRatio[f'{specie}Over{refPart}'])

cRatioToRefPart = TCanvas(f'cRatioTo{refPart}', '', 800, 800)
cRatioToRefPart.SetLogy()
hRatioToRefPart.GetYaxis().SetRangeUser(1.e-5, 1.)
hRatioToRefPart.DrawCopy('e')
cRatioToRefPart.Modified()
cRatioToRefPart.Update()

outDic = counts.copy()
outDic.update(ratio)
outUncDic = unc.copy()
outUncDic.update(uncRatio)

sumNonStrange = outDic[nonStrMes[0]] + outDic[nonStrMes[1]]
outDic['fsOverfufd'] = outDic[strangeMes] / sumNonStrange
uncSumNonStrange = np.sqrt(unc[nonStrMes[0]]**2 + unc[nonStrMes[1]]**2)
outUncDic['fsOverfufd'] = np.sqrt(uncSumNonStrange**2 / sumNonStrange**2 \
    + outUncDic[strangeMes]**2 / outDic[strangeMes]**2) * outDic['fsOverfufd']

outDic = pd.Series(outDic)
outUncDic = pd.Series(outUncDic)

outDf = pd.DataFrame(dict(value=outDic, error=outUncDic))
outDf.to_csv(args.inFileName.replace('.root', '.txt'), sep=' ', float_format='%0.6f')

input('Press enter to exit')
