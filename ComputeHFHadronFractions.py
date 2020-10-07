'''
Script for the computation of the Hc/D0 ratios from Pythia8 output trees
'''

import argparse
import pandas as pd
import numpy as np
import uproot
from ROOT import TH1F, TCanvas, TGaxis, gStyle, kBlack, kFullCircle # pylint: disable=import-error,no-name-in-module

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inFileName', help='input root file name')
args = parser.parse_args()

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

labels = {'D0': 'D^{0}', 'Dplus': 'D^{+}', 'Dstar': 'D^{*+}', 'Ds': 'D_{s}^{+}', 'Lc': '#Lambda_{c}^{+}',
          'Sigmac0': '#Sigma_{c}^{0}', 'Sigmacplus': '#Sigma_{c}^{+}', 'Sigmacplusplus': '#Sigma_{c}^{++}',
          'Xic0': '#Xi_{c}^{0}', 'Xicplus': '#Xi_{c}^{+}', 'Omegac': '#Omega_{c}^{0}', 'Jpsi': 'J/#psi'}

pdgCodes = {'D0': 421, 'Dplus': 411, 'Dstar': 413, 'Ds': 431,
            'Lc': 4122, 'Sigmac0': 4112, 'Sigmacplus': 4212, 'Sigmacplusplus': 4222,
            'Xic0': 4132, 'Xicplus': 4232, 'Omegac': 4332, 'Jpsi': 443}

hRatioToD0 = TH1F('hRatioToD0', ';;Ratio to D^{0}', len(labels)-1, 0.5, len(labels)-0.5)
hRatioToD0.SetLineColor(kBlack)
hRatioToD0.SetMarkerColor(kBlack)
hRatioToD0.SetLineWidth(2)
hRatioToD0.SetMarkerStyle(kFullCircle)

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
    if 'D0' not in specie:
        ratio[f'{specie}OverD0'] = counts[specie] / counts['D0']
        uncRatio[f'{specie}OverD0'] = np.sqrt(1./counts[specie] + 1./counts['D0']) * counts[specie] / counts['D0']
        hRatioToD0.GetXaxis().SetBinLabel(iSpecie, labels[specie])
        hRatioToD0.SetBinContent(iSpecie, ratio[f'{specie}OverD0'])
        hRatioToD0.SetBinError(iSpecie, uncRatio[f'{specie}OverD0'])

cRatioToD0 = TCanvas('cRatioToD0', '', 800, 800)
cRatioToD0.SetLogy()
hRatioToD0.GetYaxis().SetRangeUser(1.e-5, 1.)
hRatioToD0.DrawCopy('e')
cRatioToD0.Modified()
cRatioToD0.Update()

outDic = counts.copy()
outDic.update(ratio)
outUncDic = unc.copy()
outUncDic.update(uncRatio)

D0plusDplus = outDic['D0'] + outDic['Dplus']
outDic['fsOverfufd'] = outDic['Ds'] / D0plusDplus
uncD0plusDplus = np.sqrt(unc['D0']**2 + unc['Dplus']**2)
outUncDic['fsOverfufd'] = np.sqrt(
    uncD0plusDplus**2 / D0plusDplus**2 + outUncDic['Ds']**2 / outDic['Ds']**2) * outDic['fsOverfufd']

outDic = pd.Series(outDic)
outUncDic = pd.Series(outUncDic)

outDf = pd.DataFrame(dict(value=outDic, error=outUncDic))
outDf.to_csv(args.inFileName.replace('.root', '.txt'), sep=' ', float_format='%0.6f')

input('Press enter to exit')
