'''
Script for the simulation of Hb -> nuclei decays
'''

import sys
import ctypes
import argparse
import numpy as np
import pandas as pd
import yaml
from ROOT import AliPythia8, TClonesArray, TLorentzVector, TDatabasePDG, gSystem, gRandom, vector # pylint: disable=import-error,no-name-in-module
from ROOT import TH1F, TFile, TTree, TF1 # pylint: disable=import-error,no-name-in-module


def SimpleCoalescence(pc, part1, part2, part3=None):
    '''
    Method to perform simple coalescence among 2 or 3 nucleons
    '''

    Ptrack1, Ptrack2, Ptrack3, trackSum, Ptrack1CMS, Ptrack2CMS, Ptrack3CMS = (TLorentzVector() for _ in range(7))
    Ptrack1.SetXYZM(part1['px'], part1['py'], part1['pz'], TDatabasePDG.Instance().GetParticle(2212).Mass())
    Ptrack2.SetXYZM(part2['px'], part2['py'], part2['pz'], TDatabasePDG.Instance().GetParticle(2212).Mass())
    Ptrack3.SetXYZM(part3['px'], part3['py'], part3['pz'], TDatabasePDG.Instance().GetParticle(2112).Mass())
    trackSum = Ptrack1 + Ptrack2 + Ptrack3

    beta = trackSum.Beta()
    betax = beta * np.cos(trackSum.Phi()) * np.sin(trackSum.Theta())
    betay = beta * np.sin(trackSum.Phi()) * np.sin(trackSum.Theta())
    betaz = beta * np.cos(trackSum.Theta())

    Ptrack1CMS = Ptrack1
    Ptrack2CMS = Ptrack2
    Ptrack3CMS = Ptrack3

    Ptrack1CMS.Boost(-betax, -betay, -betaz)
    Ptrack2CMS.Boost(-betax, -betay, -betaz)
    Ptrack3CMS.Boost(-betax, -betay, -betaz)

    coalRadius = np.power(2, 1./6) * pc / 2

    pHe3 = {}
    if Ptrack1CMS.P() <= coalRadius and Ptrack2CMS.P() <= coalRadius and Ptrack3CMS.P() <= coalRadius:
        for el in ['px', 'py', 'pz']:
            pHe3[el] = part1[el] + part2[el] + part3[el]
        return True, pHe3

    return False, {'px': -999., 'py': -999., 'pz': -999.}


def ReadFONLL(inFileName):
    '''
    Helper function to read FONLL txt files
    '''

    dfFONLL = pd.read_csv(inFileName, sep=' ', comment='#')
    ptCent = dfFONLL['pt'].to_numpy()
    xSec = dfFONLL['central'].to_numpy()
    deltaPt = ptCent[1]-ptCent[0]
    ptMin = ptCent[0] - deltaPt / 2
    ptMax = ptCent[-1] + deltaPt / 2

    hFONLL = TH1F('hFONLL', ';#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (pb GeV^{-1} #it{c})',
                  len(ptCent), ptMin, ptMax)
    for iPt, content in enumerate(xSec):
        hFONLL.SetBinContent(iPt+1, content)

    return hFONLL


parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgFileName', help='input config yaml file')
parser.add_argument('--seed', type=int, default=42,
                    help='seed for simulation')
parser.add_argument('--nevents', type=int, default=1000,
                    help='number of events')
parser.add_argument('--outfile', metavar='text', default='AnalysisResults.root',
                    help='output file name')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlConfigFile:
    cfg = yaml.load(ymlConfigFile, yaml.FullLoader)

process = cfg['generator']['process']
tune = cfg['generator']['tune']
fileFONLL = cfg['FONLL']['filename']

gSystem.Load('liblhapdf.so')
gSystem.Load('libpythia8.so')
gSystem.Load('libAliPythia8.so')
gSystem.Setenv('PYTHIA8DATA', gSystem.ExpandPathName('$ALICE_ROOT/PYTHIA8/pythia8/xmldoc'))
gSystem.Setenv('LHAPDF', gSystem.ExpandPathName('$ALICE_ROOT/LHAPDF'))
gSystem.Setenv('LHAPATH', gSystem.ExpandPathName('$ALICE_ROOT/LHAPDF/PDFsets'))

generator = AliPythia8()

# process (SOftQCD or HardQCD)
if process == 'SoftQCD':
    generator.ReadString('SoftQCD:all = on')
elif process == 'HardQCD':
    generator.ReadString('HardQCD:hardccbar = on')
    generator.ReadString('HardQCD:hardbbbar = on')
elif process == 'HardQCDCharm':
    generator.ReadString('HardQCD:hardccbar = on')
elif process == 'HardQCDBeauty':
    generator.ReadString('HardQCD:hardbbbar = on')
else:
    print('ERROR: pythia process not defined! Plese choose among {SoftQCD, HardQCD, HardQCDCharm, HardQCDBeauty}')
    sys.exit()

# tune (Monash, CRMode0, CRMode2, CRMode3)
if tune == 'Monash':
    generator.ReadString('Tune:pp = 14')
elif tune == 'CRMode0':
    generator.ReadString('Tune:pp = 14')
    generator.ReadString('ColourReconnection:mode = 1')
    generator.ReadString('ColourReconnection:allowDoubleJunRem = off')
    generator.ReadString('ColourReconnection:m0 = 2.9')
    generator.ReadString('ColourReconnection:allowJunctions = on')
    generator.ReadString('ColourReconnection:junctionCorrection = 1.43')
    generator.ReadString('ColourReconnection:timeDilationMode = 0')
    generator.ReadString('StringPT:sigma = 0.335')
    generator.ReadString('StringZ:aLund = 0.36')
    generator.ReadString('StringZ:bLund = 0.56')
    generator.ReadString('StringFlav:probQQtoQ = 0.078')
    generator.ReadString('StringFlav:ProbStoUD = 0.2')
    generator.ReadString('StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275')
    generator.ReadString('MultiPartonInteractions:pT0Ref = 2.12')
    generator.ReadString('BeamRemnants:remnantMode = 1')
    generator.ReadString('BeamRemnants:saturation =5')
elif tune == 'CRMode2':
    generator.ReadString('Tune:pp = 14')
    generator.ReadString('ColourReconnection:mode = 1')
    generator.ReadString('ColourReconnection:allowDoubleJunRem = off')
    generator.ReadString('ColourReconnection:m0 = 0.3')
    generator.ReadString('ColourReconnection:allowJunctions = on')
    generator.ReadString('ColourReconnection:junctionCorrection = 1.20')
    generator.ReadString('ColourReconnection:timeDilationMode = 2')
    generator.ReadString('ColourReconnection:timeDilationPar = 0.18')
    generator.ReadString('StringPT:sigma = 0.335')
    generator.ReadString('StringZ:aLund = 0.36')
    generator.ReadString('StringZ:bLund = 0.56')
    generator.ReadString('StringFlav:probQQtoQ = 0.078')
    generator.ReadString('StringFlav:ProbStoUD = 0.2')
    generator.ReadString('StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275')
    generator.ReadString('MultiPartonInteractions:pT0Ref = 2.15')
    generator.ReadString('BeamRemnants:remnantMode = 1')
    generator.ReadString('BeamRemnants:saturation =5')
elif tune == 'CRMode3':
    generator.ReadString('Tune:pp = 14')
    generator.ReadString('ColourReconnection:mode = 1')
    generator.ReadString('ColourReconnection:allowDoubleJunRem = off')
    generator.ReadString('ColourReconnection:m0 = 0.3')
    generator.ReadString('ColourReconnection:allowJunctions = on')
    generator.ReadString('ColourReconnection:junctionCorrection = 1.15')
    generator.ReadString('ColourReconnection:timeDilationMode = 3')
    generator.ReadString('ColourReconnection:timeDilationPar = 0.073')
    generator.ReadString('StringPT:sigma = 0.335')
    generator.ReadString('StringZ:aLund = 0.36')
    generator.ReadString('StringZ:bLund = 0.56')
    generator.ReadString('StringFlav:probQQtoQ = 0.078')
    generator.ReadString('StringFlav:ProbStoUD = 0.2')
    generator.ReadString('StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275')
    generator.ReadString('MultiPartonInteractions:pT0Ref = 2.05')
    generator.ReadString('BeamRemnants:remnantMode = 1')
    generator.ReadString('BeamRemnants:saturation =5')
else:
    print('ERROR: pythia tune not defined! Plese choose among {Monash, CRMode0, CRMode2, CRMode3}')
    sys.exit()

generator.ReadString('Random:setSeed = on')
generator.ReadString(f'Random:seed {args.seed}')

# keep only interesting decays, to be reweighted a posteriori
generator.ReadString('5122:onMode = off')
generator.ReadString('5122:onIfMatch = 2 1 2 2101') # bRatio='0.0120000' dominant one according to https://arxiv.org/pdf/2006.16251.pdf
# generator.ReadString('5122:onIfMatch = 2 1 4 2101') // bRatio='0.4411147'
# generator.ReadString('5122:onIfMatch = 2 4 1 2101') // bRatio='0.0910000'
# generator.ReadString('5122:onIfMatch = 2 3 2 2101') // bRatio='0.0120000'
# generator.ReadString('5122:onIfMatch = 2 3 4 2101') // bRatio='0.0800000'

# define output objects
outFile = TFile(args.outfile, 'recreate')

hBR = TH1F('hBR', 'BR', 2, 0.5, 2.5)
hBR.GetXaxis().SetBinLabel(1, '#Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}')
hBR.GetXaxis().SetBinLabel(2, '#Lambda_{b}^{0} #rightarrow (>=2)p + (>=1)n / #Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}')
hBR.SetBinContent(1, 0.0120000)

hFONLLLb = ReadFONLL(fileFONLL)
hFONLLLb.Scale(0.816) # f(b -> B) from e+e- provides good normalisation for LHCb and CMS B-meson measurements
hFONLLLb.SetName('hFONLLLb')

# f(b -> Lb) / f(b -> B) from LHCb measurement https://arxiv.org/pdf/1902.06794.pdf
fFFLHCb = TF1('fracLb','([4] * ([5] + exp([6] + [7] * x))) /  (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1)  ', 0, 50)
parLbA = 1
parLbp1 = 0.0793
parLbp2 = -1.022
parLbp3 = -0.107
parBsA = 1
parBsp1 = 0.119
parBsp2 = -0.00091
parBsAvePt = 10.1
fFFLHCb.SetParameters(parBsA, parBsp1, parBsp2, parBsAvePt, parLbA, parLbp1, parLbp2, parLbp3)

for iPt in range(1, hFONLLLb.GetNbinsX()+1):
    ptCent = hFONLLLb.GetBinCenter(iPt)
    hFONLLLb.SetBinContent(iPt, hFONLLLb.GetBinContent(iPt) * fFFLHCb.Eval((ptCent if ptCent > 5 else 5)))

hHe3FromLb = hFONLLLb.Clone('hHe3FromLb')
hHe3FromLb.Reset()

treeLb = TTree('treeLb', 'treeLb')
pxDau, pyDau, pzDau, ptDau = (vector('float')(0) for _ in range(4))
pdgDau = vector('int')(0)
pxLb, pyLb, pzLb, ptLb = (np.array([0], 'f') for _ in range(4))
treeLb.Branch('ptLb', ptLb, 'ptLb/F')
treeLb.Branch('pxLb', pxLb, 'pxLb/F')
treeLb.Branch('pyLb', pyLb, 'pyLb/F')
treeLb.Branch('pzLb', pzLb, 'pzLb/F')
treeLb.Branch('ptDau', ptDau)
treeLb.Branch('pxDau', pxDau)
treeLb.Branch('pyDau', pyDau)
treeLb.Branch('pzDau', pzDau)
treeLb.Branch('pdgDau', pdgDau)

generator.Pythia8().init()

particles = TClonesArray("TParticle", 1000)
pdgHb = 5122
massLb = TDatabasePDG.Instance().GetParticle(pdgHb).Mass()
massHe3 = 2.80839160743

nEventSel = 0
for iEvent in range(args.nevents):
    ptLb[0] = hFONLLLb.GetRandom()
    phiLb = gRandom.Rndm() * 2 * np.pi
    yLb = gRandom.Rndm() * 2. - 1. # flat in -1<y<1
    pxLb[0] = ptLb[0] * np.cos(phiLb)
    pyLb[0] = ptLb[0] * np.sin(phiLb)
    mt = np.sqrt(massLb * massLb + ptLb[0] * ptLb[0])
    pzLb[0] = mt * np.sinh(yLb)
    pLb = np.sqrt(ptLb[0] * ptLb[0] + pzLb[0] * pzLb[0])
    ELb = np.sqrt(massLb * massLb + pLb * pLb)
    generator.Pythia8().event.clear()
    generator.Pythia8().event.append(pdgHb, 11, 0, 0, pxLb[0], pyLb[0], pzLb[0], ELb, massLb)
    idPart = generator.Pythia8().event[0].id()
    generator.Pythia8().particleData.mayDecay(idPart, True)
    generator.Pythia8().moreDecays()
    generator.ImportParticles(particles, "All")
    nPart = particles.GetEntriesFast()

    if iEvent%1000000 == 0:
        print(f'Lb decay number {iEvent:10d}\r')

    nProtons, nNeutrons = 0, 0
    for part in particles:
        ptDau.push_back(part.Pt())
        pxDau.push_back(part.Px())
        pyDau.push_back(part.Py())
        pzDau.push_back(part.Pz())
        pdgDau.push_back(part.GetPdgCode())
        
        if pdgDau[-1] == 2212:
            nProtons += 1
        elif pdgDau[-1] == 2112:
            nNeutrons += 1

    if nProtons >= 2 and nNeutrons >= 1:
        pCoal1, pCoal2, pCoal3 = ({} for _ in range(3))
        isSecond = False
        for iPart, pdgD in enumerate(pdgDau):
            if pdgD == 2212:
                if not isSecond:
                    pCoal1['px'] = pxDau[iPart]
                    pCoal1['py'] = pyDau[iPart]
                    pCoal1['pz'] = pzDau[iPart]
                    isSecond = True
                else:
                    pCoal2['px'] = pxDau[iPart]
                    pCoal2['py'] = pyDau[iPart]
                    pCoal2['pz'] = pzDau[iPart]
            elif pdgD == 2112:
                pCoal3['px'] = pxDau[iPart]
                pCoal3['py'] = pyDau[iPart]
                pCoal3['pz'] = pzDau[iPart]

        hasCoalesced, pCoalHe3 = SimpleCoalescence(0.2, pCoal1, pCoal2, pCoal3)
        if hasCoalesced:
            ptHe3 = np.sqrt(pCoalHe3['px']*pCoalHe3['px'] + pCoalHe3['py']*pCoalHe3['py'])
            pHe3 = np.sqrt(ptHe3*ptHe3 + pCoalHe3['pz']*pCoalHe3['pz'])
            EHe3 = np.sqrt(massHe3 * massHe3 + pLb * pLb)
            yHe3 = np.log((EHe3 + pCoalHe3['pz']) / (EHe3 - pCoalHe3['pz']))
            hHe3FromLb.Fill(ptHe3)
        nEventSel += 1
        treeLb.Fill()

    ptDau.clear()
    pxDau.clear()
    pyDau.clear()
    pzDau.clear()
    pdgDau.clear()

hBR.SetBinContent(2, nEventSel / args.nevents)
if hHe3FromLb.Integral() > 0:
    hHe3FromLb.Scale(hFONLLLb.Integral() / nEventSel * hBR.GetBinContent(1) * hBR.GetBinContent(2))

# Save histogram on file and close file.
outFile.cd()
hFONLLLb.Write()
hHe3FromLb.Write()
hBR.Write()
treeLb.Write()
outFile.Close()
