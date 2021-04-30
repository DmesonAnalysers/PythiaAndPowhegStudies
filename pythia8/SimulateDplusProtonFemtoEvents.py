'''
Script for the simulation of D+ - p CF with pythia and herwig
'''

import sys
import argparse
import numpy as np
import yaml
from ROOT import AliPythia8, THerwig6, TClonesArray, TLorentzVector, TDatabasePDG, TH1F, TFile, TDirectoryFile, TList, TObject, TNtuple # pylint: disable=import-error,no-name-in-module


def ComputeKstar(part1, part2):
    '''
    Helper function for k* computation
    
    ----------------
    Parameters
    - part1: 4 element dictionary with pdg code, px, py, pz for first particle
    - kStar: 4 element dictionary with pdg code, px, py, pz for second particle

    ----------------
    Returns
    - kStar: relative momentum between part1 and part2 in their rest frame
    '''

    mass1 = TDatabasePDG.Instance().GetParticle(abs(part1['pdg'])).Mass()
    mass2 = TDatabasePDG.Instance().GetParticle(abs(part2['pdg'])).Mass()
    SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS, trackRelK = (TLorentzVector() for _ in range(6))
    SPtrack.SetXYZM(part1['px'], part1['py'], part1['pz'], mass1)
    TPProng.SetXYZM(part2['px'], part2['py'], part2['pz'], mass2)
    trackSum = SPtrack + TPProng

    beta = trackSum.Beta()
    betax = beta * np.cos(trackSum.Phi()) * np.sin(trackSum.Theta())
    betay = beta * np.sin(trackSum.Phi()) * np.sin(trackSum.Theta())
    betaz = beta * np.cos(trackSum.Theta())

    SPtrackCMS = SPtrack
    TPProngCMS = TPProng

    SPtrackCMS.Boost(-betax, -betay, -betaz)
    TPProngCMS.Boost(-betax, -betay, -betaz)

    trackRelK = SPtrackCMS - TPProngCMS
    kStar = 0.5 * trackRelK.P()
    return kStar


def CheckOrigin(motherPdg, pdg):
    '''
    Helper function to check origin of particle

    ----------------
    Parameters
    - motherPdg: list of pdg codes of the mother particles

    ----------------
    Returns
    - origin: 0 for light, 4 for charm, 5 for beauty
    '''

    origin = 0
    absMotherPdg = abs(np.array(motherPdg))
    if any(absMotherPdg == 5):
        origin = 5
    elif any((absMotherPdg / 100) == 5):
        origin = 5
    elif any((absMotherPdg / 1000) == 5):
        origin = 5
    elif any(((absMotherPdg - 10000) / 100) == 5):
        origin = 5
    elif any(((absMotherPdg - 20000) / 100) == 5):
        origin = 5
    elif any(((absMotherPdg - 100000) / 100) == 5):
        origin = 5
    elif any(((absMotherPdg - 200000) / 100) == 5):
        origin = 5

    if origin < 5:
        if abs(pdg) == 411:
            origin = 4            
        elif any(absMotherPdg == 4):
            origin = 4
        elif any((absMotherPdg / 100) == 4):
            origin = 4
        elif any((absMotherPdg / 1000) == 4):
            origin = 4
        elif any(((absMotherPdg - 10000) / 100) == 4):
            origin = 4
        elif any(((absMotherPdg - 20000) / 100) == 4):
            origin = 4
        elif any(((absMotherPdg - 100000) / 100) == 4):
            origin = 4
        elif any(((absMotherPdg - 200000) / 100) == 4):
            origin = 4

    return origin


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

genType = cfg['generator']['type']
process = cfg['generator']['process']
tune = cfg['generator']['tune']
energy = cfg['energy']

if genType == 'Pythia8':

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
    generator.Initialize(2212, 2212, energy)

elif genType == 'Herwig6':
    generator = THerwig6()
    generator.Initialize('P                   ', 'P                   ', energy/2, energy/2, tune)
    generator.SetMODPDF(1, 19070)
    generator.SetMODPDF(2, 19070)
    generator.SetAUTPDF(1, 'HWLHAPDF            ')
    generator.SetAUTPDF(2, 'HWLHAPDF            ')
    generator.SetRMASS(4, 1.2)
    generator.SetRMASS(5, 4.75)
    generator.SetPTMIN(0.)
    generator.SetPTMAX(10000.)
    generator.SetPTRMS(0.)
    generator.SetMAXPR(10)
    generator.SetMAXER(1000)
    generator.SetENSOF(1)
    generator.PrepareRun()
else:
    print('ERROR: MC generator not defined! Plese choose among {Pythia8, Herwig6}')
    sys.exit()

# define output objects
fHistNEvents = TH1F("fHistNEvents", "Number of processed events", 3, -0.5, 2.5)
fHistNEvents.GetXaxis().SetBinLabel(1, "Read events")
fHistNEvents.GetXaxis().SetBinLabel(2, "Selected SE events")
fHistNEvents.GetXaxis().SetBinLabel(3, "ME events")
fHistNEvents.SetMinimum(0)

fMult = TH1F("fMult", ";multiplicity;counts", 1000, 0., 1000.)

fHistDplusProtonPairsDPrompt = TH1F("fHistDplusProtonPairsDPrompt", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusProtonPairsDFromB = TH1F("fHistDplusProtonPairsDFromB", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusProtonPairsDFromDstar = TH1F("fHistDplusProtonPairsDFromDstar", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusProtonPairsDpCommonMother = TH1F("fHistDplusProtonPairsDpCommonMother", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusAntiProtonPairsDPrompt = TH1F("fHistDplusAntiProtonPairsDPrompt", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusAntiProtonPairsDFromB = TH1F("fHistDplusAntiProtonPairsDFromB", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusAntiProtonPairsDFromDstar = TH1F("fHistDplusAntiProtonPairsDFromDstar", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusAntiProtonPairsDpCommonMother = TH1F("fHistDplusAntiProtonPairsDpCommonMother", ";#it{k}* (GeV/#it{c})",
                                               1000, 0., 1.)
fHistDplusProtonPairsMEDPrompt = TH1F("fHistDplusProtonPairsMEDPrompt", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusProtonPairsMEDFromB = TH1F("fHistDplusProtonPairsMEDFromB", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusProtonPairsMEDFromDstar = TH1F("fHistDplusProtonPairsMEDFromDstar", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusProtonPairsMEDpCommonMother = TH1F("fHistDplusProtonPairsMEDpCommonMother", ";#it{k}* (GeV/#it{c})",
                                             1000, 0., 1.)
fHistDplusAntiProtonPairsMEDPrompt = TH1F("fHistDplusAntiProtonPairsMEDPrompt", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusAntiProtonPairsMEDFromB = TH1F("fHistDplusAntiProtonPairsMEDFromB", ";#it{k}* (GeV/#it{c})", 1000, 0., 1.)
fHistDplusAntiProtonPairsMEDFromDstar = TH1F("fHistDplusAntiProtonPairsMEDFromDstar", ";#it{k}* (GeV/#it{c})",
                                             1000, 0., 1.)
fHistDplusAntiProtonPairsMEDpCommonMother = TH1F("fHistDplusAntiProtonPairsMEDpCommonMother", ";#it{k}* (GeV/#it{c})",
                                                 1000, 0., 1.)

fNtupleCommonMothers = TNtuple("fNtupleCommonMothers", "fNtupleCommonMothers", "pdgMother:pdgD:pdgp:kStar")
fNtupleCommonMothers.SetDirectory(0)

# buffers for mixed event
DplusME, DminusME, PrME, antiPrME = ([] for _ in range(4))
# run simulation
particles = TClonesArray('TParticle', 1000)
for iEvent in range(args.nevents):
    generator.GenerateEvent()
    generator.ImportParticles(particles, 'All')

    fHistNEvents.Fill(0)
    fMult.Fill(particles.GetEntriesFast())

    # TODO: add a multiplicity cut here

    fHistNEvents.Fill(1)

    Pr, antiPr, Dplus, Dminus = ([] for _ in range(4)) 
    for iPart, part in enumerate(particles):
        if iPart < 2: # 0 and 1 beam protons
            continue
        pdg = part.GetPdgCode()
        absPdg = abs(pdg)
        if absPdg not in [411, 2212]:
            continue

        px = part.Px()
        py = part.Py()
        pz = part.Pz()
        if abs(part.Y()) > 2:
            continue 

        # check origin
        mothIdx = part.GetFirstMother()
        mothIndices, pdgMothers = ([] for _ in range(2))
        while mothIdx > 1: # 0 and 1 beam protons
            motherPart = particles.At(mothIdx)
            mothIdx = motherPart.GetFirstMother()
            mothIndices.append(mothIdx)
            pdgMothers.append(motherPart.GetPdgCode())

        origin = CheckOrigin(pdgMothers, pdg)
        stateVec = {'pdg': pdg, 'px': px, 'py': py, 'pz': pz, 'idx': iPart,
                    'motherIdx': mothIndices, 'motherPdg': pdgMothers, 'origin': origin}

        if pdg == 411:
            Dplus.append(stateVec)
        elif pdg == -411:
            Dminus.append(stateVec)
        elif pdg == 2212:
            Pr.append(stateVec)
        elif pdg == -2212:
            antiPr.append(stateVec)

        # fill mix event buffers
        if len(DplusME) > 49:
            DplusME.pop()
            DminusME.pop()
            PrME.pop()
            antiPrME.pop()
        DplusME.append(Dplus)
        DminusME.append(Dminus)
        PrME.append(Pr)
        antiPrME.append(antiPr)

    # form pairs from same event
    for part1 in Dplus:
        for part2 in antiPr:
            if part1['idx'] not in part2['motherIdx']:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusAntiProtonPairsDPrompt.Fill(kStar)
                    else:
                        fHistDplusAntiProtonPairsDFromDstar.Fill(kStar)
                if part1['origin'] == 5:
                    fHistDplusAntiProtonPairsDFromB.Fill(kStar)
                    commonMother = next((idx for idx in part1['motherIdx'] if idx in part2['motherIdx']), None)
                    if commonMother:
                        fHistDplusAntiProtonPairsDpCommonMother.Fill(kStar)
                        idxCommonMother = part1['motherIdx'].index(commonMother)
                        fNtupleCommonMothers.Fill(part1['motherPdg'][idxCommonMother], part1['pdg'], part2['pdg'], kStar)
        for part2 in Pr:
            if part1['idx'] not in part2['motherIdx']:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusAntiProtonPairsDPrompt.Fill(kStar)
                    else:
                        fHistDplusAntiProtonPairsDFromDstar.Fill(kStar)
                if part1['origin'] == 5:
                    fHistDplusAntiProtonPairsDFromB.Fill(kStar)
                    commonMother = next((idx for idx in part1['motherIdx'] if idx in part2['motherIdx']), None)
                    if commonMother:
                        fHistDplusAntiProtonPairsDpCommonMother.Fill(kStar)
                        idxCommonMother = part1['motherIdx'].index(commonMother)
                        fNtupleCommonMothers.Fill(part1['motherPdg'][idxCommonMother], part1['pdg'], part2['pdg'], kStar)
    for part1 in Dminus:
        for part2 in Pr:
            if part1['idx'] not in part2['motherIdx']:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusProtonPairsDPrompt.Fill(kStar)
                    else:
                        fHistDplusProtonPairsDFromDstar.Fill(kStar)
                elif part1['origin'] == 5:
                    fHistDplusProtonPairsDFromB.Fill(kStar)
                    commonMother = next((idx for idx in part1['motherIdx'] if idx in part2['motherIdx']), None)
                    if commonMother:
                        fHistDplusProtonPairsDpCommonMother.Fill(kStar)
                        idxCommonMother = part1['motherIdx'].index(commonMother)
                        fNtupleCommonMothers.Fill(part1['motherPdg'][idxCommonMother], part1['pdg'], part2['pdg'], kStar)
        for part2 in antiPr:
            if part1['idx'] not in part2['motherIdx']:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusProtonPairsDPrompt.Fill(kStar)
                    else:
                        fHistDplusProtonPairsDFromDstar.Fill(kStar)
                if part1['origin'] == 5:
                    fHistDplusProtonPairsDFromB.Fill(kStar)
                    commonMother = next((idx for idx in part1['motherIdx'] if idx in part2['motherIdx']), None)
                    if commonMother:
                        fHistDplusProtonPairsDpCommonMother.Fill(kStar)
                        idxCommonMother = part1['motherIdx'].index(commonMother)
                        fNtupleCommonMothers.Fill(part1['motherPdg'][idxCommonMother], part1['pdg'], part2['pdg'], kStar)
    
    # form pairs from mixed events (mix D's from events with protons from other events)
    if len(DplusME) < 2:
        continue

    fHistNEvents.Fill(2, (len(DplusME)-1))

    for part1 in DplusME[-1]:
        for protonsFromME, antiProtonsFromME in zip(PrME[:-1], antiPrME[:-1]):
            for part2 in protonsFromME:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusAntiProtonPairsMEDPrompt.Fill(kStar)
                    else:
                        fHistDplusAntiProtonPairsMEDFromDstar.Fill(kStar)
                elif part1['origin'] == 5:
                    fHistDplusAntiProtonPairsMEDFromB.Fill(kStar)
            for part2 in antiProtonsFromME:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusAntiProtonPairsMEDPrompt.Fill(kStar)
                    else:
                        fHistDplusAntiProtonPairsMEDFromDstar.Fill(kStar)
                elif part1['origin'] == 5:
                    fHistDplusAntiProtonPairsMEDFromB.Fill(kStar)
    for part1 in DminusME[-1]:
        for protonsFromME, antiProtonsFromME in zip(PrME[:-1], antiPrME[:-1]):
            for part2 in protonsFromME:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusProtonPairsMEDPrompt.Fill(kStar)
                    else:
                        fHistDplusProtonPairsMEDFromDstar.Fill(kStar)
                elif part1['origin'] == 5:
                    fHistDplusProtonPairsMEDFromB.Fill(kStar)
            for part2 in antiProtonsFromME:
                kStar = ComputeKstar(part1, part2)
                if part1['origin'] == 4:
                    if 413 not in part1['motherPdg'] and -413 not in part1['motherPdg']:
                        fHistDplusProtonPairsMEDPrompt.Fill(kStar)
                    else:
                        fHistDplusProtonPairsMEDFromDstar.Fill(kStar)
                elif part1['origin'] == 5:
                    fHistDplusProtonPairsMEDFromB.Fill(kStar)

# output file
outFile = TFile.Open(args.outfile, 'recreate')
dirName = 'DprotonPairGen'
listName = 'coutputDprotonPairGen'
outDir = TDirectoryFile(dirName, dirName)
outDir.Write()
outDir.cd()
fNtupleCommonMothers.Write()
outList = TList()
outList.Add(fHistNEvents)
outList.Add(fMult)
outList.Add(fHistDplusProtonPairsDPrompt)
outList.Add(fHistDplusProtonPairsDFromB)
outList.Add(fHistDplusProtonPairsDFromDstar)
outList.Add(fHistDplusProtonPairsDpCommonMother)
outList.Add(fHistDplusAntiProtonPairsDPrompt)
outList.Add(fHistDplusAntiProtonPairsDFromB)
outList.Add(fHistDplusAntiProtonPairsDFromDstar)
outList.Add(fHistDplusAntiProtonPairsDpCommonMother)
outList.Add(fHistDplusProtonPairsMEDPrompt)
outList.Add(fHistDplusProtonPairsMEDFromB)
outList.Add(fHistDplusProtonPairsMEDFromDstar)
outList.Add(fHistDplusProtonPairsMEDpCommonMother)
outList.Add(fHistDplusAntiProtonPairsMEDPrompt)
outList.Add(fHistDplusAntiProtonPairsMEDFromB)
outList.Add(fHistDplusAntiProtonPairsMEDFromDstar)
outList.Add(fHistDplusAntiProtonPairsMEDpCommonMother)
outList.Write(listName, TObject.kSingleKey)
outFile.Close()
