#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <array>

#include <TPythia8.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TTree.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TLorentzVector.h>

#endif

enum DmesonOrigin{kPrompt, kFromB, kFromDstar, kFromExcitedCharm};
enum HadronisationTunes{kMonash13, kCRMode0, kCRMode2, kCRMode3};

//_____________________________________________________________________
// METHOD PROTOTYPES
void GetGenDplusProtonCF(int nevents = 100000, int optHadronisation=kMonash13, int seed = 42);
double ComputeKstar(TParticle *part1, TParticle *part2);

//_____________________________________________________________________
// METHOD IMPLEMENTATIONS
void GetGenDplusProtonCF(int nevents, int optHadronisation, int seed)
{
    // create pythia generator
    TPythia8 pythia;

    // set SoftQCD processes
    pythia.ReadString("SoftQCD:all = on ");

    // set tune
    switch(optHadronisation)
    {
        case kMonash13:
        {
            pythia.ReadString("Tune:pp = 14");
            break;
        }
        case kCRMode0:
        {
            pythia.ReadString("Tune:pp = 14");
            pythia.ReadString("ColourReconnection:mode = 1");
            pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
            pythia.ReadString("ColourReconnection:m0 = 2.9");
            pythia.ReadString("ColourReconnection:allowJunctions = on");
            pythia.ReadString("ColourReconnection:junctionCorrection = 1.43");
            pythia.ReadString("ColourReconnection:timeDilationMode = 0");
            pythia.ReadString("StringPT:sigma = 0.335");
            pythia.ReadString("StringZ:aLund = 0.36");
            pythia.ReadString("StringZ:bLund = 0.56");
            pythia.ReadString("StringFlav:probQQtoQ = 0.078");
            pythia.ReadString("StringFlav:ProbStoUD = 0.2");
            pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
            pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.12");
            pythia.ReadString("BeamRemnants:remnantMode = 1");
            pythia.ReadString("BeamRemnants:saturation =5");
            break;
        }
        case kCRMode2:
        {
            pythia.ReadString("Tune:pp = 14");
            pythia.ReadString("ColourReconnection:mode = 1");
            pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
            pythia.ReadString("ColourReconnection:m0 = 0.3");
            pythia.ReadString("ColourReconnection:allowJunctions = on");
            pythia.ReadString("ColourReconnection:junctionCorrection = 1.20");
            pythia.ReadString("ColourReconnection:timeDilationMode = 2");
            pythia.ReadString("ColourReconnection:timeDilationPar = 0.18");
            pythia.ReadString("StringPT:sigma = 0.335");
            pythia.ReadString("StringZ:aLund = 0.36");
            pythia.ReadString("StringZ:bLund = 0.56");
            pythia.ReadString("StringFlav:probQQtoQ = 0.078");
            pythia.ReadString("StringFlav:ProbStoUD = 0.2");
            pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
            pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.ReadString("BeamRemnants:remnantMode = 1");
            pythia.ReadString("BeamRemnants:saturation =5");
            break;
        }
        case kCRMode3:
        {
            pythia.ReadString("Tune:pp = 14");
            pythia.ReadString("ColourReconnection:mode = 1");
            pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
            pythia.ReadString("ColourReconnection:m0 = 0.3");
            pythia.ReadString("ColourReconnection:allowJunctions = on");
            pythia.ReadString("ColourReconnection:junctionCorrection = 1.15");
            pythia.ReadString("ColourReconnection:timeDilationMode = 3");
            pythia.ReadString("ColourReconnection:timeDilationPar = 0.073");
            pythia.ReadString("StringPT:sigma = 0.335");
            pythia.ReadString("StringZ:aLund = 0.36");
            pythia.ReadString("StringZ:bLund = 0.56");
            pythia.ReadString("StringFlav:probQQtoQ = 0.078");
            pythia.ReadString("StringFlav:ProbStoUD = 0.2");
            pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
            pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.05");
            pythia.ReadString("BeamRemnants:remnantMode = 1");
            pythia.ReadString("BeamRemnants:saturation =5");
        }
        default:
        {
            std::cerr << "Invalid hadronisation tune, exit" << std::endl;
            return;
        }
    }

    //FIXME: currently all "prompt" particles have production radius of 0. Add tunes for parton vertex information? (see: http://home.thep.lu.se/~torbjorn/pythia83html/Welcome.html)

    // set energy
    double collEnergy = 13000.;
    pythia.ReadString(Form("Beams:eCM = %f", collEnergy));

    // initialize random seed
    pythia.ReadString("Random:setSeed = on");
    pythia.ReadString(Form("Random:seed %d", seed));
    pythia.Initialize(2212 /* p */, 2212 /* p */, collEnergy /* GeV */);

    // output tree
    float kStar;
    double radiusProton, radiusD;
    short originD, pdgProton, pdgD;
    TTree treeDp("treeDp", "D-p pairs");
    treeDp.Branch("k_star", &kStar);
    treeDp.Branch("R_D", &radiusD);
    treeDp.Branch("R_p", &radiusProton);
    treeDp.Branch("origin_D", &originD);
    treeDp.Branch("pdg_D", &pdgD);
    treeDp.Branch("pdg_p", &pdgProton);

    // Array of particles
    TClonesArray* particles = new TClonesArray("TParticle", 1000);

    for (int iEvent = 0; iEvent < nevents; ++iEvent)
    {
        pythia.GenerateEvent();
       
        pythia.ImportParticles(particles, "All");
        int nPart = particles->GetEntriesFast();

        std::vector<TParticle*> DpVec{};
        std::vector<TParticle*> pVec{};
        std::vector<short> originDpVec{};
        for (int iPart = 0; iPart < nPart; iPart++)
        {
            TParticle* part = dynamic_cast<TParticle*>(particles->At(iPart));
            int pdg = part->GetPdgCode();
            int absPdg = TMath::Abs(pdg);
            if(absPdg != 411 && absPdg != 2212) // we are only interested in protons and D+/-
                continue;

            if(TMath::Abs(part->Y()) > 1.) // let's look only at mid-rapidity
                continue;

            if(absPdg == 2212)
            {
                pVec.push_back(part);
                // Add selection of primary protons only? Tag their origin?
            }
            else if(absPdg == 411)
            {
                DpVec.push_back(part);
                // check origin from B and Dstar
                int motherIdx = part->GetFirstMother();
                bool isFromB = false, isFromDstar = false, isFromExcitedCharm = false, isQuarkFound = false;
                while(motherIdx>=0) {
                    TParticle* motherPart = dynamic_cast<TParticle*>(particles->At(motherIdx));
                    motherIdx = motherPart->GetFirstMother();
                    int absPdgMom = TMath::Abs(motherPart->GetPdgCode());
                    if(absPdgMom / 100 == 5 || absPdgMom / 1000 == 5 || (absPdgMom-10000) / 100 == 5 || (absPdgMom-20000) / 100 == 5 || (absPdgMom-100000) / 100 == 5 |(absPdgMom-200000) / 100 == 5)
                        isFromB = true;
                    if(absPdgMom == 413)
                        isFromDstar = true;
                    if(absPdgMom / 100 == 4 || absPdgMom / 1000 == 4 || (absPdgMom-10000) / 100 == 4 || (absPdgMom-20000) / 100 == 4 || (absPdgMom-100000) / 100 == 4 |(absPdgMom-200000) / 100 == 4)
                        isFromExcitedCharm = true;
                }

                if(isFromB)
                    originDpVec.push_back(kFromB);
                else if(isFromDstar)
                    originDpVec.push_back(kFromDstar);
                else if(isFromExcitedCharm) // do not keep track of these, they are anyway few
                    originDpVec.push_back(kFromExcitedCharm);
                else
                    originDpVec.push_back(kPrompt);
            }
        }

        for(size_t iD=0; iD<DpVec.size(); iD++)
        {
            for(size_t ip=0; ip<pVec.size(); ip++)
            {
                kStar = ComputeKstar(DpVec[iD], pVec[ip]);
                radiusProton = pVec[ip]->Rho();
                radiusD = DpVec[iD]->Rho();
                originD = originDpVec[iD];
                pdgProton = pVec[ip]->GetPdgCode();
                pdgD = DpVec[iD]->GetPdgCode();
                if(originD < kFromExcitedCharm) // set == kPrompt to keep only prompt D
                    treeDp.Fill();
            }
        }
    }

    // Output file
    std::string outFileName = "Dplus_proton_production_radius";
    switch(optHadronisation)
    {
        case kMonash13:
        {
            outFileName += "_Monash.root";
            break;
        }
        case kCRMode0:
        {
            outFileName += "_CR_mode0.root";
            break;
        }
        case kCRMode2:
        {
            outFileName += "_CR_mode2.root";
            break;
        }
        case kCRMode3:
        {
            outFileName += "_CR_mode3.root";
            break;
        }
    }
    TFile outFile(outFileName.data(), "recreate");
    treeDp.Write();
    outFile.Close();
}

double ComputeKstar(TParticle* part1, TParticle *part2)
{
    TLorentzVector SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS;
    SPtrack.SetXYZM(part1->Px(), part1->Py(), part1->Pz(),
                    TDatabasePDG::Instance()->GetParticle(abs(part1->GetPdgCode()))->Mass());
    TPProng.SetXYZM(part2->Px(), part2->Py(), part2->Pz(),
                    TDatabasePDG::Instance()->GetParticle(abs(part2->GetPdgCode()))->Mass());
    trackSum = SPtrack + TPProng;

    double beta = trackSum.Beta();
    double betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
    double betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
    double betaz = beta * cos(trackSum.Theta());

    SPtrackCMS = SPtrack;
    TPProngCMS = TPProng;

    SPtrackCMS.Boost(-betax, -betay, -betaz);
    TPProngCMS.Boost(-betax, -betay, -betaz);

    TLorentzVector trackRelK;

    trackRelK = SPtrackCMS - TPProngCMS;
    double kStar = 0.5 * trackRelK.P();
    return kStar;
}
