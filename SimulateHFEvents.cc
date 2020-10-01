#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <array>

#include "yaml-cpp/yaml.h"

#include <TPythia8.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TTree.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#endif

void SimulateHFEvents(TString cfgFileName)
{
    //Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    if (config.IsNull())
    {
        std::cerr << "\033[31mERROR: yaml config file not found! Exit\033[0m" << std::endl;
        return;
    }
    else 
    {
        std::cout << "\n\n*******************************************" << std::endl;
        std::cout << Form("\033[32mLoading configuration from file %s\033[0m\n", cfgFileName.Data()) << std::endl;
    }

    std::string flavour = config["flavour"].as<std::string>();
    int nevents = config["nevents"].as<int>();
    double collenergy = config["collenergy"].as<double>();
    int seed = config["seed"].as<int>();

    std::array<double, 2> ptLims = config["selections"]["pt"].as<std::array<double, 2> >();
    std::array<double, 2> yLims = config["selections"]["rapidity"].as<std::array<double, 2> >();

    int tune = config["tunes"]["tune"].as<int>();
    std::string process = config["tunes"]["process"].as<std::string>();
    int paperSkandsMode = config["tunes"]["paperSkandsMode"].as<int>();;

    int CRmode = config["tunes"]["advanced"]["CRmode"].as<int>();
    int qq1toQQ0charm = config["tunes"]["advanced"]["qq1toQQ0charm"].as<int>();
    int qq1toQQ0beauty = config["tunes"]["advanced"]["qq1toQQ0beauty"].as<int>();
    int timeDilationMode = config["tunes"]["advanced"]["timeDilationMode"].as<int>();
    int doSigmaLoop = config["tunes"]["advanced"]["doSigmaLoop"].as<int>();
    int doubleJunct = config["tunes"]["advanced"]["doubleJunct"].as<int>();
    int m0 = config["tunes"]["advanced"]["m0"].as<int>();

    std::string outFileName = config["output"]["filename"].as<std::string>();

    // create pythia generator
    TPythia8 pythia;

    // Set HardQCD or SoftQCD processes
    if (process == "HardQCDCharm")
    {
        pythia.ReadString("HardQCD:hardccbar = on ");
    }
    else if (process == "HardQCDBeauty")
    {
        pythia.ReadString("HardQCD:hardbbbar = on ");
    }
    else if (process == "HardQCDCharmBeauty")
    {
        pythia.ReadString("HardQCD:hardccbar = on ");
        pythia.ReadString("HardQCD:hardbbbar = on ");
    }
    else if (process == "SoftQCD")
    {
        pythia.ReadString("SoftQCD:all = on ");
    }
    else
    {
        std::cerr << "\033[31mERROR: invalid process! Exit\033[0m" << std::endl;
        return;  
    }

    pythia.ReadString(Form("Tune:pp = %d", tune));
    pythia.ReadString(Form("Beams:eCM = %f", collenergy));

    if (paperSkandsMode == -1)
        {
            if(CRmode>-1){
                pythia.ReadString(Form("ColourReconnection:mode = %d", CRmode));     
        }
            if(!doubleJunct){
                pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        }
        else {
            pythia.ReadString("ColourReconnection:allowDoubleJunRem = on");
        }
        if(qq1toQQ0charm>0 || qq1toQQ0beauty>0){
            //      if(qq1toQQ0beauty<0)qq1toQQ0beauty=27; // to match default value in CRmode 1 (BUT THIS IS THE DEFAULT IN THE PAPER, NOT THE DEFAULT IN THE SIMULATION)
            // if(qq1toQQ0charm<0)qq1toQQ0charm=27; // to match default value in CRmode 1 (BUT THIS IS THE DEFAULT IN THE PAPER, NOT THE DEFAULT IN THE SIMULATION)
            pythia.ReadString(Form("StringFlav:probQQ1toQQ0join = 0.027,0.027,%f,%f",(Float_t)(qq1toQQ0charm/1000.),(Float_t)(qq1toQQ0beauty/1000.)));// original: 0.027, 0.027 ,0.027 , 0.027 ; noeffects
            //      pythia.ReadString("ColourReconnection:m0=0.11");// Lc reduced with 3, no change with 0.11 (default is 0.3, min allowed 0.1)
        }
        if(timeDilationMode>-1){
            pythia.ReadString(Form("ColourReconnection:timeDilationMode=%d",timeDilationMode));
        }
        if(m0>-1){
            pythia.ReadString(Form("ColourReconnection:m0=%f",(Float_t)(m0/100.)));
        }
    }
    else
    {
        if (paperSkandsMode == 0)
        {
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
            if (timeDilationMode > -1)
            { // These settings should never be used overwrting the above ones
                pythia.ReadString(Form("ColourReconnection:timeDilationMode=%d", timeDilationMode));
            }
            if (CRmode > -1)
            { // These settings should never be used overwrting the above ones
                pythia.ReadString(Form("ColourReconnection:mode = %d", CRmode));
            }
            if (m0 > -1)
            {
                pythia.ReadString(Form("ColourReconnection:m0=%f", (float)(m0 / 100.)));
            }
        }
        else if (paperSkandsMode == 2)
        {
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
            if (timeDilationMode > -1)
            { // These settings should never be used overwrting the above ones
                pythia.ReadString(Form("ColourReconnection:timeDilationMode=%d", timeDilationMode));
            }
            if (CRmode > -1)
            { // These settings should never be used overwrting the above ones
                pythia.ReadString(Form("ColourReconnection:mode = %d", CRmode));
            }
            if (m0 > -1)
            {
                pythia.ReadString(Form("ColourReconnection:m0=%f", (float)(m0 / 100.)));
            }
        }
        else if (paperSkandsMode == 3)
        {
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
            if (timeDilationMode > -1)
            { // These settings should never be used overwrting the above ones
                pythia.ReadString(Form("ColourReconnection:timeDilationMode=%d", timeDilationMode));
            }
            if (CRmode > -1)
            { // These settings should never be used overwrting the above ones
                pythia.ReadString(Form("ColourReconnection:mode = %d", CRmode));
            }
            if (m0 > -1)
            {
                pythia.ReadString(Form("ColourReconnection:m0=%f", (float)(m0 / 100.)));
            }
        }
        else
        {
            std::cerr << "\033[31mERROR: invalid mode/tune! Exit\033[0m" << std::endl;
            return;
        }
    }

    // initialize random seed
    pythia.ReadString("Random:setSeed = on");
    pythia.ReadString(Form("Random:seed %d", seed));
    pythia.Initialize(2212 /* p */, 2212 /* p */, collenergy /* GeV */);

    // *********************************************************************************************************************************

    // Array of particles
    TClonesArray* particles = new TClonesArray("TParticle", 1000);

    // Create file on which histogram(s) can be saved.
    TFile outFile(outFileName.data(), "recreate");
    TTree *treeEvents = new TTree("treeEvents", "treeEvents");
    double pT = -1., y = -1., eta = -1.;
    int evNum = -1, pdg = -1, origin = -1;
    treeEvents->Branch("pT", &pT);
    treeEvents->Branch("y", &y);
    treeEvents->Branch("eta", &eta);
    treeEvents->Branch("pdg", &pdg);
    treeEvents->Branch("origin", &origin);
    treeEvents->Branch("event_number", &evNum);

    TH1F *hNevents = new TH1F("hNevents", "hNevents", 1, -0.5, 1.5);

    // Begin event loop. Generate event; skip if generation aborted.
    for (int iEvent = 0; iEvent < nevents; ++iEvent)
    {
        pythia.GenerateEvent();
        hNevents->Fill(0);
       
        pythia.ImportParticles(particles, "All");
        int nPart = particles->GetEntriesFast();

        // TODO: add multiplicity in output tree
        // for (int iPart = 0; iPart < nPart; iPart++)
        // {
        //     TParticle* part = (TParticle*) particles->At(iPart);
        //     int ist = part->GetStatusCode();
        //     // Positive codes are final particles.
        //     if (ist <= 0) continue;
        //     float charge = TDatabasePDG::Instance()->GetParticle(part->GetPdgCode())->Charge();
        //     if (charge == 0.) continue;
          
        //     ++nCharged;
        //     if (TMath::Abs(part->Eta()) < 0.5)
        //         ++nChargedEta05;
        // }

        int nb = 0;
        int nbbar = 0;

        for (int iPart = 0; iPart < nPart; iPart++)
        {
            TParticle* part = dynamic_cast<TParticle*>(particles->At(iPart));
            pT = part->Pt();
            y = part->Y();
            eta = part->Eta();

            // apply some selection
            if(pT < ptLims[0] || pT > ptLims[1])
                continue;
            if(y < yLims[0] || y > yLims[1])
                continue;

            pdg = part->GetPdgCode();
            int absPdg = TMath::Abs(pdg);

            // select charm (beauty) "stable" meson or baryon 
            bool isCharm = false;
            bool isBeauty = false;
            if(absPdg / 100 == 4 || absPdg / 1000 == 4)
                isCharm = true;
            else if(absPdg / 100 == 5 || absPdg / 1000 == 5)
                isBeauty = true;

            origin = 4; // origin for charm (4: prompt, 5: feed-down)

            if ((isCharm && (flavour == "charm" || flavour == "charmbeauty")) || (isBeauty && (flavour == "beauty" || flavour == "charmbeauty")))
            {   
                evNum = iEvent;

                // check origin (do not go up to quark)
                if(isCharm)
                {
                    int motherIdx = part->GetFirstMother();
                    while(motherIdx>=0) {
                        TParticle* motherPart = dynamic_cast<TParticle*>(particles->At(motherIdx));
                        motherIdx = motherPart->GetFirstMother();
                        int absPdgMom = TMath::Abs(motherPart->GetPdgCode());
                        if(absPdgMom == 5 || absPdgMom / 100 == 5 || absPdgMom / 1000 == 5 || (absPdgMom-10000) / 100 == 5 || (absPdgMom-20000) / 100 == 5 || (absPdgMom-100000) / 100 == 5 || (absPdgMom-200000) / 100 == 5)
                        {
                            origin = 5;
                            break;
                        }
                    }
                }

                treeEvents->Fill();
            }
        }
    }

    // Save histogram on file and close file.
    outFile.cd();
    treeEvents->Write();
    hNevents->Write();
    outFile.Close();
}
