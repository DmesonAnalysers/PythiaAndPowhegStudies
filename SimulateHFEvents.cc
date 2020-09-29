#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <TPythia8.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TTree.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#endif

enum 
{
    kCharm,
    kBeauty,
    kCharmAndBeauty
};

int SimulateHFEvents(int nevents=1000000, int flavour=kCharm)
{
    // *********************************************************************************************************************************
    // Configure Pythia

    // Create the ROOT application environment.
    int seed = 0;
    int tune = 14; // 4=2M, 5=4C, 6=4Cx, 14=Monash
    int process = 0;
    int CRmode = -1;
    int paperSkandsMode = -1;
    int qq1toQQ0charm = -1;
    int qq1toQQ0beauty = -1;
    int timeDilationMode = -1;
    int doSigmaLoop = 0;
    int doubleJunct = 1;
    int m0 = -1;
    int collenergy = 5020;

    TPythia8 *pythia = new TPythia8();

    // Set HardQCD or SoftQCD processes
    if (process == 4)
    {
        pythia->ReadString("HardQCD:hardccbar = on ");
    }
    else if (process == 5)
    {
        pythia->ReadString("HardQCD:hardbbbar = on ");
    }
    else if (process == 45)
    {
        pythia->ReadString("HardQCD:hardccbar = on ");
        pythia->ReadString("HardQCD:hardbbbar = on ");
    }
    else if (process == 0)
    {
        pythia->ReadString("SoftQCD:all = on ");
    }

    pythia->ReadString(Form("Tune:pp = %d", tune));
    pythia->ReadString(Form("Beams:eCM = %f", (float)collenergy));

    if (paperSkandsMode == -1)
        {
            if(CRmode>-1){
                pythia->ReadString(Form("ColourReconnection:mode = %d",CRmode));     
        }
            if(!doubleJunct){
                pythia->ReadString("ColourReconnection:allowDoubleJunRem = off");
        }
        else {
            pythia->ReadString("ColourReconnection:allowDoubleJunRem = on");
        }
        if(qq1toQQ0charm>0 || qq1toQQ0beauty>0){
            //      if(qq1toQQ0beauty<0)qq1toQQ0beauty=27; // to match default value in CRmode 1 (BUT THIS IS THE DEFAULT IN THE PAPER, NOT THE DEFAULT IN THE SIMULATION)
            // if(qq1toQQ0charm<0)qq1toQQ0charm=27; // to match default value in CRmode 1 (BUT THIS IS THE DEFAULT IN THE PAPER, NOT THE DEFAULT IN THE SIMULATION)
            pythia->ReadString(Form("StringFlav:probQQ1toQQ0join = 0.027,0.027,%f,%f",(Float_t)(qq1toQQ0charm/1000.),(Float_t)(qq1toQQ0beauty/1000.)));// original: 0.027, 0.027 ,0.027 , 0.027 ; noeffects
            //      pythia->ReadString("ColourReconnection:m0=0.11");// Lc reduced with 3, no change with 0.11 (default is 0.3, min allowed 0.1)
        }
        if(timeDilationMode>-1){
            pythia->ReadString(Form("ColourReconnection:timeDilationMode=%d",timeDilationMode));
        }
        if(m0>-1){
            pythia->ReadString(Form("ColourReconnection:m0=%f",(Float_t)(m0/100.)));
        }
    }
    else
    {
        if (paperSkandsMode == 0)
        {
            pythia->ReadString("ColourReconnection:mode = 1");
            pythia->ReadString("ColourReconnection:allowDoubleJunRem = off");
            pythia->ReadString("ColourReconnection:m0 = 2.9");
            pythia->ReadString("ColourReconnection:allowJunctions = on");
            pythia->ReadString("ColourReconnection:junctionCorrection = 1.43");
            pythia->ReadString("ColourReconnection:timeDilationMode = 0");
            pythia->ReadString("StringPT:sigma = 0.335");
            pythia->ReadString("StringZ:aLund = 0.36");
            pythia->ReadString("StringZ:bLund = 0.56");
            pythia->ReadString("StringFlav:probQQtoQ = 0.078");
            pythia->ReadString("StringFlav:ProbStoUD = 0.2");
            pythia->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
            pythia->ReadString("MultiPartonInteractions:pT0Ref = 2.12");
            pythia->ReadString("BeamRemnants:remnantMode = 1");
            pythia->ReadString("BeamRemnants:saturation =5");
            if (timeDilationMode > -1)
            { // These settings should never be used overwrting the above ones
                pythia->ReadString(Form("ColourReconnection:timeDilationMode=%d", timeDilationMode));
            }
            if (CRmode > -1)
            { // These settings should never be used overwrting the above ones
                pythia->ReadString(Form("ColourReconnection:mode = %d", CRmode));
            }
            if (m0 > -1)
            {
                pythia->ReadString(Form("ColourReconnection:m0=%f", (float)(m0 / 100.)));
            }
        }
        else if (paperSkandsMode == 2)
        {
            pythia->ReadString("ColourReconnection:mode = 1");
            pythia->ReadString("ColourReconnection:allowDoubleJunRem = off");
            pythia->ReadString("ColourReconnection:m0 = 0.3");
            pythia->ReadString("ColourReconnection:allowJunctions = on");
            pythia->ReadString("ColourReconnection:junctionCorrection = 1.20");
            pythia->ReadString("ColourReconnection:timeDilationMode = 2");
            pythia->ReadString("ColourReconnection:timeDilationPar = 0.18");
            pythia->ReadString("StringPT:sigma = 0.335");
            pythia->ReadString("StringZ:aLund = 0.36");
            pythia->ReadString("StringZ:bLund = 0.56");
            pythia->ReadString("StringFlav:probQQtoQ = 0.078");
            pythia->ReadString("StringFlav:ProbStoUD = 0.2");
            pythia->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
            pythia->ReadString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia->ReadString("BeamRemnants:remnantMode = 1");
            pythia->ReadString("BeamRemnants:saturation =5");
            if (timeDilationMode > -1)
            { // These settings should never be used overwrting the above ones
                pythia->ReadString(Form("ColourReconnection:timeDilationMode=%d", timeDilationMode));
            }
            if (CRmode > -1)
            { // These settings should never be used overwrting the above ones
                pythia->ReadString(Form("ColourReconnection:mode = %d", CRmode));
            }
            if (m0 > -1)
            {
                pythia->ReadString(Form("ColourReconnection:m0=%f", (float)(m0 / 100.)));
            }
        }
        else if (paperSkandsMode == 3)
        {
            pythia->ReadString("ColourReconnection:mode = 1");
            pythia->ReadString("ColourReconnection:allowDoubleJunRem = off");
            pythia->ReadString("ColourReconnection:m0 = 0.3");
            pythia->ReadString("ColourReconnection:allowJunctions = on");
            pythia->ReadString("ColourReconnection:junctionCorrection = 1.15");
            pythia->ReadString("ColourReconnection:timeDilationMode = 3");
            pythia->ReadString("ColourReconnection:timeDilationPar = 0.073");
            pythia->ReadString("StringPT:sigma = 0.335");
            pythia->ReadString("StringZ:aLund = 0.36");
            pythia->ReadString("StringZ:bLund = 0.56");
            pythia->ReadString("StringFlav:probQQtoQ = 0.078");
            pythia->ReadString("StringFlav:ProbStoUD = 0.2");
            pythia->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
            pythia->ReadString("MultiPartonInteractions:pT0Ref = 2.05");
            pythia->ReadString("BeamRemnants:remnantMode = 1");
            pythia->ReadString("BeamRemnants:saturation =5");
            if (timeDilationMode > -1)
            { // These settings should never be used overwrting the above ones
                pythia->ReadString(Form("ColourReconnection:timeDilationMode=%d", timeDilationMode));
            }
            if (CRmode > -1)
            { // These settings should never be used overwrting the above ones
                pythia->ReadString(Form("ColourReconnection:mode = %d", CRmode));
            }
            if (m0 > -1)
            {
                pythia->ReadString(Form("ColourReconnection:m0=%f", (float)(m0 / 100.)));
            }
        }
        else
        {
            std::cerr << "Mode not yet implemented" << std::endl;
            return 0;
        }
    }

    // initialize random seed
    pythia->ReadString("Random:setSeed = on");
    pythia->ReadString(Form("Random:seed %d", seed));
    pythia->Initialize(2212 /* p */, 2212 /* p */, collenergy /* TeV */);

    // *********************************************************************************************************************************

    // Array of particles
    TClonesArray* particles = new TClonesArray("TParticle", 1000);

    // Create file on which histogram(s) can be saved.
    TString outFileName = "";
    if(flavour == kCharm)
        outFileName = "Cevents.root";
    else if(flavour == kBeauty)
        outFileName = "Bevents.root";
    else if(flavour == kCharmAndBeauty)
        outFileName = "CandBevents.root";
    TFile *outFile = new TFile(outFileName.Data(), "recreate");
    TTree *treeEvents = new TTree("treeEvents", "treeEvents");
    double pT = -1., y = -1., eta = -1.;
    int evNum = -1, pdg = -1, label = -1, mother = -1;
    treeEvents->Branch("pT", &pT);
    treeEvents->Branch("y", &y);
    treeEvents->Branch("eta", &eta);
    treeEvents->Branch("pdg", &pdg);
    treeEvents->Branch("label", &label);
    treeEvents->Branch("mother", &label);
    treeEvents->Branch("event_number", &evNum);

    TH1F *hNevents = new TH1F("hNevents", "hNevents", 1, -0.5, 1.5);

    // Begin event loop. Generate event; skip if generation aborted.
    TStopwatch t;
    t.Start();

    for (int iEvent = 0; iEvent < nevents; ++iEvent)
    {
        pythia->GenerateEvent();
        hNevents->Fill(0);
        if (iEvent % 10000 == 0)
        {
            std::cout << Form("time: CPU %f real %f", t.CpuTime(), t.RealTime()) << std::endl;
            t.Start(false);
        }
       
        pythia->ImportParticles(particles, "All");
        int nPart = particles->GetEntriesFast();

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
            TParticle* part = (TParticle*) particles->At(iPart);
            pdg = part->GetPdgCode();
            int absPdg = TMath::Abs(pdg);
            
            if ((flavour == kCharm && (absPdg == 4 || absPdg / 100 == 4 || absPdg / 1000 == 4)) || 
                (flavour == kBeauty && (absPdg == 5 || absPdg / 100 == 5 || absPdg / 1000 == 5)) ||
                (flavour == kCharmAndBeauty && (absPdg == 4 || absPdg / 100 == 4 || absPdg / 1000 == 4 || absPdg == 5 || absPdg / 100 == 5 || absPdg / 1000 == 5)))
            { // charm (beauty) quark or charm (beauty) "stable" meson or baryon                
                pT = part->Pt();
                y = part->Y();
                eta = part->Eta();
                evNum = iEvent;
                label = iPart;
                mother = part->GetFirstMother();
                treeEvents->Fill();
            }
        }
    }

    // Save histogram on file and close file.
    outFile->cd();
    treeEvents->Write();
    hNevents->Write();
    outFile->Close();

    return 0;
}
