#if !defined(__CINT__) || defined(__MAKECINT__)

#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

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
#include <TRandom3.h>
#include <TF1.h>
#include <TLorentzVector.h>

#endif

enum FONLLPred {
    kCentral, 
    kMin,
    kMax
};

const int kPdgHe = 1000020030;

//__________________________________________________________________________________________________
void ComputeLbtoHeDecays(TString cfgFileName);
bool SimpleCoalescence(double p1[3], double p2[3], double p3[3], double pHe3[3], double pc = 0.2 /* MeV/c */);
TH1F* ReadFONLLVsPt(string txtFileName, int whichPred, double ptMin, double binWidth, int nBins);

//__________________________________________________________________________________________________
void ComputeLbtoHeDecays(TString cfgFileName)
{
    //__________________________________________________________
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

    int nEvents = config["nevents"].as<int>();
    int tune = config["tunes"]["tune"].as<int>();
    int paperSkandsMode = config["tunes"]["paperSkandsMode"].as<int>();;

    int CRmode = config["tunes"]["advanced"]["CRmode"].as<int>();
    int qq1toQQ0charm = config["tunes"]["advanced"]["qq1toQQ0charm"].as<int>();
    int qq1toQQ0beauty = config["tunes"]["advanced"]["qq1toQQ0beauty"].as<int>();
    int timeDilationMode = config["tunes"]["advanced"]["timeDilationMode"].as<int>();
    int doubleJunct = config["tunes"]["advanced"]["doubleJunct"].as<int>();
    int m0 = config["tunes"]["advanced"]["m0"].as<int>();

    std::string FONLLFileName = config["FONLL"]["filename"].as<std::string>();

    std::string outFileName = config["output"]["filename"].as<std::string>();

    //__________________________________________________________
    // create and configure pythia generator
    TPythia8 pythia;
    pythia.ReadString("SoftQCD:all = on ");

    // set tune
    pythia.ReadString(Form("Tune:pp = %d", tune));
    if (paperSkandsMode == -1) // custom --> only if you are an expert user
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
            pythia.ReadString(Form("StringFlav:probQQ1toQQ0join = 0.027,0.027,%f,%f",(Float_t)(qq1toQQ0charm/1000.),(Float_t)(qq1toQQ0beauty/1000.)));
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

    // keep only interesting decays, to be reweighted a posteriori
    pythia.ReadString("5122:onMode = off");
    pythia.ReadString("5122:onIfMatch = 2 1 2 2101"); // bRatio="0.0120000" dominant one according to https://arxiv.org/pdf/2006.16251.pdf
    // pythia.ReadString("5122:onIfMatch = 2 1 4 2101"); // bRatio="0.4411147"
    // pythia.ReadString("5122:onIfMatch = 2 4 1 2101"); // bRatio="0.0910000"
    // pythia.ReadString("5122:onIfMatch = 2 3 2 2101"); // bRatio="0.0120000"
    // pythia.ReadString("5122:onIfMatch = 2 3 4 2101"); // bRatio="0.0800000"

    // init
    pythia.Pythia8()->init();

    //__________________________________________________________
    // perform the simulation

    // Output file
    TFile outFile(outFileName.data(), "recreate");

    TH1F* hBR = new TH1F("hBR", ";;BR", 2, 0.5, 2.5);
    hBR->GetXaxis()->SetBinLabel(1, "#Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}");
    hBR->GetXaxis()->SetBinLabel(2, "#Lambda_{b}^{0} #rightarrow (>=2)p + (>=1)n / #Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}");
    hBR->SetBinContent(1, 0.0120000);

    TH1F* hFONLLLb = ReadFONLLVsPt(FONLLFileName, kCentral, 0.025, 0.05, 1001);
    hFONLLLb->Scale(0.816); // f(b->B) from e+e- provides good normalisation for LHCb and CMS B-meson measurements
    hFONLLLb->SetName("hFONLLLb");

    // f(b -> Lb) / f(b -> B) from LHCb measurement https://arxiv.org/pdf/1902.06794.pdf
    TF1* fFFLHCb = new TF1("fracLb","([4] * ([5] + exp([6] + [7] * x))) /  (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1)  ", 0, 50);
    double parLbA = 1;
    double parLbp1 = 0.0793;
    double parLbp2 = -1.022;
    double parLbp3 = -0.107;
    double parBsA = 1;
    double parBsp1 = 0.119;
    double parBsp2 = -0.00091;
    double parBsAvePt = 10.1;
    fFFLHCb->SetParameters(parBsA, parBsp1, parBsp2, parBsAvePt, parLbA, parLbp1, parLbp2, parLbp3);

    for(int iPt=1; iPt<hFONLLLb->GetNbinsX()+1; iPt++)
    {
        double ptCent = hFONLLLb->GetBinCenter(iPt);
        hFONLLLb->SetBinContent(iPt, hFONLLLb->GetBinContent(iPt) * fFFLHCb->Eval(ptCent>5 ? ptCent:5));
    }

    TH1F* hHe3FromLb = (TH1F*)hFONLLLb->Clone("hHe3FromLb");
    hHe3FromLb->Reset();

    TTree *treeLb = new TTree("treeLb", "treeLb");
    std::vector<double> pxDau, pyDau, pzDau, ptDau;
    std::vector<int> pdgDau;
    double pxLb, pyLb, pzLb, ptLb;
    treeLb->Branch("ptLb", &ptLb);
    treeLb->Branch("pxLb", &pxLb);
    treeLb->Branch("pyLb", &pyLb);
    treeLb->Branch("pzLb", &pzLb);
    treeLb->Branch("ptDau", &ptDau);
    treeLb->Branch("pxDau", &pxDau);
    treeLb->Branch("pyDau", &pyDau);
    treeLb->Branch("pzDau", &pzDau);
    treeLb->Branch("pdgDau", &pdgDau);

    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    double massLb = TDatabasePDG::Instance()->GetParticle(5122)->Mass();
    int nEventSel = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        ptLb = hFONLLLb->GetRandom();
        double phiLb = gRandom->Rndm() * 2 * TMath::Pi();
        double yLb = gRandom->Rndm() * 2. - 1.; // flat in -1<y<1
        pxLb = ptLb * TMath::Cos(phiLb);
        pyLb = ptLb * TMath::Sin(phiLb);
        double mt = TMath::Sqrt(massLb * massLb + ptLb * ptLb);
        pzLb = mt * TMath::SinH(yLb);
        double pLb = TMath::Sqrt(ptLb * ptLb + pzLb * pzLb);
        double ELb = TMath::Sqrt(massLb * massLb + pLb * pLb);

        pythia.Pythia8()->event.clear();
        pythia.Pythia8()->event.append(5122, 11, 0, 0, pxLb, pyLb, pzLb, ELb, massLb);
        int idPart = pythia.Pythia8()->event[0].id();
        pythia.Pythia8()->particleData.mayDecay(idPart, true);
        pythia.Pythia8()->moreDecays();

        pythia.ImportParticles(particles, "All");
        int nPart = particles->GetEntriesFast();

        if(iEvent%1000000 == 0)
            std::cout << Form("Lb decay number %10d\r", iEvent) << std::endl;

        int nProtons = 0, nAntiProtons = 0, nNeutrons = 0, nAntiNeutrons = 0;

        for (int iPart = 1; iPart < nPart; iPart++)
        {
            TParticle* part = dynamic_cast<TParticle*>(particles->At(iPart));
            ptDau.push_back(part->Pt());
            pxDau.push_back(part->Px());
            pyDau.push_back(part->Py());
            pzDau.push_back(part->Pz());
            pdgDau.push_back(part->GetPdgCode());
            switch(pdgDau[iPart-1])
            {
                case 2212:
                    nProtons++;
                break;
                case -2212:
                    nAntiProtons++;
                break;
                case 2112:
                    nNeutrons++;
                break;
                case -2112:
                    nAntiNeutrons++;
                break;
            }
        }
        if((nProtons >= 2 and nNeutrons >= 1) || (nAntiProtons >= 2 and nAntiNeutrons >= 1))
        {
            double pCoal1[3], pCoal2[3], pCoal3[3];
            double pCoalHe3[3] = {-999., -999., -999.};
            bool isSecond = false;
            for(size_t iPart=0; iPart<pdgDau.size(); iPart++)
            {
                if(pdgDau[iPart] == 2212) {
                    if(!isSecond) {
                        pCoal1[0] = pxDau[iPart];
                        pCoal1[1] = pyDau[iPart];
                        pCoal1[2] = pzDau[iPart];
                        isSecond = true;
                    }
                    else {
                        pCoal2[0] = pxDau[iPart];
                        pCoal2[1] = pyDau[iPart];
                        pCoal2[2] = pzDau[iPart];
                    }
                }
                else if(pdgDau[iPart] == 2112) {
                    pCoal3[0] = pxDau[iPart];
                    pCoal3[1] = pyDau[iPart];
                    pCoal3[2] = pzDau[iPart];
                }
            }
            bool hasCoalesced = SimpleCoalescence(pCoal1, pCoal2, pCoal3, pCoalHe3);
            if(hasCoalesced) {
                double ptHe3 = TMath::Sqrt(pCoalHe3[0]*pCoalHe3[0] + pCoalHe3[1]*pCoalHe3[1]);
                hHe3FromLb->Fill(ptHe3);
            }
            nEventSel++;
            treeLb->Fill();
        }

        ptDau.clear();
        pxDau.clear();
        pyDau.clear();
        pzDau.clear();
        pdgDau.clear();
    }

    hBR->SetBinContent(2, static_cast<double>(nEventSel) / nEvents);
    hHe3FromLb->Scale(hFONLLLb->Integral() / hHe3FromLb->Integral() * hBR->GetBinContent(1) * hBR->GetBinContent(2));

    // Save histogram on file and close file.
    outFile.cd();
    hFONLLLb->Write();
    hHe3FromLb->Write();
    hBR->Write();
    treeLb->Write();
    outFile.Close();
}

//__________________________________________________________________________________________________
TH1F* ReadFONLLVsPt(string txtFileName, int whichPred, double ptMin, double binWidth, int nBins)
{
    if(txtFileName.find("txt") == string::npos && txtFileName.find("dat") == string::npos && txtFileName.find("csv") == string::npos) {
        std::cerr << "ERROR: Wrong file format! Exit." << std::endl;
        return nullptr;
    }

    std::ifstream inSet(txtFileName.data());
    if(!inSet) {
        std::cerr << "ERROR: Please check if "<< txtFileName.data() << " is the right path. Exit." << std::endl;
        return nullptr;
    }

    std::vector<std::string> values;
    
    TH1F* hBvsPt = new TH1F("hBvsPt", "", nBins, ptMin, ptMin + nBins*binWidth);
    double centY = -1., minY = -1., maxY = -1.;
    int iPt = 1;
    while(!inSet.eof())
    {    
        std::vector<std::string> values;
        std::string line;
        std::getline(inSet, line);
        if(line.find("#") != std::string::npos || line.find("pt") != std::string::npos)
            continue;

        size_t pos = 0;
        while((pos = line.find(' ')) != std::string::npos)
        {
            values.push_back(line.substr(0, pos));
            if(values.size()==2) {
                std::stringstream convert(values[1]);
                if( !(convert >> centY) )
                    centY = -1;
            }
            else if(values.size()==3) {
                std::stringstream convert(values[2]);
                if( !(convert >> minY) )
                    minY = -1;
            }
            else if(values.size()==4) {
                std::stringstream convert(values[3]);
                if( !(convert >> maxY) )
                    maxY = -1;
            }
            line = line.substr(pos + 1);
        }
        if(whichPred == kCentral)
            hBvsPt->SetBinContent(iPt, centY);
        else if(whichPred == kMin)
            hBvsPt->SetBinContent(iPt, minY);
        else if(whichPred == kMax)
            hBvsPt->SetBinContent(iPt, maxY);
        iPt++;
    }

    inSet.close();

    return hBvsPt;
}

//__________________________________________________________________________________________________
bool SimpleCoalescence(double p1[3], double p2[3], double p3[3], double pHe3[3], double pc) {

    // returns true if coalescence is realised, false otherwise
    TLorentzVector Ptrack1, Ptrack2, Ptrack3, trackSum, Ptrack1CMS, Ptrack2CMS, Ptrack3CMS;
    Ptrack1.SetXYZM(p1[0], p1[1], p1[2], TDatabasePDG::Instance()->GetParticle(2212)->Mass());
    Ptrack2.SetXYZM(p2[0], p2[1], p2[2], TDatabasePDG::Instance()->GetParticle(2212)->Mass());
    Ptrack3.SetXYZM(p3[0], p3[1], p3[2], TDatabasePDG::Instance()->GetParticle(2112)->Mass());
    trackSum = Ptrack1 + Ptrack2 + Ptrack3;

    double beta = trackSum.Beta();
    double betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
    double betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
    double betaz = beta * cos(trackSum.Theta());

    Ptrack1CMS = Ptrack1;
    Ptrack2CMS = Ptrack2;
    Ptrack3CMS = Ptrack3;

    Ptrack1CMS.Boost(-betax, -betay, -betaz);
    Ptrack2CMS.Boost(-betax, -betay, -betaz);
    Ptrack3CMS.Boost(-betax, -betay, -betaz);

    double coalRadius = TMath::Power(2, 1./6) * pc / 2;
    if(Ptrack1CMS.P() <= coalRadius && Ptrack2CMS.P() <= coalRadius && Ptrack3CMS.P() <= coalRadius)
    {
        for(int iEl=0; iEl<3; iEl++)
            pHe3[iEl] = p1[iEl] + p2[iEl] + p3[iEl];
        return true;
    }

    for(int iEl=0; iEl<3; iEl++)
        pHe3[0] = -1.;
    return false;
}
