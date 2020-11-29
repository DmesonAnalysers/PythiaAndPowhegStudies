#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TClonesArray.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TNtuple.h>

#include "TPythia6.h"
#include "AliTPythia8.h"
#include "AliDecayerPythia8.h"
#include "TPythia6Decayer.h"
#include "AliGenEvtGen.h"
#include "AliDecayerEvtGen.h"

#endif

enum dec
{
    kPythia6,
    kPythia8,
    kEvtGen
};

//_______________________________________________________________________________________________
// FUNCTION PROTOTYPES
void ComputeDplusFromDstarCrossSec(int nGen = 10000000, int decayer = kPythia8);
double PowLaw(double *pt, double *pars);

//_______________________________________________________________________________________________
// FUNCTION IMPLEMENTATIONS
void ComputeDplusFromDstarCrossSec(int nGen, int decayer)
{
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadTopMargin(0.035);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetTitleOffset(1.2, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleSize(0.045, "xy");
    gStyle->SetLabelSize(0.04, "xy");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kRainBow);

    // Define BRs (PDG2020)
    double BRDstarToDplus = 0.307;      // D*+ -> D+ pi0
    double BRDstar = 0.677 * 0.03951; // D*+ -> D0pi+ -> K- pi+ pi+
    double BRDplus = 0.0938;          // D+ -> K- pi+ pi+

    double mDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();

    // load D*+ and D+ cross sections
    TFile *inFileDstar = TFile::Open("crosssections/CrossSecDstar_pp5TeV_published.root");
    TH1F *hDstar = static_cast<TH1F *>(inFileDstar->Get("histoSigmaCorr"));
    hDstar->Scale(1. / BRDstar);
    hDstar->SetName("hDstar");
    hDstar->SetDirectory(0);
    hDstar->SetStats(0);
    hDstar->SetLineWidth(2);
    hDstar->SetLineColor(kAzure + 4);
    hDstar->SetMarkerColor(kAzure + 4);
    hDstar->SetMarkerStyle(kFullCircle);
    TGraphAsymmErrors* gDstar = static_cast<TGraphAsymmErrors *>(inFileDstar->Get("gSigmaCorr"));
    TGraphAsymmErrors* gDstarFD = static_cast<TGraphAsymmErrors *>(inFileDstar->Get("gSigmaCorrConservative"));
    for(int iPt=0; iPt<gDstar->GetN(); iPt++)
    {
        double pt =  -1., sigma = -1.;
        gDstar->GetPoint(iPt, pt, sigma);
        gDstar->SetPoint(iPt, pt, sigma / BRDstar);
        double systUncLow = TMath::Sqrt(gDstar->GetErrorYlow(iPt)*gDstar->GetErrorYlow(iPt) + gDstarFD->GetErrorYlow(iPt)*gDstarFD->GetErrorYlow(iPt));
        double systUncHigh = TMath::Sqrt(gDstar->GetErrorYhigh(iPt)*gDstar->GetErrorYhigh(iPt) + gDstarFD->GetErrorYhigh(iPt)*gDstarFD->GetErrorYhigh(iPt));
        gDstar->SetPointError(iPt, gDstar->GetErrorXlow(iPt), gDstar->GetErrorXhigh(iPt), systUncLow / BRDstar, systUncHigh / BRDstar);
    }
    gDstar->SetName("gDstar");
    gDstar->SetLineWidth(2);
    gDstar->SetLineColor(kAzure + 4);
    gDstar->SetMarkerColor(kAzure + 4);
    gDstar->SetFillColorAlpha(kAzure + 4, 0.2);
    inFileDstar->Close();
    TFile *inFileDplus = TFile::Open("crosssections/CrossSecDplusML_pp5TeV.root");
    TH1F *hDplus = static_cast<TH1F *>(inFileDplus->Get("histoSigmaCorr"));
    hDplus->Scale(1. / BRDplus);
    hDplus->SetName("hDplus");
    hDplus->SetStats(0);
    hDplus->SetDirectory(0);
    hDplus->SetLineWidth(2);
    hDplus->SetLineColor(kGreen + 2);
    hDplus->SetMarkerColor(kGreen + 2);
    hDplus->SetMarkerStyle(kFullCircle);
    TGraphAsymmErrors* gDplus = static_cast<TGraphAsymmErrors *>(inFileDplus->Get("gSigmaCorr"));
    TGraphAsymmErrors* gDplusFD = static_cast<TGraphAsymmErrors *>(inFileDplus->Get("gSigmaCorrConservative"));
    for(int iPt=0; iPt<gDplus->GetN(); iPt++)
    {
        double pt =  -1., sigma = -1.;
        gDplus->GetPoint(iPt, pt, sigma);
        gDplus->SetPoint(iPt, pt, sigma / BRDplus);
        double systUncLow = TMath::Sqrt(gDplus->GetErrorYlow(iPt)*gDplus->GetErrorYlow(iPt) + gDplusFD->GetErrorYlow(iPt)*gDplusFD->GetErrorYlow(iPt));
        double systUncHigh = TMath::Sqrt(gDplus->GetErrorYhigh(iPt)*gDplus->GetErrorYhigh(iPt) + gDplusFD->GetErrorYhigh(iPt)*gDplusFD->GetErrorYhigh(iPt));
        gDplus->SetPointError(iPt, gDplus->GetErrorXlow(iPt), gDplus->GetErrorXhigh(iPt), systUncLow / BRDplus, systUncHigh / BRDplus);
    }
    gDplus->SetName("gDplus");
    gDplus->SetLineWidth(2);
    gDplus->SetLineColor(kGreen + 2);
    gDplus->SetMarkerColor(kGreen + 2);
    gDplus->SetFillColorAlpha(kGreen + 2, 0.2);
    inFileDplus->Close();
    TH1F* hDplusFromDstar = static_cast<TH1F *>(hDplus->Clone("hDplusFromDstar"));
    hDplusFromDstar->Reset();
    hDplusFromDstar->SetLineColor(kRed+1);
    hDplusFromDstar->SetMarkerColor(kRed+1);
    TH1F* hDplusFromDstarSys[2] = {static_cast<TH1F *>(hDplus->Clone("hDplusFromDstarSysLow")), static_cast<TH1F *>(hDplus->Clone("hDplusFromDstarSysHigh"))};
    TGraphAsymmErrors* gDplusFromDstar = static_cast<TGraphAsymmErrors *>(gDplus->Clone("gDplusFromDstar"));
    gDplusFromDstar->SetName("gDplusFromDstar");
    gDplusFromDstar->SetLineColor(kRed + 1);
    gDplusFromDstar->SetMarkerColor(kRed + 1);
    gDplusFromDstar->SetFillColorAlpha(kRed + 1, 0.2);

    TVirtualMCDecayer *pdec = nullptr;
    if (decayer == kPythia6)
    {
        gSystem->Load("liblhapdf.so");
        gSystem->Load("libEGPythia6.so");
        gSystem->Load("libpythia6.so");
        pdec = new TPythia6Decayer();
    }
    else if (decayer == kPythia8)
    {
        gSystem->Load("liblhapdf.so");
        gSystem->Load("libpythia8.so");
        gSystem->Load("libAliPythia8.so");
        gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
        gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
        gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
        pdec = new AliDecayerPythia8();
    }
    else if (decayer == kEvtGen)
    {
        gSystem->Load("liblhapdf.so");
        gSystem->Load("libpythia8.so");
        gSystem->Load("libAliPythia8.so");
        gSystem->Load("libPhotos.so");
        gSystem->Load("libEvtGen.so");
        gSystem->Load("libEvtGenExternal.so");
        gSystem->Load("libTEvtGen.so");
        pdec = new AliDecayerEvtGen();
    }
    pdec->Init();

    TClonesArray array = TClonesArray("TParticle", 100);
    TClonesArray arraySys[2] = {TClonesArray("TParticle", 100), TClonesArray("TParticle", 100)};
    TLorentzVector vec = TLorentzVector();
    TLorentzVector vecSys[2] = {TLorentzVector(), TLorentzVector()};

    TH2F* hPtCorr = new TH2F("hPtCorr", ";#it{p}_{T} (GeV/#it{c}) (D*^{+});#it{p}_{T} (GeV/#it{c}) (D^{+})", 500, 0., 50., 500, 0., 50.);

    // simulate the decays
    TString decName = "";
    if(decayer == kPythia6)
        decName = "Pythia6";
    else if(decayer == kPythia8)
        decName = "Pythia8";
    else if(decayer == kEvtGen)
        decName = "EvtGen";

    TFile outFile(Form("DplusFromDstar_pp5TeV_decayer_%s.root", decName.Data()), "recreate");
    TDirectoryFile dir = TDirectoryFile("SmearDistr", "SmearDistr");
    dir.Write();

    TF1 *fPowLawDstar = new TF1("fPowLawDstar", PowLaw, 0, 50, 4);
    TF1 *fPowLawDstarSys[2] = {new TF1("fPowLawDstarSysLow", PowLaw, 0, 50, 4), new TF1("fPowLawDstarSysHigh", PowLaw, 0, 50, 4)};
    TH1F* hDstarSmeared = static_cast<TH1F*>(hDstar->Clone("hDstarSmeared"));
    TH1F* hDstarSys[2] = {static_cast<TH1F*>(hDstar->Clone("hDstarSysLow")), static_cast<TH1F*>(hDstar->Clone("hDstarSysHigh"))};
    TNtuple* nDplusFromDstar = new TNtuple("nDplusFromDstar", "nDplusFromDstar", "pT:DplusFromDstarCrossSec");
    // repeat procedure 100 times smearing the input cross section to consider statistical uncertainty
    const int nTrials = 100;
    for(int iTrial = 0; iTrial < nTrials; iTrial++)
    {
        TString trialString = Form("%03d", iTrial);
        std::cout << "Smearing of input Dstar cross section no " << trialString << "\r" << std::flush;
        hDstarSmeared->Reset();
        for(int iPt = 1; iPt < hDstarSmeared->GetNbinsX()+1; iPt++)
        {
            if(iTrial < nTrials-1)
            {
                hDstarSmeared->SetBinContent(iPt, gRandom->Gaus(hDstar->GetBinContent(iPt), hDstar->GetBinError(iPt)));
                hDstarSmeared->SetBinError(iPt, hDstar->GetBinError(iPt));
            }
            else
            {
                hDstarSmeared->SetBinContent(iPt, hDstar->GetBinContent(iPt));
                hDstarSmeared->SetBinError(iPt, hDstar->GetBinError(iPt));
                hDstarSys[0]->SetBinContent(iPt, hDstar->GetBinContent(iPt)-gDstar->GetErrorYlow(iPt));
                hDstarSys[0]->SetBinError(iPt, hDstar->GetBinError(iPt));
                hDstarSys[1]->SetBinContent(iPt, hDstar->GetBinContent(iPt)+gDstar->GetErrorYhigh(iPt));
                hDstarSys[1]->SetBinError(iPt, hDstar->GetBinError(iPt));
            }
        }
        // parametrise cross section with power-law function
        fPowLawDstar->SetParameters(hDstarSmeared->Integral(), 2.6, 2.8, 2.);
        fPowLawDstar->SetLineColor(kAzure + 4);
        hDstarSmeared->Fit("fPowLawDstar", "LIME0Q");
        if(iTrial == nTrials-1) // evaluate syst unc
        {
            for(int iSys = 0; iSys < 2; iSys++)
            {
                fPowLawDstarSys[iSys]->SetParameters(hDstarSys[iSys]->Integral(), 2.6, 2.8, 2.);
                hDstarSys[iSys]->Fit(fPowLawDstarSys[iSys], "LIME0Q");
            }
        }

        hDplusFromDstar->Reset();

        // simulate D* -> D+ pi0 decay
        double massDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
        for (int iGen = 0; iGen < nGen; iGen++)
        {
            double ptDstar = fPowLawDstar->GetRandom();
            double phiDstar = gRandom->Rndm() * 2 * TMath::Pi();
            double yDstar = gRandom->Rndm() * 2. - 1.; // flat in -1<y<1
            double px = ptDstar * TMath::Cos(phiDstar);
            double py = ptDstar * TMath::Sin(phiDstar);
            double mt = TMath::Sqrt(massDstar * massDstar + ptDstar * ptDstar);
            double pz = mt * TMath::SinH(yDstar);
            double pDstar = TMath::Sqrt(ptDstar * ptDstar + pz * pz);
            double E = TMath::Sqrt(massDstar * massDstar + pDstar * pDstar);
            vec.SetPxPyPzE(px, py, pz, E);

            pdec->Decay(413, &vec);
            int nEntries = pdec->ImportParticles(&array);

            for (int iPart = 0; iPart < nEntries; iPart++)
            {
                TParticle *part = static_cast<TParticle*>(array.At(iPart));
                int pdgDau = TMath::Abs(part->GetPdgCode());
                if (pdgDau == 411)
                {
                    double ptDplus = part->Pt();
                    double yDplus = part->Y();
                    if(TMath::Abs(yDplus) < 0.5) {
                        hDplusFromDstar->Fill(ptDplus);
                        hPtCorr->Fill(ptDstar, ptDplus);
                    }
                }
            }
        }
        if(iTrial == nTrials-1)
        {
            for(int iSys = 0; iSys < 2; iSys++)
            {
                hDplusFromDstarSys[iSys]->Reset();
                for (int iGen = 0; iGen < nGen; iGen++)
                {
                    double ptDstar = fPowLawDstarSys[iSys]->GetRandom();
                    double phiDstar = gRandom->Rndm() * 2 * TMath::Pi();
                    double yDstar = gRandom->Rndm() * 2. - 1.; // flat in -1<y<1
                    double px = ptDstar * TMath::Cos(phiDstar);
                    double py = ptDstar * TMath::Sin(phiDstar);
                    double mt = TMath::Sqrt(massDstar * massDstar + ptDstar * ptDstar);
                    double pz = mt * TMath::SinH(yDstar);
                    double pDstar = TMath::Sqrt(ptDstar * ptDstar + pz * pz);
                    double E = TMath::Sqrt(massDstar * massDstar + pDstar * pDstar);
                    vecSys[iSys].SetPxPyPzE(px, py, pz, E);
                    pdec->Decay(413, &vecSys[iSys]);
                    int nEntries = pdec->ImportParticles(&arraySys[iSys]);
                    for (int iPart = 0; iPart < nEntries; iPart++)
                    {
                        TParticle *part = static_cast<TParticle*>(arraySys[iSys].At(iPart));
                        int pdgDau = TMath::Abs(part->GetPdgCode());
                        if (pdgDau == 411)
                        {
                            double ptDplus = part->Pt();
                            double yDplus = part->Y();
                            if(TMath::Abs(yDplus) < 0.5) {
                                hDplusFromDstarSys[iSys]->Fill(ptDplus);
                            }
                        }
                    }
                }
            }
        }

        hDplusFromDstar->Scale(fPowLawDstar->Integral(0, 100) / hDplusFromDstar->Integral() * BRDstarToDplus, "width");
        if(iTrial == nTrials-1)
        {
            for(int iSys = 0; iSys < 2; iSys++)
                hDplusFromDstarSys[iSys]->Scale(fPowLawDstarSys[iSys]->Integral(0, 100) / hDplusFromDstarSys[iSys]->Integral() * BRDstarToDplus, "width");
            for(int iPt = 1; iPt < hDplusFromDstar->GetNbinsX()+1; iPt++)
            {
                gDplusFromDstar->SetPoint(iPt, hDplusFromDstar->GetBinCenter(iPt), hDplusFromDstar->GetBinContent(iPt));
                double sysUncLow = hDplusFromDstar->GetBinContent(iPt) - hDplusFromDstarSys[0]->GetBinContent(iPt);
                double sysUncHigh = hDplusFromDstarSys[1]->GetBinContent(iPt) - hDplusFromDstar->GetBinContent(iPt);
                if(sysUncLow < 0)
                    sysUncLow = sysUncHigh;
                if(sysUncHigh < 0)
                    sysUncHigh = sysUncLow;
                gDplusFromDstar->SetPointError(iPt, hDplusFromDstar->GetBinWidth(iPt)/2, hDplusFromDstar->GetBinWidth(iPt)/2, sysUncLow, sysUncHigh);
            }
        }

        if(iTrial < nTrials-1)
        {
            for(int iPt = 1; iPt < hDplusFromDstar->GetNbinsX()+1; iPt++)
            {
                float array4Ntuple[2] = {static_cast<float>(hDplusFromDstar->GetBinCenter(iPt)), static_cast<float>(hDplusFromDstar->GetBinContent(iPt))};
                nDplusFromDstar->Fill(array4Ntuple);
            }
        }
    }

    std::cout << endl;

    for(int iPt = 1; iPt < hDplusFromDstar->GetNbinsX()+1; iPt++)
    {
        TH1F* hDistr = new TH1F(Form("hDistrPt%d", iPt), "", 100, hDplusFromDstar->GetBinContent(iPt)*0.1, hDplusFromDstar->GetBinContent(iPt)*3);
        nDplusFromDstar->Project(Form("hDistrPt%d", iPt), "DplusFromDstarCrossSec", Form("pT > %f && pT < %f", hDplusFromDstar->GetBinLowEdge(iPt), hDplusFromDstar->GetXaxis()->GetBinUpEdge(iPt)));
        dir.cd();
        hDistr->Write();
        hDplusFromDstar->SetBinError(iPt, hDistr->GetRMS());
    }

    TH1F* hFrac = static_cast<TH1F *>(hDplusFromDstar->Clone("hFrac"));
    hFrac->Divide(hDplus);
    TGraphAsymmErrors* gFrac = static_cast<TGraphAsymmErrors *>(gDplusFromDstar->Clone("gFrac"));
    for(int iPt = 1; iPt < hFrac->GetNbinsX()+1; iPt++)
    {
        gFrac->SetPoint(iPt, hFrac->GetBinCenter(iPt), hFrac->GetBinContent(iPt));
        double relSysUncDplusLow = gDplus->GetErrorYlow(iPt) / hDplus->GetBinContent(iPt);
        double relSysUncDplusHigh = gDplus->GetErrorYhigh(iPt) / hDplus->GetBinContent(iPt);
        double relSysUncDplusFromDstarLow = gDplusFromDstar->GetErrorYlow(iPt) / hDplusFromDstar->GetBinContent(iPt);
        double relSysUncDplusFromDstarHigh = gDplusFromDstar->GetErrorYhigh(iPt) / hDplusFromDstar->GetBinContent(iPt);
        double sysUncLow = TMath::Sqrt(relSysUncDplusHigh*relSysUncDplusHigh + relSysUncDplusFromDstarLow*relSysUncDplusFromDstarLow) * hFrac->GetBinContent(iPt);
        double sysUncHigh = TMath::Sqrt(relSysUncDplusLow*relSysUncDplusLow + relSysUncDplusFromDstarHigh*relSysUncDplusFromDstarHigh) * hFrac->GetBinContent(iPt);
        gFrac->SetPointError(iPt, hFrac->GetBinWidth(iPt)/2, hFrac->GetBinWidth(iPt)/2, sysUncLow, sysUncHigh);
    }

    TLegend* leg = new TLegend(0.4, 0.7, 0.75, 0.9);
    leg->SetTextSize(0.045);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);    
    leg->AddEntry(hDstar, "D*^{+}", "p");    
    leg->AddEntry(fPowLawDstar, "power law fit to D*^{+}", "l");    
    leg->AddEntry(hDplus, "D^{+}", "p");    
    leg->AddEntry(hDplusFromDstar, "D^{+} #leftarrow D*^{+}", "p");    

    TCanvas *cPtCorr = new TCanvas("cPtCorr", "", 500, 500);
    cPtCorr->SetLogz();
    hPtCorr->Draw("colz");

    TCanvas *cResult = new TCanvas("cResult", "", 1000, 500);
    cResult->Divide(2, 1);
    cResult->cd(1)->DrawFrame(0., hDstar->GetMinimum() * 0.05, 36., hDstar->GetMaximum() * 5, ";#it{p}_{T} (GeV/#it{c}); d#sigma/(d#it{p}_{T}d#it{y}) (pb GeV^{-1} #it{c})");
    cResult->cd(1)->SetLogy();
    gDstar->Draw("2");
    hDstar->DrawCopy("same");
    fPowLawDstar->DrawCopy("same");
    gDplus->Draw("2");
    hDplus->DrawCopy("same");
    gDplusFromDstar->Draw("2");
    hDplusFromDstar->DrawCopy("same");
    leg->Draw();
    cResult->cd(2)->DrawFrame(0., 0., 36., 1., ";#it{p}_{T} (GeV/#it{c}); D^{+} #leftarrow D*^{+} / D^{+}");
    gFrac->Draw("2");
    hFrac->DrawCopy("same");
    cResult->Modified();
    cResult->Update();

    outFile.cd();
    nDplusFromDstar->Write();
    cPtCorr->Write();
    hPtCorr->Write();
    cResult->Write();
    gDstar->Write();
    hDstar->Write();
    fPowLawDstar->Write();
    gDplus->Write();
    hDplus->Write();
    gDplusFromDstar->Write();
    hDplusFromDstar->Write();
    gFrac->Write();
    hFrac->Write();
    outFile.Close();
}

//_______________________________________________________________________________________________
double PowLaw(double *pt, double *pars)
{
    return pars[0] * pt[0] / TMath::Power((1 + TMath::Power(pt[0] / pars[1], pars[3])), pars[2]);
}
