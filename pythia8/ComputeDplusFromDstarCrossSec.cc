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
void ComputeDplusFromDstarCrossSec(int nGen = 1000000, int decayer = kPythia8);
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
    inFileDplus->Close();
    TH1F* hDplusFromDstar = static_cast<TH1F *>(hDplus->Clone("hDplusFromDstar"));
    hDplusFromDstar->Reset();
    hDplusFromDstar->SetLineColor(kRed+1);
    hDplusFromDstar->SetMarkerColor(kRed+1);

    // parametrise cross section with power-law function
    TF1 *fPowLawDstar = new TF1("fPowLawDstar", PowLaw, 0, 50, 4);
    fPowLawDstar->SetParameters(hDstar->Integral(), 2.6, 2.8, 2.);
    fPowLawDstar->SetLineColor(kAzure + 4);
    hDstar->Fit("fPowLawDstar", "LIME0");

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
    TLorentzVector vec = TLorentzVector();

    TH2F* hPtCorr = new TH2F("hPtCorr", ";#it{p}_{T} (GeV/#it{c}) (D*^{+});#it{p}_{T} (GeV/#it{c}) (D^{+})", 500, 0., 50., 500, 0., 50.);

    // simulate D* -> D+ pi0 decay
    // TODO: take into account all uncertainties
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

    hDplusFromDstar->Scale(fPowLawDstar->Integral(0, 100) / hDplusFromDstar->Integral() * BRDstarToDplus, "width");
    TH1F* hFrac = static_cast<TH1F *>(hDplusFromDstar->Clone("hFrac"));
    hFrac->Divide(hDplus);

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
    hDstar->Draw("same");
    fPowLawDstar->Draw("same");
    hDplus->Draw("same");
    hDplusFromDstar->Draw("same");
    leg->Draw();
    cResult->cd(2)->DrawFrame(0., 0., 36., 1., ";#it{p}_{T} (GeV/#it{c}); D^{+} #leftarrow D*^{+} / D^{+}");
    hFrac->Draw("same");

    TString decName = "";
    if(decayer == kPythia6)
        decName = "Pythia6";
    else if(decayer == kPythia8)
        decName = "Pythia8";
    else if(decayer == kEvtGen)
        decName = "EvtGen";

    TFile outFile(Form("DplusFromDstar_pp5TeV_decayer_%s.root", decName.Data()), "recreate");
    cPtCorr->Write();
    hPtCorr->Write();
    cResult->Write();
    hDstar->Write();
    fPowLawDstar->Write();
    hDplus->Write();
    hDplusFromDstar->Write();
    hFrac->Write();
    outFile.Close();
}

//_______________________________________________________________________________________________
double PowLaw(double *pt, double *pars)
{
    return pars[0] * pt[0] / TMath::Power((1 + TMath::Power(pt[0] / pars[1], pars[3])), pars[2]);
}
