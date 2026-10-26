#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

// ================================================================
// Fit_LambdaInvMass_SN_3p85.C
//
// hLambda_InvMass を pol2+gaus でフィットし、
// S/N および Purity を計算して PDF に出力する。
// 出力PDF: pdf/Lambda_InvMass_SN_fit_3p85.pdf
// ================================================================

void Fit_LambdaInvMass_SN_3p85() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(1);

  TString inputFile =
    "/star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/"
    "auau3p85_anaLambdaNuclearId/auau3p85_all1.root";
  TFile *fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Error: Cannot open " << inputFile << std::endl; return;
  }

  // ---------- histogram search ----------
  TH1 *h = (TH1*)fIn->Get("hLambda_InvMass");
  if (!h) {
    TDirectory *d = (TDirectory*)fIn->Get("lambda");
    if (d) h = (TH1*)d->Get("hLambda_InvMass");
  }
  if (!h) {
    std::cerr << "Error: hLambda_InvMass not found" << std::endl;
    fIn->Close(); return;
  }

  // ---------- fit ----------
  double fitMin = 1.08, fitMax = 1.15;
  TF1 *fTotal = new TF1("fTotal_SN","gaus(0)+pol2(3)", fitMin, fitMax);
  fTotal->SetParameters(h->GetMaximum(), 1.11596, 0.003, 0, 0, 0);
  fTotal->SetParLimits(1, 1.110, 1.122);
  fTotal->SetParLimits(2, 0.001, 0.008);
  fTotal->SetParNames("Amp","Mean","Sigma","p0","p1","p2");
  fTotal->SetLineColor(kBlue);
  fTotal->SetLineWidth(2);

  h->Fit(fTotal, "S R");

  double mu    = fTotal->GetParameter(1);
  double sigma = fTotal->GetParameter(2);
  double amu   = fTotal->GetParError(1);
  double asig  = fTotal->GetParError(2);

  double nsig_range = 3.0;
  double sigMin = mu - nsig_range * sigma;
  double sigMax = mu + nsig_range * sigma;
  double binW   = h->GetBinWidth(1);

  // Gaussian component
  TF1 *fGaus = new TF1("fGaus_SN","gaus", fitMin, fitMax);
  fGaus->SetParameters(fTotal->GetParameter(0), mu, sigma);
  fGaus->SetLineColor(kRed);
  fGaus->SetLineStyle(2);
  fGaus->SetLineWidth(2);

  // Background component (pol2)
  TF1 *fBg = new TF1("fBg_SN","pol2", fitMin, fitMax);
  fBg->SetParameters(fTotal->GetParameter(3),
                     fTotal->GetParameter(4),
                     fTotal->GetParameter(5));
  fBg->SetLineColor(kGreen+2);
  fBg->SetLineStyle(2);
  fBg->SetLineWidth(2);

  double N_sig = fGaus->Integral(sigMin, sigMax) / binW;
  double N_bg  = fBg->Integral(sigMin, sigMax) / binW;
  double SN    = (N_bg > 0) ? N_sig / N_bg : -1;
  double purity= N_sig / (N_sig + N_bg);

  std::cout << std::endl;
  std::cout << "=== Lambda S/N for AuAu 3.85 GeV ===" << std::endl;
  std::cout << "  Fit range  : " << fitMin << " -- " << fitMax << " GeV/c^2" << std::endl;
  std::cout << "  Mean       : " << mu    << " +/- " << amu  << " GeV/c^2" << std::endl;
  std::cout << "  Sigma      : " << sigma << " +/- " << asig << " GeV/c^2" << std::endl;
  std::cout << "  Sig window : " << sigMin << " -- " << sigMax << " (+-3sigma)" << std::endl;
  std::cout << "  N_sig      : " << N_sig  << std::endl;
  std::cout << "  N_bg       : " << N_bg   << std::endl;
  std::cout << "  S/N        : " << SN     << std::endl;
  std::cout << "  Purity     : " << purity << std::endl;
  std::cout << "======================================" << std::endl;

  // ---------- draw ----------
  TCanvas *c1 = new TCanvas("c1_SN","Lambda Invariant Mass Fit",900,700);
  c1->SetLeftMargin(0.13);
  c1->SetBottomMargin(0.12);

  h->SetLineColor(kBlack);
  h->SetLineWidth(1);
  h->GetXaxis()->SetTitle("M_{p#pi^{-}} [GeV/c^{2}]");
  h->GetYaxis()->SetTitle("Counts");
  h->GetXaxis()->SetRangeUser(1.08, 1.15);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->SetTitle("Lambda Invariant Mass  (AuAu 3.85 GeV, all centralities)");
  h->Draw("HIST");

  fTotal->Draw("SAME");
  fGaus->Draw("SAME");
  fBg->Draw("SAME");

  // Signal window lines
  double ymax = h->GetMaximum();
  TLine *lMin = new TLine(sigMin, 0, sigMin, ymax * 0.8);
  lMin->SetLineColor(kRed); lMin->SetLineStyle(2); lMin->SetLineWidth(2); lMin->Draw("SAME");
  TLine *lMax = new TLine(sigMax, 0, sigMax, ymax * 0.8);
  lMax->SetLineColor(kRed); lMax->SetLineStyle(2); lMax->SetLineWidth(2); lMax->Draw("SAME");

  // Legend
  TLegend *leg = new TLegend(0.14, 0.62, 0.55, 0.88);
  leg->AddEntry(h,      "Data (all CentBin)", "l");
  leg->AddEntry(fTotal, "Fit: Gaus + Pol2",  "l");
  leg->AddEntry(fGaus,  "Signal (Gaussian)", "l");
  leg->AddEntry(fBg,    "Background (Pol2)", "l");
  leg->Draw();

  // Stat box (fit result + S/N)
  TPaveText *pt = new TPaveText(0.57, 0.55, 0.95, 0.88, "NDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.030);
  pt->AddText(Form("Mean  = %.5f GeV/c^{2}", mu));
  pt->AddText(Form("#sigma = %.4f GeV/c^{2}", sigma));
  pt->AddText("Signal window: #pm3#sigma");
  pt->AddText(Form("N_{sig} = %.3e", N_sig));
  pt->AddText(Form("N_{bg}  = %.3e", N_bg));
  pt->AddText(Form("S/N  = %.3f", SN));
  pt->AddText(Form("Purity = %.4f", purity));
  pt->Draw();

  // ---------- output ----------
  TString outDir =
    "/star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/"
    "auau3p85_anaLambdaNuclearId/pdf";
  gSystem->mkdir(outDir, kTRUE);
  TString outFile = Form("%s/Lambda_InvMass_SN_fit_3p85.pdf", outDir.Data());
  c1->Print(outFile);
  std::cout << "Saved: " << outFile << std::endl;

  fIn->Close();
}
