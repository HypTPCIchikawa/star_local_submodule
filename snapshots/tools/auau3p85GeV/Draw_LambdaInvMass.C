#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

void Draw_LambdaInvMass() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  TString inputFile = "rootfile/auau3p85_anaLambdaNuclearId/auau3p85_all1.root";
  TFile *fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
    return;
  }

  TString outDir = "rootfile/auau3p85_anaLambdaNuclearId/pdf";
  gSystem->mkdir(outDir, kTRUE);

  TString histName = "hLambda_InvMass";
  TH1 *h = (TH1 *)fIn->Get(histName);

  if (!h) {
    TDirectory *dir = (TDirectory *)fIn->Get("lambda");
    if (dir) h = (TH1 *)dir->Get(histName);
  }

  if (!h) {
    std::cerr << "Error: histogram not found" << std::endl;
    fIn->Close();
    return;
  }

  TCanvas *c1 = new TCanvas("c1", "Lambda Invariant Mass", 800, 600);
  c1->SetLeftMargin(0.15);

  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle("M_{p#pi^{-}} [GeV/c^{2}]");
  h->GetYaxis()->SetTitle("Counts");
  h->GetYaxis()->SetTitleOffset(1.5);
  h->Draw("HIST");

  double mean = 1.11596;
  double sigma = 2.96e-03;
  double nsigma = 3.0;

  double massMin = mean - nsigma * sigma;
  double massMax = mean + nsigma * sigma;

  double massSB1Min = mean - 4.0 * nsigma * sigma;
  double massSB1Max = massMin;

  double massSB2Min = massMax;
  double massSB2Max = mean + 4.0 * nsigma * sigma;

  double yMax = h->GetMaximum();

  TLine *lineMin = new TLine(massMin, 0, massMin, yMax);
  lineMin->SetLineColor(kRed);
  lineMin->SetLineStyle(2);
  lineMin->SetLineWidth(2);
  lineMin->Draw("SAME");

  TLine *lineMax = new TLine(massMax, 0, massMax, yMax);
  lineMax->SetLineColor(kRed);
  lineMax->SetLineStyle(2);
  lineMax->SetLineWidth(2);
  lineMax->Draw("SAME");

  TLine *lineSB1Min = new TLine(massSB1Min, 0, massSB1Min, yMax);
  lineSB1Min->SetLineColor(kBlue);
  lineSB1Min->SetLineStyle(2);
  lineSB1Min->SetLineWidth(2);
  lineSB1Min->Draw("SAME");

  TLine *lineSB2Max = new TLine(massSB2Max, 0, massSB2Max, yMax);
  lineSB2Max->SetLineColor(kBlue);
  lineSB2Max->SetLineStyle(2);
  lineSB2Max->SetLineWidth(2);
  lineSB2Max->Draw("SAME");

  TString outFile = Form("%s/Lambda_InvMass.pdf", outDir.Data());
  c1->Print(outFile);
  std::cout << "Saved " << outFile << std::endl;
  fIn->Close();
}