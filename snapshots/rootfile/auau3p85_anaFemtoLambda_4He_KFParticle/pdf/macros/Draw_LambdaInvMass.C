// Generated adaptation; legacy original preserved.
#include "KFSummarySupport.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

void Draw_LambdaInvMass() {
 const double massFitSettings[9] = {1.0800000000000001,1.1499999999999999,1.1100000000000001,1.1200000000000001,0.001,0.0080000000000000002,0.0030000000000000001,10,2};
 KFReport::Init("/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/logs/Draw_LambdaInvMass_fit_audit.tsv",1.1159600000000001,0.00296,massFitSettings);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  TString inputFile = "/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root";

  TFile *fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
    return;
  }

  TString outDir = "/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf";

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

  double mean = 1.1159600000000001;
  double sigma = 0.00296;
  double nsigma = 3;

  double massMin = mean - nsigma * sigma;
  double massMax = mean + nsigma * sigma;

  double massSB1Min = mean - 4*nsigma * sigma;
  double massSB1Max = massMin;

  double massSB2Min = massMax;
  double massSB2Max = mean + 4*nsigma * sigma;

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

  TString outFile = Form("%s/Lambda_InvMass_legacy.pdf", outDir.Data());
  c1->Print(outFile);
  std::cout << "Saved " << outFile << std::endl;
  fIn->Close();
}