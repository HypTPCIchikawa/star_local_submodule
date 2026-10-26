// Generated adaptation; legacy original preserved.
#include "KFSummarySupport.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

void getPurity(TH1D* hY, double minM, double maxM, double &P, double &err_P,
  double &N_sig, double &err_N_sig, double &N_tot, double &err_N_tot) {
  KFReport::Purity(hY,minM,maxM,P,err_P,N_sig,err_N_sig,N_tot,err_N_tot);
}

void draw_LambdaNuclearId_relamom_sum_CFsub_method5() {
 const double massFitSettings[9] = {1.0800000000000001,1.1499999999999999,1.1100000000000001,1.1200000000000001,0.001,0.0080000000000000002,0.0030000000000000001,10,2};
 KFReport::Init("/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/logs/draw_LambdaNuclearId_relamom_sum_CFsub_method5_fit_audit.tsv",1.1159600000000001,0.00296,massFitSettings);
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

  const char *species[] = {"4He"};
  int nSpecies = 1;

  double meanLambda = 1.1159600000000001;
  double sigmaLambda = 0.00296;
  double minM = meanLambda - 3*sigmaLambda;
  double maxM = meanLambda + 3*sigmaLambda;

  TDirectory *dirNuc = (TDirectory *)fIn->Get("true");
  TDirectory *dirMix = (TDirectory *)fIn->Get("mix");

  for (int i = 0; i < nSpecies; ++i) {
    TString sp = species[i];

    int nRebin = 1;
    if (sp == "t" || sp == "3He")
      nRebin = 4;
    else if (sp == "4He")
      nRebin = 4;

    TString outFile =
        Form("%s/kstar_sum_CFsub_method5_%s.pdf", outDir.Data(), sp.Data());

    TCanvas *c_start = new TCanvas("c_start");
    c_start->Print(outFile + "[");

    for (int groupSet = 0; groupSet < 2; ++groupSet) {
      int nGroups = (groupSet == 0) ? 2 : 3;

      TCanvas *cCF = new TCanvas(Form("cCF_%s_groupSet%d", sp.Data(), groupSet),
                                 Form("CF(k*) - %s", sp.Data()), 1200, 600);
      cCF->Divide(nGroups, 1);

      TCanvas *cPurity =
          new TCanvas(Form("cPurity_%s_groupSet%d", sp.Data(), groupSet),
                      Form("Purity(k*) - %s", sp.Data()), 1200, 600);
      cPurity->Divide(nGroups, 1);

      TCanvas *cNSig =
          new TCanvas(Form("cNSig_%s_groupSet%d", sp.Data(), groupSet),
                      Form("NSig(k*) - %s", sp.Data()), 1200, 600);
      cNSig->Divide(nGroups, 1);

      for (int group = 0; group < nGroups; ++group) {
        TString centLabel;
        int jStart, jEnd;
        if (groupSet == 0) {
if (group == 0) { centLabel = "20-80%"; jStart = 0; jEnd = 5; }
else if (group == 1) { centLabel = "0-20%"; jStart = 6; jEnd = 8; }
} else {
if (group == 0) { centLabel = "0-10%"; jStart = 7; jEnd = 8; }
else if (group == 1) { centLabel = "10-20%"; jStart = 6; jEnd = 6; }
else if (group == 2) { centLabel = "20-60%"; jStart = 2; jEnd = 5; }
}

        TH2 *h2True_sum = 0;
        TH2 *h2Mix_sum = 0;

        for (int j = jStart; j <= jEnd; ++j) {
          TString hNameTrue = Form("hKstarMass_%s_CentBin%d", sp.Data(), j);
          TString hNameMix =
              Form("hKstarMass_Mixed_%s_CentBin%d", sp.Data(), j);

          TH2 *h2True = 0;
          TH2 *h2Mix = 0;

          if (dirNuc)
            h2True = (TH2 *)dirNuc->Get(hNameTrue);
          if (dirMix)
            h2Mix = (TH2 *)dirMix->Get(hNameMix);

          if (!h2True)
            h2True = (TH2 *)fIn->Get(hNameTrue);
          if (!h2Mix)
            h2Mix = (TH2 *)fIn->Get(hNameMix);

          if (h2True) {
            if (!h2True_sum) {
              h2True_sum = (TH2 *)h2True->Clone(
                  Form("h2True_sum_%s_GS%d_G%d", sp.Data(), groupSet, group));
            } else {
              h2True_sum->Add(h2True);
            }
          }
          if (h2Mix) {
            if (!h2Mix_sum) {
              h2Mix_sum = (TH2 *)h2Mix->Clone(
                  Form("h2Mix_sum_%s_GS%d_G%d", sp.Data(), groupSet, group));
            } else {
              h2Mix_sum->Add(h2Mix);
            }
          }
        }

        if (!h2True_sum || !h2Mix_sum)
          continue;

        h2True_sum->RebinX(nRebin);
        h2Mix_sum->RebinX(nRebin);

        int nBinsX = h2True_sum->GetNbinsX();

        TH1D *hCF = new TH1D(
            Form("hCF_%s_GS%d_G%d", sp.Data(), groupSet, group),
            Form("%s %s CF;k* (GeV/c);CF", sp.Data(), centLabel.Data()), nBinsX,
            h2True_sum->GetXaxis()->GetXmin(),
            h2True_sum->GetXaxis()->GetXmax());
        TH1D *hPTrue = new TH1D(
            Form("hPTrue_%s_GS%d_G%d", sp.Data(), groupSet, group),
            Form("%s %s Purity;k* (GeV/c);Purity", sp.Data(), centLabel.Data()),
            nBinsX, h2True_sum->GetXaxis()->GetXmin(),
            h2True_sum->GetXaxis()->GetXmax());
        TH1D *hPMix = new TH1D(
            Form("hPMix_%s_GS%d_G%d", sp.Data(), groupSet, group),
            Form("%s %s Purity;k* (GeV/c);Purity", sp.Data(), centLabel.Data()),
            nBinsX, h2True_sum->GetXaxis()->GetXmin(),
            h2True_sum->GetXaxis()->GetXmax());

        TH1D *hSigTrue =
            new TH1D(Form("hSigTrue_%s_GS%d_G%d", sp.Data(), groupSet, group),
                     Form("%s %s N_{sig};k* (GeV/c);Counts", sp.Data(),
                          centLabel.Data()),
                     nBinsX, h2True_sum->GetXaxis()->GetXmin(),
                     h2True_sum->GetXaxis()->GetXmax());
        TH1D *hSigMix =
            new TH1D(Form("hSigMix_%s_GS%d_G%d", sp.Data(), groupSet, group),
                     Form("%s %s N_{sig};k* (GeV/c);Counts", sp.Data(),
                          centLabel.Data()),
                     nBinsX, h2True_sum->GetXaxis()->GetXmin(),
                     h2True_sum->GetXaxis()->GetXmax());

        int W = 2; // Sliding window half-width (merges 2*W + 1 bins)

        for (int bx = 1; bx <= nBinsX; ++bx) {
          int bx_start = TMath::Max(1, bx - W);
          int bx_end = TMath::Min(nBinsX, bx + W);

          // Sliding window projection for Purity
          TH1D *hYTrue_window = h2True_sum->ProjectionY(
              Form("hYTrue_w_GS%d_G%d_B%d", groupSet, group, bx), bx_start,
              bx_end);
          TH1D *hYMix_window = h2Mix_sum->ProjectionY(
              Form("hYMix_w_GS%d_G%d_B%d", groupSet, group, bx), bx_start,
              bx_end);

          double Pt, ePt, Nsigt_w, eNsigt_w, Ntott_w, eNtott_w;
          double Pm, ePm, Nsigm_w, eNsigm_w, Ntotm_w, eNtotm_w;

          getPurity(hYTrue_window, minM, maxM, Pt, ePt, Nsigt_w, eNsigt_w,
                    Ntott_w, eNtott_w);
          getPurity(hYMix_window, minM, maxM, Pm, ePm, Nsigm_w, eNsigm_w,
                    Ntotm_w, eNtotm_w);

          // Fallback to 0 if fit failed horribly
          if (Pt < 0 || Pt > 1) {
            Pt = 0;
            ePt = 0;
          }
          if (Pm < 0 || Pm > 1) {
            Pm = 0;
            ePm = 0;
          }

          hPTrue->SetBinContent(bx, Pt);
          hPTrue->SetBinError(bx, ePt);
          hPMix->SetBinContent(bx, Pm);
          hPMix->SetBinError(bx, ePm);

          delete hYTrue_window;
          delete hYMix_window;

          // Single fine bin projection for N_tot
          TH1D *hYTrue = h2True_sum->ProjectionY(
              Form("hYTrue_GS%d_G%d_B%d", groupSet, group, bx), bx, bx);
          TH1D *hYMix = h2Mix_sum->ProjectionY(
              Form("hYMix_GS%d_G%d_B%d", groupSet, group, bx), bx, bx);

          int binMin = hYTrue->FindBin(minM);
          int binMax = hYTrue->FindBin(maxM);

          double err_Ntott;
          double Ntott = hYTrue->IntegralAndError(binMin, binMax, err_Ntott);

          double err_Ntotm;
          double Ntotm = hYMix->IntegralAndError(binMin, binMax, err_Ntotm);

          double N_sig_t = Ntott * Pt;
          double err_N_sig_t = TMath::Sqrt(TMath::Power(err_Ntott * Pt, 2) +
                                           TMath::Power(Ntott * ePt, 2));

          double N_sig_m = Ntotm * Pm;
          double err_N_sig_m = TMath::Sqrt(TMath::Power(err_Ntotm * Pm, 2) +
                                           TMath::Power(Ntotm * ePm, 2));

          hSigTrue->SetBinContent(bx, N_sig_t);
          hSigTrue->SetBinError(bx, err_N_sig_t);
          hSigMix->SetBinContent(bx, N_sig_m);
          hSigMix->SetBinError(bx, err_N_sig_m);

          delete hYTrue;
          delete hYMix;
        }

        int bin1 = hSigTrue->FindBin(0.40001000000000003);
        int bin2 = hSigTrue->FindBin(0.59999000000000002);
        double intTrue = hSigTrue->Integral(bin1, bin2);
        double intMix = hSigMix->Integral(bin1, bin2);

        if (intMix > 0) {
          hSigMix->Scale(intTrue / intMix);
        }

        for (int bx = 1; bx <= nBinsX; ++bx) {
          double nt = hSigTrue->GetBinContent(bx);
          double et = hSigTrue->GetBinError(bx);
          double nm = hSigMix->GetBinContent(bx);
          double em = hSigMix->GetBinError(bx);

          if (nt > 0 && nm > 0 && TMath::Finite(nt) && TMath::Finite(nm) && TMath::Finite(et) && TMath::Finite(em)) {
            double cf = nt / nm;
            double err_cf = TMath::Sqrt(TMath::Power(et/nm,2)+TMath::Power(nt*em/(nm*nm),2));
            hCF->SetBinContent(bx, cf);
            hCF->SetBinError(bx, err_cf);
          }
        }

        cNSig->cd(group + 1);
        hSigTrue->SetTitle(
            Form("%s (%s) N_{sig}", sp.Data(), centLabel.Data()));
        hSigTrue->SetMarkerStyle(20);
        hSigTrue->SetMarkerColor(kRed);
        hSigTrue->SetLineColor(kRed);
        hSigTrue->GetXaxis()->SetRangeUser(0,0.40000000000000002);
        double maxH = hSigTrue->GetMaximum();
        double maxMh = hSigMix->GetMaximum();
        hSigTrue->SetMaximum(TMath::Max(maxH, maxMh) * 1.2);
        hSigTrue->SetMinimum(0);
        hSigTrue->Draw("E1");

        hSigMix->SetMarkerStyle(24);
        hSigMix->SetMarkerColor(kBlue);
        hSigMix->SetLineColor(kBlue);
        hSigMix->Draw("E1 SAME");

        if (group == 0) {
          TLegend *legSig = new TLegend(0.4, 0.75, 0.88, 0.88);
          legSig->AddEntry(hSigTrue, "True N_{sig}", "pl");
          legSig->AddEntry(hSigMix, "Mixed N_{sig} (Scaled)", "pl");
          legSig->Draw();
        }

        cPurity->cd(group + 1);
        hPTrue->SetTitle(Form("%s (%s) Purity", sp.Data(), centLabel.Data()));
        hPTrue->SetMarkerStyle(20);
        hPTrue->SetMarkerColor(kRed);
        hPTrue->SetLineColor(kRed);
        hPTrue->SetMinimum(0.0);
        hPTrue->SetMaximum(1.2);
        hPTrue->GetXaxis()->SetRangeUser(0,0.40000000000000002);
        KFReport::DrawMeasured(hPTrue,"E1");

        hPMix->SetMarkerStyle(24);
        hPMix->SetMarkerColor(kBlue);
        hPMix->SetLineColor(kBlue);
        KFReport::DrawMeasured(hPMix,"E1 SAME");

        if (group == 0) {
          TLegend *legP = new TLegend(0.4, 0.75, 0.88, 0.88);
          legP->AddEntry(hPTrue, "True Purity", "pl");
          legP->AddEntry(hPMix, "Mixed Purity", "pl");
          legP->Draw();
        }

        cCF->cd(group + 1);
        hCF->SetTitle(Form("%s (%s) CF", sp.Data(), centLabel.Data()));
        hCF->SetMarkerStyle(20);
        hCF->SetMarkerColor(kBlack);
        hCF->SetLineColor(kBlack);
        hCF->SetMinimum(0.0);
        hCF->SetMaximum(2.0);
        hCF->GetXaxis()->SetRangeUser(0,0.40000000000000002);
        KFReport::DrawMeasured(hCF,"E1");

        TLine *l1 = new TLine(0,1.0,0.40000000000000002,1.0);
        l1->SetLineStyle(2);
        l1->Draw("SAME");
      } // end group loop

      cNSig->Print(outFile);
      cPurity->Print(outFile);
      cCF->Print(outFile);
    } // end groupSet loop

    TCanvas *c_end = new TCanvas("c_end");
    c_end->Print(outFile + "]");

    std::cout << "Successfully created " << outFile << std::endl;
  } // end species loop
}
