#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <math.h>

// ============================================================
// draw_LambdaNuclearId_relamom_sum_CFsub_method5_compare1s2s.C
//
// Method5 (sliding-window purity) の CF を
//   2sigma cut (auau3p85_all1.root, Black) vs
//   1sigma cut (auau3p85_all1_Nuclear1sigma.root, Red)
// で重ね描き比較するマクロ。
// 出力PDF: pdf/kstar_sum_CFsub_method5_compare1s2s_<species>.pdf
// ============================================================

void getPurity5(TH1D *hY, double minM, double maxM, double &P, double &err_P,
                double &N_sig, double &err_N_sig, double &N_tot,
                double &err_N_tot) {
  if (!hY || hY->GetEntries() < 10) {
    P = 0;
    err_P = 0;
    N_sig = 0;
    err_N_sig = 0;
    N_tot = 0;
    err_N_tot = 0;
    return;
  }
  TF1 *fTotal = new TF1("fTotal5", "gaus(0)+pol2(3)", 1.08, 1.15);
  fTotal->SetParameters(hY->GetMaximum(), 1.11596, 0.003, 0, 0, 0);
  fTotal->SetParLimits(1, 1.11, 1.12);
  fTotal->SetParLimits(2, 0.001, 0.008);
  TFitResultPtr r = hY->Fit(fTotal, "S Q N", "", 1.08, 1.15);

  int binMin = hY->FindBin(minM);
  int binMax = hY->FindBin(maxM);
  double exactMinM = hY->GetBinLowEdge(binMin);
  double exactMaxM = hY->GetBinLowEdge(binMax) + hY->GetBinWidth(binMax);
  N_tot = hY->IntegralAndError(binMin, binMax, err_N_tot);

  if (r->IsValid() && r->Status() == 0) {
    TF1 *fBg = new TF1("fBg5", "pol2", 1.08, 1.15);
    fBg->SetParameters(fTotal->GetParameter(3), fTotal->GetParameter(4),
                       fTotal->GetParameter(5));
    double covBg[9];
    TMatrixDSym cov = r->GetCovarianceMatrix();
    covBg[0] = cov(3, 3);
    covBg[1] = cov(3, 4);
    covBg[2] = cov(3, 5);
    covBg[3] = cov(4, 3);
    covBg[4] = cov(4, 4);
    covBg[5] = cov(4, 5);
    covBg[6] = cov(5, 3);
    covBg[7] = cov(5, 4);
    covBg[8] = cov(5, 5);

    double binW = hY->GetBinWidth(binMin);
    double N_bg = fBg->Integral(exactMinM, exactMaxM) / binW;
    double err_N_bg =
        fBg->IntegralError(exactMinM, exactMaxM, fBg->GetParameters(), covBg) /
        binW;

    N_sig = N_tot - N_bg;
    err_N_sig = sqrt(err_N_tot * err_N_tot + err_N_bg * err_N_bg);
    if (N_tot > 0) {
      P = N_sig / N_tot;
      err_P = (N_bg > 0 && N_tot > 0)
                  ? (N_bg / N_tot) *
                        sqrt(((err_N_bg / N_bg) * (err_N_bg / N_bg)) +
                             ((err_N_tot / N_tot) * (err_N_tot / N_tot)))
                  : 0;
    } else {
      P = 0;
      err_P = 0;
    }
    delete fBg;
  } else {
    P = 0;
    err_P = 0;
    N_sig = 0;
    err_N_sig = 0;
  }
  delete fTotal;
}

void getCentInfo(int groupSet, int group, TString &centLabel, int &jStart,
                 int &jEnd) {
  if (groupSet == 0) {
    centLabel = (group == 0) ? "20-80%" : "0-20%";
    jStart = (group == 0) ? 0 : 6;
    jEnd = (group == 0) ? 5 : 8;
  } else {
    if (group == 0) {
      centLabel = "0-10%";
      jStart = 7;
      jEnd = 8;
    } else if (group == 1) {
      centLabel = "10-20%";
      jStart = 6;
      jEnd = 6;
    } else {
      centLabel = "20-60%";
      jStart = 2;
      jEnd = 5;
    }
  }
}

TH1D *buildCF_method5(TFile *fIn, const TString &sp, int nRebin, int groupSet,
                      int group, int jStart, int jEnd, double minM, double maxM,
                      const TString &sfx) {
  TDirectory *dirNuc = (TDirectory *)fIn->Get("true");
  TDirectory *dirMix = (TDirectory *)fIn->Get("mix");

  TH2D *h2True_sum = 0;
  TH2D *h2Mix_sum = 0;

  for (int j = jStart; j <= jEnd; ++j) {
    TString hNameTrue = Form("hKstarMass_%s_CentBin%d", sp.Data(), j);
    TString hNameMix = Form("hKstarMass_Mixed_%s_CentBin%d", sp.Data(), j);
    TH2D *h2T =
        dirNuc ? (TH2D *)dirNuc->Get(hNameTrue) : (TH2D *)fIn->Get(hNameTrue);
    TH2D *h2M =
        dirMix ? (TH2D *)dirMix->Get(hNameMix) : (TH2D *)fIn->Get(hNameMix);

    if (h2T) {
      if (!h2True_sum)
        h2True_sum = (TH2D *)h2T->Clone(
            Form("h2T_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
      else
        h2True_sum->Add(h2T);
    }
    if (h2M) {
      if (!h2Mix_sum)
        h2Mix_sum = (TH2D *)h2M->Clone(
            Form("h2M_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
      else
        h2Mix_sum->Add(h2M);
    }
  }

  if (!h2True_sum || !h2Mix_sum) {
    if (h2True_sum)
      delete h2True_sum;
    if (h2Mix_sum)
      delete h2Mix_sum;
    return 0;
  }

  h2True_sum->RebinX(nRebin);
  h2Mix_sum->RebinX(nRebin);
  int nBinsX = h2True_sum->GetNbinsX();
  int W = 2; // Sliding window half-width

  TH1D *hSigTrue = new TH1D(
      Form("hST_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()), ";;N",
      nBinsX, h2True_sum->GetXaxis()->GetXmin(),
      h2True_sum->GetXaxis()->GetXmax());
  TH1D *hSigMix = new TH1D(
      Form("hSM_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()), ";;N",
      nBinsX, h2True_sum->GetXaxis()->GetXmin(),
      h2True_sum->GetXaxis()->GetXmax());
  TH1D *hCF = new TH1D(
      Form("hCF5_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()),
      Form("%s CF (method5);k* (GeV/c);CF", sp.Data()), nBinsX,
      h2True_sum->GetXaxis()->GetXmin(), h2True_sum->GetXaxis()->GetXmax());

  for (int bx = 1; bx <= nBinsX; ++bx) {
    int bx_s = TMath::Max(1, bx - W);
    int bx_e = TMath::Min(nBinsX, bx + W);

    TH1D *hYTw = h2True_sum->ProjectionY(
        Form("hYTw_%s_%d%s", sp.Data(), bx, sfx.Data()), bx_s, bx_e);
    TH1D *hYMw = h2Mix_sum->ProjectionY(
        Form("hYMw_%s_%d%s", sp.Data(), bx, sfx.Data()), bx_s, bx_e);
    double Pt, ePt, Nsigt_w, eNsigt_w, Ntott_w, eNtott_w;
    double Pm, ePm, Nsigm_w, eNsigm_w, Ntotm_w, eNtotm_w;
    getPurity5(hYTw, minM, maxM, Pt, ePt, Nsigt_w, eNsigt_w, Ntott_w, eNtott_w);
    getPurity5(hYMw, minM, maxM, Pm, ePm, Nsigm_w, eNsigm_w, Ntotm_w, eNtotm_w);
    if (Pt < 0 || Pt > 1) {
      Pt = 0;
      ePt = 0;
    }
    if (Pm < 0 || Pm > 1) {
      Pm = 0;
      ePm = 0;
    }
    delete hYTw;
    delete hYMw;

    TH1D *hYT = h2True_sum->ProjectionY(
        Form("hYT_%s_%d%s", sp.Data(), bx, sfx.Data()), bx, bx);
    TH1D *hYM = h2Mix_sum->ProjectionY(
        Form("hYM_%s_%d%s", sp.Data(), bx, sfx.Data()), bx, bx);
    int bMin = hYT->FindBin(minM), bMax = hYT->FindBin(maxM);
    double eNt, Nt = hYT->IntegralAndError(bMin, bMax, eNt);
    double eNm, Nm = hYM->IntegralAndError(bMin, bMax, eNm);
    double N_sig_t = Nt * Pt;
    double err_N_sig_t =
        sqrt(((eNt * Pt) * (eNt * Pt)) + ((Nt * ePt) * (Nt * ePt)));
    double N_sig_m = Nm * Pm;
    double err_N_sig_m =
        sqrt(((eNm * Pm) * (eNm * Pm)) + ((Nm * ePm) * (Nm * ePm)));
    hSigTrue->SetBinContent(bx, N_sig_t);
    hSigTrue->SetBinError(bx, err_N_sig_t);
    hSigMix->SetBinContent(bx, N_sig_m);
    hSigMix->SetBinError(bx, err_N_sig_m);
    delete hYT;
    delete hYM;
  }

  int bin1 = hSigTrue->FindBin(0.40001);
  int bin2 = hSigTrue->FindBin(0.59999);
  double intT = hSigTrue->Integral(bin1, bin2);
  double intM = hSigMix->Integral(bin1, bin2);
  if (intM > 0)
    hSigMix->Scale(intT / intM);

  for (int bx = 1; bx <= nBinsX; ++bx) {
    double nt = hSigTrue->GetBinContent(bx), et = hSigTrue->GetBinError(bx);
    double nm = hSigMix->GetBinContent(bx), em = hSigMix->GetBinError(bx);
    if (nm > 0 && nt > 0) {
      double cf = nt / nm;
      hCF->SetBinContent(bx, cf);
      hCF->SetBinError(
          bx, cf * sqrt(((et / nt) * (et / nt)) + ((em / nm) * (em / nm))));
    }
  }
  hCF->GetXaxis()->SetRangeUser(0.0, 0.6);

  delete h2True_sum;
  delete h2Mix_sum;
  delete hSigTrue;
  delete hSigMix;
  return hCF;
}

void draw_LambdaNuclearId_relamom_sum_CFsub_method5_compare1s2s() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  TString inputFile2s =
      "/star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/"
      "auau3p85_anaLambdaNuclearId/auau3p85_all1.root";
  TString inputFile1s =
      "/star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/"
      "auau3p85_anaLambdaNuclearId/auau3p85_all1_Nuclear1sigma.root";

  TFile *fIn2s = TFile::Open(inputFile2s, "READ");
  if (!fIn2s || fIn2s->IsZombie()) {
    std::cerr << "Error: Cannot open " << inputFile2s << std::endl;
    return;
  }
  TFile *fIn1s = TFile::Open(inputFile1s, "READ");
  if (!fIn1s || fIn1s->IsZombie()) {
    std::cerr << "Error: Cannot open " << inputFile1s << std::endl;
    return;
  }

  TString outDir = "/star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/"
                   "auau3p85_anaLambdaNuclearId/pdf";
  gSystem->mkdir(outDir, kTRUE);

  double meanLambda = 1.11596;
  double sigmaLambda = 2.96e-03;
  double minM = meanLambda - 3.0 * sigmaLambda;
  double maxM = meanLambda + 3.0 * sigmaLambda;

  const char *species[] = {"d", "t", "3He", "4He"};
  int nSpecies = 4;

  for (int i = 0; i < nSpecies; ++i) {
    TString sp = species[i];
    int nRebin = 1;
    if (sp == "t" || sp == "3He")
      nRebin = 4;
    else if (sp == "4He")
      nRebin = 10;

    TString outFile = Form("%s/kstar_sum_CFsub_method5_compare1s2s_%s.pdf",
                           outDir.Data(), sp.Data());

    for (int groupSet = 0; groupSet < 2; ++groupSet) {
      int nGroups = (groupSet == 0) ? 2 : 3;

      TCanvas *cComp = new TCanvas(
          Form("cComp5_%s_GS%d", sp.Data(), groupSet),
          Form("CF_method5 Compare (2#sigma vs 1#sigma) Lambda-%s", sp.Data()),
          1200, 600);
      cComp->Divide(nGroups, 1);

      for (int group = 0; group < nGroups; ++group) {
        TString centLabel;
        int jStart, jEnd;
        getCentInfo(groupSet, group, centLabel, jStart, jEnd);

        TH1D *hCF_2s = buildCF_method5(fIn2s, sp, nRebin, groupSet, group,
                                       jStart, jEnd, minM, maxM, "_2s");
        TH1D *hCF_1s = buildCF_method5(fIn1s, sp, nRebin, groupSet, group,
                                       jStart, jEnd, minM, maxM, "_1s");

        cComp->cd(group + 1);
        double yMax_full = (sp == "d") ? 25.0 : 2.0;
        if (hCF_2s) {
          hCF_2s->SetTitle(
              Form("CF_m5 Compare %s (%s)", sp.Data(), centLabel.Data()));
          hCF_2s->SetLineColor(kBlack);
          hCF_2s->SetMarkerColor(kBlack);
          hCF_2s->SetMarkerStyle(20);
          hCF_2s->SetMarkerSize(0.8);
          hCF_2s->GetYaxis()->SetRangeUser(0.0, yMax_full);
          hCF_2s->Draw("E1");
        }
        if (hCF_1s) {
          hCF_1s->SetLineColor(kRed);
          hCF_1s->SetMarkerColor(kRed);
          hCF_1s->SetMarkerStyle(24);
          hCF_1s->SetMarkerSize(0.8);
          hCF_1s->GetYaxis()->SetRangeUser(0.0, yMax_full);
          if (!hCF_2s) {
            hCF_1s->SetTitle(
                Form("CF_m5 Compare %s (%s)", sp.Data(), centLabel.Data()));
            hCF_1s->Draw("E1");
          } else {
            hCF_1s->Draw("E1 same");
          }
        }

        TLine *l1 = new TLine(0.0, 1.0, 0.6, 1.0);
        l1->SetLineStyle(2);
        l1->Draw("SAME");

        if (group == 0) {
          TLegend *leg = new TLegend(0.40, 0.72, 0.88, 0.88);
          if (hCF_2s)
            leg->AddEntry(hCF_2s, "p vs dE/dx cut: #pm2#sigma", "pe");
          if (hCF_1s)
            leg->AddEntry(hCF_1s, "p vs dE/dx cut: #pm1#sigma", "pe");
          leg->Draw();
        }
      } // group

      if (groupSet == 0) {
        cComp->Print(outFile + "(");
      } else {
        cComp->Print(outFile + ")");
      }

      delete cComp;
    } // groupSet

    std::cout << "Saved: " << outFile << std::endl;
  } // species

  fIn2s->Close();
  fIn1s->Close();
}
