#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <math.h>

// ============================================================
// draw_LambdaNuclearId_relamom_sum_compare1s2s.C
//
// 2sigma cut (auau3p85_all1.root, Black) vs
// 1sigma cut (auau3p85_all1_Nuclear1sigma.root, Red)
// の Correlation Function を重ね描き比較するマクロ。
// 出力PDF: pdf/kstar_sum_compare1s2s_<species>.pdf
// ============================================================

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

// Helper: accumulate histograms and compute CF (SB-subtracted)
TH1 *buildCF_sb(TFile *fIn, const TString &sp, int nRebin, int groupSet,
                int group, int jStart, int jEnd, const TString &sfx) {
  TDirectory *dirNuc = (TDirectory *)fIn->Get("true");
  TDirectory *dirMix = (TDirectory *)fIn->Get("mix");

  TH1 *h = 0, *hmix = 0;
  TH1 *hSBNeg_orig = 0, *hSBPos_orig = 0, *hmixSBNeg_orig = 0,
      *hmixSBPos_orig = 0;

  for (int j = jStart; j <= jEnd; ++j) {
    // True kstar
    TString histName = Form("hKstar_%s_CentBin%d", sp.Data(), j);
    TH1 *h_tmp =
        dirNuc ? (TH1 *)dirNuc->Get(histName) : (TH1 *)fIn->Get(histName);
    if (h_tmp) {
      if (!h) {
        h = (TH1 *)h_tmp->Clone(
            Form("h_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!h->GetSumw2N())
          h->Sumw2();
      } else {
        TH1 *h_add = (TH1 *)h_tmp->Clone();
        if (!h_add->GetSumw2N())
          h_add->Sumw2();
        h->Add(h_add);
        delete h_add;
      }
    }

    // Mixed kstar
    TString mixName = Form("hKstar_Mixed_%s_CentBin%d", sp.Data(), j);
    TH1 *hmix_tmp =
        dirMix ? (TH1 *)dirMix->Get(mixName) : (TH1 *)fIn->Get(mixName);
    if (hmix_tmp) {
      if (!hmix) {
        hmix = (TH1 *)hmix_tmp->Clone(
            Form("hmix_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hmix->GetSumw2N())
          hmix->Sumw2();
      } else {
        TH1 *hmix_add = (TH1 *)hmix_tmp->Clone();
        if (!hmix_add->GetSumw2N())
          hmix_add->Sumw2();
        hmix->Add(hmix_add);
        delete hmix_add;
      }
    }

    // SBNeg
    TString nameSBNeg = Form("hKstar_%s_SBNeg_CentBin%d", sp.Data(), j);
    TH1 *hSBNeg_tmp =
        dirNuc ? (TH1 *)dirNuc->Get(nameSBNeg) : (TH1 *)fIn->Get(nameSBNeg);
    if (hSBNeg_tmp) {
      if (!hSBNeg_orig) {
        hSBNeg_orig = (TH1 *)hSBNeg_tmp->Clone(Form(
            "hSBNeg_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hSBNeg_orig->GetSumw2N())
          hSBNeg_orig->Sumw2();
      } else {
        TH1 *hSBNeg_add = (TH1 *)hSBNeg_tmp->Clone();
        if (!hSBNeg_add->GetSumw2N())
          hSBNeg_add->Sumw2();
        hSBNeg_orig->Add(hSBNeg_add);
        delete hSBNeg_add;
      }
    }

    // SBPos
    TString nameSBPos = Form("hKstar_%s_SBPos_CentBin%d", sp.Data(), j);
    TH1 *hSBPos_tmp =
        dirNuc ? (TH1 *)dirNuc->Get(nameSBPos) : (TH1 *)fIn->Get(nameSBPos);
    if (hSBPos_tmp) {
      if (!hSBPos_orig) {
        hSBPos_orig = (TH1 *)hSBPos_tmp->Clone(Form(
            "hSBPos_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hSBPos_orig->GetSumw2N())
          hSBPos_orig->Sumw2();
      } else {
        TH1 *hSBPos_add = (TH1 *)hSBPos_tmp->Clone();
        if (!hSBPos_add->GetSumw2N())
          hSBPos_add->Sumw2();
        hSBPos_orig->Add(hSBPos_add);
        delete hSBPos_add;
      }
    }

    // Mixed SBNeg
    TString mixNameSBNeg =
        Form("hKstar_Mixed_%s_SBNeg_CentBin%d", sp.Data(), j);
    TH1 *hmixSBNeg_tmp = dirMix ? (TH1 *)dirMix->Get(mixNameSBNeg)
                                : (TH1 *)fIn->Get(mixNameSBNeg);
    if (hmixSBNeg_tmp) {
      if (!hmixSBNeg_orig) {
        hmixSBNeg_orig = (TH1 *)hmixSBNeg_tmp->Clone(Form(
            "hmixSBNeg_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hmixSBNeg_orig->GetSumw2N())
          hmixSBNeg_orig->Sumw2();
      } else {
        TH1 *hmixSBNeg_add = (TH1 *)hmixSBNeg_tmp->Clone();
        if (!hmixSBNeg_add->GetSumw2N())
          hmixSBNeg_add->Sumw2();
        hmixSBNeg_orig->Add(hmixSBNeg_add);
        delete hmixSBNeg_add;
      }
    }

    // Mixed SBPos
    TString mixNameSBPos =
        Form("hKstar_Mixed_%s_SBPos_CentBin%d", sp.Data(), j);
    TH1 *hmixSBPos_tmp = dirMix ? (TH1 *)dirMix->Get(mixNameSBPos)
                                : (TH1 *)fIn->Get(mixNameSBPos);
    if (hmixSBPos_tmp) {
      if (!hmixSBPos_orig) {
        hmixSBPos_orig = (TH1 *)hmixSBPos_tmp->Clone(Form(
            "hmixSBPos_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hmixSBPos_orig->GetSumw2N())
          hmixSBPos_orig->Sumw2();
      } else {
        TH1 *hmixSBPos_add = (TH1 *)hmixSBPos_tmp->Clone();
        if (!hmixSBPos_add->GetSumw2N())
          hmixSBPos_add->Sumw2();
        hmixSBPos_orig->Add(hmixSBPos_add);
        delete hmixSBPos_add;
      }
    }
  }

  if (!h || !hmix) {
    if (h)
      delete h;
    if (hmix)
      delete hmix;
    if (hSBNeg_orig)
      delete hSBNeg_orig;
    if (hSBPos_orig)
      delete hSBPos_orig;
    if (hmixSBNeg_orig)
      delete hmixSBNeg_orig;
    if (hmixSBPos_orig)
      delete hmixSBPos_orig;
    return 0;
  }

  if (hSBNeg_orig)
    hSBNeg_orig->Scale(1.0 / 3.0);
  if (hSBPos_orig)
    hSBPos_orig->Scale(1.0 / 3.0);
  if (hmixSBNeg_orig)
    hmixSBNeg_orig->Scale(1.0 / 3.0);
  if (hmixSBPos_orig)
    hmixSBPos_orig->Scale(1.0 / 3.0);

  h->Rebin(nRebin);
  hmix->Rebin(nRebin);
  TH1 *hmix_unscaled = (TH1 *)hmix->Clone(Form(
      "hmix_unscaled_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));

  // 5ページ目相当: Sideband Subtraction
  TH1 *h5_true = (TH1 *)h->Clone(
      Form("h5_true_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
  TH1 *hSB_tot_true = 0;
  if (hSBNeg_orig && hSBPos_orig) {
    hSB_tot_true = (TH1 *)hSBNeg_orig->Clone(Form(
        "hSB_tot_true_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
    if (!hSB_tot_true->GetSumw2N())
      hSB_tot_true->Sumw2();
    hSB_tot_true->Add(hSBPos_orig);
    hSB_tot_true->Rebin(nRebin);
  } else if (hSBNeg_orig) {
    hSB_tot_true = (TH1 *)hSBNeg_orig->Clone(Form(
        "hSB_tot_true_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
    if (!hSB_tot_true->GetSumw2N())
      hSB_tot_true->Sumw2();
    hSB_tot_true->Rebin(nRebin);
  } else if (hSBPos_orig) {
    hSB_tot_true = (TH1 *)hSBPos_orig->Clone(Form(
        "hSB_tot_true_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
    if (!hSB_tot_true->GetSumw2N())
      hSB_tot_true->Sumw2();
    hSB_tot_true->Rebin(nRebin);
  }
  if (hSB_tot_true) {
    h5_true->Add(hSB_tot_true, -1.0);
  }

  TH1 *h5_mix = 0;
  if (hmix_unscaled) {
    h5_mix = (TH1 *)hmix_unscaled->Clone(
        Form("h5_mix_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
    TH1 *hmixSB_tot = 0;
    if (hmixSBNeg_orig && hmixSBPos_orig) {
      hmixSB_tot = (TH1 *)hmixSBNeg_orig->Clone(Form(
          "hmixSB_tot_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
      if (!hmixSB_tot->GetSumw2N())
        hmixSB_tot->Sumw2();
      hmixSB_tot->Add(hmixSBPos_orig);
      hmixSB_tot->Rebin(nRebin);
    } else if (hmixSBNeg_orig) {
      hmixSB_tot = (TH1 *)hmixSBNeg_orig->Clone(Form(
          "hmixSB_tot_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
      if (!hmixSB_tot->GetSumw2N())
        hmixSB_tot->Sumw2();
      hmixSB_tot->Rebin(nRebin);
    } else if (hmixSBPos_orig) {
      hmixSB_tot = (TH1 *)hmixSBPos_orig->Clone(Form(
          "hmixSB_tot_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
      if (!hmixSB_tot->GetSumw2N())
        hmixSB_tot->Sumw2();
      hmixSB_tot->Rebin(nRebin);
    }
    if (hmixSB_tot) {
      h5_mix->Add(hmixSB_tot, -1.0);
      delete hmixSB_tot;
    }

    // 0.4 - 0.6 GeV/c で積分してスケールを合わせる
    int bin1 = h5_true->FindBin(0.40001);
    int bin2 = h5_true->FindBin(0.59999);
    double int5_true = h5_true->Integral(bin1, bin2);
    double int5_mix = h5_mix->Integral(bin1, bin2);
    if (int5_mix > 0) {
      h5_mix->Scale(int5_true / int5_mix);
    }
  }

  // 6ページ目相当: CF (SB Subtracted)
  TH1 *hCF_sb = 0;
  if (h5_mix) {
    hCF_sb = (TH1 *)h5_true->Clone(
        Form("hCF_sb_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
    hCF_sb->Divide(h5_mix);
    hCF_sb->GetYaxis()->SetTitle("C(k*) = True_sub / Mixed_sub");
    hCF_sb->GetXaxis()->SetRangeUser(0.0, 0.6);
  }

  // Cleanup
  delete h;
  delete hmix;
  delete hmix_unscaled;
  if (hSBNeg_orig)
    delete hSBNeg_orig;
  if (hSBPos_orig)
    delete hSBPos_orig;
  if (hmixSBNeg_orig)
    delete hmixSBNeg_orig;
  if (hmixSBPos_orig)
    delete hmixSBPos_orig;
  if (hSB_tot_true)
    delete hSB_tot_true;
  if (h5_true)
    delete h5_true;
  if (h5_mix)
    delete h5_mix;

  return hCF_sb;
}

void draw_LambdaNuclearId_relamom_sum_compare1s2s() {
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

  const char *species[] = {"d", "t", "3He", "4He"};
  int nSpecies = 4;

  for (int i = 0; i < nSpecies; ++i) {
    TString sp = species[i];

    int nRebin = 1;
    if (sp == "t" || sp == "3He")
      nRebin = 4;
    else if (sp == "4He")
      nRebin = 10;

    TString outFile =
        Form("%s/kstar_sum_compare1s2s_%s.pdf", outDir.Data(), sp.Data());

    for (int groupSet = 0; groupSet < 2; ++groupSet) {
      int nGroups = (groupSet == 0) ? 2 : 3;

      // 6ページ目と同じ2行構成のキャンバスを作成 (上段: 全域, 下段: Zoom)
      TCanvas *c6 = new TCanvas(
          Form("c6_comp_%s_groupSet%d", sp.Data(), groupSet),
          Form("CF (SB Subtracted) Compare for Lambda - %s", sp.Data()), 1200,
          1200);
      c6->Divide(nGroups, 2);

      for (int group = 0; group < nGroups; ++group) {
        TString centLabel;
        int jStart, jEnd;
        getCentInfo(groupSet, group, centLabel, jStart, jEnd);

        TH1 *hCF_2s =
            buildCF_sb(fIn2s, sp, nRebin, groupSet, group, jStart, jEnd, "_2s");
        TH1 *hCF_1s =
            buildCF_sb(fIn1s, sp, nRebin, groupSet, group, jStart, jEnd, "_1s");

        // --- 上段: 全域 [0.0, 0.6] ---
        c6->cd(group + 1);
        double yMax_full = (sp == "d") ? 25.0 : 2.;
        if (hCF_2s) {
          hCF_2s->SetTitle(
              Form("CF (SB Subtracted) %s (%s)", sp.Data(), centLabel.Data()));
          hCF_2s->SetLineColor(kBlack);
          hCF_2s->SetMarkerColor(kBlack);
          hCF_2s->SetMarkerStyle(20);
          hCF_2s->SetMarkerSize(0.8);
          hCF_2s->GetYaxis()->SetTitle("C(k*) = True_sub / Mixed_sub");
          hCF_2s->GetXaxis()->SetRangeUser(0.0, 0.6);
          hCF_2s->GetYaxis()->SetRangeUser(0.0, yMax_full);
          hCF_2s->Draw("E1");
        }
        if (hCF_1s) {
          hCF_1s->SetLineColor(kRed);
          hCF_1s->SetMarkerColor(kRed);
          hCF_1s->SetMarkerStyle(24);
          hCF_1s->SetMarkerSize(0.8);
          hCF_1s->GetXaxis()->SetRangeUser(0.0, 0.6);
          hCF_1s->GetYaxis()->SetRangeUser(0.0, yMax_full);
          if (!hCF_2s) {
            hCF_1s->SetTitle(Form("CF (SB Subtracted) %s (%s)", sp.Data(),
                                  centLabel.Data()));
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

        // --- 下段: Zoom [0.0, 0.15] ---
        c6->cd(group + 1 + nGroups);
        TH1 *hCF_2s_zoom =
            hCF_2s ? (TH1 *)hCF_2s->Clone(Form("hCF_2s_zoom_%s_GS%d_G%d",
                                               sp.Data(), groupSet, group))
                   : 0;
        TH1 *hCF_1s_zoom =
            hCF_1s ? (TH1 *)hCF_1s->Clone(Form("hCF_1s_zoom_%s_GS%d_G%d",
                                               sp.Data(), groupSet, group))
                   : 0;

        double yMax_zoom = (sp == "d") ? 30.0 : ((sp == "4He") ? 5.0 : 3.0);
        if (hCF_2s_zoom) {
          hCF_2s_zoom->SetTitle(Form("CF Zoom (SB Subtracted) %s (%s)",
                                     sp.Data(), centLabel.Data()));
          hCF_2s_zoom->GetXaxis()->SetRangeUser(0.0, 0.15);
          hCF_2s_zoom->GetYaxis()->SetRangeUser(0.0, yMax_zoom);
          hCF_2s_zoom->SetLineColor(kBlack);
          hCF_2s_zoom->SetMarkerColor(kBlack);
          hCF_2s_zoom->SetMarkerStyle(20);
          hCF_2s_zoom->SetMarkerSize(0.8);
          hCF_2s_zoom->Draw("E1");
        }
        if (hCF_1s_zoom) {
          hCF_1s_zoom->GetXaxis()->SetRangeUser(0.0, 0.15);
          hCF_1s_zoom->GetYaxis()->SetRangeUser(0.0, yMax_zoom);
          hCF_1s_zoom->SetLineColor(kRed);
          hCF_1s_zoom->SetMarkerColor(kRed);
          hCF_1s_zoom->SetMarkerStyle(24);
          hCF_1s_zoom->SetMarkerSize(0.8);
          if (!hCF_2s_zoom) {
            hCF_1s_zoom->SetTitle(Form("CF Zoom (SB Subtracted) %s (%s)",
                                       sp.Data(), centLabel.Data()));
            hCF_1s_zoom->Draw("E1");
          } else {
            hCF_1s_zoom->Draw("E1 same");
          }
        }

        TLine *lZ = new TLine(0.0, 1.0, 0.15, 1.0);
        lZ->SetLineStyle(2);
        lZ->Draw("SAME");

        if (group == 0) {
          TLegend *legZ = new TLegend(0.40, 0.72, 0.88, 0.88);
          if (hCF_2s_zoom)
            legZ->AddEntry(hCF_2s_zoom, "p vs dE/dx cut: #pm2#sigma", "pe");
          if (hCF_1s_zoom)
            legZ->AddEntry(hCF_1s_zoom, "p vs dE/dx cut: #pm1#sigma", "pe");
          legZ->Draw();
        }
      } // end group loop

      if (groupSet == 0) {
        c6->Print(outFile + "(");
      } else {
        c6->Print(outFile + ")");
      }

      delete c6;
    } // end groupSet loop

    std::cout << "Successfully saved: " << outFile << std::endl;
  } // end species loop

  fIn2s->Close();
  fIn1s->Close();
}
