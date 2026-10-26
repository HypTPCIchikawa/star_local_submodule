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
// draw_LambdaNuclearId_relamom_sum_CFsub_compare1s2s.C
//
// 2sigma cut (auau3p85_all1.root, Black) vs
// 1sigma cut (auau3p85_all1_Nuclear1sigma.root, Red)
// の CF (Constant Purity Corrected) を重ね描き比較するマクロ。
// 出力PDF: pdf/kstar_sum_CFsub_compare1s2s_<species>.pdf
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

TH1 *buildCF_CFsub(TFile *fIn, const TString &sp, int nRebin, int groupSet,
                   int group, int jStart, int jEnd, const TString &sfx,
                   double sn_for_purity = 2.172) {
  TDirectory *dirNuc = (TDirectory *)fIn->Get("true");
  TDirectory *dirMix = (TDirectory *)fIn->Get("mix");

  TH1 *h = 0, *hmix = 0;
  TH1 *hSBNeg_orig = 0, *hSBPos_orig = 0, *hmixSBNeg_orig = 0,
      *hmixSBPos_orig = 0;

  for (int j = jStart; j <= jEnd; ++j) {
    TString histName = Form("hKstar_%s_CentBin%d", sp.Data(), j);
    TH1 *t = dirNuc ? (TH1 *)dirNuc->Get(histName) : (TH1 *)fIn->Get(histName);
    if (t) {
      if (!h) {
        h = (TH1 *)t->Clone(
            Form("h_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!h->GetSumw2N())
          h->Sumw2();
      } else {
        TH1 *a = (TH1 *)t->Clone();
        if (!a->GetSumw2N())
          a->Sumw2();
        h->Add(a);
        delete a;
      }
    }

    TString mixName = Form("hKstar_Mixed_%s_CentBin%d", sp.Data(), j);
    TH1 *m = dirMix ? (TH1 *)dirMix->Get(mixName) : (TH1 *)fIn->Get(mixName);
    if (m) {
      if (!hmix) {
        hmix = (TH1 *)m->Clone(
            Form("hmix_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hmix->GetSumw2N())
          hmix->Sumw2();
      } else {
        TH1 *a = (TH1 *)m->Clone();
        if (!a->GetSumw2N())
          a->Sumw2();
        hmix->Add(a);
        delete a;
      }
    }

    TString nameSBNeg = Form("hKstar_%s_SBNeg_CentBin%d", sp.Data(), j);
    TH1 *sn_in =
        dirNuc ? (TH1 *)dirNuc->Get(nameSBNeg) : (TH1 *)fIn->Get(nameSBNeg);
    if (sn_in) {
      if (!hSBNeg_orig) {
        hSBNeg_orig = (TH1 *)sn_in->Clone(Form(
            "hSBNeg_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hSBNeg_orig->GetSumw2N())
          hSBNeg_orig->Sumw2();
      } else {
        TH1 *a = (TH1 *)sn_in->Clone();
        if (!a->GetSumw2N())
          a->Sumw2();
        hSBNeg_orig->Add(a);
        delete a;
      }
    }

    TString nameSBPos = Form("hKstar_%s_SBPos_CentBin%d", sp.Data(), j);
    TH1 *sp_in =
        dirNuc ? (TH1 *)dirNuc->Get(nameSBPos) : (TH1 *)fIn->Get(nameSBPos);
    if (sp_in) {
      if (!hSBPos_orig) {
        hSBPos_orig = (TH1 *)sp_in->Clone(Form(
            "hSBPos_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hSBPos_orig->GetSumw2N())
          hSBPos_orig->Sumw2();
      } else {
        TH1 *a = (TH1 *)sp_in->Clone();
        if (!a->GetSumw2N())
          a->Sumw2();
        hSBPos_orig->Add(a);
        delete a;
      }
    }

    TString mixNameSBNeg =
        Form("hKstar_Mixed_%s_SBNeg_CentBin%d", sp.Data(), j);
    TH1 *msn_in = dirMix ? (TH1 *)dirMix->Get(mixNameSBNeg)
                         : (TH1 *)fIn->Get(mixNameSBNeg);
    if (msn_in) {
      if (!hmixSBNeg_orig) {
        hmixSBNeg_orig = (TH1 *)msn_in->Clone(Form(
            "hmixSBNeg_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hmixSBNeg_orig->GetSumw2N())
          hmixSBNeg_orig->Sumw2();
      } else {
        TH1 *a = (TH1 *)msn_in->Clone();
        if (!a->GetSumw2N())
          a->Sumw2();
        hmixSBNeg_orig->Add(a);
        delete a;
      }
    }

    TString mixNameSBPos =
        Form("hKstar_Mixed_%s_SBPos_CentBin%d", sp.Data(), j);
    TH1 *msp_in = dirMix ? (TH1 *)dirMix->Get(mixNameSBPos)
                         : (TH1 *)fIn->Get(mixNameSBPos);
    if (msp_in) {
      if (!hmixSBPos_orig) {
        hmixSBPos_orig = (TH1 *)msp_in->Clone(Form(
            "hmixSBPos_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
        if (!hmixSBPos_orig->GetSumw2N())
          hmixSBPos_orig->Sumw2();
      } else {
        TH1 *a = (TH1 *)msp_in->Clone();
        if (!a->GetSumw2N())
          a->Sumw2();
        hmixSBPos_orig->Add(a);
        delete a;
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

  int bin1 = h->FindBin(0.40001);
  int bin2 = h->FindBin(0.59999);
  double intTrue = h->Integral(bin1, bin2);
  double intMix = hmix->Integral(bin1, bin2);
  if (intMix > 0) {
    hmix->Scale(intTrue / intMix);
  }

  TH1 *hCF_sig = (TH1 *)h->Clone(
      Form("hCF_sig_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
  hCF_sig->Divide(hmix);

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

  TH1 *hCF_bg = 0;
  if (hSB_tot_true && hmixSB_tot) {
    double intSB_true = hSB_tot_true->Integral(bin1, bin2);
    double intSB_mix = hmixSB_tot->Integral(bin1, bin2);
    if (intSB_mix > 0) {
      hmixSB_tot->Scale(intSB_true / intSB_mix);
    }
    hCF_bg = (TH1 *)hSB_tot_true->Clone(
        Form("hCF_bg_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
    hCF_bg->Divide(hmixSB_tot);
  }

  TH1 *hCF_sub = 0;
  if (hCF_sig && hCF_bg) {
    double purity = sn_for_purity / (sn_for_purity + 1.0);
    hCF_sub = (TH1 *)hCF_sig->Clone(
        Form("hCF_sub_%s_GS%d_G%d%s", sp.Data(), groupSet, group, sfx.Data()));
    if (purity > 0) {
      TH1 *hTerm2 = (TH1 *)hCF_bg->Clone("hTerm2_tmp");
      hTerm2->Scale(1.0 - purity);
      hCF_sub->Add(hTerm2, -1.0);
      hCF_sub->Scale(1.0 / purity);
      delete hTerm2;
    }
  }

  delete h;
  delete hmix;
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
  if (hmixSB_tot)
    delete hmixSB_tot;
  if (hCF_sig)
    delete hCF_sig;
  if (hCF_bg)
    delete hCF_bg;

  return hCF_sub;
}

void draw_LambdaNuclearId_relamom_sum_CFsub_compare1s2s() {
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

  double sn_2s = 2.172;
  double sn_1s = 2.172;

  for (int i = 0; i < nSpecies; ++i) {
    TString sp = species[i];

    int nRebin = 1;
    if (sp == "t" || sp == "3He")
      nRebin = 4;
    else if (sp == "4He")
      nRebin = 10;

    TString outFile =
        Form("%s/kstar_sum_CFsub_compare1s2s_%s.pdf", outDir.Data(), sp.Data());

    for (int groupSet = 0; groupSet < 2; ++groupSet) {
      int nGroups = (groupSet == 0) ? 2 : 3;

      TCanvas *c6 = new TCanvas(
          Form("c6_cfsub_comp_%s_groupSet%d", sp.Data(), groupSet),
          Form("CF (Constant Purity Corrected) Compare for Lambda - %s",
               sp.Data()),
          1200, 1200);
      c6->Divide(nGroups, 2);

      for (int group = 0; group < nGroups; ++group) {
        TString centLabel;
        int jStart, jEnd;
        getCentInfo(groupSet, group, centLabel, jStart, jEnd);

        TH1 *hCF_2s = buildCF_CFsub(fIn2s, sp, nRebin, groupSet, group, jStart,
                                    jEnd, "_2s", sn_2s);
        TH1 *hCF_1s = buildCF_CFsub(fIn1s, sp, nRebin, groupSet, group, jStart,
                                    jEnd, "_1s", sn_1s);

        // --- 上段: 全域 [0.0, 0.6] ---
        c6->cd(group + 1);
        double yMax_full = (sp == "d") ? 25.0 : 2.;
        if (hCF_2s) {
          hCF_2s->SetTitle(Form("CF (Purity Corrected) %s (%s)", sp.Data(),
                                centLabel.Data()));
          hCF_2s->SetLineColor(kBlack);
          hCF_2s->SetMarkerColor(kBlack);
          hCF_2s->SetMarkerStyle(20);
          hCF_2s->SetMarkerSize(0.8);
          hCF_2s->GetYaxis()->SetTitle("CF (Constant Purity Corrected)");
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
            hCF_1s->SetTitle(Form("CF (Purity Corrected) %s (%s)", sp.Data(),
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
          hCF_2s_zoom->SetTitle(Form("CF Zoom (Purity Corrected) %s (%s)",
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
            hCF_1s_zoom->SetTitle(Form("CF Zoom (Purity Corrected) %s (%s)",
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
