// Generated adaptation; legacy original preserved.
#include "KFSummarySupport.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include <iostream>


const double hbarc = 197.3269804; // MeV fm

// Dawson integral: D(x) = exp(-x^2) * integral_0^x exp(t^2) dt
// Computed using a stable series expansion to avoid overflow
double Dawson(double x) {
  // For small x: use Taylor series D(x) = x - 2x^3/3 + 4x^5/15 - ...
  // For large x: use asymptotic D(x) ~ 1/(2x) * (1 + 1/(2x^2) + ...)
  // For intermediate x: numerical integration with substitution t = x*u
  //   D(x) = x * integral_0^1 exp(x^2*(u^2-1)) du
  //   integrand is always <= 1, numerically stable!
  if (x < 0)
    return -Dawson(-x);

  if (x < 1e-6)
    return x;

  if (x > 20.0) {
    // Asymptotic: D(x) = 1/(2x) + 1/(4x^3) + ...
    return 0.5 / x * (1.0 + 0.5 / (x * x));
  }

  // Use substitution t = x*u, so exp(t^2 - x^2) = exp(x^2*(u^2-1))
  // which is always in (-inf, 0] -> exp value in [0,1]
  int n = 2000;
  double h = 1.0 / n;
  double sum = 0.0;
  // Simpson's rule on integrand exp(x^2*(u^2-1))
  for (int i = 0; i <= n; i++) {
    double u = i * h;
    double y = exp(x * x * (u * u - 1.0));
    double weight = (i == 0 || i == n) ? 1.0 : ((i % 2 != 0) ? 4.0 : 2.0);
    sum += weight * y;
  }
  sum *= h / 3.0;
  return x * sum;
}

// F1(x) = D(x)/x  where D(x) is the Dawson integral
// = e^{-x^2} * integral_0^x e^{t^2} dt / x
double F1(double x) {
  if (x < 1e-6)
    return 1.0;
  return Dawson(x) / x;
}

// F2(x) = (1 - e^{-x^2}) / x
double F2(double x) {
  if (x < 1e-6)
    return x;
  return (1.0 - exp(-x * x)) / x;
}

// F3 correction: 1 - r_eff / (2*sqrt(pi)*R)
// Called as F3(reff/R) = 1 - (reff/R)/(2*sqrt(pi))
double F3(double x) { return 1.0 - x / (2.0 * sqrt(TMath::Pi())); }

// Lednicky-Lyuboshits Formula
double LL_Formula(double q_MeV, double R, double a0, double reff) {
  double q = q_MeV / hbarc; // convert to fm^-1

  double qcotdelta = -1.0 / a0 + 0.5 * reff * q * q;
  double denom = qcotdelta * qcotdelta + q * q;

  double mag_f_sq = (denom > 0) ? 1.0 / denom : a0 * a0;
  double re_f = (denom > 0) ? qcotdelta / denom : -a0;
  double im_f = (denom > 0) ? q / denom : 0.0;

  double term1 = mag_f_sq / (2.0 * R * R) * F3(reff / R);
  double term2 = 2.0 * re_f / (sqrt(TMath::Pi()) * R) * F1(2.0 * q * R);
  double term3 = im_f / R * F2(2.0 * q * R);

  return 1.0 + term1 + term2 - term3;
}



// Wrapper for TF1: a0 and reff are fixed
double fitLL_fixed(double *x, double *p) {
    double q_GeV = x[0];
    double q_MeV = q_GeV * 1000.0;
    if (q_MeV < 1e-3) q_MeV = 1e-3;
    
    double R = p[0];
    double norm = p[1];
    
    double a0 = 4.0099999999999998;   // fixed
    double reff = 1.8400000000000001; // fixed
    
    return norm * LL_Formula(q_MeV, R, a0, reff);
}

// Wrapper for TF1: a0, reff, R are free
double fitLL_free(double *x, double *p) {
    double q_GeV = x[0];
    double q_MeV = q_GeV * 1000.0;
    if (q_MeV < 1e-3) q_MeV = 1e-3;
    
    double R = p[0];
    double a0 = p[1];
    double reff = p[2];
    double norm = p[3];
    
    return norm * LL_Formula(q_MeV, R, a0, reff);
}

#include "TLatex.h"



void getPurity(TH1D* hY, double minM, double maxM, double &P, double &err_P,
  double &N_sig, double &err_N_sig, double &N_tot, double &err_N_tot) {
  KFReport::Purity(hY,minM,maxM,P,err_P,N_sig,err_N_sig,N_tot,err_N_tot);
}

void Fit_LambdaNuclearId_relamom_sum_CFsub_method3() {
 const double massFitSettings[9] = {1.0800000000000001,1.1499999999999999,1.1100000000000001,1.1200000000000001,0.001,0.0080000000000000002,0.0030000000000000001,10,2};
 KFReport::Init("/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/logs/Fit_LambdaNuclearId_relamom_sum_CFsub_method3_fit_audit.tsv",1.1159600000000001,0.00296,massFitSettings);
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
        if (sp == "t" || sp == "3He") nRebin = 4;
        else if (sp == "4He") nRebin = 4;
        
        TString outFile = Form("%s/Fit_LL_method3_%s.pdf", outDir.Data(), sp.Data());
        
        TCanvas *c_start = new TCanvas("c_start");
        c_start->Print(outFile + "[");

        for (int groupSet = 0; groupSet < 2; ++groupSet) {
            int nGroups = (groupSet == 0) ? 2 : 3;

            TCanvas *cFit = new TCanvas(Form("cCF_%s_groupSet%d", sp.Data(), groupSet), Form("CF(k*) - %s", sp.Data()), 1200, 600);
            cFit->Divide(nGroups, 1);
            
            TCanvas *cPurity = new TCanvas(Form("cPurity_%s_groupSet%d", sp.Data(), groupSet), Form("Purity(k*) - %s", sp.Data()), 1200, 600);
            cPurity->Divide(nGroups, 1);
            
            TCanvas *cNSig = new TCanvas(Form("cNSig_%s_groupSet%d", sp.Data(), groupSet), Form("NSig(k*) - %s", sp.Data()), 1200, 600);
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
                    TString hNameMix = Form("hKstarMass_Mixed_%s_CentBin%d", sp.Data(), j);
                    
                    TH2 *h2True = 0;
                    TH2 *h2Mix = 0;
                    
                    if (dirNuc) h2True = (TH2*)dirNuc->Get(hNameTrue);
                    if (dirMix) h2Mix = (TH2*)dirMix->Get(hNameMix);
                    
                    if (!h2True) h2True = (TH2*)fIn->Get(hNameTrue);
                    if (!h2Mix) h2Mix = (TH2*)fIn->Get(hNameMix);

                    if (h2True) {
                        if (!h2True_sum) {
                            h2True_sum = (TH2*)h2True->Clone(Form("h2True_sum_%s_GS%d_G%d", sp.Data(), groupSet, group));
                        } else {
                            h2True_sum->Add(h2True);
                        }
                    }
                    if (h2Mix) {
                        if (!h2Mix_sum) {
                            h2Mix_sum = (TH2*)h2Mix->Clone(Form("h2Mix_sum_%s_GS%d_G%d", sp.Data(), groupSet, group));
                        } else {
                            h2Mix_sum->Add(h2Mix);
                        }
                    }
                }

                if (!h2True_sum || !h2Mix_sum) continue;
                
                h2True_sum->RebinX(nRebin);
                h2Mix_sum->RebinX(nRebin);
                
                int nBinsX = h2True_sum->GetNbinsX();
                
                TH1D *hCF = new TH1D(Form("hCF_%s_GS%d_G%d", sp.Data(), groupSet, group), Form("%s %s CF;k* (GeV/c);CF", sp.Data(), centLabel.Data()), nBinsX, h2True_sum->GetXaxis()->GetXmin(), h2True_sum->GetXaxis()->GetXmax());
                TH1D *hPTrue = new TH1D(Form("hPTrue_%s_GS%d_G%d", sp.Data(), groupSet, group), Form("%s %s Purity;k* (GeV/c);Purity", sp.Data(), centLabel.Data()), nBinsX, h2True_sum->GetXaxis()->GetXmin(), h2True_sum->GetXaxis()->GetXmax());
                TH1D *hPMix = new TH1D(Form("hPMix_%s_GS%d_G%d", sp.Data(), groupSet, group), Form("%s %s Purity;k* (GeV/c);Purity", sp.Data(), centLabel.Data()), nBinsX, h2True_sum->GetXaxis()->GetXmin(), h2True_sum->GetXaxis()->GetXmax());
                
                TH1D *hSigTrue = new TH1D(Form("hSigTrue_%s_GS%d_G%d", sp.Data(), groupSet, group), Form("%s %s N_{sig};k* (GeV/c);Counts", sp.Data(), centLabel.Data()), nBinsX, h2True_sum->GetXaxis()->GetXmin(), h2True_sum->GetXaxis()->GetXmax());
                TH1D *hSigMix = new TH1D(Form("hSigMix_%s_GS%d_G%d", sp.Data(), groupSet, group), Form("%s %s N_{sig};k* (GeV/c);Counts", sp.Data(), centLabel.Data()), nBinsX, h2True_sum->GetXaxis()->GetXmin(), h2True_sum->GetXaxis()->GetXmax());
                
                for (int bx = 1; bx <= nBinsX; ++bx) {
                    TH1D *hYTrue = h2True_sum->ProjectionY(Form("hYTrue_%s_GS%d_G%d_B%d", sp.Data(), groupSet, group, bx), bx, bx);
                    TH1D *hYMix = h2Mix_sum->ProjectionY(Form("hYMix_%s_GS%d_G%d_B%d", sp.Data(), groupSet, group, bx), bx, bx);
                    
                    double P_t, err_P_t, N_sig_t, err_N_sig_t, N_tot_t, err_N_tot_t;
                    double P_m, err_P_m, N_sig_m, err_N_sig_m, N_tot_m, err_N_tot_m;
                    
                    getPurity(hYTrue, minM, maxM, P_t, err_P_t, N_sig_t, err_N_sig_t, N_tot_t, err_N_tot_t);
                    getPurity(hYMix, minM, maxM, P_m, err_P_m, N_sig_m, err_N_sig_m, N_tot_m, err_N_tot_m);
                    
                    hPTrue->SetBinContent(bx, P_t);
                    hPTrue->SetBinError(bx, err_P_t);
                    hPMix->SetBinContent(bx, P_m);
                    hPMix->SetBinError(bx, err_P_m);
                    
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
                hSigTrue->SetTitle(Form("%s (%s) N_{sig}", sp.Data(), centLabel.Data()));
                hSigTrue->SetMarkerStyle(20);
                hSigTrue->SetMarkerColor(kRed);
                hSigTrue->SetLineColor(kRed);
                hSigTrue->GetXaxis()->SetRangeUser(0,0.29999999999999999);
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
                hPTrue->GetXaxis()->SetRangeUser(0,0.29999999999999999);
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
                
                cFit->cd(group + 1);
                hCF->SetTitle(Form("%s (%s) CF", sp.Data(), centLabel.Data()));
                hCF->SetMarkerStyle(20);
                hCF->SetMarkerColor(kBlack);
                hCF->SetLineColor(kBlack);
                hCF->SetMinimum(0.0);
                hCF->SetMaximum(3.0);
                hCF->GetXaxis()->SetRangeUser(0,0.29999999999999999);
                KFReport::DrawMeasured(hCF,"E1");
                // --- FITTING ---
                double fit_max = 0.40000000000000002;
                
                TF1 *fit1 = new TF1(Form("fit1_%d_%d", groupSet, group), fitLL_fixed, 0, fit_max, 2);
                fit1->SetParameters(3,1);
                fit1->SetParNames("R", "Norm");
                fit1->SetParLimits(0,0.10000000000000001,6);
                fit1->SetLineColor(kBlue);
                KFReport::FitLL(hCF,fit1,"R");
                KFReport::FitLL(hCF,fit1,"R");
                
                TF1 *fit2 = new TF1(Form("fit2_%d_%d", groupSet, group), fitLL_free, 0, fit_max, 4);
                fit2->SetParameters(3,4.0099999999999998,1.8400000000000001,1);
                fit2->SetParNames("R", "a_{0}", "r_{eff}", "Norm");
                fit2->SetParLimits(0,0.10000000000000001,6);
                fit2->SetLineColor(kRed);
                KFReport::FitLL(hCF,fit2,"R+");
                KFReport::FitLL(hCF,fit2,"R+");
                
                TLegend *leg = new TLegend(0.52,0.76,0.88,0.89);
                leg->SetBorderSize(0);
                leg->AddEntry(hCF, "Data", "pe");
                leg->AddEntry(fit1, "LL: fixed", "l");
                leg->AddEntry(fit2, "LL: free", "l");
                leg->Draw();

                KFReport::DrawLLWarnings(fit1,fit2);
    
                
                TLine *l1 = new TLine(0,1.0,0.29999999999999999,1.0);
                l1->SetLineStyle(2);
                l1->Draw("SAME");
            } // end group loop
            
            //cNSig->Print(outFile);
            //cPurity->Print(outFile);
            cFit->Print(outFile);
        } // end groupSet loop
        
        TCanvas *c_end = new TCanvas("c_end");
        c_end->Print(outFile + "]");
        
        std::cout << "Successfully created " << outFile << std::endl;
    } // end species loop
}
