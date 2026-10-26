// Generated adaptation; legacy original preserved.
#include "KFSummarySupport.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"
#include "TLine.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include "TPaveText.h"
#include <iostream>

void getPurity(TH1D* hY, double minM, double maxM, double &P, double &err_P,
  double &N_sig, double &err_N_sig, double &N_tot, double &err_N_tot) {
  KFReport::Purity(hY,minM,maxM,P,err_P,N_sig,err_N_sig,N_tot,err_N_tot);
}

void draw_LambdaNuclearId_relamom_CFsub_method3() {
 const double massFitSettings[9] = {1.0800000000000001,1.1499999999999999,1.1100000000000001,1.1200000000000001,0.001,0.0080000000000000002,0.0030000000000000001,10,2};
 KFReport::Init("/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/logs/draw_LambdaNuclearId_relamom_CFsub_method3_fit_audit.tsv",1.1159600000000001,0.00296,massFitSettings);
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
    TString outFile = outDir + "/kstar_cent_CFsub_method3_4He.pdf";
    
    TCanvas *c_start = new TCanvas("c_start");
    c_start->Print(outFile + "[");

    const char *species[] = {"4He"};
    int nSpecies = 1;
    int nCentBins = 9;
    const char *centLabels[] = {"70-80%", "60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%", "5-10%", "0-5%"};

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
        
        TCanvas *cCF = new TCanvas(Form("cCF_%s", sp.Data()), Form("CF(k*) - %s", sp.Data()), 1200, 900);
        cCF->Divide(3, 3);
        
        TCanvas *cPurity = new TCanvas(Form("cPurity_%s", sp.Data()), Form("Purity(k*) - %s", sp.Data()), 1200, 900);
        cPurity->Divide(3, 3);
        
        TCanvas *cNSig = new TCanvas(Form("cNSig_%s", sp.Data()), Form("NSig(k*) - %s", sp.Data()), 1200, 900);
        cNSig->Divide(3, 3);

        for (int j = 0; j < nCentBins; ++j) {
            TH2 *h2True = 0;
            TH2 *h2Mix = 0;
            
            TString hNameTrue = Form("hKstarMass_%s_CentBin%d", sp.Data(), j);
            TString hNameMix = Form("hKstarMass_Mixed_%s_CentBin%d", sp.Data(), j);
            
            if (dirNuc) h2True = (TH2*)dirNuc->Get(hNameTrue);
            if (dirMix) h2Mix = (TH2*)dirMix->Get(hNameMix);
            
            if (!h2True) h2True = (TH2*)fIn->Get(hNameTrue);
            if (!h2Mix) h2Mix = (TH2*)fIn->Get(hNameMix);

            if (!h2True || !h2Mix) continue;
            
            h2True->RebinX(nRebin);
            h2Mix->RebinX(nRebin);
            
            int nBinsX = h2True->GetNbinsX();
            
            TH1D *hCF = new TH1D(Form("hCF_%s_Cent%d", sp.Data(), j), Form("%s CentBin%d CF;k* (GeV/c);CF", sp.Data(), j), nBinsX, h2True->GetXaxis()->GetXmin(), h2True->GetXaxis()->GetXmax());
            TH1D *hPTrue = new TH1D(Form("hPTrue_%s_Cent%d", sp.Data(), j), Form("%s CentBin%d Purity;k* (GeV/c);Purity", sp.Data(), j), nBinsX, h2True->GetXaxis()->GetXmin(), h2True->GetXaxis()->GetXmax());
            TH1D *hPMix = new TH1D(Form("hPMix_%s_Cent%d", sp.Data(), j), Form("%s CentBin%d Purity;k* (GeV/c);Purity", sp.Data(), j), nBinsX, h2True->GetXaxis()->GetXmin(), h2True->GetXaxis()->GetXmax());
            
            TH1D *hSigTrue = new TH1D(Form("hSigTrue_%s_Cent%d", sp.Data(), j), Form("%s CentBin%d N_{sig};k* (GeV/c);Counts", sp.Data(), j), nBinsX, h2True->GetXaxis()->GetXmin(), h2True->GetXaxis()->GetXmax());
            TH1D *hSigMix = new TH1D(Form("hSigMix_%s_Cent%d", sp.Data(), j), Form("%s CentBin%d N_{sig};k* (GeV/c);Counts", sp.Data(), j), nBinsX, h2True->GetXaxis()->GetXmin(), h2True->GetXaxis()->GetXmax());
            
            for (int bx = 1; bx <= nBinsX; ++bx) {
                TH1D *hYTrue = h2True->ProjectionY(Form("hYTrue_%s_C%d_B%d", sp.Data(), j, bx), bx, bx);
                TH1D *hYMix = h2Mix->ProjectionY(Form("hYMix_%s_C%d_B%d", sp.Data(), j, bx), bx, bx);
                
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
            
            int bin1 = hSigTrue->FindBin(0.30001);
            int bin2 = hSigTrue->FindBin(0.39999000000000001);
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
            
            cNSig->cd(j + 1);
            hSigTrue->SetTitle(Form("%s CentBin%d (%s) N_{sig}", sp.Data(), j, centLabels[j]));
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
            
            if (j == 0) {
                TLegend *legSig = new TLegend(0.4, 0.75, 0.88, 0.88);
                legSig->AddEntry(hSigTrue, "True N_{sig}", "pl");
                legSig->AddEntry(hSigMix, "Mixed N_{sig} (Scaled)", "pl");
                legSig->Draw();
            }
            
            cPurity->cd(j + 1);
            hPTrue->SetTitle(Form("%s CentBin%d (%s) Purity", sp.Data(), j, centLabels[j]));
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
            
            if (j == 0) {
                TLegend *legP = new TLegend(0.4, 0.75, 0.88, 0.88);
                legP->AddEntry(hPTrue, "True Purity", "pl");
                legP->AddEntry(hPMix, "Mixed Purity", "pl");
                legP->Draw();
            }
            
            cCF->cd(j + 1);
            hCF->SetTitle(Form("%s CentBin%d (%s) CF", sp.Data(), j, centLabels[j]));
            hCF->SetMarkerStyle(20);
            hCF->SetMarkerColor(kBlack);
            hCF->SetLineColor(kBlack);
            hCF->SetMinimum(0.0);
            hCF->SetMaximum(3.0);
            hCF->GetXaxis()->SetRangeUser(0,0.29999999999999999);
            KFReport::DrawMeasured(hCF,"E1");
            
            TLine *l1 = new TLine(0,1.0,0.29999999999999999,1.0);
            l1->SetLineStyle(2);
            l1->Draw("SAME");
        }
        
        cNSig->Print(outFile);
        cPurity->Print(outFile);
        cCF->Print(outFile);
    }

    TCanvas *c_end = new TCanvas("c_end");
    c_end->Print(outFile + "]");
    
    std::cout << "Successfully created " << outFile << std::endl;
}
