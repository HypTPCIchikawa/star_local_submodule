// Display-only redraw: no analysis, cut changes, rebinning or smoothing.
// Read complete existing count histograms; normalize by FULL ordinary-bin sum.
// Arguments control the visible range and the unchanged counting window.
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TColor.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

namespace {
Bool_t LroAligned(TH1D* h, Double_t low, Double_t high, Int_t& first, Int_t& last) {
  const Double_t width=h->GetBinWidth(1);
  if (!(high>low) || low<h->GetXaxis()->GetXmin()-1.e-9 ||
      high>h->GetXaxis()->GetXmax()+1.e-9) return kFALSE;
  first=h->FindBin(low+width*1.e-6);
  last=h->FindBin(high-width*1.e-6);
  return first>=1 && last<=h->GetNbinsX() && first<=last &&
    TMath::Abs(h->GetBinLowEdge(first)-low)<1.e-9 &&
    TMath::Abs(h->GetBinLowEdge(last)+h->GetBinWidth(last)-high)<1.e-9;
}
Int_t LroFail(const char* message) {
  std::cerr<<"[drawLambdaOverlayReadable] ERROR: "<<message<<std::endl;
  return 1;
}
}

Int_t drawLambdaOverlayReadable(const char* inputRoot, const char* baselineKey,
    const char* selectedKey, const char* selectedLabel, const char* outputStem,
    Double_t displayLow, Double_t displayHigh,
    Double_t windowLow, Double_t windowHigh) {
  if (!inputRoot || !baselineKey || !selectedKey || !selectedLabel ||
      !outputStem || !*outputStem) return LroFail("missing argument");
  const char* extensions[4]={".pdf",".png",".root","_display.txt"};
  for (Int_t j=0;j<4;++j)
    if (!gSystem->AccessPathName(TString(outputStem)+extensions[j]))
      return LroFail("refusing to overwrite an existing output");
  TFile input(inputRoot,"READ");
  if (input.IsZombie() || input.TestBit(TFile::kRecovered))
    return LroFail("input must be a complete, non-recovered ROOT file");

  const char* keys[3]={"HelixLambdaCounts",baselineKey,selectedKey};
  const char* aliases[3]={"HelixCounts","KfBaselineCounts","KfSelectedCounts"};
  TH1D* sources[3]={0,0,0};
  TH1D* counts[3]={0,0,0};
  TH1D* draw[2][3];
  Double_t norm[3],nwin[3],maxima[2]={0.,0.};
  Int_t first=0,last=0,wfirst=0,wlast=0;
  for (Int_t j=0;j<3;++j) {
    sources[j]=dynamic_cast<TH1D*>(input.Get(keys[j]));
    if (!sources[j]) return LroFail("missing TH1D input key");
    counts[j]=dynamic_cast<TH1D*>(sources[j]->Clone(aliases[j]));
    counts[j]->SetDirectory(0);
    counts[j]->GetXaxis()->SetRange(1,counts[j]->GetNbinsX());
    if (!LroAligned(counts[j],displayLow,displayHigh,first,last) ||
        !LroAligned(counts[j],windowLow,windowHigh,wfirst,wlast))
      return LroFail("display/count range must align with complete input bins");
    if (j && counts[j]->GetNbinsX()!=counts[0]->GetNbinsX())
      return LroFail("different bin counts");
    for (Int_t i=0;i<=counts[j]->GetNbinsX()+1;++i) {
      if (!TMath::Finite(counts[j]->GetBinContent(i)) ||
          !TMath::Finite(counts[j]->GetBinError(i)) ||
          counts[j]->GetBinContent(i)<0.) return LroFail("invalid count bin");
      if (j && (counts[j]->GetBinLowEdge(i)!=counts[0]->GetBinLowEdge(i) ||
                counts[j]->GetBinWidth(i)!=counts[0]->GetBinWidth(i)))
        return LroFail("different bin edges");
    }
    norm[j]=counts[j]->Integral(1,counts[j]->GetNbinsX());
    nwin[j]=counts[j]->Integral(wfirst,wlast);
    if (!(norm[j]>0.)) return LroFail("empty normalization domain");
    for (Int_t p=0;p<2;++p) {
      draw[p][j]=dynamic_cast<TH1D*>(counts[j]->Clone(Form("%s_%s",
                  p ? "Normalized" : "Absolute",aliases[j])));
      draw[p][j]->SetDirectory(0);
      if (p) draw[p][j]->Scale(1./norm[j]);
      // The explicit full-range denominator above is computed BEFORE zooming.
      draw[p][j]->GetXaxis()->SetRange(first,last);
      for (Int_t i=first;i<=last;++i)
        maxima[p]=TMath::Max(maxima[p],draw[p][j]->GetBinContent(i));
    }
  }
  const Double_t width=counts[0]->GetBinWidth(1);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0); gStyle->SetTextFont(42);
  gStyle->SetGridColor(TColor::GetColor("#e0e3e7"));
  gStyle->SetGridStyle(3); gStyle->SetGridWidth(1);
  const Int_t colors[3]={TColor::GetColor("#30343b"),
    TColor::GetColor("#0072b2"),TColor::GetColor("#d55e00")};
  const Int_t styles[3]={1,7,1};
  TCanvas* canvas=new TCanvas("LambdaReadableOverlay","Lambda comparison",1800,1000);
  canvas->SetFillColor(kWhite);
  TGraph* markers[2][3]={{0,0,0},{0,0,0}};
  for (Int_t p=0;p<2;++p) {
    canvas->cd();
    TPad* pad=new TPad(Form("ReadablePanel%d",p),"",p*0.5,0.06,(p+1)*0.5,0.755);
    pad->SetFillColor(kWhite); pad->SetLeftMargin(0.145); pad->SetRightMargin(0.045);
    pad->SetTopMargin(0.09); pad->SetBottomMargin(0.15);
    pad->SetTicks(1,1); pad->SetGridy(); pad->Draw(); pad->cd();
    TH1D* frame=dynamic_cast<TH1D*>(draw[p][0]->Clone(Form("Frame%d",p)));
    frame->SetDirectory(0); frame->Reset("ICES"); frame->SetTitle("");
    frame->SetMinimum(0.); frame->SetMaximum(maxima[p]*1.16);
    frame->GetXaxis()->SetRange(first,last);
    frame->GetXaxis()->SetTitle("M_{p#pi^{-}} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(Form(p ? "Fraction per %.3g MeV/c^{2} bin" :
                                      "Candidates per %.3g MeV/c^{2} bin",width*1000.));
    frame->GetXaxis()->SetNdivisions(505,kFALSE);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetTitleSize(0.05); frame->GetXaxis()->SetLabelSize(0.041);
    frame->GetYaxis()->SetTitleSize(0.046); frame->GetYaxis()->SetLabelSize(0.039);
    frame->GetXaxis()->SetTitleOffset(1.15); frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelFont(42); frame->GetXaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetLabelFont(42); frame->GetYaxis()->SetTitleFont(42);
    frame->Draw("AXIS");
    const Int_t order[3]={1,2,0};
    for (Int_t q=0;q<3;++q) {
      const Int_t j=order[q];
      draw[p][j]->SetLineColor(colors[j]); draw[p][j]->SetLineStyle(styles[j]);
      draw[p][j]->SetLineWidth(j==0 ? 1 : 2);
      draw[p][j]->SetFillStyle(0); draw[p][j]->SetMarkerColor(colors[j]);
      draw[p][j]->SetMarkerStyle(j==0 ? 20 : 25);
      draw[p][j]->SetMarkerSize(j==0 ? 0.65 : 1.);
      draw[p][j]->Draw("HIST SAME");
    }
    // Sparse markers help distinguish coincident steps without shifting masses.
    // Every original bin remains visible in each histogram line.
    for (Int_t j=0;j<3;++j) {
      if (j==1) continue;
      markers[p][j]=new TGraph();
      markers[p][j]->SetName(Form("Markers%d_%d",p,j));
      for (Int_t i=first;i<=last;i+=2) {
        if (draw[p][j]->GetBinContent(i)>0.)
          markers[p][j]->SetPoint(markers[p][j]->GetN(),
            draw[p][j]->GetBinCenter(i),draw[p][j]->GetBinContent(i));
      }
      markers[p][j]->SetMarkerColor(colors[j]); markers[p][j]->SetLineColor(colors[j]);
      markers[p][j]->SetMarkerStyle(j==0 ? 20 : 25);
      markers[p][j]->SetMarkerSize(j==0 ? 0.65 : 1.);
      markers[p][j]->Draw("P SAME");
    }
    TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.041);
    text.DrawLatex(0.145,0.95,p ? "Normalized shape" : "Absolute counts");
    pad->RedrawAxis();
  }
  canvas->cd();
  TLatex heading; heading.SetNDC(); heading.SetTextFont(42);
  heading.SetTextSize(0.031);
  heading.DrawLatex(0.045,0.951,"Au+Au 13.5 GeV  |  #Lambda #rightarrow p#pi^{-}  |  10,000 input events");
  TLegend* legend=new TLegend(0.045,0.815,0.975,0.925);
  legend->SetNColumns(1); legend->SetBorderSize(0); legend->SetFillStyle(0);
  legend->SetTextFont(42); legend->SetTextSize(0.026); legend->SetMargin(0.055);
  legend->AddEntry(draw[0][0],Form("Helix Imp5: Nwin = %.0f",nwin[0]),"lp");
  legend->AddEntry(draw[0][1],Form("KF baseline: Nwin = %.0f",nwin[1]),"l");
  legend->AddEntry(draw[0][2],Form("%s: Nwin = %.0f",selectedLabel,nwin[2]),"lp");
  legend->Draw();
  heading.SetTextSize(0.022);
  heading.DrawLatex(0.045,0.776,Form("Nwin: %.3f #leq M < %.3f GeV/c^{2} (including background).  Normalization unchanged: %.2f #leq M < %.2f.",
    windowLow,windowHigh,counts[0]->GetXaxis()->GetXmin(),counts[0]->GetXaxis()->GetXmax()));
  heading.SetTextSize(0.019);
  heading.DrawLatex(0.075,0.025,Form("Display: %.2f-%.2f GeV/c^{2}.  Original %.3g MeV bins; no rebinning, smoothing or offsets.  Markers every 2 bins; lines show all bins.",
    displayLow,displayHigh,width*1000.));
  canvas->Modified(); canvas->Update();

  std::ostringstream report; report<<std::setprecision(17);
  report<<"Display-only redraw; sourceRoot="<<inputRoot<<"\n"
    <<"display=["<<displayLow<<","<<displayHigh<<")\n"
    <<"fullNormalization=["<<counts[0]->GetXaxis()->GetXmin()<<","
    <<counts[0]->GetXaxis()->GetXmax()<<"); excludes flow; no bin-width factor\n"
    <<"countWindow=["<<windowLow<<","<<windowHigh<<")\n";
  for (Int_t j=0;j<3;++j) {
    for (Int_t i=0;i<=counts[j]->GetNbinsX()+1;++i) {
      if (counts[j]->GetBinContent(i)!=sources[j]->GetBinContent(i) ||
          counts[j]->GetBinError(i)!=sources[j]->GetBinError(i))
        return LroFail("source counts/errors changed");
      if (draw[0][j]->GetBinContent(i)!=counts[j]->GetBinContent(i) ||
          TMath::Abs(draw[1][j]->GetBinContent(i)-counts[j]->GetBinContent(i)/norm[j])>1.e-12)
        return LroFail("display bin values changed unexpectedly");
    }
    if (TMath::Abs(draw[1][j]->Integral(1,draw[1][j]->GetNbinsX())-1.)>1.e-12)
      return LroFail("full-domain normalized integral is not one");
    report<<keys[j]<<": Nwin="<<nwin[j]<<"; fullIntegral="<<norm[j]
      <<"; displayedRawIntegral="<<counts[j]->Integral(first,last)<<"; bins="<<counts[j]->GetNbinsX()<<"\n";
  }
  report<<"PASS: all source and absolute count bins/errors unchanged; normalized bins equal original/fullIntegral; unit sum over original full range. No cut or statistical metric recomputation.\n";
  TFile output(TString(outputStem)+".root","CREATE");
  if (output.IsZombie()) return LroFail("cannot create output ROOT");
  Int_t bytes=0;
  for (Int_t j=0;j<3;++j) bytes+=counts[j]->Write();
  for (Int_t p=0;p<2;++p) for (Int_t j=0;j<3;++j) bytes+=draw[p][j]->Write();
  TObjString metadata(report.str().c_str()); bytes+=metadata.Write("DisplayAudit")+canvas->Write();
  output.Flush(); Bool_t failed=output.TestBit(TFile::kWriteError); output.Close();
  if (bytes<=0 || failed || output.TestBit(TFile::kWriteError)) return LroFail("ROOT write failure");
  std::ofstream audit((TString(outputStem)+"_display.txt").Data());
  audit<<report.str(); audit.close();
  if (!audit) return LroFail("cannot write audit");
  canvas->Print(TString(outputStem)+".pdf");
  canvas->Print(TString(outputStem)+".png");
  for (Int_t j=0;j<4;++j) {
    FileStat_t st;
    if (gSystem->GetPathInfo(TString(outputStem)+extensions[j],st)!=0 || st.fSize<=0)
      return LroFail("missing output product");
  }
  std::cout<<report.str();
  return 0;
}
