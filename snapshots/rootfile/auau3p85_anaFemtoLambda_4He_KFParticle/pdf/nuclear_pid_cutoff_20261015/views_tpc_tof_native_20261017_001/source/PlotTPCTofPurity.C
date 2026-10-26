// Presentation only. Input purities are never fitted or averaged here.
void PlotTPCTofPurity(const char*,const char*,const char*,const char*,double,double,double,double);
#ifndef __CINT__
#include "PlotNuclearCutoff.C"

namespace NuclearTPCTofPlot {
using NuclearCutoffPlot::Row;
using NuclearCutoffPlot::Text;
using NuclearCutoffPlot::Number;

bool Finite(double x) { return std::isfinite(x); }
void Caption(double x,double y,const std::string& label,double size=.025) {
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(size);
  text.DrawLatex(x,y,label.c_str());
}
const Row& Meta(const std::vector<Row>& rows,const char* species) {
  const Row* found=0;
  for(unsigned i=0;i<rows.size();++i)if(Text(rows[i],"species")==species) {
    if(found)throw std::runtime_error("Duplicate species metadata");found=&rows[i];
  }
  if(!found)throw std::runtime_error("Missing species metadata");return *found;
}
void DrawCurve(const std::vector<Row>& rows,const char* species,const char* regime,
               int color,int marker,TLegend* legend,const char* label) {
  TGraphAsymmErrors* band=new TGraphAsymmErrors;
  TGraph* valid=new TGraph;TGraph* unassessed=new TGraph;
  int n=0,v=0,u=0;
  for(unsigned i=0;i<rows.size();++i) {
    const Row& row=rows[i];
    if(Text(row,"species")!=species||Text(row,"regime")!=regime)continue;
    const double lo=Number(row,"rLo"),hi=Number(row,"rHi");
    const double y=Number(row,"purityNominal"),pLo=Number(row,"purityLow"),pHi=Number(row,"purityHigh");
    if(!(hi>lo))throw std::runtime_error("Invalid purity interval");
    if(!Finite(y))continue;
    if(y<0||y>1)throw std::runtime_error("Purity outside physical range");
    const double center=.5*(lo+hi);
    if(Finite(pLo)&&Finite(pHi)) {
      if(pLo<0||pHi>1||pHi<pLo)throw std::runtime_error("Invalid model envelope");
      // Cover both the envelope and the displayed reference, without moving data.
      band->SetPoint(n,center,y);
      band->SetPointError(n,.5*(hi-lo),.5*(hi-lo),std::max(0.,y-pLo),std::max(0.,pHi-y));++n;
    }
    if(Number(row,"quality")>.5)valid->SetPoint(v++,center,y);
    else unassessed->SetPoint(u++,center,y);
  }
  band->SetFillColor(color-9);band->SetFillStyle(3001);band->SetLineColor(color);
  if(band->GetN())band->Draw("2 SAME");
  valid->SetMarkerColor(color);valid->SetLineColor(color);valid->SetMarkerStyle(marker);valid->SetMarkerSize(.58);
  if(valid->GetN())valid->Draw("P SAME");
  unassessed->SetMarkerColor(kRed+1);unassessed->SetMarkerStyle(24);unassessed->SetMarkerSize(.75);
  if(unassessed->GetN())unassessed->Draw("P SAME");
  if(legend) {
    legend->AddEntry(valid,label,"p");
    if(unassessed->GetN())legend->AddEntry(unassessed,"Open red: quality unassessed","p");
  }
}
void PlotPanel(const std::vector<Row>& rows,const std::vector<Row>& metadata,
               const std::vector<Row>& ranges,const char* species,const char* label,
               double rMin,double split,double rMax,double threshold,bool highOnly) {
  const Row& meta=Meta(metadata,species);
  const double displayMax=Number(meta,"displayMax");
  if(!(displayMax>split&&displayMax<=rMax))throw std::runtime_error("Invalid species display range");
  rMax=displayMax;
  const bool hasHigh=Number(meta,"hasHigh")>.5;
  const bool exact=Number(meta,"exactStudyWindow")>.5;
  if(hasHigh&&!exact)throw std::runtime_error("TOF sample differs from requested study window");
  gPad->SetTopMargin(.17);gPad->SetBottomMargin(.13);gPad->SetLeftMargin(.13);gPad->SetRightMargin(.04);gPad->SetGridy();
  TH1* frame=gPad->DrawFrame(highOnly?split:rMin,0,rMax,1.06);
  frame->SetTitle(";R = |p/q| [GeV/c];Selected-track purity");
  frame->GetXaxis()->SetTitleSize(.043);frame->GetYaxis()->SetTitleSize(.043);
  frame->GetYaxis()->SetTitleOffset(1.35);
  if(!highOnly)NuclearCutoffPlot::DrawRangeVeto(ranges,species,rMin,split);
  TLegend* legend=new TLegend(.18,.16,.88,.34);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(.025);
  if(!highOnly)DrawCurve(rows,species,"TPC",kBlue+1,20,legend,"TPC only (20 fit conditions)");
  if(hasHigh)DrawCurve(rows,species,"TOF",exact?kOrange+7:kMagenta+1,exact?21:25,legend,Text(meta,"curveLabel").c_str());
  else {
    TBox* missing=new TBox(split,0,rMax,1.06);missing->SetFillColor(kGray+1);missing->SetFillStyle(3005);missing->Draw();
    Caption(.18,.54,Text(meta,"missingLabel"),.030);
  }
  TLine* passing=new TLine(highOnly?split:rMin,threshold,rMax,threshold);passing->SetLineStyle(2);passing->Draw();
  if(!highOnly) {
    TLine* boundary=new TLine(split,0,split,1.06);boundary->SetLineStyle(7);boundary->SetLineWidth(2);boundary->Draw();
  }
  Caption(.13,.965,label,.037);
  Caption(.13,.923,Text(meta,"conditionLabel"),.024);
  Caption(.13,.882,"Bands: model sensitivity, not a confidence interval",.022);
  legend->Draw();
}
}

void PlotTPCTofPurity(const char* tablePath,const char* metadataPath,const char* rangesPath,
                     const char* outputPdf,double rMin,double split,double rMax,double threshold) {
  using namespace NuclearTPCTofPlot;
  try {
    if(!tablePath||!*tablePath||!metadataPath||!*metadataPath||!rangesPath||!*rangesPath||!outputPdf||!*outputPdf||
       !(rMin<split&&split<rMax)||!(threshold>0&&threshold<1))throw std::runtime_error("Explicit tables, ranges and threshold required");
    if(!gSystem->AccessPathName(outputPdf))throw std::runtime_error("Refusing to overwrite PDF");
    const std::vector<Row> rows=NuclearCutoffPlot::Read(tablePath),meta=NuclearCutoffPlot::Read(metadataPath),ranges=NuclearCutoffPlot::Read(rangesPath);
    const char* species[]={"d","t","he3","he4"};const char* labels[]={"d","t","^{3}He","^{4}He"};
    for(unsigned i=0;i<rows.size();++i) {
      const std::string regime=Text(rows[i],"regime");
      const double lo=Number(rows[i],"rLo"),hi=Number(rows[i],"rHi");
      if(regime=="TPC") {if(lo<rMin-1.e-8||hi>split+1.e-8)throw std::runtime_error("TPC point crosses regime boundary");}
      else if(regime=="TOF") {if(lo<split-1.e-8||hi>rMax+1.e-8)throw std::runtime_error("TOF point outside explicit range");}
      else throw std::runtime_error("Unknown selection regime");
    }
    gStyle->SetOptStat(0);gStyle->SetPaperSize(26,20);
    TCanvas canvas("tpc_tof_purity","TPC and TPC+TOF nuclear purity",1400,1080);
    canvas.Print((std::string(outputPdf)+"[").c_str());
    for(int page=0;page<2;++page) {
      canvas.Clear();canvas.Divide(2,2);
      for(int s=0;s<4;++s) {canvas.cd(s+1);PlotPanel(rows,meta,ranges,species[s],labels[s],rMin,split,rMax,threshold,page==1);}
      canvas.Print(outputPdf);
    }
    canvas.Clear();canvas.cd();
    Caption(.08,.94,"Selection and interpretation",.040);
    double y=.86;
    for(int s=0;s<4;++s) {
      const Row& row=Meta(meta,species[s]);
      Caption(.08,y,std::string(labels[s])+": "+Text(row,"conditionLabel"),.029);y-=.055;
      Caption(.10,y,Text(row,"detailLabel"),.024);y-=.09;
    }
    Caption(.08,.23,"Low R: unchanged 80-run study; inclusive TPC best-species selection.",.026);
    Caption(.08,.18,"High R: purity conditional on the stated TOF window; not TOF efficiency.",.026);
    Caption(.08,.13,"A curve is not evidence of valid fit: open red points remain unassessed.",.026);
    Caption(.08,.08,"Neither result is a centrality/m_{T}/CF-pair-conditioned purity measurement.",.026);
    canvas.Print(outputPdf);canvas.Print((std::string(outputPdf)+"]").c_str());
  } catch(const std::exception& error) {Error("PlotTPCTofPurity","%s",error.what());gSystem->Exit(1);}
}
#endif
