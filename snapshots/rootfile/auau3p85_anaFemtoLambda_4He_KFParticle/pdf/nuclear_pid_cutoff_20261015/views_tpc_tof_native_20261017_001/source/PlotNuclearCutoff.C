// Presentation only: all thresholds and input tables are explicit arguments.
void PlotNuclearCutoff(const char*,const char*,const char*,const char*,double,double,double,const char* = 0);
#ifndef __CINT__
#include "TCanvas.h"
#include "TFile.h"
#include "TError.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace NuclearCutoffPlot {
typedef std::map<std::string,std::string> Row;
std::vector<std::string> Fields(const std::string& text) {
  std::vector<std::string> values;std::istringstream in(text);std::string v;
  while(std::getline(in,v,'\t'))values.push_back(v);return values;
}
std::vector<Row> Read(const char* path) {
  std::ifstream in(path);if(!in)throw std::runtime_error("Missing explicit plot table");
  std::string line;std::getline(in,line);std::vector<std::string> header=Fields(line);
  std::vector<Row> rows;
  while(std::getline(in,line)) {
    if(line.empty())continue;std::vector<std::string> fields=Fields(line);
    if(fields.size()!=header.size())throw std::runtime_error("Malformed plot table");
    Row row;for(unsigned i=0;i<header.size();++i)row[header[i]]=fields[i];rows.push_back(row);
  }
  return rows;
}
std::string Text(const Row& row,const char* key) {
  Row::const_iterator i=row.find(key);if(i==row.end())throw std::runtime_error(std::string("Missing plot column: ")+key);return i->second;
}
double Number(const Row& row,const char* key) {
  const std::string text=Text(row,key);
  if(text=="True"||text=="true")return 1;
  if(text=="False"||text=="false")return 0;
  if(text=="nan"||text=="None"||text=="null"||text=="")return std::numeric_limits<double>::quiet_NaN();
  std::istringstream in(text);double value;
  if(!(in>>value))throw std::runtime_error(std::string("Bad plot number: ")+key);return value;
}
void Label(const char* title,const char* subtitle) {
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.034);
  text.DrawLatex(.13,.96,title);text.SetTextSize(.024);text.DrawLatex(.13,.925,subtitle);
}
bool InRanges(const std::vector<Row>& ranges,const char* species,double r) {
  for(unsigned i=0;i<ranges.size();++i)if(Text(ranges[i],"species")==species&&r>=Number(ranges[i],"min")&&r<Number(ranges[i],"max"))return true;
  return false;
}
double DrawRangeVeto(const std::vector<Row>& ranges,const char* species,double rMin,double rMax) {
  double previous=rMin,first=std::numeric_limits<double>::quiet_NaN();
  for(unsigned i=0;i<ranges.size();++i)if(Text(ranges[i],"species")==species) {
    const double lo=Number(ranges[i],"min"),hi=Number(ranges[i],"max");
    if(!(hi>lo)||lo<previous-1.e-8)throw std::runtime_error("Unsorted/overlapping adopted plot ranges");
    if(!std::isfinite(first))first=lo;
    if(lo>previous){TBox* box=new TBox(previous,0,lo,1.06);box->SetFillColor(kGray+1);box->SetFillStyle(3004);box->SetLineColor(kGray+1);box->Draw();}
    TLine* accepted=new TLine(lo,.025,hi,.025);accepted->SetLineColor(kGreen+2);accepted->SetLineWidth(4);accepted->Draw();
    previous=hi;
  }
  if(previous<rMax){TBox* box=new TBox(previous,0,rMax,1.06);box->SetFillColor(kGray+1);box->SetFillStyle(3004);box->SetLineColor(kGray+1);box->Draw();}
  return first;
}
}

void PlotNuclearCutoff(const char* cachePath,const char* tablePath,const char* gatePath,
                      const char* outputPdf,double rMin,double rMax,double purityThreshold,const char* adoptedRangePath) {
  using namespace NuclearCutoffPlot;
  try {
    if(!cachePath||!tablePath||!gatePath||!outputPdf||!(rMax>rMin)||!(purityThreshold>0&&purityThreshold<1))
      throw std::runtime_error("Explicit plot inputs/range/threshold required");
    if(!gSystem->AccessPathName(outputPdf))throw std::runtime_error("Refusing to overwrite prior PDF");
    const std::vector<Row> rows=Read(tablePath),gates=Read(gatePath);
    const bool useRanges=adoptedRangePath&&*adoptedRangePath;
    const std::vector<Row> ranges=useRanges?Read(adoptedRangePath):std::vector<Row>();
    TFile cache(cachePath,"READ");if(cache.IsZombie())throw std::runtime_error("Cannot open cache");
    TH2* input=dynamic_cast<TH2*>(cache.Get("ALL"));if(!input)throw std::runtime_error("Missing ALL histogram");
    const char* names[]={"d","t","he3","he4"};
    const char* labels[]={"d","t","^{3}He","^{4}He"};
    gStyle->SetOptStat(0);gStyle->SetPaperSize(26,20);
    TCanvas canvas("nuclear_cutoff_report","Nuclear PID cutoff",1300,1000);
    canvas.Print((std::string(outputPdf)+"[").c_str());
    canvas.Divide(2,2);
    for(int s=0;s<4;++s) {
      canvas.cd(s+1);gPad->SetGridy();gPad->SetTopMargin(.13);
      TH1* frame=gPad->DrawFrame(rMin,0,rMax,1.06);
      frame->SetTitle(";R = |p/q| [GeV/c];TPC-only selection purity");
      TGraphAsymmErrors* band=new TGraphAsymmErrors;TGraph* nominal=new TGraph;TGraph* unsupported=new TGraph;
      int n=0,u=0;double cut=std::numeric_limits<double>::quiet_NaN();
      const double adoptedMin=useRanges?DrawRangeVeto(ranges,names[s],rMin,rMax):cut;
      for(unsigned i=0;i<rows.size();++i)if(Text(rows[i],"species")==names[s]) {
        const double lo=Number(rows[i],"rLo"),hi=Number(rows[i],"rHi"),x=.5*(lo+hi);
        const double y=Number(rows[i],"purityNominal"),low=Number(rows[i],"purityLow"),high=Number(rows[i],"purityHigh");
        if(std::isfinite(y)) {
          nominal->SetPoint(n,x,y);
          if(std::isfinite(low)&&std::isfinite(high)) {
            band->SetPoint(n,x,y);band->SetPointError(n,.5*(hi-lo),.5*(hi-lo),std::max(0.,y-low),std::max(0.,high-y));
          } else {band->SetPoint(n,x,y);band->SetPointError(n,0,0,0,0);}
          ++n;
          if(Number(rows[i],"quality")<.5)unsupported->SetPoint(u++,x,y);
        }
        const double value=Number(rows[i],"cutoff");if(std::isfinite(value))cut=value;
      }
      if(useRanges)cut=adoptedMin;
      band->SetFillColor(kAzure-9);band->SetLineColor(kAzure+1);if(band->GetN())band->Draw("2 SAME");
      nominal->SetMarkerStyle(20);nominal->SetMarkerSize(.45);nominal->SetMarkerColor(kBlue+1);nominal->SetLineColor(kBlue+1);if(nominal->GetN())nominal->Draw("P SAME");
      unsupported->SetMarkerStyle(24);unsupported->SetMarkerSize(.7);unsupported->SetMarkerColor(kRed+1);if(unsupported->GetN())unsupported->Draw("P SAME");
      TLine* line=new TLine(rMin,purityThreshold,rMax,purityThreshold);line->SetLineStyle(2);line->Draw();
      if(std::isfinite(cut)){TLine* lower=new TLine(cut,0,cut,1);lower->SetLineColor(kGreen+2);lower->SetLineWidth(2);lower->Draw();}
      Label(labels[s],useRanges?"Band: model sensitivity; hatched R excluded (may be unassessed)":"Band: model/transfer sensitivity; red: quality not established");
      if(s==0) {
        TLegend* legend=new TLegend(.17,.16,.61,.32);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.026);
        legend->AddEntry(nominal,"Gaussian + linear reference","p");legend->AddEntry(band,"Evaluated model envelope","f");legend->AddEntry(unsupported,"Unassessed / quality failed","p");legend->Draw();
      }
    }
    canvas.Print(outputPdf);canvas.Clear();canvas.Divide(2,2);
    for(int s=0;s<4;++s) {
      canvas.cd(s+1);gPad->SetLogz();gPad->SetRightMargin(.16);gPad->SetTopMargin(.13);
      TH2* background=dynamic_cast<TH2*>(input->Clone((std::string("cutoff_colz_")+names[s]).c_str()));background->SetDirectory(0);
      background->GetXaxis()->SetRangeUser(rMin,rMax);background->SetTitle(";R = |p/q| [GeV/c];dE/dx [keV/cm]");background->Draw("COLZ");
      TGraph* center=new TGraph;TGraph* below=new TGraph;TGraph* above=new TGraph;int n=0;
      for(unsigned i=0;i<gates.size();++i)if(Text(gates[i],"species")==names[s]) {
        const double r=Number(gates[i],"r"),mean=Number(gates[i],"mean"),sigma=Number(gates[i],"sigma"),nsigma=Number(gates[i],"nSigma");
        if(r<rMin||r>rMax||!std::isfinite(mean)||!(sigma>0))continue;
        center->SetPoint(n,r,mean);below->SetPoint(n,r,mean-nsigma*sigma);above->SetPoint(n,r,mean+nsigma*sigma);++n;
      }
      center->SetLineColor(kRed+1);center->SetLineWidth(2);if(center->GetN())center->Draw("L SAME");
      below->SetLineColor(kMagenta+1);above->SetLineColor(kMagenta+1);below->SetLineStyle(2);above->SetLineStyle(2);if(below->GetN())below->Draw("L SAME");if(above->GetN())above->Draw("L SAME");
      Label(labels[s],"ALL; current center and +/-nSigma guides (best-species also applied)");
    }
    canvas.Print(outputPdf);canvas.Clear();canvas.Divide(2,2);
    for(int s=0;s<4;++s) {
      canvas.cd(s+1);gPad->SetLogy();gPad->SetTopMargin(.13);
      double largest=1;for(unsigned i=0;i<rows.size();++i)if(Text(rows[i],"species")==names[s])largest=std::max(largest,Number(rows[i],"rawSelected"));
      TH1* frame=gPad->DrawFrame(rMin,.5,rMax,largest*1.5);frame->SetTitle(";R = |p/q| [GeV/c];Selected tracks / reported R bin");
      TGraph* count=new TGraph;TGraph* accepted=new TGraph;int n=0,a=0;
      for(unsigned i=0;i<rows.size();++i)if(Text(rows[i],"species")==names[s]) {
        const double x=.5*(Number(rows[i],"rLo")+Number(rows[i],"rHi")),y=Number(rows[i],"rawSelected");
        if(y>0){count->SetPoint(n++,x,y);if(useRanges?InRanges(ranges,names[s],x):(Number(rows[i],"quality")>.5&&Number(rows[i],"purityLow")>=purityThreshold))accepted->SetPoint(a++,x,y);}
      }
      count->SetLineColor(kGray+2);count->SetLineWidth(2);if(count->GetN())count->Draw("L SAME");
      accepted->SetMarkerColor(kGreen+2);accepted->SetMarkerStyle(20);accepted->SetMarkerSize(.6);if(accepted->GetN())accepted->Draw("P SAME");
      Label(labels[s],useRanges?"Gray: current selected counts; green: adopted R ranges":"Gray: current selected counts; green: evaluated low envelope passes");
    }
    canvas.Print(outputPdf);canvas.Print((std::string(outputPdf)+"]").c_str());
  } catch(const std::exception& error) {
    Error("PlotNuclearCutoff","%s",error.what());gSystem->Exit(1);
  }
}
#endif
