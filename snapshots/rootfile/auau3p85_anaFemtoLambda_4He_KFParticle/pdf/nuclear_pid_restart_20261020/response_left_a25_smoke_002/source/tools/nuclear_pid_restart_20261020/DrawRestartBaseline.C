// ROOT5 ACLiC. Raw ALL/TOF-tag spectra plus the frozen Git PID curves: no fits.
int DrawRestartBaseline(const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct RestartBaselinePoint { double r,mean,sigma; };
void BaselineRequire(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}
std::vector<std::string> BaselineFields(const std::string& line) {
  std::vector<std::string> result;std::istringstream input(line);std::string field;
  while(std::getline(input,field,'\t'))result.push_back(field);
  return result;
}
double BaselineNumber(const std::string& text) {
  std::istringstream input(text);double value=0.;char extra;
  BaselineRequire(bool(input>>value) && !(input>>extra) && TMath::Finite(value),
                  "Invalid finite numeric field: "+text);
  return value;
}
std::map<std::string,double> BaselineSettings(const char* path) {
  std::ifstream input(path);BaselineRequire(input.good(),"Cannot read mandatory settings TSV");
  std::string line;
  BaselineRequire(bool(std::getline(input,line)) && line=="key\tvalue",
                  "Settings TSV header must be key<TAB>value");
  std::map<std::string,double> result;
  while(std::getline(input,line)) {
    const std::vector<std::string> row=BaselineFields(line);
    BaselineRequire(row.size()==2 && !row[0].empty(),"Settings TSV needs two nonempty columns");
    BaselineRequire(result.find(row[0])==result.end(),"Duplicate settings key: "+row[0]);
    result[row[0]]=BaselineNumber(row[1]);
  }
  BaselineRequire(input.eof(),"Settings read did not finish normally");
  const char* required[]={"rLowMin","rLowMax","rHighMin","rHighMax",
    "dedxMin","dedxMax","nsigma","logMinimum","canvasW","canvasH",
    "tagDMass2Min","tagDMass2Max","tagTMass2Min","tagTMass2Max",
    "tagHe3Mass2Min","tagHe3Mass2Max"};
  BaselineRequire(result.size()==sizeof(required)/sizeof(required[0]),
                  "Unexpected or missing settings fields");
  for(unsigned i=0;i<sizeof(required)/sizeof(required[0]);++i)
    BaselineRequire(result.find(required[i])!=result.end(),
                    std::string("Missing mandatory setting: ")+required[i]);
  BaselineRequire(result["rLowMin"]>0. && result["rLowMax"]>result["rLowMin"] &&
                  result["rHighMin"]>0. && result["rHighMax"]>result["rHighMin"],
                  "Invalid explicit rigidity display range");
  BaselineRequire(result["dedxMax"]>result["dedxMin"] && result["dedxMin"]>=0. &&
                  result["nsigma"]>0. && result["logMinimum"]>0.,
                  "Invalid dE/dx, band or log-Z display setting");
  BaselineRequire(result["canvasW"]>=1. && result["canvasH"]>=1. &&
                  std::floor(result["canvasW"])==result["canvasW"] &&
                  std::floor(result["canvasH"])==result["canvasH"],
                  "Canvas dimensions must be explicit positive integers");
  const char* prefixes[]={"tagD","tagT","tagHe3"};
  for(unsigned i=0;i<sizeof(prefixes)/sizeof(prefixes[0]);++i) {
    const std::string prefix(prefixes[i]);
    BaselineRequire(result[prefix+"Mass2Max"]>result[prefix+"Mass2Min"],
                    "Invalid explicit mass-squared window: "+prefix);
  }
  return result;
}
std::map<std::string,std::vector<RestartBaselinePoint> > BaselineCurves(const char* path) {
  std::ifstream input(path);BaselineRequire(input.good(),"Cannot read mandatory Git-curve TSV");
  std::string line;
  BaselineRequire(bool(std::getline(input,line)) && line=="species\tr\tmean\tsigma",
                  "Curve TSV header must be species<TAB>r<TAB>mean<TAB>sigma");
  const char* names[]={"d","t","he3","he4"};std::set<std::string> accepted;
  for(unsigned i=0;i<sizeof(names)/sizeof(names[0]);++i)accepted.insert(names[i]);
  std::map<std::string,std::vector<RestartBaselinePoint> > result;
  while(std::getline(input,line)) {
    const std::vector<std::string> row=BaselineFields(line);
    BaselineRequire(row.size()==4 && accepted.find(row[0])!=accepted.end(),
                    "Unknown species or malformed Git-curve TSV row");
    RestartBaselinePoint point;
    point.r=BaselineNumber(row[1]);point.mean=BaselineNumber(row[2]);point.sigma=BaselineNumber(row[3]);
    BaselineRequire(point.r>0. && point.mean>0. && point.sigma>0.,
                    "Curve inputs must be supported positive Git values; no fallback accepted");
    std::vector<RestartBaselinePoint>& points=result[row[0]];
    BaselineRequire(points.empty() || point.r>points.back().r,
                    "Curve rigidity samples must strictly increase within each species");
    points.push_back(point);
  }
  BaselineRequire(input.eof(),"Git-curve TSV read did not finish normally");
  for(unsigned i=0;i<sizeof(names)/sizeof(names[0]);++i)
    BaselineRequire(result[names[i]].size()>=2,std::string("Missing Git curve: ")+names[i]);
  return result;
}
TString BaselineTagLabel(const char* key,std::map<std::string,double>& settings) {
  if(std::string(key)=="ALL")return "ALL: no TOF m^{2} selection";
  const char* prefix=std::string(key)=="TAG_D" ? "tagD" :
                     (std::string(key)=="TAG_T" ? "tagT" : "tagHe3");
  TString window;
  window.Form("%.6g < m^{2} < %.6g",settings[std::string(prefix)+"Mass2Min"],
               settings[std::string(prefix)+"Mass2Max"]);
  TString result;
  if(std::string(key)=="TAG_D")result="TAG_D: d + ^{4}He shared m^{2} window; NOT pure ^{4}He; ";
  else if(std::string(key)=="TAG_T")result="TAG_T: t m^{2} window; ";
  else result="TAG_HE3: ^{3}He m^{2} window; ";
  return result+window;
}
}

int DrawRestartBaseline(const char* cachePath,const char* curvesTSV,
                        const char* settingsTSV,const char* outPDF) {
  TFile* cache=0;TCanvas* canvas=0;bool pdfOpen=false;
  try {
    BaselineRequire(cachePath && cachePath[0] && curvesTSV && curvesTSV[0] &&
                    settingsTSV && settingsTSV[0] && outPDF && outPDF[0],
                    "Explicit cache, curves, settings and output paths are required");
    BaselineRequire(gSystem->AccessPathName(outPDF),"Output PDF must be new");
    std::map<std::string,double> settings=BaselineSettings(settingsTSV);
    std::map<std::string,std::vector<RestartBaselinePoint> > curves=BaselineCurves(curvesTSV);
    cache=TFile::Open(cachePath,"READ");
    BaselineRequire(cache && !cache->IsZombie(),"Cannot open the immutable raw cache");
    const char* keys[]={"ALL","TAG_D","TAG_T","TAG_HE3"};
    for(unsigned i=0;i<sizeof(keys)/sizeof(keys[0]);++i) {
      TH2* raw=dynamic_cast<TH2*>(cache->Get(keys[i]));
      BaselineRequire(raw && raw->GetDimension()==2,std::string("Missing raw TH2 key: ")+keys[i]);
      BaselineRequire(raw->GetXaxis()->GetXmin()<=settings["rLowMin"] &&
                      raw->GetXaxis()->GetXmax()>=settings["rLowMax"] &&
                      raw->GetXaxis()->GetXmin()<=settings["rHighMin"] &&
                      raw->GetXaxis()->GetXmax()>=settings["rHighMax"] &&
                      raw->GetYaxis()->GetXmin()<=settings["dedxMin"] &&
                      raw->GetYaxis()->GetXmax()>=settings["dedxMax"],
                      "Requested baseline display range extends beyond the raw axes");
    }
    gStyle->SetOptStat(0);
    canvas=new TCanvas("restart_git_pid_baseline","Original Git PID baseline",
                       static_cast<int>(settings["canvasW"]),static_cast<int>(settings["canvasH"]));
    canvas->Print((TString(outPDF)+"[").Data());pdfOpen=true;
    unsigned pages=0;
    const char* species[]={"d","t","he3","he4"};
    const char* labels[]={"d","t","^{3}He","^{4}He"};
    const int colors[]={kBlue+1,kOrange+7,kGreen+2,kMagenta+1};
    for(unsigned key=0;key<sizeof(keys)/sizeof(keys[0]);++key)for(unsigned range=0;range<2;++range) {
      const double lo=settings[range==0 ? "rLowMin" : "rHighMin"];
      const double hi=settings[range==0 ? "rLowMax" : "rHighMax"];
      canvas->Clear();canvas->SetLogz(1);canvas->SetLeftMargin(.10);canvas->SetRightMargin(.15);
      canvas->SetTopMargin(.15);canvas->SetBottomMargin(.11);
      TString name;name.Form("baseline_%s_%u",keys[key],range);
      TH2* histogram=static_cast<TH2*>(cache->Get(keys[key])->Clone(name.Data()));
      histogram->SetDirectory(0);
      TString title;title.Form("Original Git PID baseline (NO FIT): %s;R = |p/q| [GeV/c];dE/dx [keV/cm]",keys[key]);
      histogram->SetTitle(title.Data());
      histogram->GetXaxis()->SetRangeUser(lo,hi);
      histogram->GetYaxis()->SetRangeUser(settings["dedxMin"],settings["dedxMax"]);
      histogram->SetMinimum(settings["logMinimum"]);histogram->Draw("COLZ");
      std::vector<TGraph*> graphs;
      TLegend legend(.65,.62,.83,.83);legend.SetBorderSize(0);legend.SetFillStyle(0);
      for(unsigned sp=0;sp<sizeof(species)/sizeof(species[0]);++sp) {
        TGraph* center=new TGraph();TGraph* lower=new TGraph();TGraph* upper=new TGraph();
        const std::vector<RestartBaselinePoint>& points=curves[species[sp]];
        for(size_t i=0;i<points.size();++i)if(points[i].r>=lo && points[i].r<=hi) {
          center->SetPoint(center->GetN(),points[i].r,points[i].mean);
          lower->SetPoint(lower->GetN(),points[i].r,points[i].mean-settings["nsigma"]*points[i].sigma);
          upper->SetPoint(upper->GetN(),points[i].r,points[i].mean+settings["nsigma"]*points[i].sigma);
        }
        center->SetLineColor(colors[sp]);center->SetLineWidth(3);
        lower->SetLineColor(colors[sp]);upper->SetLineColor(colors[sp]);
        lower->SetLineWidth(1);upper->SetLineWidth(1);lower->SetLineStyle(2);upper->SetLineStyle(2);
        if(center->GetN()>=2){center->Draw("L SAME");lower->Draw("L SAME");upper->Draw("L SAME");}
        legend.AddEntry(center,labels[sp],"l");graphs.push_back(center);graphs.push_back(lower);graphs.push_back(upper);
      }
      legend.Draw();
      TLatex text;text.SetNDC();text.SetTextSize(.025);
      text.DrawLatex(.10,.915,BaselineTagLabel(keys[key],settings).Data());
      TString note;note.Form("Raw counts, log Z; solid: Git #mu(R); dashed: #mu #pm %.6g#sigma; no parameter update",settings["nsigma"]);
      text.DrawLatex(.10,.875,note.Data());
      canvas->Modified();canvas->Update();canvas->Print(outPDF);++pages;
      std::cout<<"RESTART_BASELINE_PAGE page="<<pages<<" cacheKey="<<keys[key]
               <<" rMin="<<lo<<" rMax="<<hi<<" noFit=1"<<std::endl;
      canvas->Clear();
      for(size_t i=0;i<graphs.size();++i)delete graphs[i];
      delete histogram;
    }
    canvas->Print((TString(outPDF)+"]").Data());pdfOpen=false;
    BaselineRequire(pages==8,"Baseline PDF did not receive exactly eight pages");
    delete canvas;canvas=0;cache->Close();delete cache;cache=0;
    std::cout<<"RESTART_BASELINE_PDF_CLOSED pages="<<pages<<" noFit=1 noPidUpdate=1"<<std::endl;
    return 0;
  } catch(const std::exception& error) {
    if(canvas && pdfOpen)canvas->Print((TString(outPDF)+"]").Data());
    delete canvas;if(cache){cache->Close();delete cache;}
    std::cerr<<"RESTART_BASELINE_FAIL "<<error.what()<<std::endl;
    return 1;
  }
}
#endif
