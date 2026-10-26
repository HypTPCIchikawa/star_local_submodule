// Presentation only: saved native data/model counts are never refitted or normalized.
void PlotCutoffSpectra(const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace CutoffSpectra {
void Require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
typedef std::map<std::string,std::string> Settings;
double Number(const Settings& c,const char* key){Settings::const_iterator it=c.find(key);Require(it!=c.end(),std::string("Missing explicit plot setting: ")+key);std::istringstream in(it->second);double v;Require(bool(in>>v),std::string("Invalid plot number: ")+key);in>>std::ws;Require(in.eof()&&std::isfinite(v),"Invalid finite plot setting");return v;}
Settings ReadSettings(const char* path){std::ifstream in(path);Require(bool(in),"Missing plot settings TSV");Settings out;std::string line;while(std::getline(in,line)){size_t tab=line.find('\t');Require(tab!=std::string::npos&&line.find('\t',tab+1)==std::string::npos,"Malformed settings row");Require(out.insert(std::make_pair(line.substr(0,tab),line.substr(tab+1))).second,"Duplicate plot setting");}Require(in.eof(),"Unreadable plot settings");return out;}
struct Model {std::string label,path;TFile* file;};
struct Selection {int index;double rigidity;};
std::vector<Model> ReadModels(const char* path){std::ifstream in(path);Require(bool(in),"Missing model paths TSV");std::string line;Require(bool(std::getline(in,line))&&line=="label\tpath","Wrong model TSV header");std::vector<Model> out;
  while(std::getline(in,line)){size_t p=line.find('\t');Require(p!=std::string::npos&&line.find('\t',p+1)==std::string::npos,"Malformed model row");Model m;m.label=line.substr(0,p);m.path=line.substr(p+1);Require(!m.label.empty()&&!m.path.empty(),"Empty model label/path");m.file=TFile::Open(m.path.c_str(),"READ");Require(m.file&&!m.file->IsZombie()&&!m.file->TestBit(TFile::kRecovered),"Unusable saved model ROOT");out.push_back(m);}
  Require(in.eof()&&!out.empty(),"At least one explicit comparison model file required");return out;
}
std::vector<Selection> ReadSelections(const char* path){std::ifstream in(path);Require(bool(in),"Missing explicit spectrum selections");std::string line;Require(bool(std::getline(in,line))&&line=="index\trigidity","Wrong spectrum-selection header");std::vector<Selection> out;while(std::getline(in,line)){std::istringstream row(line);Selection s;std::string extra;Require(bool(row>>s.index>>s.rigidity)&&!(row>>extra)&&s.index>=0&&std::isfinite(s.rigidity)&&s.rigidity>0,"Malformed spectrum selection");out.push_back(s);}Require(in.eof()&&!out.empty(),"Empty spectrum selection");return out;}
TH1* Load(TFile& file,const char* prefix,int index){TH1* h=dynamic_cast<TH1*>(file.Get(Form("%s_%03d",prefix,index)));Require(h!=0,std::string("Missing saved histogram: ")+prefix);return h;}
void VerifyRaw(TH1& a,TH1& b){Require(a.GetNbinsX()==b.GetNbinsX(),"Compared model native axes differ");for(int bin=0;bin<=a.GetNbinsX()+1;++bin){Require(a.GetBinContent(bin)==b.GetBinContent(bin)&&a.GetBinError(bin)==b.GetBinError(bin),"Compared models have different raw count/error/flow");if(bin>=1&&bin<=a.GetNbinsX())Require(a.GetBinLowEdge(bin)==b.GetBinLowEdge(bin)&&a.GetBinWidth(bin)==b.GetBinWidth(bin),"Compared native bin edges differ");}}
TH1* Display(TH1& native,int rebin,const std::string& name){Require(native.GetNbinsX()%rebin==0,"Display grouping does not divide native bins");TH1* h=static_cast<TH1*>(native.Clone(name.c_str()));h->SetDirectory(0);double before=h->Integral(0,h->GetNbinsX()+1);h->Rebin(rebin);Require(std::fabs(h->Integral(0,h->GetNbinsX()+1)-before)<1e-10*std::max(1.,std::fabs(before)),"Display rebin changed saved count sum");return h;}
void DrawPanel(TFile& file,const std::string& label,const Selection& selection,int rebin,double low,double high,bool log,double logMin,int serial){
  const char* keys[]={"raw_native","prediction_total","prediction_smooth","prediction_d","prediction_t","prediction_he3","prediction_he4"};const char* labels[]={"ALL native counts","Total fit","Smooth BG","d","t","^{3}He","^{4}He"};int colors[]={kBlack,kRed+1,kGray+2,kBlue+1,kOrange+7,kGreen+2,kMagenta+1};
  std::vector<TH1*> hist;for(int i=0;i<7;++i){std::ostringstream name;name<<"cutoff_view_"<<serial<<'_'<<i;hist.push_back(Display(*Load(file,keys[i],selection.index),rebin,name.str()));}
  gPad->SetTopMargin(.15);gPad->SetBottomMargin(.13);gPad->SetLeftMargin(.13);gPad->SetRightMargin(.04);gPad->SetLogy(log);
  double maximum=1;int first=hist[0]->FindBin(low+1e-9),last=hist[0]->FindBin(high-1e-9);for(int i=0;i<7;++i)for(int b=first;b<=last;++b)maximum=std::max(maximum,hist[i]->GetBinContent(b));
  TH1* raw=hist[0];raw->GetXaxis()->SetRangeUser(low,high);raw->SetTitle("");raw->GetXaxis()->SetTitle("dE/dx [keV/cm]");raw->GetYaxis()->SetTitle(Form("Counts / %.3g keV/cm",raw->GetBinWidth(1)));raw->GetXaxis()->SetTitleSize(.047);raw->GetYaxis()->SetTitleSize(.045);raw->GetYaxis()->SetTitleOffset(1.35);raw->SetMarkerStyle(20);raw->SetMarkerSize(.3);raw->SetLineColor(kBlack);raw->SetMinimum(log?logMin:0);raw->SetMaximum(maximum*(log?10:1.3));raw->Draw("E1");
  for(int i=1;i<7;++i){hist[i]->SetLineColor(colors[i]);hist[i]->SetLineWidth(i==1?2:1);hist[i]->SetLineStyle(i==2?2:1);hist[i]->Draw("HIST SAME");}raw->Draw("E1 SAME");
  TLegend* legend=new TLegend(.56,.49,.95,.84);legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(.029);for(int i=0;i<7;++i)legend->AddEntry(hist[i],labels[i],i==0?"pe":"l");legend->Draw();
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.039);text.DrawLatex(.13,.94,Form("R = %.3f GeV/c; %s",selection.rigidity,label.c_str()));text.SetTextSize(.026);text.DrawLatex(.13,.89,log?"Log scale; fixed projected TAG response":"Linear scale; no display normalization");
  gPad->Modified();gPad->Update();
  // Canvas owns drawing references until it is cleared by the caller.
}
}
#endif
void PlotCutoffSpectra(const char* modelsPath,const char* selectionsPath,const char* settingsPath,const char* outputPdf){
#ifndef __CINT__
  using namespace CutoffSpectra;
  try{Require(modelsPath&&selectionsPath&&settingsPath&&outputPdf,"Explicit spectrum plot inputs required");Require(gSystem->AccessPathName(outputPdf),"Existing spectrum PDF will not be overwritten");Settings cfg=ReadSettings(settingsPath);double low=Number(cfg,"dedxMin"),high=Number(cfg,"dedxMax"),logMin=Number(cfg,"logMinimum"),group=Number(cfg,"displayDedxRebinFactor");Require(high>low&&logMin>0&&group>=1&&group==std::floor(group),"Invalid explicit spectrum view settings");int rebin=int(group);std::vector<Model> models=ReadModels(modelsPath);std::vector<Selection> selections=ReadSelections(selectionsPath);
    gStyle->SetOptStat(0);gStyle->SetPaperSize(26,20);TCanvas canvas("nuclear_cutoff_spectra","Nuclear cutoff spectra",1300,1000);canvas.Print((std::string(outputPdf)+"[").c_str());int serial=0;
    for(unsigned i=0;i<selections.size();++i){const Selection& s=selections[i];double fitLo=high,fitHi=low;for(unsigned m=0;m<models.size();++m){VerifyRaw(*Load(*models[0].file,"raw_native",s.index),*Load(*models[m].file,"raw_native",s.index));TH1* mask=Load(*models[m].file,"prediction_fit_mask",s.index);for(int b=1;b<=mask->GetNbinsX();++b)if(mask->GetBinContent(b)>0){fitLo=std::min(fitLo,mask->GetBinLowEdge(b));fitHi=std::max(fitHi,mask->GetXaxis()->GetBinUpEdge(b));}}
      for(int view=0;view<2;++view)for(unsigned pair=0;pair<models.size();pair+=2){double a=view==0?low:std::max(low,fitLo),b=view==0?high:std::min(high,fitHi);if(!(b>a)){a=low;b=high;}canvas.Clear();canvas.Divide(2,2);for(unsigned m=0;m<2&&pair+m<models.size();++m)for(int log=0;log<2;++log){canvas.cd(2*m+log+1);DrawPanel(*models[pair+m].file,models[pair+m].label,s,rebin,a,b,log!=0,logMin,serial++);}canvas.Print(outputPdf);}
    }
    canvas.Print((std::string(outputPdf)+"]").c_str());for(unsigned m=0;m<models.size();++m){models[m].file->Close();delete models[m].file;}std::cout<<"NUCLEAR_CUTOFF_SPECTRA_COMPLETE pages="<<2*((models.size()+1)/2)*selections.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"PlotCutoffSpectra: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
