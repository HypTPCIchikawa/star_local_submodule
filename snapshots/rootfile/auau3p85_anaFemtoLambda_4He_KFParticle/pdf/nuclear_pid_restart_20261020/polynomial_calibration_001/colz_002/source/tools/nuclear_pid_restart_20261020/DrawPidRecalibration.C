// Raw COLZ with original and candidate PID curves. No refit or rebin.
void DrawPidRecalibration(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include <TFile.h>
#include <TSystem.h>
#include <TH2.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>
#include <TAxis.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <iostream>
namespace PidRecalibrationView {
struct Point {double r,mean,sigma;};typedef std::map<std::string,std::vector<Point> > Curves;
std::string Trim(const std::string& s){size_t a=s.find_first_not_of(" \t\r\n");return a==std::string::npos?"":s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);}
std::map<std::string,std::string> Settings(const char* path){std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read settings");std::map<std::string,std::string> out;std::string line;
  while(std::getline(in,line)){line=Trim(line);if(line.empty()||line[0]=='#')continue;size_t t=line.find('\t');if(t==std::string::npos)throw std::runtime_error("Settings must be key TAB value");std::string k=Trim(line.substr(0,t)),v=Trim(line.substr(t+1));if(k=="key"&&v=="value")continue;if(k.empty()||v.empty()||out.count(k))throw std::runtime_error("Invalid/duplicate setting");out[k]=v;}return out;}
std::string Get(const std::map<std::string,std::string>& c,const std::string& k){std::map<std::string,std::string>::const_iterator i=c.find(k);if(i==c.end())throw std::runtime_error("Missing setting "+k);return i->second;}
double Number(const std::map<std::string,std::string>& c,const std::string& k){std::istringstream in(Get(c,k));double x;std::string extra;if(!(in>>x)||(in>>extra)||!TMath::Finite(x))throw std::runtime_error("Invalid plot number");return x;}
bool Less(const Point& a,const Point& b){return a.r<b.r;}
Curves ReadCurves(const char* path){std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read curves");Curves out;std::string line;
  while(std::getline(in,line)){line=Trim(line);if(line.empty()||line[0]=='#')continue;std::istringstream row(line);std::string s,extra;Point p;row>>s;if(s=="species")continue;
    if(!(row>>p.r>>p.mean>>p.sigma)||(row>>extra)||!TMath::Finite(p.r)||!TMath::Finite(p.mean)||!TMath::Finite(p.sigma)||p.r<=0||p.sigma<=0)throw std::runtime_error("Invalid curve row");
    if(s!="d"&&s!="t"&&s!="he3"&&s!="he4")throw std::runtime_error("Unknown species");out[s].push_back(p);}
  const char* s[]={"d","t","he3","he4"};for(int j=0;j<4;++j){std::vector<Point>& p=out[s[j]];if(p.size()<2)throw std::runtime_error("Missing curve points");std::sort(p.begin(),p.end(),Less);for(unsigned i=1;i<p.size();++i)if(p[i].r<=p[i-1].r)throw std::runtime_error("Duplicate R");}return out;}
TGraph* Graph(const std::vector<Point>& p,double n,int color,int style,int width){TGraph* g=new TGraph(p.size());for(unsigned i=0;i<p.size();++i)g->SetPoint(i,p[i].r,p[i].mean+n*p[i].sigma);g->SetLineColor(color);g->SetLineStyle(style);g->SetLineWidth(width);return g;}
void Text(double x,double y,const char* s,double size){TLatex t;t.SetNDC();t.SetTextSize(size);t.DrawLatex(x,y,s);}
void Run(const char* cache,const char* oldPath,const char* newPath,const char* settings,const char* pdf){
  std::map<std::string,std::string> c=Settings(settings);Curves oldC=ReadCurves(oldPath),newC=ReadCurves(newPath);double n=Number(c,"selectionNSigma");if(n<=0)throw std::runtime_error("Invalid gate");
  double xl[]={Number(c,"rLowLo"),Number(c,"rHighLo")},xh[]={Number(c,"rLowHi"),Number(c,"rHighHi")},yl[]={Number(c,"dedxLowLo"),Number(c,"dedxHighLo")},yh[]={Number(c,"dedxLowHi"),Number(c,"dedxHighHi")};
  TFile file(cache,"READ");if(file.IsZombie())throw std::runtime_error("Cannot open cache");const char* modes[]={"ALL: no m^{2} requirement","TAG_D: d / ^{4}He-compatible m^{2}","TAG_T: t-compatible m^{2}","TAG_HE3: ^{3}He-compatible m^{2}"};const char* keys[]={"allKey","tagDKey","tagTKey","tagHe3Key"};const char* s[]={"d","t","he3","he4"};const char* labels[]={"d","t","^{3}He","^{4}He"};int colors[]={kBlue+1,kOrange+7,kGreen+2,kMagenta+1};
  TCanvas canvas("pidRecalibrationCanvas","PID polynomial calibration",1400,850);gStyle->SetOptStat(0);canvas.Print((std::string(pdf)+"[").c_str());unsigned pages=0;
  for(unsigned m=0;m<4;++m){TH2* source=dynamic_cast<TH2*>(file.Get(Get(c,keys[m]).c_str()));if(!source)throw std::runtime_error("Missing raw COLZ");
    for(unsigned a=0;a<2;++a){if(!(xl[a]<xh[a]&&yl[a]<yh[a]))throw std::runtime_error("Invalid view range");canvas.Clear();canvas.SetLeftMargin(.085);canvas.SetRightMargin(.15);canvas.SetTopMargin(.10);canvas.SetBottomMargin(.12);canvas.SetLogz(1);
      TH2* h=dynamic_cast<TH2*>(source->Clone(Form("pidRecalibration_%u_%u",m,a)));h->SetDirectory(0);h->SetTitle("");h->GetXaxis()->SetTitle("p/q [GeV/c]");h->GetYaxis()->SetTitle("dE/dx [keV/cm]");h->GetXaxis()->SetRangeUser(xl[a],xh[a]);h->GetYaxis()->SetRangeUser(yl[a],yh[a]);h->Draw("COLZ");
      TLegend legend(.57,.63,.83,.88);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.026);std::vector<TGraph*> graphs;
      for(unsigned j=0;j<4;++j)for(int b=-1;b<=1;++b){TGraph* g=Graph(oldC[s[j]],b*n,kGray+2,2,b==0?2:1);g->Draw("L SAME");graphs.push_back(g);if(j==0&&b==0)legend.AddEntry(g,"Git original (center / gates)","l");}
      for(unsigned j=0;j<4;++j)for(int b=-1;b<=1;++b){TGraph* g=Graph(newC[s[j]],b*n,colors[j],b==0?1:7,b==0?3:1);g->Draw("L SAME");graphs.push_back(g);if(b==0)legend.AddEntry(g,Form("%s candidate",labels[j]),"l");}
      legend.Draw();Text(.085,.947,modes[m],.033);Text(.085,.915,Form("%s R; centers and #pm%.3g#sigma gates (candidate, not production)",a==0?"Low":"High",n),.026);Text(.085,.035,"Raw counts unchanged. Supplied curve support only; no extrapolation.",.022);canvas.RedrawAxis();canvas.Print(pdf);++pages;canvas.Clear();for(unsigned i=0;i<graphs.size();++i)delete graphs[i];delete h;
    }}canvas.Print((std::string(pdf)+"]").c_str());std::cout<<"PID_RECALIBRATION_DRAW_COMPLETE pages="<<pages<<" rawCountsUnchanged=1 noExtrapolation=1"<<std::endl;
}
}
#endif
void DrawPidRecalibration(const char* cache,const char* oldPath,const char* newPath,const char* settings,const char* pdf){
#ifndef __CINT__
  try{PidRecalibrationView::Run(cache,oldPath,newPath,settings,pdf);}catch(const std::exception& e){std::cerr<<"PID_RECALIBRATION_DRAW_FATAL: "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
