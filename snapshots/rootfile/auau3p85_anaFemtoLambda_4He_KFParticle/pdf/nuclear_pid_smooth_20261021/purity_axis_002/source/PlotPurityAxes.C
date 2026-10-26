void PlotPurityAxes(const char*);
#ifndef __CINT__
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cmath>
namespace SmoothPlot {
struct Point {std::string band,sp;double lo,hi,p,err,low,high;bool good;};
std::vector<std::string> split(const std::string& s){std::vector<std::string> v;std::istringstream f(s);std::string x;while(std::getline(f,x,'\t'))v.push_back(x);return v;}
std::vector<Point> Read(const std::string& p){std::ifstream f(p.c_str());if(!f)throw std::runtime_error("Cannot read "+p);std::string line;std::getline(f,line);std::vector<std::string> h=split(line);std::vector<Point> rows;while(std::getline(f,line)){std::vector<std::string> v=split(line);std::map<std::string,std::string> d;for(unsigned k=0;k<h.size();++k)d[h[k]]=v.at(k);Point r;r.band=d["band"];r.sp=d["species"];r.lo=atof(d["rLo"].c_str());r.hi=atof(d["rHi"].c_str());r.p=atof(d["purity"].c_str());r.err=atof(d["statConditionalBound"].c_str());r.low=atof(d["purityLow"].c_str());r.high=atof(d["purityHigh"].c_str());r.good=d["quality"]=="1";rows.push_back(r);}return rows;}
const char* sp[4]={"d","t","he3","he4"};const char* label[4]={"d","t","^{3}He","^{4}He"};int colors[4]={kBlue+1,kOrange+7,kGreen+2,kMagenta+1};
void Text(double x,double y,const std::string& s,double size=.035){TLatex t;t.SetNDC();t.SetTextSize(size);t.DrawLatex(x,y,s.c_str());}
void Purity(const std::string& dir){std::vector<Point> now=Read(dir+"/candidate_summary.tsv"),old=Read(dir+"/historical_summary.tsv"),git=Read(dir+"/original_summary.tsv");double axisLo[4],axisHi[4];std::ifstream ranges((dir+"/axes.tsv").c_str());if(!ranges)throw std::runtime_error("Missing display axes");std::string row;std::getline(ranges,row);for(int k=0;k<4;++k){std::string species;if(!(ranges>>species>>axisLo[k]>>axisHi[k])||species!=sp[k]||axisHi[k]<=axisLo[k])throw std::runtime_error("Invalid display axes");}TCanvas c("purity","Purity",1300,950);std::string pdf=dir+"/purity_vs_p_over_q.pdf";c.Print((pdf+"[").c_str());
for(int page=0;page<2;++page){c.Clear();c.Divide(2,2);for(int s=0;s<4;++s){c.cd(s+1);gPad->SetGrid();gPad->SetBottomMargin(.14);gPad->SetLeftMargin(.13);double xmin=page?2.:axisLo[s],xmax=page?5.:axisHi[s],ymin=page?.9:-.07;if(page)for(unsigned k=0;k<now.size();++k)if(now[k].good&&now[k].band=="high"&&now[k].sp==sp[s])ymin=std::max(0.,std::min(ymin,std::floor((now[k].low-.02)*10)/10));TH1* frame=gPad->DrawFrame(xmin,ymin,xmax,1.16);frame->SetTitle(Form("%s;R = |p/q| [GeV/c];Purity S/(S+B)",label[s]));frame->GetXaxis()->SetTitleSize(.045);frame->GetYaxis()->SetTitleSize(.045);TLine line(2.,ymin,2.,1.16);line.SetLineStyle(2);line.DrawClone();TLine half(xmin,.5,xmax,.5);half.SetLineStyle(2);half.SetLineColor(kRed+1);if(!page)half.DrawClone();
TLegend legend(.16,.16,.61,.36);legend.SetTextSize(.027);legend.SetFillColor(kWhite);legend.SetBorderSize(0);
for(int collection=0;collection<3;++collection){if(collection==1)continue;std::vector<Point>& data=collection==0?old:collection==1?git:now;for(int band=0;band<2;++band){TGraphAsymmErrors env;TGraphErrors stat;TGraph missing;int n=0,m=0;for(unsigned j=0;j<data.size();++j){const Point& r=data[j];if(r.sp!=sp[s]||r.band!=(band?"high":"low")||r.hi<=xmin)continue;if(!r.good){if(collection==2&&!page)missing.SetPoint(m++,.5*(r.lo+r.hi),-.035);continue;}double x=.5*(r.lo+r.hi)+(collection==1?-.007:0);env.SetPoint(n,x,r.p);env.SetPointError(n,.5*(r.hi-r.lo),.5*(r.hi-r.lo),std::max(0.,r.p-r.low),std::max(0.,r.high-r.p));stat.SetPoint(n,x,r.p);stat.SetPointError(n,0,std::isfinite(r.err)?r.err:0);++n;}
int color=collection==0?kGray+2:collection==1?kGreen+2:(band?kRed+1:kBlue+1);env.SetLineColor(color);env.SetMarkerColor(color);env.SetMarkerStyle(collection==0?24:collection==1?25:20);env.SetMarkerSize(collection==2?.65:.5);env.SetLineWidth(collection==2?2:1);env.DrawClone("PZ SAME");stat.SetLineColor(color);stat.SetLineWidth(1);if(collection==2)stat.DrawClone("PZ SAME");missing.SetMarkerStyle(5);missing.SetMarkerColor(kGray+1);missing.SetMarkerSize(.5);if(m)missing.DrawClone("P SAME");if(!band)legend.AddEntry(env.DrawClone("PZ SAME"),collection==0?"17 Oct reference":collection==1?"Git original (4 nominal)":"Adjusted: ALL below 2", "p");else if(collection==2)legend.AddEntry(env.DrawClone("PZ SAME"),"Adjusted: TOF at/above 2","p");}}
legend.DrawClone();Text(.16,.90,"2021 P24iy, inclusive track QA",.030);Text(.16,.855,"Model range + conditional yield error; x = unassessed",.025);
}c.Print(pdf.c_str());}c.Print((pdf+"]").c_str());c.SaveAs((dir+"/purity_vs_p_over_q_highR.png").c_str());
}
}
#endif
void PlotPurityAxes(const char* directory){
#ifndef __CINT__
try{gStyle->SetOptStat(0);gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);SmoothPlot::Purity(directory);std::cout<<"PURITY_AXES_COMPLETE"<<std::endl;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
