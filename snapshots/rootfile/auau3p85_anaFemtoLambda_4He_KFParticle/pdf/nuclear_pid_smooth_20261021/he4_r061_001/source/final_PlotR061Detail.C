void PlotR061Detail(const char*);
#ifndef __CINT__
#include "TCanvas.h"
#include "TH1.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>
namespace R061Detail {
std::vector<std::string> Split(const std::string& s){std::vector<std::string> v;std::istringstream f(s);std::string x;while(std::getline(f,x,'\t'))v.push_back(x);return v;}
struct Row{std::string response,bg;double r,p,m3,m4;bool good,numerical;};
void Text(double x,double y,const char* s){TLatex t;t.SetNDC();t.SetTextSize(.027);t.DrawLatex(x,y,s);}
void Draw(const std::string& directory){
std::ifstream settings((directory+"/detail_settings.tsv").c_str());std::map<std::string,double> c;std::string key;double value;while(settings>>key>>value)c[key]=value;if(c.size()!=4)throw std::runtime_error("Missing explicit detail settings");
std::ifstream input((directory+"/metrics.tsv").c_str());if(!input)throw std::runtime_error("Missing metrics");std::string line;std::getline(input,line);std::vector<std::string> fields=Split(line);std::vector<Row> rows;while(std::getline(input,line)){std::vector<std::string> v=Split(line);std::map<std::string,std::string> d;for(unsigned i=0;i<fields.size();++i)d[fields[i]]=v.at(i);Row r;r.response=d["response"];r.bg=d["model"];r.r=atof(d["r"].c_str());r.p=atof(d["purity"].c_str());r.m3=atof(d["he3Mean"].c_str());r.m4=atof(d["he4Mean"].c_str());r.good=d["supported"]=="1";r.numerical=d["numericalValid"]=="1";rows.push_back(r);}
std::ifstream anchorFile((directory+"/anchors.tsv").c_str());if(!anchorFile)throw std::runtime_error("Missing anchors");std::getline(anchorFile,line);std::vector<double> ar,a3,a4;double x,y,z;while(anchorFile>>x>>y>>z){ar.push_back(x);a3.push_back(y);a4.push_back(z);}if(ar.empty())throw std::runtime_error("Empty anchors");
const char* responses[2]={"gaussian","lognormal"};const char* backgrounds[2]={"linear","exponential"};int colors[4]={kBlue+1,kRed+1,kGreen+2,kMagenta+1};TCanvas canvas("r061detail","Low R He4",1400,900);std::string pdf=directory+"/he4_low_r_detail.pdf";canvas.Print((pdf+"[").c_str());canvas.SetGrid();canvas.SetBottomMargin(.13);TH1* frame=canvas.DrawFrame(c["minimumR"],-.06,c["detailMaxR"],1.1);frame->SetTitle("^{4}He: ALL, no TOF;R = |p/q| [GeV/c];Purity S/(S+B)");TLine half(c["minimumR"],c["purityThreshold"],c["detailMaxR"],c["purityThreshold"]);half.SetLineStyle(2);half.DrawClone();TLegend legend(.58,.16,.88,.37);legend.SetBorderSize(0);legend.SetTextSize(.026);
for(int response=0;response<2;++response)for(int bg=0;bg<2;++bg){int color=colors[2*response+bg];TGraph good,bad;int ng=0,nb=0;for(unsigned i=0;i<rows.size();++i){Row r=rows[i];if(r.response!=responses[response]||r.bg!=backgrounds[bg]||r.r>c["detailMaxR"]||!r.numerical||!std::isfinite(r.p))continue;if(r.good)good.SetPoint(ng++,r.r,r.p);else bad.SetPoint(nb++,r.r,r.p);}good.SetMarkerStyle(20+response);good.SetMarkerColor(color);good.SetLineColor(color);good.SetMarkerSize(.9);TObject* obj=good.DrawClone("LP SAME");bad.SetMarkerStyle(24+response);bad.SetMarkerColor(color);bad.SetMarkerSize(1.);bad.DrawClone("P SAME");legend.AddEntry(obj,Form("%s + %s",responses[response],backgrounds[bg]),"lp");}legend.DrawClone();Text(.14,.91,"Strict R > 0.61; native 0.01 GeV/c bins; filled: historical quality, open: unassessed");Text(.14,.87,"Model diagnostics: backgrounds retained; mean/width sensitivity envelope not included");canvas.Print(pdf.c_str());
canvas.Clear();canvas.Divide(1,2);for(int species=0;species<2;++species){canvas.cd(species+1);gPad->SetGrid();gPad->SetBottomMargin(.14);TH1* f=gPad->DrawFrame(c["minimumR"],species?48:34,c["peakMaxR"],species?72:55);f->SetTitle(Form("%s peak positions;R = |p/q| [GeV/c];dE/dx [keV/cm]",species?"^{4}He":"^{3}He"));TLegend leg(.57,.48,.88,.84);leg.SetBorderSize(0);leg.SetTextSize(.032);for(int response=0;response<2;++response)for(int bg=0;bg<2;++bg){TGraph g;int n=0;for(unsigned i=0;i<rows.size();++i){Row r=rows[i];if(r.response!=responses[response]||r.bg!=backgrounds[bg]||r.r>c["peakMaxR"]||!r.numerical)continue;g.SetPoint(n++,r.r,species?r.m4:r.m3);}g.SetLineColor(colors[2*response+bg]);g.SetMarkerColor(colors[2*response+bg]);g.SetMarkerStyle(20+response);g.SetMarkerSize(.65);TObject* obj=g.DrawClone("LP SAME");leg.AddEntry(obj,Form("%s + %s",responses[response],backgrounds[bg]),"lp");}TGraph anchors;for(unsigned i=0;i<ar.size();++i)anchors.SetPoint(i,ar[i],species?a4[i]:a3[i]);anchors.SetMarkerStyle(29);anchors.SetMarkerColor(kBlack);anchors.SetMarkerSize(2);leg.AddEntry(anchors.DrawClone("P SAME"),"User approximate peaks","p");leg.DrawClone();Text(.14,.88,species?"He4 mean/mode and width fitted to ALL":"He3 low-R center from user anchors; width fitted to ALL");}canvas.Print(pdf.c_str());canvas.Print((pdf+"]").c_str());std::cout<<"R061_DETAIL_COMPLETE"<<std::endl;
}
}
#endif
void PlotR061Detail(const char* directory){
#ifndef __CINT__
try{gStyle->SetOptStat(0);R061Detail::Draw(directory);}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
