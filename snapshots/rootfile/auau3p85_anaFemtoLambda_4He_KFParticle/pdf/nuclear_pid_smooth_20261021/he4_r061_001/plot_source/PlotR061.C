void PlotR061(const char*,const char*,int=0);
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
std::vector<std::string> split(const std::string& s){std::vector<std::string> v;std::istringstream f(s);std::string x;while(std::getline(f,x,'\t'))v.push_back(x);return v;}
void Text(double x,double y,const std::string& s,double size=.035){TLatex t;t.SetNDC();t.SetTextSize(size);t.DrawLatex(x,y,s.c_str());}
struct FitRow {std::string band,sp,model,root;int index;double lo,hi,fitlo,fithi,mean,sigma,p,linearLo,linearHi,logLo,logHi;double he3Mean,he3Sigma;double chi2,ndf,deviance;bool identity;bool numerical,supported;};
void Atlas(const std::string& dir){std::ifstream input((dir+"/atlas.tsv").c_str());std::string line;std::getline(input,line);std::map<std::string,std::map<int,std::vector<FitRow> > > groups;while(std::getline(input,line)){std::vector<std::string> v=split(line);FitRow r;r.band=v[0];r.sp=v[1];r.model=v[2];r.index=atoi(v[3].c_str());r.lo=atof(v[4].c_str());r.hi=atof(v[5].c_str());r.fitlo=atof(v[6].c_str());r.fithi=atof(v[7].c_str());r.mean=atof(v[8].c_str());r.sigma=atof(v[9].c_str());r.root=v[10];r.p=atof(v[11].c_str());r.numerical=v[12]=="1";r.supported=v[13]=="1";r.linearLo=atof(v[14].c_str());r.linearHi=atof(v[15].c_str());r.logLo=atof(v[16].c_str());r.logHi=atof(v[17].c_str());r.chi2=atof(v[18].c_str());r.ndf=atof(v[19].c_str());r.deviance=atof(v[20].c_str());r.identity=v[21]=="1";r.he3Mean=atof(v[22].c_str());r.he3Sigma=atof(v[23].c_str());groups[r.band+"_"+r.sp][int(std::floor(r.lo*10000+.5))].push_back(r);}
std::map<std::string,TFile*> files;TCanvas c("atlas","All R fits",1500,1500);for(std::map<std::string,std::map<int,std::vector<FitRow> > >::iterator group=groups.begin();group!=groups.end();++group){std::string pdf=dir+"/fit_atlas_"+group->first+".pdf";c.Print((pdf+"[").c_str());int pages=0;for(std::map<int,std::vector<FitRow> >::iterator bin=group->second.begin();bin!=group->second.end();++bin){std::vector<FitRow>& rows=bin->second;if(rows.size()!=3)throw std::runtime_error("Atlas needs three models");c.Clear();c.Divide(2,3,.003,.003);for(int model=0;model<3;++model){FitRow r=rows[model];if(!files.count(r.root))files[r.root]=new TFile(r.root.c_str(),"READ");TFile* file=files[r.root];for(int log=0;log<2;++log){c.cd(2*model+log+1);gPad->SetLogy(log);gPad->SetLeftMargin(.12);gPad->SetBottomMargin(.14);TH1* raw0=dynamic_cast<TH1*>(file->Get(Form("raw_native_%03d",r.index)));if(!raw0)throw std::runtime_error("Missing raw spectrum");TH1* raw=dynamic_cast<TH1*>(raw0->Clone("raw_view"));raw->SetDirectory(0);raw->Rebin(5);double lo=log?r.logLo:r.linearLo,hi=log?r.logHi:r.linearHi;raw->GetXaxis()->SetRangeUser(lo,hi);raw->SetTitle(Form("%s %s R=[%.3f,%.3f) %s;dE/dx [keV/cm];Counts / 0.5 keV/cm",r.sp.c_str(),r.band.c_str(),r.lo,r.hi,r.model.c_str()));raw->GetXaxis()->SetTitleSize(.045);raw->GetYaxis()->SetTitleSize(.045);double max=0;for(int b=1;b<=raw->GetNbinsX();++b)if(raw->GetBinCenter(b)>=lo&&raw->GetBinCenter(b)<=hi)max=std::max(max,raw->GetBinContent(b)+raw->GetBinError(b));raw->SetMinimum(log?.5:0);raw->SetMaximum(std::max(1.,max)*(log?4.:1.22));raw->SetMarkerStyle(20);raw->SetMarkerSize(.25);raw->SetLineColor(kBlack);TH1* drawnRaw=raw->DrawCopy("E");TLegend leg(.56,.49,.89,.84);leg.SetBorderSize(0);leg.SetTextSize(.031);leg.AddEntry(drawnRaw,"Data","lep");const char* parts[6]={"total","smooth","d","t","he3","he4"};int partcolors[6]={kRed+1,kGray+2,kBlue+1,kOrange+7,kGreen+2,kMagenta+1};for(int part=0;part<6;++part){if(part==1&&r.model.find("none_")==0)continue;TH1* h0=dynamic_cast<TH1*>(file->Get(Form("prediction_%s_%03d",parts[part],r.index)));if(!h0)throw std::runtime_error("Missing fit component");TH1* h=dynamic_cast<TH1*>(h0->Clone("part_view"));h->SetDirectory(0);h->Rebin(5);TGraph g;int n=0;for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=std::max(lo,r.fitlo)&&h->GetBinCenter(b)<=std::min(hi,r.fithi)){g.SetPoint(n++,h->GetBinCenter(b),h->GetBinContent(b));}g.SetLineColor(partcolors[part]);g.SetLineWidth(part==0?2:1);g.SetLineStyle(part==1?2:1);TObject* drawn=g.DrawClone("L SAME");leg.AddEntry(drawn,parts[part],"l");delete h;}leg.DrawClone();Text(.15,.80,Form("#chi^{2}/NDF=%.1f/%.0f = %.2f; D=%.1f",r.chi2,r.ndf,r.ndf>0?r.chi2/r.ndf:0,r.deviance),.027);Text(.15,.74,Form("Fit=[%.1f,%.1f]; identity=%d",r.fitlo,r.fithi,r.identity),.026);Text(.15,.68,Form("Peaks: ^{3}He=%.2f, ^{4}He=%.2f",r.he3Mean,r.mean),.026);TLine edge(r.fitlo,log?.5:0,r.fitlo,std::max(1.,max));edge.SetLineStyle(3);edge.SetLineColor(kGray+1);edge.DrawClone();Text(.15,.86,Form("P=%.4f; numerical=%d; historical quality=%d",r.p,r.numerical,r.supported),.031);delete raw;}}c.Print(pdf.c_str());++pages;}c.Print((pdf+"]").c_str());std::cout<<"ATLAS_COMPLETE "<<group->first<<" pages="<<pages<<std::endl;}for(std::map<std::string,TFile*>::iterator it=files.begin();it!=files.end();++it){it->second->Close();delete it->second;}}
}
#endif
void PlotR061(const char* directory,const char* cache,int atlas){
#ifndef __CINT__
try{gStyle->SetOptStat(0);gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);SmoothPlot::Atlas(directory);std::cout<<"SMOOTH_PLOT_COMPLETE"<<std::endl;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
