void DrawPidPurity(const char*,const char*);
#ifndef __CINT__
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
namespace PidPurityView {
struct Point {std::string species,mode;double r,p,e,low;bool assessed;};
void Text(double x,double y,const char* t,double size){TLatex l;l.SetNDC();l.SetTextSize(size);l.DrawLatex(x,y,t);}
void Run(const char* table,const char* pdf){std::ifstream in(table);if(!in)throw std::runtime_error("Cannot read purity table");std::string line;std::getline(in,line);std::vector<Point> points;
  while(std::getline(in,line)){std::istringstream row(line);Point q;int models,assessed;std::string reason;double signal,count,chi;if(!(row>>q.species>>q.r>>q.mode>>q.p>>q.e>>q.low>>models>>assessed>>reason>>signal>>count>>chi))throw std::runtime_error("Malformed purity row");q.assessed=assessed!=0;points.push_back(q);}
  const char* names[]={"d","t","he3","he4"};const char* labels[]={"d","t","^{3}He","^{4}He"};gStyle->SetOptStat(0);gStyle->SetPaperSize(28.3,20.);TCanvas c("puritySummaryCanvas","Nuclear purity with quality masking",1600,1100);c.SetCanvasSize(1600,1100);c.Print((std::string(pdf)+"[").c_str());
  for(unsigned page=0;page<2;++page){c.Clear();c.Divide(2,2,.005,.005);std::vector<TH1F*> frames;std::vector<TGraph*> graphs;
    for(unsigned s=0;s<4;++s){c.cd(s+1);gPad->SetLeftMargin(.15);gPad->SetRightMargin(.035);gPad->SetTopMargin(.22);gPad->SetBottomMargin(.16);double xmax=page==0?5.:2.;TH1F* f=new TH1F(Form("purityFrame_%u_%u",page,s),"",100,.4,xmax);f->SetDirectory(0);frames.push_back(f);f->SetMinimum(-.10);f->SetMaximum(1.08);f->GetXaxis()->SetTitle("p/q [GeV/c]");f->GetYaxis()->SetTitle("Purity S/(S+B)");f->GetYaxis()->SetTitleOffset(1.05);f->Draw("AXIS");TLine threshold(.4,.5,xmax,.5);threshold.SetLineStyle(2);threshold.SetLineColor(kRed+1);threshold.DrawClone();TLine transition(2.,-.1,2.,1.08);transition.SetLineStyle(3);transition.DrawClone();
      TGraphErrors *all=new TGraphErrors(),*tof=new TGraphErrors();TGraph *lower=new TGraph(),*bad=new TGraph();all->SetMarkerColor(kBlue+1);all->SetLineColor(kBlue+1);all->SetMarkerStyle(20);tof->SetMarkerColor(kGreen+2);tof->SetLineColor(kGreen+2);tof->SetMarkerStyle(21);lower->SetMarkerColor(kBlack);lower->SetMarkerStyle(25);bad->SetMarkerColor(kGray+1);bad->SetMarkerStyle(5);bad->SetMarkerSize(.65);int valid=0,total=0;
      for(unsigned j=0;j<points.size();++j){const Point& q=points[j];if(q.species!=names[s]||q.r>xmax)continue;++total;if(q.assessed){++valid;TGraphErrors* g=q.mode=="ALL"?all:tof;int n=g->GetN();g->SetPoint(n,q.r,q.p);g->SetPointError(n,0.,q.e);lower->SetPoint(lower->GetN(),q.r,q.low);}else bad->SetPoint(bad->GetN(),q.r,-.05);}
      all->Draw("P SAME");tof->Draw("P SAME");lower->Draw("P SAME");bad->Draw("P SAME");graphs.push_back(all);graphs.push_back(tof);graphs.push_back(lower);graphs.push_back(bad);
      TLegend legend(.50,.80,.97,.98);legend.SetNColumns(2);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.028);legend.AddEntry(all,"ALL, R<2","pe");legend.AddEntry(tof,"TOF, R#geq2","pe");legend.AddEntry(lower,"Valid-model minimum","p");legend.AddEntry(bad,"Unassessed (not P=0)","p");legend.DrawClone();Text(.12,.955,Form("%s; accepted %d/%d",labels[s],valid,total),.040);Text(.12,.89,"Common polynomial candidates; NOT production-adopted",.028);Text(.12,.045,"Bad fits excluded. Errors conditional on fixed curves; no isolated veto.",.025);gPad->RedrawAxis();
    }c.Print(pdf);c.Clear();for(unsigned j=0;j<graphs.size();++j)delete graphs[j];for(unsigned j=0;j<frames.size();++j)delete frames[j];}
  c.Print((std::string(pdf)+"]").c_str());std::cout<<"PID_PURITY_DRAW_COMPLETE pages=2"<<std::endl;
}
}
#endif
void DrawPidPurity(const char* table,const char* pdf){
#ifndef __CINT__
  try{PidPurityView::Run(table,pdf);}catch(const std::exception& e){std::cerr<<"PID_PURITY_DRAW_FATAL: "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
