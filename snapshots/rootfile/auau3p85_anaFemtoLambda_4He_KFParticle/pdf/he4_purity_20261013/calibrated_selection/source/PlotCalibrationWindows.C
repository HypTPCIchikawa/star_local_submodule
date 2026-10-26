// Plot serialized, fixed-selection model integrals. No refitting or ROOT input writes.
void PlotCalibrationWindows(const char* points,const char* curves,const char* output);
#ifndef __CINT__
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TPad.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
namespace CalibrationWindows {
struct Point {double r,p,low,high,ratio;int selection,model,chosen,quality;};
struct Curve {double r,mean,sigma,oldMean,oldSigma;int species,active;};
const char* Labels[]={"Original calibration |n#sigma_{4}|<2", "Applied calibration |n#sigma_{4}|<2", "Applied calibration + TPC best species"};
int Colors[]={kBlack,kRed+1,kBlue+1};
TGraph* Graph(const std::vector<Point>& points,int selection,int model,bool chosen,bool ratio=false,int quality=-1) {
  TGraph* g=new TGraph();
  for(unsigned i=0;i<points.size();++i){const Point& p=points[i];if(p.selection!=selection||(model>=0&&p.model!=model)||(chosen&&!p.chosen)||(quality>=0&&p.quality!=quality))continue;g->SetPoint(g->GetN(),p.r,ratio?p.ratio:p.p);}
  g->SetLineColor(model<0?Colors[selection]:(model==0?kRed+1:kBlue+1));g->SetMarkerColor(g->GetLineColor());
  g->SetLineWidth(2);g->SetMarkerStyle(quality==1?20:24);g->SetMarkerSize(.55);return g;
}
void Caption(const char* extra="") {TLatex t;t.SetNDC();t.SetTextSize(.025);t.DrawLatex(.12,.94,"Au+Au FXT 3.85 GeV beam energy | 2021 P24iy | centrality-integrated, charge-inclusive");t.DrawLatex(.12,.90,"Model-dependent TPC proxy: open markers = reference fits; no TOF / pair-conditioned purity");if(*extra)t.DrawLatex(.12,.86,extra);}
void Frame(const char* name,const char* title,double maximum=1.08){TH1F* frame=new TH1F(name,title,1,.4,2.5);frame->SetMinimum(0);frame->SetMaximum(maximum);frame->SetDirectory(0);frame->GetXaxis()->SetTitleOffset(1.1);frame->GetYaxis()->SetTitleOffset(1.2);frame->Draw();}
}
#endif
void PlotCalibrationWindows(const char* points,const char* curves,const char* output) {
#ifndef __CINT__
  using namespace CalibrationWindows;
  try {
    std::ifstream input(points);if(!input)throw std::runtime_error("missing serialized window points");
    std::vector<Point> data;Point p;while(input>>p.r>>p.selection>>p.model>>p.chosen>>p.quality>>p.p>>p.low>>p.high>>p.ratio)data.push_back(p);
    if(data.empty())throw std::runtime_error("empty window points");
    std::ifstream curveInput(curves);if(!curveInput)throw std::runtime_error("missing calibration curves");
    std::vector<Curve> curveData;Curve c;while(curveInput>>c.r>>c.species>>c.mean>>c.sigma>>c.oldMean>>c.oldSigma>>c.active)curveData.push_back(c);
    gStyle->SetOptStat(0);gStyle->SetPaperSize(30,20);TCanvas canvas("calibration_windows","Calibration windows",1500,1000);
    const TString book=TString(output)+"/purity_calibrated_windows.pdf";canvas.Print(book+"[");
    canvas.SetLeftMargin(.12);canvas.SetBottomMargin(.12);canvas.SetTopMargin(.20);canvas.SetRightMargin(.04);
    Frame("purity_frame",";R=|p/q| [GeV/c];Model S/(S+B)");
    TLegend legend(.45,.17,.92,.35);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.027);
    for(int selection=0;selection<3;++selection){
      TGraphAsymmErrors* band=new TGraphAsymmErrors();
      for(unsigned j=0;j<data.size();++j){const Point& q=data[j];if(q.selection!=selection||!q.chosen)continue;int k=band->GetN();band->SetPoint(k,q.r,q.p);band->SetPointError(k,0,0,q.p-q.low,q.high-q.p);}
      band->SetFillColor(Colors[selection]);band->SetFillStyle(3004+selection);band->Draw("3 same");
      TGraph* graph=Graph(data,selection,-1,true);graph->Draw("LP same");Graph(data,selection,-1,true,false,1)->Draw("P same");legend.AddEntry(graph,Labels[selection],"lp");
    }
    legend.Draw();Caption("Hatched range = linear/exponential spread at the SAME selection; not a statistical error");canvas.Print(book);canvas.Print(TString(output)+"/purity_calibrated_windows.png");
    canvas.Clear();canvas.Divide(3,1);
    for(int selection=0;selection<3;++selection){canvas.cd(selection+1);gPad->SetLeftMargin(.16);gPad->SetBottomMargin(.13);gPad->SetTopMargin(.18);Frame(Form("model_%d",selection),Form("%s;R [GeV/c];Model S/(S+B)",Labels[selection]));TLegend* l=new TLegend(.2,.2,.85,.31);l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(.035);for(int model=0;model<2;++model){TGraph* g=Graph(data,selection,model,false);g->Draw("LP same");Graph(data,selection,model,false,false,1)->Draw("P same");l->AddEntry(g,model==0?"linear background":"exponential background","lp");}l->Draw();TLatex t;t.SetNDC();t.SetTextSize(.033);t.DrawLatex(.18,.88,"TPC proxy; open = reference");}
    canvas.Print(book);
    canvas.Clear();canvas.cd();canvas.SetTopMargin(.20);double maximum=1.;for(unsigned i=0;i<data.size();++i)if(data[i].chosen)maximum=std::max(maximum,data[i].ratio);
    Frame("signal_ratio",";R=|p/q| [GeV/c];Model signal count / original-gate model signal count",maximum*1.2);
    TLegend ratioLegend(.45,.68,.92,.8);ratioLegend.SetBorderSize(0);ratioLegend.SetFillStyle(0);ratioLegend.SetTextSize(.027);
    for(int selection=1;selection<3;++selection){TGraph* g=Graph(data,selection,-1,true,true);g->Draw("LP same");Graph(data,selection,-1,true,true,1)->Draw("P same");ratioLegend.AddEntry(g,Labels[selection],"lp");}
    ratioLegend.Draw();Caption("Model window acceptance only: this ratio is NOT an independently measured reconstruction efficiency");canvas.Print(book);
    canvas.Clear();canvas.Divide(2,2);
    for(int species=2;species<4;++species)for(int coordinate=0;coordinate<2;++coordinate){canvas.cd((species-2)*2+coordinate+1);gPad->SetLeftMargin(.15);gPad->SetBottomMargin(.13);gPad->SetTopMargin(.15);TGraph* active=new TGraph();TGraph* legacy=new TGraph();double top=0;for(unsigned j=0;j<curveData.size();++j){const Curve& q=curveData[j];if(q.species!=species)continue;double a=coordinate?q.sigma:q.mean,b=coordinate?q.oldSigma:q.oldMean;active->SetPoint(active->GetN(),q.r,a);legacy->SetPoint(legacy->GetN(),q.r,b);top=std::max(top,std::max(a,b));}Frame(Form("cal_%d_%d",species,coordinate),Form("%s %s;R [GeV/c];dE/dx [keV/cm]",species==2?"^{3}He":"^{4}He",coordinate?"Gaussian width":"center"),top*1.25);active->SetLineColor(kBlue+1);active->SetLineWidth(2);legacy->SetLineColor(kBlack);legacy->SetLineStyle(2);legacy->SetLineWidth(2);legacy->Draw("L same");active->Draw("L same");TLegend* l=new TLegend(.48,.69,.93,.82);l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(.035);l->AddEntry(active,"Applied function (+ legacy gaps)","l");l->AddEntry(legacy,"Original function","l");l->Draw();}
    canvas.Print(book);canvas.Print(book+"]");std::cout<<"CALIBRATION_WINDOWS_COMPLETE"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"PlotCalibrationWindows: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
