void DrawFittedPIDComparison(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "ComparePIDSelections.C"
#include "TGraphErrors.h"
#include <set>
namespace FittedPid {
using namespace PidCompare;
const char* Models[]={"gaussian_linear_nominal","gaussian_exponential_nominal","lognormal_linear_nominal","lognormal_exponential_nominal"};
const char* ModelLabels[]={"Gaussian + linear BG","Gaussian + exp BG","Lognormal + linear BG","Lognormal + exp BG"};
const int FitColors[]={kBlue+1,kCyan+2,kMagenta+2,kGreen+2};
const int FitMarkers[]={20,24,22,26};
bool Matches(const Row& row,const char* region,const char* species){return row.find("region")->second==region&&row.find("species")->second==species;}
double Y(const Row& row,int panel){return N(row,panel==0?"mean":panel==1?"sigma":panel==2?"centerPull":"widthRatio");}
double E(const Row& row,int panel){return N(row,panel==0?"meanError":panel==1?"sigmaError":panel==2?"centerPullDiagonalError":"widthRatioError");}
bool Shown(const Row& row){return N(row,"plotted")>.5;}
void Page(TCanvas& canvas,const Rows& points,const Rows& curves,const char* region,int species,const TString& pdf,const TString& output){
  double xmin=1.e300,xmax=-1.e300;int total=0,shown=0,gray=0;
  for(unsigned i=0;i<points.size();++i)if(Matches(points[i],region,Species[species])){xmin=std::min(xmin,N(points[i],"rLo"));xmax=std::max(xmax,N(points[i],"rHi"));++total;if(Shown(points[i])){++shown;if(N(points[i],"targetFreeSupported")<.5)++gray;}}
  Require(Finite(xmin)&&Finite(xmax)&&xmax>xmin,"No fit interval");canvas.Clear();canvas.Divide(2,2);
  for(int panel=0;panel<4;++panel){canvas.cd(panel+1);Pad();double ymin=1.e300,ymax=-1.e300;
    for(unsigned i=0;i<points.size();++i)if(Matches(points[i],region,Species[species])&&Shown(points[i])){double y=Y(points[i],panel),e=E(points[i],panel);if(Finite(y)&&Finite(e)&&e>=0){ymin=std::min(ymin,y-e);ymax=std::max(ymax,y+e);}}
    if(panel<2)for(unsigned i=0;i<curves.size();++i)if(curves[i].find("species")->second==Species[species]&&N(curves[i],"R")>=xmin&&N(curves[i],"R")<=xmax)for(int p=0;p<3;++p){double y=Value(curves[i],p,panel,2);if(Finite(y)){ymin=std::min(ymin,y);ymax=std::max(ymax,y);}}
    if(!Finite(ymin)||!Finite(ymax)||ymin>ymax){ymin=0;ymax=1;}double span=ymax-ymin;if(span<1.e-10)span=std::max(.05,std::fabs(ymax)*.1);ymin-=.15*span;ymax+=.32*span;if(panel<2)ymin=std::max(0.,ymin);
    TH1F* frame=gPad->DrawFrame(xmin,ymin,xmax,ymax);const char* ytitles[]={"Fit #mu / mode (keV/cm)","Fit #sigma / curvature width (keV/cm)","(#mu_{config}-#mu_{fit}) / #sigma_{fit}","#sigma_{config} / #sigma_{fit}"};frame->SetTitle(Form(";R=|p/q| (GeV/c);%s",ytitles[panel]));
    TLegend* legend=new TLegend(.15,.69,.93,.87);legend->SetBit(TObject::kCanDelete);legend->SetNColumns(2);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.021);
    if(panel<2)for(int p=0;p<3;++p){TGraph* g=Curve(curves,Species[species],p,panel,xmin,xmax,2);g->Draw("L SAME");legend->AddEntry(g,ProfileLabels[p],"l");}
    else{double reference=panel==2?0:1;if(reference>ymin&&reference<ymax)Line(xmin,reference,xmax,reference,kGray+1,3);}
    for(int model=0;model<4;++model){TGraphErrors* graph=new TGraphErrors();graph->SetBit(TObject::kCanDelete);graph->SetLineColor(FitColors[model]);graph->SetMarkerColor(FitColors[model]);graph->SetMarkerStyle(FitMarkers[model]);graph->SetMarkerSize(.45);TGraphErrors* diagnostic=new TGraphErrors();diagnostic->SetBit(TObject::kCanDelete);diagnostic->SetLineColor(kGray+1);diagnostic->SetMarkerColor(kGray+1);diagnostic->SetMarkerStyle(FitMarkers[model]);diagnostic->SetMarkerSize(.32);
      for(unsigned i=0;i<points.size();++i)if(Matches(points[i],region,Species[species])&&points[i].find("label")->second==Models[model]&&Shown(points[i])){double x=N(points[i],"r"),y=Y(points[i],panel),e=E(points[i],panel);if(!Finite(y)||!Finite(e)||e<0)continue;TGraphErrors* use=N(points[i],"targetFreeSupported")>.5?graph:diagnostic;int j=use->GetN();use->SetPoint(j,x,y);use->SetPointError(j,.5*(N(points[i],"rHi")-N(points[i],"rLo")),e);}
      diagnostic->Draw("P SAME");graph->Draw("P SAME");legend->AddEntry(graph,ModelLabels[model],"p");}
    legend->Draw();Text(.13,.95,Form("%s | %s | %s",Labels[species],region,panel==0?"center":panel==1?"width":panel==2?"center residual":"width ratio"),.029);Text(.14,.035,"G: mean/std; LN: mode/curvature. Fixed config has no fitted uncertainty.",.020);Text(.14,.01,Form("Response-quality shown %d/%d; gray identity/reference-unverified %d; not purity validation",shown,total,gray),.019);
  }
  canvas.cd();canvas.Print(pdf);canvas.Print(output+Form("/fitted_%s_%s.png",Species[species],region));
}
}
#endif
void DrawFittedPIDComparison(const char* pointPath,const char* curvePath,const char* coveragePath,const char* settingsPath,const char* outputPath){
#ifndef __CINT__
  using namespace FittedPid;
  try{Rows points=Read(pointPath),curves=Read(curvePath),coverage=Read(coveragePath);Row cfg=Settings(settingsPath);TString output(outputPath),pdf=output+"/fitted_parameter_comparison.pdf";Require(gSystem->AccessPathName(pdf),"Refusing existing fit-parameter PDF");gStyle->SetOptStat(0);TCanvas canvas("fitted_pid_parameters","Fitted response vs configured selection",int(N(cfg,"canvasWidth")),int(N(cfg,"canvasHeight")));canvas.Print(pdf+"[");
    Text(.07,.94,"Re-estimated projected fit response vs immutable configured PID curves",.034);Text(.07,.86,"Low R: target mean/width fitted to ALL. High R: fitted after the specified TOF tag.",.027);Text(.07,.80,"Gaussian: arithmetic mean and standard deviation. Lognormal: mode and local curvature width.",.025);Text(.07,.74,"Only nominal 4 response/background models are plotted; no smoothing or calibration adoption.",.025);Text(.07,.67,"Colored points pass fitted-response precision and reference/identity screens.",.026);Text(.07,.61,"Gray points have stable numerical response but fail reference/identity support; not validated isotope peaks.",.023);Text(.07,.55,"Numerical failures / unstable response errors remain in the TSV and coverage counts, not zero-valued fits.",.023);Text(.07,.47,"The old baseline is the previous report's immutable source snapshot, never the live mainconf.",.025);Text(.07,.40,"Center residual and width-ratio errors use fit-error diagonal terms only (no config error or mu-sigma covariance).",.022);Text(.07,.33,"Finite-R-bin projected responses, NOT intrinsic calibration; adopted low-R views repeat the same fits.",.025);Text(.07,.26,"He4 TOF study uses the saved D m2 window, not the historical broader He4 [3,4] window.",.026);Text(.07,.19,"Native and optional coarse-refit results occupy separate pages; never favorable-point splicing.",.025);Text(.07,.11,Form("Coverage rows: %u. Detailed fit counts and omission reasons: coverage.tsv / fit_points.tsv",unsigned(coverage.size())),.024);Text(.07,.05,"No original ROOT / m2 cache read, no new fit, no PID coefficient change in this rendering.",.024);canvas.Print(pdf);canvas.Print(output+"/fitted_parameter_definitions.png");
    std::set<std::string> regions;for(unsigned i=0;i<points.size();++i)regions.insert(points[i].find("region")->second);
    const char* order[]={"low_ALL","low_adopted_ALL","high_native_TOF","high_coarse_TOF"};for(unsigned region=0;region<4;++region)if(regions.count(order[region]))for(int species=0;species<4;++species){bool exists=false;for(unsigned i=0;i<points.size();++i)exists=exists||Matches(points[i],order[region],Species[species]);if(exists)Page(canvas,points,curves,order[region],species,pdf,output);}
    canvas.Print(pdf+"]");std::cout<<"FITTED_PID_COMPARISON_PASS"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"DrawFittedPIDComparison: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
