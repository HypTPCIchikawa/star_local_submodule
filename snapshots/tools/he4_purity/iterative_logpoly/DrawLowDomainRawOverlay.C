// DISPLAY ONLY: original raw TH2 + saved log10-polynomial coefficients, no fit.
void DrawLowDomainRawOverlay(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
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
namespace LowDomainOverlay {
typedef std::map<std::string,std::string> Settings;
struct Curve {std::string label;std::vector<double> c;};
void Check(bool valid,const char* what){if(!valid)throw std::runtime_error(what);}
double Number(const Settings& s,const std::string& key){
  Settings::const_iterator it=s.find(key);Check(it!=s.end(),"missing explicit display setting");
  std::istringstream in(it->second);double x=0;Check(bool(in>>x),"invalid numeric setting");in>>std::ws;
  Check(in.eof()&&TMath::Finite(x),"invalid display setting");return x;
}
Settings Read(const char* path){std::ifstream in(path);Check(bool(in),"missing settings TSV");Settings s;std::string line;
  while(std::getline(in,line)){size_t k=line.find('\t');Check(k!=std::string::npos,"malformed settings TSV");
    Check(s.insert(std::make_pair(line.substr(0,k),line.substr(k+1))).second,"duplicate display key");}return s;}
std::vector<Curve> Curves(const char* path){std::ifstream in(path);Check(bool(in),"missing coefficient TSV");std::vector<Curve> out;std::string line;
  while(std::getline(in,line)){std::istringstream row(line);Curve c;row>>c.label;c.c.resize(14);
    for(int i=0;i<14;++i)Check(bool(row>>c.c[i])&&TMath::Finite(c.c[i]),"invalid coefficient");
    std::string extra;Check(!(row>>extra),"extra coefficient field");out.push_back(c);}
  Check(out.size()==3&&out[0].label=="latest"&&out[1].label=="outer12"&&out[2].label=="initialQR","unexpected curve sets");return out;}
double Mean(double* x,double* p){double z=std::log10(x[0]);return p[0]+z*(p[1]+z*(p[2]+z*p[3]))+p[7]*(p[4]+z*(p[5]+z*p[6]));}
double Sigma(double* x,double* p){double z=std::log10(x[0]);return p[0]+z*(p[1]+z*p[2]);}
TF1* Center(const Curve& c,int sp,int band,double lo,double hi,const std::string& prefix){
  TF1* f=new TF1(Form("%s_%s_%d_%d",prefix.c_str(),c.label.c_str(),sp,band),Mean,lo,hi,8);
  for(int i=0;i<7;++i)f->SetParameter(i,c.c[7*sp+i]);f->SetParameter(7,band);f->SetNpx(1000);
  f->SetLineColor(sp==0?kBlue+2:kRed+1);f->SetLineWidth(band==0?3:2);f->SetLineStyle(band==0?1:(std::abs(band)==1?2:3));return f;
}
void Caption(const Settings& s,const char* title){
  TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.026);t.DrawLatex(.13,.962,title);
  t.SetTextSize(.019);t.DrawLatex(.13,.060,"DISPLAY ONLY: pre-PID charge-inclusive raw tracks; model-dependent components; PID NOT adopted.");
  t.DrawLatex(.13,.038,Form("Latest numerical convergence=%s; shared Minuit status=%d; curves only for 0.40 #leq p/q #leq 0.70 GeV/c.",
      Number(s,"latestConverged")?"YES":"NO",int(Number(s,"latestGlobalStatus"))));
}
void Map(TH2D* h,const Settings& s,double ymin,double ymax){
  gPad->SetLeftMargin(.13);gPad->SetRightMargin(.14);gPad->SetTopMargin(.12);gPad->SetBottomMargin(.15);
  gPad->SetLogz(int(Number(s,"logZ")));h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->SetMinimum(Number(s,"minimumZ"));h->SetTitle("");h->GetXaxis()->SetTitle("p/q (GeV/c)");h->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  h->GetZaxis()->SetTitle("Raw tracks / displayed bin");h->GetXaxis()->SetTitleSize(.036);h->GetYaxis()->SetTitleSize(.036);
  h->GetXaxis()->SetLabelSize(.031);h->GetYaxis()->SetLabelSize(.031);h->GetZaxis()->SetLabelSize(.025);
  h->Draw("COLZ");for(int k=0;k<2;++k){double x=Number(s,k?"domainMax":"domainMin");TLine* l=new TLine(x,ymin,x,ymax);l->SetLineStyle(7);l->SetLineColor(kGray+2);l->Draw("same");}
}
void SavePage(TCanvas& canvas,const TString& pdf,const TString& output,const char* name){canvas.Print(pdf);canvas.SaveAs(output+"/"+name+".png");}
}
#endif
void DrawLowDomainRawOverlay(const char* sourceRoot,const char* output,const char* settings,const char* coefficients,const char* expected){
#ifndef __CINT__
  using namespace LowDomainOverlay;
  try{
    Settings s=Read(settings);Check(Number(s,"schemaVersion")==1,"unsupported overlay schema");
    TString out(output);const char* names[]={"raw_colz_low_domain_overlay.pdf","overlay_audit.json","raw_colz_full.png","raw_colz_nuclear.png","width_comparison.png","center_comparison.png"};
    for(unsigned i=0;i<sizeof(names)/sizeof(names[0]);++i)Check(gSystem->AccessPathName(out+"/"+names[i]),"existing overlay will not be replaced");
    Check(!gSystem->AccessPathName(out+"/context.json"),"overlay output must be provenance-prepared");
    std::vector<Curve> curves=Curves(coefficients);double lo=Number(s,"domainMin"),hi=Number(s,"domainMax");
    Check(lo==.4&&hi==.7,"only authorized low-domain curves are accepted");
    // Reproduce independently exported grid using the same original formula.
    double maxError=0;unsigned points=0;std::ifstream grid(expected);Check(bool(grid),"missing expected curve grid");
    int mode,sp;double r,mu,sd;
    while(grid>>mode>>sp>>r>>mu>>sd){Check(mode>=0&&mode<3&&sp>=0&&sp<2&&r>=lo&&r<=hi,"invalid expected grid point");
      double p[8];for(int j=0;j<7;++j)p[j]=curves[mode].c[7*sp+j];p[7]=0;double x[]={r};
      double computed=Mean(x,p),w=Sigma(x,p+4);Check(w>0&&TMath::Finite(computed)&&TMath::Finite(w),"nonpositive or nonfinite width/center");
      double err=std::max(std::fabs(computed-mu)/std::max(1.,std::fabs(mu)),std::fabs(w-sd)/std::max(1.,std::fabs(sd)));
      maxError=std::max(maxError,err);++points;}
    Check(grid.eof()&&points>=1800&&maxError<Number(s,"curveRelativeTolerance"),"saved formula/grid mismatch");
    TFile input(sourceRoot,"READ");Check(!input.IsZombie(),"cannot read original ROOT");
    TH2* raw=dynamic_cast<TH2*>(input.Get(s["rawHistogram"].c_str()));Check(raw!=0,"missing raw input TH2");
    double xmin=Number(s,"displayRigidityMin"),xmax=Number(s,"displayRigidityMax"),tol=Number(s,"axisTolerance");
    int first=raw->GetXaxis()->FindBin(xmin+tol),last=raw->GetXaxis()->FindBin(xmax-tol);
    Check(first>=1&&last>=first&&last<=raw->GetNbinsX(),"raw R-axis lacks display domain");
    Check(std::fabs(raw->GetXaxis()->GetBinLowEdge(first)-xmin)<tol&&std::fabs(raw->GetXaxis()->GetBinUpEdge(last)-xmax)<tol,"display R range must match native edges");
    std::vector<double> xe(last-first+2),ye(raw->GetNbinsY()+1);
    for(unsigned i=0;i<xe.size();++i)xe[i]=raw->GetXaxis()->GetBinLowEdge(first+i);
    for(unsigned i=0;i<ye.size();++i)ye[i]=raw->GetYaxis()->GetBinLowEdge(i+1);
    Check(std::fabs(ye.front()-Number(s,"displayDedxMin"))<tol&&std::fabs(ye.back()-Number(s,"displayDedxMax"))<tol,"raw dEdx axis differs from declared display");
    int rebin=int(Number(s,"dedxRebinFactor"));Check(rebin>0&&raw->GetNbinsY()%rebin==0,"invalid exact display rebin");
    TH2D map("raw_display","",last-first+1,&xe[0],raw->GetNbinsY(),&ye[0]);map.SetDirectory(0);map.Sumw2();
    double before=0,after=0;unsigned cells=0;
    for(int bx=first;bx<=last;++bx)for(int by=0;by<=raw->GetNbinsY()+1;++by){double count=raw->GetBinContent(bx,by),error=raw->GetBinError(bx,by);
      Check(TMath::Finite(count)&&count>=0&&TMath::Finite(error)&&error>=0,"invalid raw count/error");
      map.SetBinContent(bx-first+1,by,count);map.SetBinError(bx-first+1,by,error);before+=count;++cells;}
    map.RebinY(rebin);
    for(int bx=1;bx<=map.GetNbinsX();++bx)for(int by=0;by<=map.GetNbinsY()+1;++by)after+=map.GetBinContent(bx,by);
    Check(before==after,"display rebin changed raw count sum");
    gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);gStyle->SetPalette(1);
    TCanvas canvas("low_domain_colz","",int(Number(s,"canvasWidth")),int(Number(s,"canvasHeight")));
    TString pdf=out+"/raw_colz_low_domain_overlay.pdf";canvas.Print(pdf+"[");
    for(int page=0;page<2;++page){canvas.Clear();Map(&map,s,page?Number(s,"nuclearDedxMin"):Number(s,"displayDedxMin"),Number(s,"displayDedxMax"));
      TLegend legend(.16,.78,.62,.88);legend.SetFillColor(kWhite);legend.SetBorderSize(0);legend.SetTextSize(.025);
      for(int isotope=0;isotope<2;++isotope)for(int band=-2;band<=2;++band){TF1* f=Center(curves[0],isotope,band,lo,hi,Form("page%d",page));f->Draw("same");if(!band)legend.AddEntry(f,isotope?"^{4}He assigned center (latest)":"^{3}He assigned center (latest)","l");}
      legend.Draw();TLatex key;key.SetNDC();key.SetTextSize(.022);key.DrawLatex(.17,.752,"Solid: center; dashed: #pm#sigma; dotted: #pm2#sigma (not purity/identity).");
      Caption(s,page?"Latest saved low-domain functions on raw COLZ: nuclear region":"Latest saved low-domain functions on raw COLZ: full recorded dE/dx");
      SavePage(canvas,pdf,out,page?"raw_colz_nuclear":"raw_colz_full");}
    canvas.Clear();canvas.SetLogz(0);canvas.SetLeftMargin(.13);canvas.SetRightMargin(.08);canvas.SetTopMargin(.12);canvas.SetBottomMargin(.15);
    double maxWidth=0;for(unsigned c=0;c<curves.size();++c)for(int isotope=0;isotope<2;++isotope)for(int i=0;i<=300;++i){double x=lo+(hi-lo)*i/300.;maxWidth=std::max(maxWidth,Sigma(&x,&curves[c].c[7*isotope+4]));}
    TH1D frame("width_frame",";p/q (GeV/c);Gaussian #sigma (keV/cm)",100,xmin,xmax);frame.SetMinimum(0);frame.SetMaximum(maxWidth*1.15);frame.Draw();
    TLegend widthLegend(.16,.66,.64,.88);widthLegend.SetFillColor(kWhite);widthLegend.SetBorderSize(0);widthLegend.SetTextSize(.024);
    for(unsigned c=0;c<curves.size();++c)for(int isotope=0;isotope<2;++isotope){TF1* f=new TF1(Form("width_%d_%d",int(c),isotope),Sigma,lo,hi,3);for(int j=0;j<3;++j)f->SetParameter(j,curves[c].c[7*isotope+4+j]);
      f->SetLineColor(isotope?kRed+1:kBlue+2);f->SetLineStyle(c+1);f->SetLineWidth(c?2:3);f->SetNpx(1000);f->Draw("same");
      widthLegend.AddEntry(f,Form("^{%d}He %s",isotope?4:3,curves[c].label.c_str()),"l");}
    widthLegend.Draw();Caption(s,"Width comparison: latest vs outer12 vs initial FREE-derived QR");
    TLatex note;note.SetNDC();note.SetTextSize(.023);note.DrawLatex(.17,.60,"Latest / outer12: intrinsic shared kernels; initial QR: folded-width seed only.");
    SavePage(canvas,pdf,out,"width_comparison");
    canvas.Clear();Map(&map,s,Number(s,"nuclearDedxMin"),Number(s,"displayDedxMax"));
    TLegend centers(.16,.65,.60,.88);centers.SetFillColor(kWhite);centers.SetBorderSize(0);centers.SetTextSize(.024);
    for(unsigned c=0;c<curves.size();++c)for(int isotope=0;isotope<2;++isotope){TF1* f=Center(curves[c],isotope,0,lo,hi,"reference");f->SetLineStyle(c+1);f->SetLineWidth(c?2:3);f->Draw("same");centers.AddEntry(f,Form("^{%d}He %s",isotope?4:3,curves[c].label.c_str()),"l");}
    centers.Draw();Caption(s,"Center comparison on the identical raw COLZ: study references only");SavePage(canvas,pdf,out,"center_comparison");canvas.Print(pdf+"]");
    std::ofstream audit((out+"/overlay_audit.json").Data());Check(bool(audit),"cannot write display audit");
    audit<<std::setprecision(17)<<"{\"displayOnly\":true,\"pages\":4,\"rawCellsCopied\":"<<cells<<",\"rawCountsBeforeDisplayRebin\":"<<before<<",\"rawCountsAfterDisplayRebin\":"<<after
      <<",\"nativeRBinWidth\":"<<raw->GetXaxis()->GetBinWidth(first)<<",\"nativeDedxBinWidth\":"<<raw->GetYaxis()->GetBinWidth(1)<<",\"displayDedxBinWidth\":"<<map.GetYaxis()->GetBinWidth(1)
      <<",\"curveGridPoints\":"<<points<<",\"maximumFormulaRelativeError\":"<<maxError<<",\"curveDomain\":[0.4,0.7],\"outsideDomainDrawn\":false,\"fitPerformed\":false,\"pidAdopted\":false}\n";
    std::cout<<"LOW_DOMAIN_RAW_OVERLAY_COMPLETE pages=4 maximumFormulaRelativeError="<<maxError<<std::endl;
  }catch(const std::exception& e){std::cerr<<"DrawLowDomainRawOverlay: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
