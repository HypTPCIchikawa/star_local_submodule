// ROOT5 postprocessing of an existing pre-PID TH2; never updates the input file.
#ifndef __CINT__
#include "../../StMaker/common/NuclearIdDeDxVsMom.h"
#endif
#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
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

#ifndef __CINT__
namespace He4Purity {
std::string Quote(const std::string& s) {
  std::ostringstream o; o << '"';
  for(unsigned i=0;i<s.size();++i) { if(s[i]=='"'||s[i]=='\\') o << '\\'; if(s[i]=='\n') o << "\\n"; else o << s[i]; }
  o << '"'; return o.str();
}
double Required(const std::map<std::string,std::string>& c,const std::string& key) {
  std::map<std::string,std::string>::const_iterator i=c.find(key);
  if(i==c.end()) throw std::runtime_error("missing explicit setting: "+key);
  std::istringstream v(i->second); double d; v>>d;
  if(!v || !TMath::Finite(d)) throw std::runtime_error("invalid numeric setting: "+key);
  return d;
}
std::map<std::string,std::string> Settings(const char* path) {
  std::ifstream f(path); if(!f) throw std::runtime_error("cannot read explicit settings TSV");
  std::map<std::string,std::string> out; std::string line;
  while(std::getline(f,line)) { const size_t split=line.find('\t'); if(split==std::string::npos) throw std::runtime_error("invalid settings TSV"); out[line.substr(0,split)]=line.substr(split+1); }
  return out;
}
double GaussianIntegral(const std::vector<double>& p,int index,double lo,double hi) {
  if(!(p[index+2]>0)) return 0;
  const double k=std::sqrt(2.)*p[index+2];
  return p[index]*p[index+2]*std::sqrt(TMath::Pi()/2.)*(TMath::Erf((hi-p[index+1])/k)-TMath::Erf((lo-p[index+1])/k));
}
double SmoothIntegral(const std::vector<double>& p,bool exponential,double fitLo,double fitHi,double lo,double hi) {
  if(exponential) {
    const double slope=p[4];
    if(std::fabs(slope)<1.e-8) return p[3]*(hi-lo);
    return p[3]*(std::exp(slope*hi)-std::exp(slope*lo))/slope;
  }
  return p[3]*(hi-lo)+(p[4]-p[3])/((fitHi-fitLo)*2.)*((hi-fitLo)*(hi-fitLo)-(lo-fitLo)*(lo-fitLo));
}
double Purity(const std::vector<double>& p,bool exp,double fitLo,double fitHi,double gate,double* signal=0,double* background=0) {
  const double s=GaussianIntegral(p,0,-gate,gate);
  double b=SmoothIntegral(p,exp,fitLo,fitHi,-gate,gate);
  for(unsigned k=5;k+2<p.size();k+=3) b+=GaussianIntegral(p,k,-gate,gate);
  if(signal) *signal=s; if(background) *background=b;
  return s+b>0?s/(s+b):-1;
}
struct Fit {
  std::string model; TF1* total; TF1* signal; TF1* background;
  double purity,error,chi2,deviance,aic,gateSignal,gateBackground;
  int status,cov,ndf; bool valid,quality,atBoundary;
  std::vector<double> par,err;
  Fit():total(0),signal(0),background(0),purity(-1),error(0),chi2(0),deviance(0),aic(0),gateSignal(0),gateBackground(0),status(-1),cov(-1),ndf(0),valid(false),quality(false),atBoundary(false) {}
};
struct Slice {
  double lo,hi,r,mean,sigma,counts,gateCounts,drift,gateLoDrift,gateHiDrift;
  TH1D* h; std::vector<Fit> fits; int chosen;
  Slice():lo(0),hi(0),r(0),mean(0),sigma(0),counts(0),gateCounts(0),drift(0),gateLoDrift(0),gateHiDrift(0),h(0),chosen(-1) {}
};
std::string GaussianFormula(int i) { std::ostringstream f;f<<"["<<i<<"]*exp(-0.5*((x-["<<i+1<<"])/["<<i+2<<"])*((x-["<<i+1<<"])/["<<i+2<<"]))";return f.str(); }
Fit FitOne(TH1D* h,double r,double mu,double sig,bool exponential,const std::map<std::string,std::string>& cfg,int number) {
  Fit out; out.model=exponential?"exponential":"linear";
  const double low=std::max(Required(cfg,"fitNSigmaMin"),h->GetXaxis()->GetXmin());
  const double high=std::min(Required(cfg,"fitNSigmaMax"),h->GetXaxis()->GetXmax());
  const double gate=Required(cfg,"gateNSigma");
  if(!(low < -gate && high > gate)) return out;
  std::vector<double> neighborMu,neighborSig;
  for(int s=0;s<3;++s) {
    const NuclearIdDeDxVsMom::SpeciesIndex sp=static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s);
    const double nm=(NuclearIdDeDxVsMom::GetMean(sp,r)-mu)/sig, ns=NuclearIdDeDxVsMom::GetSigma(sp,r)/sig;
    if(ns>0 && nm>low && nm<high) { neighborMu.push_back(nm);neighborSig.push_back(ns); }
  }
  std::ostringstream smooth;
  if(exponential) smooth<<"[3]*exp([4]*x)";
  else smooth<<"[3]+([4]-[3])*(x-("<<std::setprecision(17)<<low<<"))/("<<high-low<<")";
  std::string formula=GaussianFormula(0)+"+"+smooth.str();
  for(unsigned j=0;j<neighborMu.size();++j) formula+="+"+GaussianFormula(5+3*j);
  const TString name=Form("he4_%s_%d",out.model.c_str(),number);
  out.total=new TF1(name,formula.c_str(),low,high); out.total->SetNpx(400);
  const double maximum=std::max(1.,h->GetMaximum());
  const double peak=std::max(1.,h->GetBinContent(h->FindBin(0.)));
  const double bg=std::max(1.,h->GetBinContent(h->FindBin(high-.5))*.2);
  out.total->SetParameter(0,peak);out.total->SetParLimits(0,0,20.*maximum);
  out.total->SetParameter(1,0);out.total->SetParLimits(1,-Required(cfg,"signalMeanShiftMax"),Required(cfg,"signalMeanShiftMax"));
  out.total->SetParameter(2,1.);out.total->SetParLimits(2,Required(cfg,"signalSigmaMin"),Required(cfg,"signalSigmaMax"));
  out.total->SetParameter(3,bg);out.total->SetParLimits(3,0.,20.*maximum);
  out.total->SetParameter(4,exponential?-.2:bg);
  out.total->SetParLimits(4,exponential?-3.:0.,exponential?3.:20.*maximum);
  for(unsigned j=0;j<neighborMu.size();++j) {
    const int index=5+3*j; const double nm=neighborMu[j],ns=neighborSig[j];
    out.total->SetParameter(index,std::max(1.,h->GetBinContent(h->FindBin(nm))));out.total->SetParLimits(index,0,20.*maximum);
    out.total->SetParameter(index+1,nm);out.total->SetParLimits(index+1,nm-Required(cfg,"neighborMeanShiftMax")*ns,nm+Required(cfg,"neighborMeanShiftMax")*ns);
    out.total->SetParameter(index+2,ns);out.total->SetParLimits(index+2,Required(cfg,"neighborSigmaFactorMin")*ns,Required(cfg,"neighborSigmaFactorMax")*ns);
  }
  TFitResultPtr result=h->Fit(out.total,"SQLNR0I","",low,high);
  out.status=int(result);out.cov=result.Get()?result->CovMatrixStatus():-1;
  const int np=out.total->GetNpar();
  out.par.resize(np);out.err.resize(np);
  for(int p=0;p<np;++p) { out.par[p]=out.total->GetParameter(p);out.err[p]=out.total->GetParError(p); }
  out.purity=Purity(out.par,exponential,low,high,gate,&out.gateSignal,&out.gateBackground);
  // The fit functions predict bin heights. Integrals / original u-bin width are counts.
  const double width=h->GetBinWidth(1);out.gateSignal/=width;out.gateBackground/=width;
  if(result.Get()) {
    std::vector<double> gradient(np);
    for(int p=0;p<np;++p) {
      const double step=std::max(1.e-7,std::fabs(out.par[p])*1.e-5);
      std::vector<double> plus=out.par,minus=out.par;plus[p]+=step;minus[p]-=step;
      gradient[p]=(Purity(plus,exponential,low,high,gate)-Purity(minus,exponential,low,high,gate))/(2*step);
    }
    double variance=0;
    for(int a=0;a<np;++a) for(int b=0;b<np;++b) variance+=gradient[a]*result->CovMatrix(a,b)*gradient[b];
    out.error=variance>0?std::sqrt(variance):0;
  }
  int bins=0;bool positive=true;
  for(int b=1;b<=h->GetNbinsX();++b) {
    const double x=h->GetBinCenter(b);if(x<low || x>=high)continue;
    const double n=h->GetBinContent(b),expected=out.total->Integral(h->GetXaxis()->GetBinLowEdge(b),h->GetXaxis()->GetBinUpEdge(b))/width;
    if(!(expected>0 && TMath::Finite(expected))) {positive=false;continue;}
    ++bins;out.chi2+=(n-expected)*(n-expected)/expected;
    out.deviance+=2.*(expected-n+(n>0?n*std::log(n/expected):0));
  }
  out.ndf=bins-np;out.aic=out.deviance+2*np;
  for(int p=0;p<np;++p) {
    double l=0,u=0;out.total->GetParLimits(p,l,u);
    if(u>l && (std::fabs(out.par[p]-l)<1.e-4*(u-l) || std::fabs(out.par[p]-u)<1.e-4*(u-l))) {
      // A vanishing nuisance amplitude is allowed. Signal shape limits are not.
      if(p==0 || p==1 || p==2 || (p>=5 && (p-5)%3!=0))out.atBoundary=true;
    }
  }
  out.valid=result.Get()&&result->IsValid()&&out.status==0&&out.cov>=2&&positive&&!out.atBoundary&&out.ndf>0&&TMath::Finite(out.purity)&&out.purity>0&&out.purity<=1&&TMath::Finite(out.error);
  out.quality=out.valid&&out.chi2/out.ndf<=Required(cfg,"maximumChi2Ndf")&&out.error/out.purity<=Required(cfg,"maximumRelativePurityError");
  out.signal=new TF1(name+"_signal",GaussianFormula(0).c_str(),low,high);out.signal->SetParameters(&out.par[0]);
  std::string bgformula=smooth.str();for(unsigned j=0;j<neighborMu.size();++j)bgformula+="+"+GaussianFormula(5+3*j);
  out.background=new TF1(name+"_bg",bgformula.c_str(),low,high);out.background->SetParameters(&out.par[0]);
  out.total->SetLineColor(exponential?kBlue+1:kRed+1);out.total->SetLineWidth(2);
  out.signal->SetLineColor(kGreen+2);out.signal->SetLineWidth(2);
  out.background->SetLineColor(kMagenta+1);out.background->SetLineStyle(2);out.background->SetLineWidth(2);
  return out;
}
void DrawSlice(const Slice& s,double gate,bool logY) {
  gPad->SetLeftMargin(.12);gPad->SetBottomMargin(.12);gPad->SetTopMargin(.12);gPad->SetRightMargin(.03);
  if(logY)gPad->SetLogy();
  s.h->GetXaxis()->SetRangeUser(-7.,4.);s.h->SetMinimum(.5);s.h->SetMaximum(std::max(10.,s.h->GetMaximum()*8));
  s.h->SetTitle(Form("R=%.2f-%.2f GeV/c;calibrated n#sigma_{^{4}He};tracks / saved dE/dx bin",s.lo,s.hi));
  s.h->SetLineColor(kBlack);s.h->SetMarkerStyle(20);s.h->SetMarkerSize(.25);s.h->Draw("E");
  for(unsigned j=0;j<s.fits.size();++j)if(s.fits[j].total)s.fits[j].total->Draw("same");
  if(s.chosen>=0){s.fits[s.chosen].signal->Draw("same");s.fits[s.chosen].background->Draw("same");}
  for(int sign=-1;sign<=1;sign+=2){TLine* l=new TLine(sign*gate,.5,sign*gate,std::max(10.,s.h->GetMaximum()));l->SetLineStyle(3);l->Draw();}
  TLatex t;t.SetNDC();t.SetTextSize(.035);
  if(s.chosen>=0){const Fit& f=s.fits[s.chosen];t.DrawLatex(.13,.84,Form("P=%.3f #pm %.3f (%s)%s",f.purity,f.error,f.model.c_str(),f.quality?"":" [poor fit]"));t.DrawLatex(.13,.79,Form("#chi^{2}/ndf=%.1f | status=%d | covariance=%d",f.ndf>0?f.chi2/f.ndf:-1,f.status,f.cov));}
  else t.DrawLatex(.13,.84,"No valid peak fit / insufficient counts");
}
void WriteFit(std::ostream& o,const Fit& f) {
  o<<"{\"model\":"<<Quote(f.model)<<",\"status\":"<<f.status<<",\"covarianceStatus\":"<<f.cov<<",\"valid\":"<<(f.valid?"true":"false")<<",\"quality\":"<<(f.quality?"true":"false")<<",\"atBoundary\":"<<(f.atBoundary?"true":"false")<<",\"purity\":"<<f.purity<<",\"purityError\":"<<f.error<<",\"chi2\":"<<f.chi2<<",\"ndf\":"<<f.ndf<<",\"deviance\":"<<f.deviance<<",\"aic\":"<<f.aic<<",\"gateSignal\":"<<f.gateSignal<<",\"gateBackground\":"<<f.gateBackground<<",\"parameters\":[";
  for(unsigned j=0;j<f.par.size();++j){if(j)o<<',';o<<f.par[j];}o<<"],\"parameterErrors\":[";
  for(unsigned j=0;j<f.err.size();++j){if(j)o<<',';o<<f.err[j];}o<<"]}";
}
}
#endif

void FitHe4Purity(const char* input,const char* output,const char* settings,const char* label) {
#ifndef __CINT__
  using namespace He4Purity;
  try {
    const std::map<std::string,std::string> cfg=Settings(settings);
    const double rmin=Required(cfg,"rigidityMin"),rmax=Required(cfg,"rigidityMax"),step=Required(cfg,"rigidityStep"),gate=Required(cfg,"gateNSigma");
    if(!(rmin>=.4&&rmax<=5.&&rmax>rmin&&step>0&&gate>0))throw std::runtime_error("invalid calibration range / gate");
    const TString out(output);const char* products[]={"purity.json","purity.csv","he4_purity.pdf","he4_fits.pdf","he4_fits.root"};
    for(unsigned i=0;i<5;++i)if(!gSystem->AccessPathName(out+"/"+products[i]))throw std::runtime_error("refusing to overwrite existing purity product");
    TFile* f=TFile::Open(input,"READ");if(!f||f->IsZombie())throw std::runtime_error("cannot open input ROOT");
    TH2* source=dynamic_cast<TH2*>(f->Get(cfg.find("histogram")->second.c_str()));if(!source)throw std::runtime_error("missing preselection TH2");
    TH2* selected=dynamic_cast<TH2*>(f->Get(cfg.find("selectedHistogram")->second.c_str()));
    const double binWidth=source->GetXaxis()->GetBinWidth(1);
    if(std::fabs(step/binWidth-std::floor(step/binWidth+.5))>1.e-6)throw std::runtime_error("rigidity step not aligned to original TH2");
    std::vector<He4Purity::Slice> slices;
    std::ofstream csv((out+"/purity.csv").Data());csv<<std::setprecision(12)<<"rigidity_low,rigidity_high,rigidity_center,gate_count,model,purity,purity_error,model_difference,fit_chi2_ndf,quality,status,covariance,calibration_mean,calibration_sigma,calibration_mean_drift,gate_low_drift,gate_high_drift\n";
    TFile fitRoot(out+"/he4_fits.root","NEW");if(fitRoot.IsZombie())throw std::runtime_error("cannot create fit ROOT");
    gStyle->SetOptStat(0);gStyle->SetOptFit(0);
    for(double lo=rmin;lo<rmax-1.e-8;lo+=step){
      Slice s;s.lo=lo;s.hi=std::min(rmax,lo+step);s.r=(s.lo+s.hi)/2.;
      const int first=source->GetXaxis()->FindBin(s.lo+binWidth*1.e-6),last=source->GetXaxis()->FindBin(s.hi-binWidth*1.e-6);
      if(std::fabs(source->GetXaxis()->GetBinLowEdge(first)-lo)>1.e-6 || std::fabs(source->GetXaxis()->GetBinUpEdge(last)-s.hi)>1.e-6)throw std::runtime_error("unaligned rigidity interval");
      s.mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,s.r);s.sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,s.r);
      s.drift=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,s.hi)-NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,s.lo);
      s.gateLoDrift=s.drift-gate*(NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,s.hi)-NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,s.lo));
      s.gateHiDrift=s.drift+gate*(NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,s.hi)-NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,s.lo));
      TH1D* raw=source->ProjectionY(Form("raw_%03d",int(slices.size())),first,last,"e");raw->SetDirectory(0);
      s.h=new TH1D(Form("projection_%03d",int(slices.size())),"",raw->GetNbinsX(),(raw->GetXaxis()->GetXmin()-s.mean)/s.sigma,(raw->GetXaxis()->GetXmax()-s.mean)/s.sigma);s.h->SetDirectory(0);s.h->Sumw2();
      for(int b=0;b<=raw->GetNbinsX()+1;++b){s.h->SetBinContent(b,raw->GetBinContent(b));s.h->SetBinError(b,raw->GetBinError(b));}
      for(int b=1;b<=s.h->GetNbinsX();++b){const double x=s.h->GetBinCenter(b);if(x>=Required(cfg,"fitNSigmaMin")&&x<Required(cfg,"fitNSigmaMax"))s.counts+=s.h->GetBinContent(b);if(x>=-gate&&x<gate)s.gateCounts+=s.h->GetBinContent(b);}
      if(s.counts>=Required(cfg,"minimumFitCounts") && s.gateCounts>=Required(cfg,"minimumGateCounts"))for(int model=0;model<2;++model){s.fits.push_back(FitOne(s.h,s.r,s.mean,s.sigma,model==1,cfg,int(slices.size())));const int idx=int(s.fits.size())-1;if(s.fits[idx].valid&&(s.chosen<0||s.fits[idx].aic<s.fits[s.chosen].aic))s.chosen=idx;}
      double diff=-1;if(s.fits.size()==2&&s.fits[0].valid&&s.fits[1].valid)diff=std::fabs(s.fits[0].purity-s.fits[1].purity);
      csv<<s.lo<<','<<s.hi<<','<<s.r<<','<<s.gateCounts<<',';
      if(s.chosen>=0){const Fit& best=s.fits[s.chosen];csv<<best.model<<','<<best.purity<<','<<best.error<<','<<diff<<','<<best.chi2/best.ndf<<','<<int(best.quality)<<','<<best.status<<','<<best.cov;}
      else csv<<"unavailable,,,,,0,,";
      csv<<','<<s.mean<<','<<s.sigma<<','<<s.drift<<','<<s.gateLoDrift<<','<<s.gateHiDrift<<'\n';
      fitRoot.cd();raw->Write();s.h->Write();for(unsigned j=0;j<s.fits.size();++j)if(s.fits[j].total){s.fits[j].total->Write();s.fits[j].signal->Write();s.fits[j].background->Write();}delete raw;
      slices.push_back(s);if(slices.size()%20==0)std::cout<<"He4 purity: fitted "<<slices.size()<<" rigidity intervals"<<std::endl;
    }
    csv.close();
    std::ofstream json((out+"/purity.json").Data());json<<std::setprecision(17)<<"{\"schemaVersion\":1,\"inputRoot\":"<<Quote(input)<<",\"histogram\":"<<Quote(cfg.find("histogram")->second)<<",\"sampleLabel\":"<<Quote(label)<<",\"scope\":\"charge-inclusive prePID tracks; inclusive centrality; TPC calibration gate only; no TOF or best-species conditioning\",\"rigidityBinWidth\":"<<binWidth<<",\"dedxBinWidth\":"<<source->GetYaxis()->GetBinWidth(1)<<",\"inputEntries\":"<<source->GetEntries()<<",\"inputRegularCount\":"<<source->Integral()<<",\"selectedTPCEntries\":"<<(selected?selected->GetEntries():-1)<<",\"gateNSigma\":"<<gate<<",\"gatePolicy\":\"continuous calibration nSigma at saved rigidity-bin center; observed counts use dEdx bin centers; no fitted-sigma gate\",\"slices\":[";
    for(unsigned i=0;i<slices.size();++i){const Slice& s=slices[i];if(i)json<<',';json<<"{\"lo\":"<<s.lo<<",\"hi\":"<<s.hi<<",\"r\":"<<s.r<<",\"calibrationMean\":"<<s.mean<<",\"calibrationSigma\":"<<s.sigma<<",\"calibrationMeanDrift\":"<<s.drift<<",\"gateLowDrift\":"<<s.gateLoDrift<<",\"gateHighDrift\":"<<s.gateHiDrift<<",\"fitCounts\":"<<s.counts<<",\"gateCounts\":"<<s.gateCounts<<",\"chosen\":"<<s.chosen<<",\"fits\":[";for(unsigned j=0;j<s.fits.size();++j){if(j)json<<',';WriteFit(json,s.fits[j]);}json<<"]}";}json<<"]}\n";json.close();
    const TString book=out+"/he4_fits.pdf";TCanvas pages("he4_fit_pages","He4 fit audit",1500,1100);pages.Print(book+"[");
    for(unsigned start=0;start<slices.size();start+=4){pages.Clear();pages.Divide(2,2);for(unsigned j=0;j<4&&start+j<slices.size();++j){pages.cd(j+1);DrawSlice(slices[start+j],gate,Required(cfg,"plotLogY")!=0);}pages.Print(book);}pages.Print(book+"]");
    TCanvas plot("he4_purity_plot","He4 purity",1400,1000);plot.SetLeftMargin(.12);plot.SetBottomMargin(.12);plot.SetTopMargin(.16);plot.SetGridy();
    TH1D frame("he4_purity_frame",";R=|p/q| (GeV/c);TPC gate purity S/(S+B)",100,rmin,rmax);frame.SetMinimum(0);frame.SetMaximum(1.08);frame.Draw();
    TGraphErrors good,bad;good.SetName("purity_validated");bad.SetName("purity_poor_fit");good.SetMarkerStyle(20);good.SetMarkerColor(kBlue+1);bad.SetMarkerStyle(24);bad.SetMarkerColor(kRed+1);
    int nGood=0,nBad=0;
    for(unsigned j=0;j<slices.size();++j){const Slice& s=slices[j];if(s.chosen<0)continue;const Fit& fit=s.fits[s.chosen];TGraphErrors* graph=fit.quality?&good:&bad;int& n=fit.quality?nGood:nBad;graph->SetPoint(n,s.r,fit.purity);graph->SetPointError(n,(s.hi-s.lo)/2,fit.error);++n;}
    good.Draw("P SAME");bad.Draw("P SAME");TLegend legend(.14,.18,.6,.29);legend.SetBorderSize(0);legend.AddEntry(&good,"Fit passes status, covariance and residual checks","p");legend.AddEntry(&bad,"Fit converges but model residuals are poor","p");legend.Draw();
    TLatex title;title.SetNDC();title.SetTextSize(.03);title.DrawLatex(.12,.95,label);title.DrawLatex(.12,.905,Form("^{4}He calibration #pm%.1f n#sigma | Gaussian + smooth BG + nearby-nucleus peaks",gate));title.SetTextSize(.026);title.DrawLatex(.12,.06,Form("%d quality fits / %lu bins; charge-inclusive TPC estimate, before best-species / TOF",nGood,static_cast<unsigned long>(slices.size())));
    plot.Print(out+"/he4_purity.pdf");plot.Print(out+"/he4_purity.png");fitRoot.cd();good.Write();bad.Write();plot.Write();fitRoot.Close();f->Close();
    std::cout<<"HE4_PURITY_COMPLETE quality="<<nGood<<" exploratory="<<nBad<<" bins="<<slices.size()<<std::endl;
  } catch(const std::exception& e){std::cerr<<"FitHe4Purity: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
