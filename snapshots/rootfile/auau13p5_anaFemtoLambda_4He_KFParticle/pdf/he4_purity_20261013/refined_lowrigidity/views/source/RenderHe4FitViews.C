// Postprocess stored fits only. Open source ROOT READ; never refit or change cuts.
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TString.h"
#include "TMath.h"
#include <algorithm>
#include <cctype>
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
namespace He4Views {
struct Model {
  int present,valid,quality,ndf,status,cov; double purity,error,chi2;
  TF1 *total,*signal,*background;
  Model():present(0),valid(0),quality(0),ndf(0),status(-1),cov(-1),purity(-1),error(0),chi2(0),total(0),signal(0),background(0) {}
};
struct Slice {
  int index,chosen,diagnostic; double lo,hi,mean,sigma; TH1D *h,*raw; Model fits[2];
  Slice():index(-1),chosen(-1),diagnostic(0),lo(0),hi(0),mean(0),sigma(0),h(0),raw(0) {}
};
double Required(const std::map<std::string,double>& c,const char* key) {
  std::map<std::string,double>::const_iterator i=c.find(key);
  if(i==c.end()||!TMath::Finite(i->second)) throw std::runtime_error(std::string("missing setting: ")+key);
  return i->second;
}
std::map<std::string,double> Settings(const char* file) {
  std::ifstream f(file); if(!f)throw std::runtime_error("cannot open explicit view settings");
  std::map<std::string,double> out;std::string key;double value;
  while(f>>key>>value)out[key]=value;
  if(!f.eof())throw std::runtime_error("invalid view settings");return out;
}
std::string PhysicalFormula(TF1* source,double mean,double sigma) {
  std::ostringstream axis; axis<<std::setprecision(17)<<"((x-("<<mean<<"))/("<<sigma<<"))";
  // Keep parameter tokens and copy their exact doubles below. ROOT5's "p"
  // expansion rounds parameter values in its textual output.
  const std::string original=source->GetExpFormula().Data();
  if(original.empty())throw std::runtime_error("stored TF1 is not a self-contained formula");
  std::string result;
  // Replace only the variable token x, preserving identifiers such as exp.
  for(size_t i=0;i<original.size();) {
    if(std::isalpha(static_cast<unsigned char>(original[i]))||original[i]=='_') {
      size_t j=i+1;while(j<original.size()&&(std::isalnum(static_cast<unsigned char>(original[j]))||original[j]=='_'))++j;
      const std::string token=original.substr(i,j-i);result+=token=="x"?axis.str():token;i=j;
    } else result+=original[i++];
  }
  return result;
}
TF1* Physical(TF1* source,const Slice& s,const char* suffix) {
  if(!source)return 0;
  const std::string formula=PhysicalFormula(source,s.mean,s.sigma);
  TF1* f=new TF1(Form("physical_%s_%d",suffix,s.index),formula.c_str(),
      s.mean+s.sigma*source->GetXmin(),s.mean+s.sigma*source->GetXmax());f->SetNpx(600);
  for(int p=0;p<f->GetNpar();++p)f->SetParameter(p,source->GetParameter(p));
  for(int i=0;i<=20;++i) {
    const double u=source->GetXmin()+(source->GetXmax()-source->GetXmin())*i/20.;
    const double old=source->Eval(u),now=f->Eval(s.mean+s.sigma*u);
    if(!TMath::Finite(now)||std::fabs(now-old)>1.e-8*std::max(1.,std::fabs(old)))
      throw std::runtime_error("physical-axis TF1 transformation failed equivalence check");
  }
  return f;
}
double VisibleMax(TH1D* h,double lo,double hi) {
  double maximum=1.;for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=lo&&h->GetBinCenter(b)<=hi)
    maximum=std::max(maximum,h->GetBinContent(b)+h->GetBinError(b));return maximum;
}
void Draw(const Slice& s,const std::map<std::string,double>& cfg,bool physical,bool logY,const char* label) {
  gPad->SetLeftMargin(.16);gPad->SetRightMargin(.035);gPad->SetBottomMargin(.14);gPad->SetTopMargin(.27);gPad->SetLogy(logY);
  TH1D* original=physical?s.raw:s.h;
  TH1D* h=static_cast<TH1D*>(original->Clone(Form("display_%d_%d_%d",s.index,int(physical),int(logY))));h->SetDirectory(0);
  const double ulo=Required(cfg,physical?"physicalDisplayNSigmaMin":"fullDisplayNSigmaMin");
  const double uhi=Required(cfg,physical?"physicalDisplayNSigmaMax":"fullDisplayNSigmaMax");
  const double lo=physical?s.mean+s.sigma*ulo:ulo,hi=physical?s.mean+s.sigma*uhi:uhi;
  const double ymin=logY?Required(cfg,"logMinimum"):0.,ymax=VisibleMax(h,lo,hi)*Required(cfg,logY?"logMaximumFactor":"linearMaximumFactor");
  h->GetXaxis()->SetRangeUser(lo,hi);h->SetMinimum(ymin);h->SetMaximum(ymax);
  h->SetTitle(Form(";%s;tracks / %.3g keV/cm fit bin",physical?"dE/dx (keV/cm)":"calibrated n#sigma_{^{4}He}",s.raw->GetBinWidth(1)));
  h->GetXaxis()->SetTitleSize(.04);h->GetYaxis()->SetTitleSize(.04);h->GetYaxis()->SetTitleOffset(1.65);
  h->GetXaxis()->SetLabelSize(.033);h->GetYaxis()->SetLabelSize(.033);
  h->SetMarkerStyle(20);h->SetMarkerSize(physical?.35:.24);h->SetLineColor(kBlack);h->Draw("E1");
  TLegend* legend=new TLegend(.57,.765,.955,.905);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.023);
  legend->AddEntry(h,"Saved pre-PID tracks","lep");
  for(int j=0;j<2;++j)if(s.fits[j].total) {
    TF1* total=physical?Physical(s.fits[j].total,s,j?"exp_total":"lin_total"):static_cast<TF1*>(s.fits[j].total->Clone(Form("total_view_%d_%d",s.index,j)));
    total->SetLineColor(j?kBlue+1:kRed+1);total->SetLineWidth(2);total->Draw("same");
    legend->AddEntry(total,j?"Exponential BG model":"Linear BG model","l");
  }
  TF1* signal=0;
  if(s.chosen>=0&&s.chosen<2&&s.fits[s.chosen].valid) {
    const Model& fit=s.fits[s.chosen];
    signal=physical?Physical(fit.signal,s,"chosen_signal"):static_cast<TF1*>(fit.signal->Clone(Form("signal_view_%d",s.index)));
    TF1* bg=physical?Physical(fit.background,s,"chosen_bg"):static_cast<TF1*>(fit.background->Clone(Form("bg_view_%d",s.index)));
    signal->SetLineColor(kGreen+2);signal->SetLineWidth(2);signal->Draw("same");
    bg->SetLineColor(kMagenta+1);bg->SetLineStyle(2);bg->SetLineWidth(2);bg->Draw("same");
    legend->AddEntry(signal,"Selected ^{4}He Gaussian","l");legend->AddEntry(bg,"Selected background","l");
    const double center=physical?s.mean+s.sigma*fit.signal->GetParameter(1):fit.signal->GetParameter(1);
    TLine* c=new TLine(center,ymin,center,ymax);c->SetLineColor(kGreen+2);c->SetLineStyle(2);c->Draw();
  }
  for(int sign=-1;sign<=1;sign+=2) {
    const double gate=sign*Required(cfg,"gateNSigma"),x=physical?s.mean+s.sigma*gate:gate;
    TLine* line=new TLine(x,ymin,x,ymax);line->SetLineColor(kGray+2);line->SetLineStyle(3);line->Draw();
  }
  legend->Draw();TLatex text;text.SetNDC();text.SetTextSize(.035);
  text.DrawLatex(.16,.935,Form("R=%.2f-%.2f GeV/c [%s]",s.lo,s.hi,logY?"log":"linear"));
  text.SetTextSize(.028);
  if(s.chosen>=0&&s.chosen<2) {
    const Model& fit=s.fits[s.chosen];
    text.DrawLatex(.16,.865,Form("P=%.3f #pm %.3f%s",fit.purity,fit.error,fit.quality?"":" [poor fit]"));
    text.DrawLatex(.16,.815,Form("#chi^{2}/ndf=%.1f",fit.ndf>0?fit.chi2/fit.ndf:-1));
  } else {
    text.SetTextColor(kRed+1);text.DrawLatex(.16,.865,"No accepted fit");text.SetTextColor(kBlack);
    text.DrawLatex(.16,.815,"Does not imply no peak");
  }
  if(physical) {
    text.SetTextSize(.019);text.DrawLatex(.16,.975,label);
    text.DrawLatex(.16,.027,Form("Gray: original calibration gate; dE/dx fit bin = %.3g keV/cm",s.raw->GetBinWidth(1)));
  }
}
}
#endif

void RenderHe4FitViews(const char* input,const char* output,const char* settings,const char* manifest,const char* label) {
#ifndef __CINT__
  using namespace He4Views;
  try {
    const std::map<std::string,double> cfg=Settings(settings);
    const TString out(output),linear=out+"/he4_fits_linear.pdf",physical=out+"/he4_low_rigidity_dedx.pdf";
    if(!gSystem->AccessPathName(linear)||!gSystem->AccessPathName(physical))throw std::runtime_error("refusing to overwrite view PDFs");
    TFile* source=TFile::Open(input,"READ");if(!source||source->IsZombie())throw std::runtime_error("cannot read stored fits");
    std::ifstream rows(manifest);if(!rows)throw std::runtime_error("cannot read explicit slice manifest");
    std::vector<He4Views::Slice> slices;std::string line;
    while(std::getline(rows,line)) {
      Slice s;std::istringstream row(line);row>>s.index>>s.lo>>s.hi>>s.mean>>s.sigma>>s.chosen>>s.diagnostic;
      for(int j=0;j<2;++j) {
        Model& model=s.fits[j];row>>model.present>>model.valid>>model.quality>>model.purity>>model.error>>model.chi2>>model.ndf>>model.status>>model.cov;
        if(model.present) {
          const TString name=Form("he4_%s_%d",j?"exponential":"linear",s.index);
          model.total=dynamic_cast<TF1*>(source->Get(name));model.signal=dynamic_cast<TF1*>(source->Get(name+"_signal"));model.background=dynamic_cast<TF1*>(source->Get(name+"_bg"));
          if(!model.total||!model.signal||!model.background)throw std::runtime_error("missing stored fitted functions");
        }
      }
      if(!row||s.index!=int(slices.size())||!(s.sigma>0))throw std::runtime_error("invalid slice metadata");
      s.h=dynamic_cast<TH1D*>(source->Get(Form("projection_%03d",s.index)));
      s.raw=dynamic_cast<TH1D*>(source->Get(Form("raw_%03d",s.index)));
      if(!s.h||!s.raw||s.h->GetNbinsX()!=s.raw->GetNbinsX())throw std::runtime_error("missing or inconsistent saved projection");
      for(int b=0;b<=s.h->GetNbinsX()+1;++b)if(s.h->GetBinContent(b)!=s.raw->GetBinContent(b)||s.h->GetBinError(b)!=s.raw->GetBinError(b))throw std::runtime_error("coordinate transform changed counts/errors");
      slices.push_back(s);
    }
    if(slices.empty())throw std::runtime_error("empty slice manifest");
    gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetPaperSize(30,20);
    TCanvas full("he4_stored_linear","Stored He4 fits linear",1500,1000);full.Print(linear+"[");
    for(unsigned start=0;start<slices.size();start+=4) {
      full.Clear();full.Divide(2,2);
      for(unsigned j=0;j<4&&start+j<slices.size();++j){full.cd(j+1);Draw(slices[start+j],cfg,false,false,label);}
      full.Print(linear);if(start==0)full.Print(out+"/he4_fits_linear_first.png");
    }
    full.Print(linear+"]");
    TCanvas detail("he4_physical_views","Physical dEdx log and linear",1600,800);detail.Print(physical+"[");
    std::vector<int> indices;
    for(unsigned i=0;i<slices.size();++i)if(slices[i].diagnostic) {
      detail.Clear();detail.Divide(2,1);detail.cd(1);Draw(slices[i],cfg,true,true,label);detail.cd(2);Draw(slices[i],cfg,true,false,label);
      detail.Print(physical);detail.Print(out+Form("/he4_dedx_R_%03d.png",int(1000*(slices[i].lo+slices[i].hi)/2.+.5)));indices.push_back(i);
    }
    detail.Print(physical+"]");
    std::ofstream audit((out+"/render_audit.json").Data());audit<<"{\"no_refit\":true,\"slices\":"<<slices.size()<<",\"diagnostic_indices\":[";
    for(unsigned i=0;i<indices.size();++i){if(i)audit<<',';audit<<indices[i];}
    audit<<"],\"fit_dedx_bin_widths\":[";
    for(unsigned i=0;i<slices.size();++i){if(i)audit<<',';audit<<slices[i].raw->GetBinWidth(1);}
    audit<<"],\"physical_curve_coordinate_equivalence_checked\":true,\"counts_errors_unchanged\":true}\n";audit.close();
    source->Close();delete source;std::cout<<"HE4_FIT_VIEWS_COMPLETE"<<std::endl;gSystem->Exit(0);
  } catch(const std::exception& error) {std::cerr<<"HE4_FIT_VIEWS_FAILED: "<<error.what()<<std::endl;gSystem->Exit(1);}
#endif
}
