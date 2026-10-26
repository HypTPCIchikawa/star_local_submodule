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
  int present,valid,quality,ndf,status,cov; double purity,error,chi2,aic;
  std::vector<int> neighborSpecies;
  TF1 *total,*signal,*background;
  Model():present(0),valid(0),quality(0),ndf(0),status(-1),cov(-1),purity(-1),error(0),chi2(0),aic(0),total(0),signal(0),background(0) {}
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
int ComponentModel(const Slice& s,bool showNeighbors) {
  if(s.chosen>=0&&s.chosen<2&&s.fits[s.chosen].valid)return s.chosen;
  if(!showNeighbors)return -1;
  int best=-1;
  for(int j=0;j<2;++j)if(s.fits[j].total&&TMath::Finite(s.fits[j].aic)
      &&(best<0||s.fits[j].aic<s.fits[best].aic))best=j;
  return best;
}
TF1* Neighbor(TF1* source,int component,int species,int slice) {
  const int p=5+3*component;
  if(p+2>=source->GetNpar())throw std::runtime_error("neighbor component parameter index out of range");
  TF1* f=new TF1(Form("neighbor_%d_%d_%s",slice,species,source->GetName()),"[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",source->GetXmin(),source->GetXmax());
  f->SetParameters(source->GetParameter(p),source->GetParameter(p+1),source->GetParameter(p+2));f->SetNpx(600);return f;
}
void CheckComponents(const Model& model,int kind) {
  if(!model.total)return;
  for(int i=0;i<=20;++i) {
    const double x=model.total->GetXmin()+(model.total->GetXmax()-model.total->GetXmin())*i/20.;
    const double total=model.total->Eval(x),sum=model.signal->Eval(x)+model.background->Eval(x);
    if(!TMath::Finite(sum)||std::fabs(sum-total)>1.e-8*std::max(1.,std::fabs(total)))throw std::runtime_error("signal+background component closure failed");
    double smooth=0.;
    if(kind==1)smooth=model.total->GetParameter(3)*std::exp(model.total->GetParameter(4)*x);
    else smooth=model.total->GetParameter(3)+(model.total->GetParameter(4)-model.total->GetParameter(3))*(x-model.total->GetXmin())/(model.total->GetXmax()-model.total->GetXmin());
    for(unsigned j=0;j<model.neighborSpecies.size();++j) {
      const int p=5+3*j;const double sigma=model.total->GetParameter(p+2);
      if(!(sigma>0))throw std::runtime_error("nonpositive neighbor width");
      const double u=(x-model.total->GetParameter(p+1))/sigma;
      smooth+=model.total->GetParameter(p)*std::exp(-.5*u*u);
    }
    if(std::fabs(smooth-model.background->Eval(x))>1.e-8*std::max(1.,std::fabs(total)))throw std::runtime_error("individual neighbor component closure failed");
  }
}
void Draw(const Slice& s,const std::map<std::string,double>& cfg,bool physical,bool logY,const char* label) {
  const bool showNeighbors=Required(cfg,"showNeighborComponents")!=0;
  gPad->SetLeftMargin(.16);gPad->SetRightMargin(.035);gPad->SetBottomMargin(.14);gPad->SetTopMargin(showNeighbors?.32:.27);gPad->SetLogy(logY);
  TH1D* original=physical?s.raw:s.h;
  TH1D* h=static_cast<TH1D*>(original->Clone(Form("display_%d_%d_%d",s.index,int(physical),int(logY))));h->SetDirectory(0);
  const double ulo=Required(cfg,physical?"physicalDisplayNSigmaMin":"fullDisplayNSigmaMin");
  const double uhi=Required(cfg,physical?"physicalDisplayNSigmaMax":"fullDisplayNSigmaMax");
  const double requestedLo=physical?s.mean+s.sigma*ulo:ulo,requestedHi=physical?s.mean+s.sigma*uhi:uhi;
  const double lo=std::max(requestedLo,h->GetXaxis()->GetXmin()),hi=std::min(requestedHi,h->GetXaxis()->GetXmax());
  if(!(hi>lo))throw std::runtime_error("display interval outside stored axis");
  const double ymin=logY?Required(cfg,"logMinimum"):0.,ymax=VisibleMax(h,lo,hi)*Required(cfg,logY?"logMaximumFactor":"linearMaximumFactor");
  h->GetXaxis()->SetRangeUser(lo,hi);h->SetMinimum(ymin);h->SetMaximum(ymax);
  h->SetTitle(Form(";%s;tracks / %.3g keV/cm fit bin",physical?"dE/dx (keV/cm)":"calibrated n#sigma_{^{4}He}",s.raw->GetBinWidth(1)));
  h->GetXaxis()->SetTitleSize(showNeighbors?.036:.04);h->GetYaxis()->SetTitleSize(showNeighbors?.036:.04);h->GetYaxis()->SetTitleOffset(1.65);
  h->GetXaxis()->SetLabelSize(.033);h->GetYaxis()->SetLabelSize(.033);
  h->SetMarkerStyle(20);h->SetMarkerSize(physical?.35:.24);h->SetLineColor(kBlack);h->Draw("E1");
  TLegend* legend=new TLegend(.57,showNeighbors?.705:.765,.955,.905);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(showNeighbors?.021:.023);
  legend->AddEntry(h,"Saved pre-PID tracks","lep");
  for(int j=0;j<2;++j)if(s.fits[j].total) {
    TF1* total=physical?Physical(s.fits[j].total,s,j?"exp_total":"lin_total"):static_cast<TF1*>(s.fits[j].total->Clone(Form("total_view_%d_%d",s.index,j)));
    total->SetLineColor(j?kBlue+1:kRed+1);total->SetLineWidth(2);total->Draw("same");
    legend->AddEntry(total,j?"Exponential BG model":"Linear BG model","l");
  }
  TF1* signal=0;const int componentModel=ComponentModel(s,showNeighbors);
  if(componentModel>=0) {
    const Model& fit=s.fits[componentModel];
    signal=physical?Physical(fit.signal,s,"chosen_signal"):static_cast<TF1*>(fit.signal->Clone(Form("signal_view_%d",s.index)));
    TF1* bg=physical?Physical(fit.background,s,"chosen_bg"):static_cast<TF1*>(fit.background->Clone(Form("bg_view_%d",s.index)));
    signal->SetLineColor(kGreen+2);signal->SetLineWidth(2);signal->Draw("same");
    bg->SetLineColor(kMagenta+1);bg->SetLineStyle(2);bg->SetLineWidth(2);bg->Draw("same");
    legend->AddEntry(signal,componentModel==s.chosen?"Selected ^{4}He Gaussian":"^{4}He: unaccepted model","l");legend->AddEntry(bg,"Summed background","l");
    if(showNeighbors)for(unsigned j=0;j<fit.neighborSpecies.size();++j) {
      const int species=fit.neighborSpecies[j];
      TF1* rawNeighbor=Neighbor(fit.total,j,species,s.index);
      TF1* neighbor=physical?Physical(rawNeighbor,s,Form("neighbor_%d",species)):rawNeighbor;
      const char* names[]={"D Gaussian (tail)","T (triton) Gaussian","^{3}He Gaussian"};
      const int colors[]={kViolet+2,kOrange+7,kCyan+2};
      neighbor->SetLineColor(colors[species]);neighbor->SetLineStyle(species==1?4:7);neighbor->SetLineWidth(2);neighbor->Draw("same");
      legend->AddEntry(neighbor,names[species],"l");
    }
    const double center=physical?s.mean+s.sigma*fit.signal->GetParameter(1):fit.signal->GetParameter(1);
    if(center>=lo&&center<=hi){TLine* c=new TLine(center,ymin,center,ymax);c->SetLineColor(kGreen+2);c->SetLineStyle(2);c->Draw();}
  }
  for(int sign=-1;sign<=1;sign+=2) {
    const double gate=sign*Required(cfg,"gateNSigma"),x=physical?s.mean+s.sigma*gate:gate;
    if(x>=lo&&x<=hi){TLine* line=new TLine(x,ymin,x,ymax);line->SetLineColor(kGray+2);line->SetLineStyle(3);line->Draw();}
  }
  legend->Draw();TLatex text;text.SetNDC();text.SetTextSize(.035);
  text.DrawLatex(.16,.935,Form("R=%.2f-%.2f GeV/c [%s]",s.lo,s.hi,logY?"log":"linear"));
  text.SetTextSize(showNeighbors?.025:.028);
  if(s.chosen>=0&&s.chosen<2) {
    const Model& fit=s.fits[s.chosen];
    text.DrawLatex(.16,.865,Form("P=%.3f #pm %.3f%s",fit.purity,fit.error,fit.quality?"":" [poor fit]"));
    text.DrawLatex(.16,.815,Form("#chi^{2}/ndf=%.1f",fit.ndf>0?fit.chi2/fit.ndf:-1));
  } else {
    text.SetTextColor(kRed+1);text.DrawLatex(.16,.865,"No accepted fit");text.SetTextColor(kBlack);
    text.DrawLatex(.16,.815,"Does not imply no peak");
    if(componentModel>=0){text.SetTextSize(.02);text.DrawLatex(.16,.765,"Components from unaccepted model");}
  }
  if(physical||showNeighbors) {
    const bool sameCoverage=Required(cfg,"physicalDisplayNSigmaMin")==Required(cfg,"fullDisplayNSigmaMin")&&Required(cfg,"physicalDisplayNSigmaMax")==Required(cfg,"fullDisplayNSigmaMax");
    text.SetTextSize(.019);text.DrawLatex(.16,.975,physical&&sameCoverage?"Physical dE/dx: same full coverage as upper panels (no refit)":label);
    text.SetTextSize(.015);
    const TF1* support=componentModel>=0?s.fits[componentModel].total:(s.fits[0].total?s.fits[0].total:s.fits[1].total);
    if(support) {
      const double fitLo=physical?s.mean+s.sigma*support->GetXmin():support->GetXmin();
      const double fitHi=physical?s.mean+s.sigma*support->GetXmax():support->GetXmax();
      text.DrawLatex(.16,.056,Form("Stored fit range: %.2f to %.2f %s (curves only within fitted support)",fitLo,fitHi,physical?"keV/cm":"calibrated n#sigma"));
    } else text.DrawLatex(.16,.056,"No stored fit functions for this slice");
    text.DrawLatex(.16,.025,Form("Gray: original calibration gate; dE/dx fit bin = %.3g keV/cm%s",s.raw->GetBinWidth(1),(requestedHi>h->GetXaxis()->GetXmax()||requestedLo<h->GetXaxis()->GetXmin())?"; display clipped to stored axis":""));
  }
}
}
#endif

void RenderHe4FitViews(const char* input,const char* output,const char* settings,const char* manifest,const char* label,const char* diagnosticPdf) {
#ifndef __CINT__
  using namespace He4Views;
  try {
    const std::map<std::string,double> cfg=Settings(settings);
    const TString out(output),linear=out+"/he4_fits_linear.pdf",physical=out+"/"+diagnosticPdf;
    const bool showNeighbors=Required(cfg,"showNeighborComponents")!=0;
    const int diagnosticPanels=int(Required(cfg,"diagnosticPanelsPerPage"));
    if(diagnosticPanels!=2&&diagnosticPanels!=4)throw std::runtime_error("invalid explicit panel count");
    if(!gSystem->AccessPathName(linear)||!gSystem->AccessPathName(physical))throw std::runtime_error("refusing to overwrite view PDFs");
    TFile* source=TFile::Open(input,"READ");if(!source||source->IsZombie())throw std::runtime_error("cannot read stored fits");
    std::ifstream rows(manifest);if(!rows)throw std::runtime_error("cannot read explicit slice manifest");
    std::vector<He4Views::Slice> slices;std::string line;
    while(std::getline(rows,line)) {
      Slice s;std::istringstream row(line);row>>s.index>>s.lo>>s.hi>>s.mean>>s.sigma>>s.chosen>>s.diagnostic;
      for(int j=0;j<2;++j) {
        Model& model=s.fits[j];row>>model.present>>model.valid>>model.quality>>model.purity>>model.error>>model.chi2>>model.ndf>>model.status>>model.cov>>model.aic;
        int neighbors=0;row>>neighbors;if(neighbors<0||neighbors>3)throw std::runtime_error("invalid neighbor count");
        for(int n=0;n<neighbors;++n){int species=-1;row>>species;if(species<0||species>2)throw std::runtime_error("invalid neighbor species");model.neighborSpecies.push_back(species);}
        if(model.present) {
          const TString name=Form("he4_%s_%d",j?"exponential":"linear",s.index);
          model.total=dynamic_cast<TF1*>(source->Get(name));model.signal=dynamic_cast<TF1*>(source->Get(name+"_signal"));model.background=dynamic_cast<TF1*>(source->Get(name+"_bg"));
          if(!model.total||!model.signal||!model.background)throw std::runtime_error("missing stored fitted functions");
          if(showNeighbors){if(model.total->GetNpar()!=5+3*neighbors)throw std::runtime_error("neighbor metadata disagrees with stored function");CheckComponents(model,j);}
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
    TCanvas detail("he4_physical_views","Physical dEdx log and linear",1600,diagnosticPanels==4?1400:800);detail.Print(physical+"[");
    std::vector<int> indices;
    for(unsigned i=0;i<slices.size();++i)if(slices[i].diagnostic) {
      detail.Clear();detail.Divide(2,diagnosticPanels==4?2:1);
      if(diagnosticPanels==4){detail.cd(1);Draw(slices[i],cfg,false,true,label);detail.cd(2);Draw(slices[i],cfg,false,false,label);}
      const int offset=diagnosticPanels==4?2:0;
      detail.cd(offset+1);Draw(slices[i],cfg,true,true,label);detail.cd(offset+2);Draw(slices[i],cfg,true,false,label);
      detail.Print(physical);detail.Print(out+Form("/he4_dedx_R_%03d.png",int(1000*(slices[i].lo+slices[i].hi)/2.+.5)));indices.push_back(i);
    }
    detail.Print(physical+"]");
    std::ofstream audit((out+"/render_audit.json").Data());audit<<"{\"no_refit\":true,\"slices\":"<<slices.size()<<",\"diagnostic_indices\":[";
    for(unsigned i=0;i<indices.size();++i){if(i)audit<<',';audit<<indices[i];}
    audit<<"],\"diagnostic_panels_per_page\":"<<diagnosticPanels<<",\"individual_neighbor_components\":"<<(showNeighbors?"true":"false")<<",\"component_sum_equivalence_checked\":"<<(showNeighbors?"true":"false")<<",\"diagnostic_components\":[";
    for(unsigned n=0;n<indices.size();++n){
      if(n)audit<<',';const Slice& s=slices[indices[n]];const int model=ComponentModel(s,showNeighbors);
      const double requestedLo=s.mean+s.sigma*Required(cfg,"physicalDisplayNSigmaMin"),requestedHi=s.mean+s.sigma*Required(cfg,"physicalDisplayNSigmaMax");
      const double displayLo=std::max(requestedLo,s.raw->GetXaxis()->GetXmin()),displayHi=std::min(requestedHi,s.raw->GetXaxis()->GetXmax());
      const double fullLo=std::max(Required(cfg,"fullDisplayNSigmaMin"),s.h->GetXaxis()->GetXmin()),fullHi=std::min(Required(cfg,"fullDisplayNSigmaMax"),s.h->GetXaxis()->GetXmax());
      const double fullDedxLo=s.mean+s.sigma*fullLo,fullDedxHi=s.mean+s.sigma*fullHi;
      const double tolerance=1.e-10*std::max(1.,std::max(std::fabs(displayLo),std::fabs(displayHi)));
      const bool sameCoverage=std::fabs(fullDedxLo-displayLo)<tolerance&&std::fabs(fullDedxHi-displayHi)<tolerance;
      TAxis fullAxis(*s.h->GetXaxis()),physicalAxis(*s.raw->GetXaxis());fullAxis.SetRangeUser(fullLo,fullHi);physicalAxis.SetRangeUser(displayLo,displayHi);
      const bool sameBins=fullAxis.GetFirst()==physicalAxis.GetFirst()&&fullAxis.GetLast()==physicalAxis.GetLast();
      const double upperCounts=s.h->Integral(fullAxis.GetFirst(),fullAxis.GetLast()),lowerCounts=s.raw->Integral(physicalAxis.GetFirst(),physicalAxis.GetLast());
      audit<<std::setprecision(17)<<"{\"index\":"<<s.index<<",\"model\":"<<model<<",\"accepted\":"<<(model>=0&&model==s.chosen?"true":"false")<<",\"neighborSpecies\":[";
      if(model>=0)for(unsigned j=0;j<s.fits[model].neighborSpecies.size();++j){if(j)audit<<',';audit<<s.fits[model].neighborSpecies[j];}
      audit<<"],\"stored_dedx_range\":["<<s.raw->GetXaxis()->GetXmin()<<','<<s.raw->GetXaxis()->GetXmax()<<"],\"physical_requested_dedx\":["<<requestedLo<<','<<requestedHi<<"],\"physical_display_dedx\":["<<displayLo<<','<<displayHi<<"],\"physical_display_clipped\":"<<((requestedLo<displayLo||requestedHi>displayHi)?"true":"false")
           <<",\"full_display_dedx\":["<<fullDedxLo<<','<<fullDedxHi<<"],\"same_physical_coverage\":"<<(sameCoverage?"true":"false")<<",\"same_visible_bins\":"<<(sameBins?"true":"false")<<",\"same_visible_counts\":"<<(upperCounts==lowerCounts?"true":"false")<<",\"upper_visible_counts\":"<<upperCounts<<",\"lower_visible_counts\":"<<lowerCounts<<",\"upper_visible_bins\":["<<fullAxis.GetFirst()<<','<<fullAxis.GetLast()<<"],\"lower_visible_bins\":["<<physicalAxis.GetFirst()<<','<<physicalAxis.GetLast()<<"],\"neighbor_centers\":[";
      if(model>=0)for(unsigned j=0;j<s.fits[model].neighborSpecies.size();++j){
        if(j)audit<<',';const double u=s.fits[model].total->GetParameter(6+3*j),dedx=s.mean+s.sigma*u;
        audit<<"{\"species\":"<<s.fits[model].neighborSpecies[j]<<",\"nsigma\":"<<u<<",\"dedx\":"<<dedx<<",\"visible_physical\":"<<((dedx>=displayLo&&dedx<=displayHi)?"true":"false")<<"}";
      }
      audit<<"]}";
    }
    audit<<"],\"fit_dedx_bin_widths\":[";
    for(unsigned i=0;i<slices.size();++i){if(i)audit<<',';audit<<slices[i].raw->GetBinWidth(1);}
    audit<<"],\"physical_curve_coordinate_equivalence_checked\":true,\"counts_errors_unchanged\":true,\"physical_y_convention\":\"counts_per_saved_fit_bin_no_jacobian_rescaling\"}\n";audit.close();
    source->Close();delete source;std::cout<<"HE4_FIT_VIEWS_COMPLETE"<<std::endl;gSystem->Exit(0);
  } catch(const std::exception& error) {std::cerr<<"HE4_FIT_VIEWS_FAILED: "<<error.what()<<std::endl;gSystem->Exit(1);}
#endif
}
