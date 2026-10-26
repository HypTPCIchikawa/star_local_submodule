// ROOT 5 QA for saved, unweighted Lambda purity-study CSVs.
// No PicoDst event loop, selection, smoothing, rebinning, or mass constraints.
// Usage (ACLiC): drawLambdaPurityStudy(scanDirectory, outputPrefix).
// Required inputs: histograms.csv, chosen_selections.csv, chosen_metrics.csv,
// plot_manifest.tsv. Manifest keys are documented in ReadManifest below.
// Fits are independent Gaussian+pol2 binned-Poisson cross-checks, not rankings.
// Failed fits retain diagnostics but their physics estimates are blank in CSV.

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TStyle.h>
#include <TObjString.h>
#include <TMath.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

Int_t drawLambdaPurityStudy(const char*, const char*);
#ifndef __CINT__
namespace LambdaPurityQA {
typedef std::map<std::string,std::string> Row;
struct Range { double lo,hi; Range(double a=0,double b=0):lo(a),hi(b){} };
struct Settings {
  double lo,hi,bw,displayLo,displayHi,windowLo,windowHi;
  double sigmaLo,sigmaHi,meanLo,meanHi,minEntries;
  std::string title,rawManifest;
  std::vector<Range> fits;
};
struct Choice { std::string id,label,role; };
struct Metric { double n,s,b,purity,purityError,retention; };
struct Fit {
  TF1 *model,*background;
  int status,covariance,ndf;
  double s,b,purity,purityError,sError,bError,mean,sigma,deviance,devianceP;
  std::string flags;
  Fit():model(0),background(0),status(-999),covariance(-1),ndf(0),
    s(0),b(0),purity(0),purityError(0),sError(0),bError(0),mean(0),sigma(0),deviance(0),devianceP(0){}
  bool valid()const{return flags.empty();}
};
static void Require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
static std::string Join(const std::string&a,const std::string&b){return a+"/"+b;}
static std::string Key(const std::string&a,const std::string&b){return a+"/"+b;}
static std::string CsvQuote(const std::string&s){std::string q="\"";for(size_t i=0;i<s.size();++i){if(s[i]=='\"')q+='\"';q+=s[i];}return q+'\"';}
static std::vector<std::string> Csv(const std::string&line){
  std::vector<std::string> fields; std::string field; bool quoted=false;
  for(size_t i=0;i<line.size();++i){char c=line[i];
    if(c=='\"'){if(quoted && i+1<line.size() && line[i+1]=='\"'){field+='\"';++i;}else quoted=!quoted;}
    else if(c==','&&!quoted){fields.push_back(field);field.clear();}
    else if(c!='\r')field+=c;
  }
  Require(!quoted,"unterminated CSV quoted field");fields.push_back(field);return fields;
}
static std::vector<Row> ReadCsv(const std::string&path){
  std::ifstream in(path.c_str());Require(bool(in),"cannot open "+path);
  std::string line;Require(bool(std::getline(in,line)),"empty CSV "+path);
  std::vector<std::string> header=Csv(line);std::vector<Row> rows;
  while(std::getline(in,line)){if(line.empty())continue;std::vector<std::string> values=Csv(line);
    Require(values.size()==header.size(),"wrong CSV field count in "+path);
    Row row;for(size_t i=0;i<header.size();++i)row[header[i]]=values[i];rows.push_back(row);
  }Require(in.eof(),"read failure "+path);return rows;
}
static std::string Get(const Row&r,const std::string&key){Row::const_iterator p=r.find(key);Require(p!=r.end(),"missing field "+key);return p->second;}
static double Number(const std::string&s){char*end=0;double v=std::strtod(s.c_str(),&end);Require(end&&end!=s.c_str()&&!*end&&TMath::Finite(v),"invalid finite number: "+s);return v;}
static double N(const Row&r,const std::string&key){return Number(Get(r,key));}
static Settings ReadManifest(const std::string&path){
  std::ifstream in(path.c_str());Require(bool(in),"cannot open manifest");
  Row values;Settings c;std::string line;
  while(std::getline(in,line)){c.rawManifest+=line+"\n";if(line.empty()||line[0]=='#')continue;
    size_t tab=line.find('\t');Require(tab!=std::string::npos,"manifest must be key TAB value");
    std::string value=line.substr(tab+1);if(!value.empty()&&value[value.size()-1]=='\r')value.erase(value.size()-1);
    Require(values.insert(std::make_pair(line.substr(0,tab),value)).second,"duplicate manifest key");
  }
  c.lo=N(values,"massMin");c.hi=N(values,"massMax");c.bw=N(values,"massBinWidth");
  c.displayLo=N(values,"displayMin");c.displayHi=N(values,"displayMax");
  c.windowLo=N(values,"windowLow");c.windowHi=N(values,"windowHigh");
  c.sigmaLo=N(values,"fitSigmaMin");c.sigmaHi=N(values,"fitSigmaMax");
  c.meanLo=N(values,"fitMeanMin");c.meanHi=N(values,"fitMeanMax");
  c.minEntries=N(values,"fitMinEntries");c.title=Get(values,"title");
  Require(c.hi>c.lo&&c.bw>0&&c.displayLo>=c.lo&&c.displayHi<=c.hi&&c.displayLo<c.displayHi,
          "invalid mass/display range");
  Require(c.windowLo>c.lo&&c.windowHi<c.hi&&c.windowHi>c.windowLo&&c.meanLo>=c.windowLo&&
    c.meanHi<=c.windowHi&&c.meanHi>c.meanLo&&c.sigmaLo>0&&c.sigmaHi>c.sigmaLo&&c.minEntries>=0,"invalid fit/window bounds");
  std::stringstream ranges(Get(values,"fitRanges"));std::string token;
  while(std::getline(ranges,token,';')){size_t colon=token.find(':');Require(colon!=std::string::npos,"invalid fitRanges");
    Range r(Number(token.substr(0,colon)),Number(token.substr(colon+1)));
    Require(r.lo>=c.lo&&r.hi<=c.hi&&r.lo<c.windowLo&&r.hi>c.windowHi,"fit range must include signal and sidebands");c.fits.push_back(r);
  }Require(!c.fits.empty(),"no fit ranges");return c;
}
static int Edge(const TH1D*h,double x){
  for(int i=1;i<=h->GetNbinsX()+1;++i)if(std::fabs(h->GetBinLowEdge(i)-x)<1.e-9)return i;
  throw std::runtime_error("range not aligned with histogram bin edges");
}
static int Color(size_t i){const char*hex[]={"#30343b","#0072b2","#d55e00","#009e73","#a64d9e","#c78d00"};return TColor::GetColor(hex[i%6]);}
static int Marker(size_t i){const int m[]={20,24,25,26,27,28};return m[i%6];}
static void HistStyle(TH1D*h,size_t i){h->SetStats(false);h->SetTitle("");h->SetLineColor(Color(i));h->SetMarkerColor(Color(i));h->SetLineWidth(i?2:1);h->SetLineStyle(i==1?7:(i==3?9:1));h->SetMarkerStyle(Marker(i));h->SetMarkerSize(i?.9:.65);}
static void Axis(TH1D*h,const char*y){
  h->GetXaxis()->SetTitle("M_{p#pi^{-}} [GeV/c^{2}]");h->GetYaxis()->SetTitle(y);
  h->GetXaxis()->SetTitleSize(.05);h->GetYaxis()->SetTitleSize(.047);
  h->GetXaxis()->SetLabelSize(.041);h->GetYaxis()->SetLabelSize(.039);
  h->GetXaxis()->SetTitleOffset(1.2);h->GetYaxis()->SetTitleOffset(1.55);h->GetXaxis()->SetNdivisions(504,false);
}
static void Save(TCanvas*c,const std::string&stem,TFile&out){
  c->Modified();c->Update();c->Print((stem+".pdf").c_str());c->Print((stem+".png").c_str());
  out.cd();Require(c->Write()>0,"canvas write failure");
  const char*ext[]={".pdf",".png"};for(int j=0;j<2;++j){FileStat_t st;Require(gSystem->GetPathInfo((stem+ext[j]).c_str(),st)==0&&st.fSize>0,"missing plot product");}
}
static void Overlay(const std::vector<Choice>&choices,const std::vector<size_t>&which,
  const std::map<std::string,TH1D*>&hist,const std::map<std::string,Metric>&metrics,
  const Settings&c,const std::string&stem,int index,TFile&out){
  TCanvas*canvas=new TCanvas(Form("PurityOverlay%d",index),"KF Lambda purity comparison",1800,1060);
  std::vector<TH1D*> raw;
  for(int panel=0;panel<2;++panel){
    canvas->cd();TPad*pad=new TPad(Form("PurityOverlay%dPad%d",index,panel),"",panel*.5,.06,(panel+1)*.5,.72);
    pad->SetLeftMargin(.15);pad->SetRightMargin(.04);pad->SetBottomMargin(.15);pad->SetTopMargin(.09);pad->SetTicks(1,1);pad->SetGridy();pad->Draw();pad->cd();
    std::vector<TH1D*> drawn;double max=0;
    for(size_t k=0;k<which.size();++k){size_t j=which[k];TH1D*source=hist.find(Key(choices[j].id,"full"))->second;
      TH1D*h=(TH1D*)source->Clone(Form("Overlay%dP%dC%u",index,panel,(unsigned)j));h->SetDirectory(0);
      double norm=h->Integral(1,h->GetNbinsX());Require(norm>0,"empty normalization histogram");
      if(panel)h->Scale(1./norm);HistStyle(h,j);h->GetXaxis()->SetRange(Edge(h,c.displayLo),Edge(h,c.displayHi)-1);
      max=TMath::Max(max,h->GetMaximum());drawn.push_back(h);if(!panel)raw.push_back(h);
      for(int b=1;b<=h->GetNbinsX();++b)Require(std::fabs(h->GetBinContent(b)-source->GetBinContent(b)/(panel?norm:1.))<1.e-10,"normalization audit failed");
      if(panel)Require(std::fabs(h->Integral(1,h->GetNbinsX())-1.)<1.e-10,"normalized full range not unity");
      out.cd();h->Write();pad->cd();
    }
    TH1D*frame=(TH1D*)drawn[0]->Clone(Form("Overlay%dFrame%d",index,panel));frame->SetDirectory(0);frame->Reset("ICES");frame->SetMinimum(0);frame->SetMaximum(max*1.16);
    Axis(frame,Form(panel?"Fraction per %.3g MeV/c^{2} bin":"Candidates per %.3g MeV/c^{2} bin",1000*c.bw));frame->Draw("AXIS");
    for(size_t k=drawn.size();k>0;--k)drawn[k-1]->Draw("HIST SAME");
    for(size_t k=0;k<drawn.size();++k){TGraph*marks=new TGraph;for(int b=Edge(drawn[k],c.displayLo);b<Edge(drawn[k],c.displayHi);b+=2)
      if(drawn[k]->GetBinContent(b)>0)marks->SetPoint(marks->GetN(),drawn[k]->GetBinCenter(b),drawn[k]->GetBinContent(b));
      marks->SetMarkerColor(Color(which[k]));marks->SetMarkerStyle(Marker(which[k]));marks->SetMarkerSize(which[k]?.9:.65);marks->Draw("P SAME");}
    TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.041);text.DrawLatex(.15,.95,panel?"Unit-area normalized shape":"Raw candidate counts");pad->RedrawAxis();
  }
  canvas->cd();TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.029);t.DrawLatex(.045,.96,c.title.c_str());
  TLegend*legend=new TLegend(.045,.795,.97,.938);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextFont(42);legend->SetTextSize(.025);legend->SetMargin(.055);
  for(size_t k=0;k<which.size();++k){size_t j=which[k];const Metric&m=metrics.find(Key(choices[j].id,"full"))->second;
    legend->AddEntry(raw[k],Form("%s: Nwin = %.0f; purity = %.1f%%; S/S_{baseline} = %.1f%%",choices[j].label.c_str(),m.n,100*m.purity,100*m.retention),"lp");}
  legend->Draw();t.SetTextSize(.021);t.DrawLatex(.045,.75,"Selections chosen on training events; this overlay uses all 10,000 inputs (descriptive, not held-out validation).");
  t.SetTextSize(.019);t.DrawLatex(.07,.025,Form("Fixed Nwin: [%.3f, %.3f); normalization: [%.2f, %.2f). No mass constraint or smoothing. Purity = S/(S+B), sideband estimate.",c.windowLo,c.windowHi,c.lo,c.hi));
  Save(canvas,stem,out);
}
static void Frontier(const std::vector<Choice>&choices,const std::map<std::string,Metric>&metrics,
  const Settings&c,const std::string&stem,TFile&out){
  TCanvas*canvas=new TCanvas("PurityRetention","Train-selected cuts and held-out validation",1400,1000);
  canvas->SetLeftMargin(.12);canvas->SetRightMargin(.05);canvas->SetBottomMargin(.12);canvas->SetTopMargin(.30);canvas->SetGrid();
  double xmin=1,xmax=1,ymin=1,ymax=0;
  for(size_t j=0;j<choices.size();++j)for(int q=0;q<2;++q){const Metric&m=metrics.find(Key(choices[j].id,q?"validation":"train"))->second;
    xmin=std::min(xmin,m.retention);xmax=std::max(xmax,m.retention);ymin=std::min(ymin,m.purity-m.purityError);ymax=std::max(ymax,m.purity+m.purityError);}
  double dx=std::max(.1,xmax-xmin),dy=std::max(.05,ymax-ymin);
  TH1D*frame=new TH1D("PurityRetentionFrame","",1,std::max(0.,xmin-.12*dx),xmax+.15*dx);frame->SetDirectory(0);frame->SetStats(false);
  frame->SetMinimum(std::max(0.,ymin-.2*dy));frame->SetMaximum(std::min(1.,ymax+.2*dy));
  frame->GetXaxis()->SetTitle("Background-subtracted signal retention, S/S_{baseline}");frame->GetYaxis()->SetTitle("Purity, S/(S+B)");frame->GetYaxis()->SetTitleOffset(1.4);frame->Draw("AXIS");
  TLegend*legend=new TLegend(.13,.76,.94,.91);legend->SetNColumns(2);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.024);
  for(size_t j=0;j<choices.size();++j){const Metric&a=metrics.find(Key(choices[j].id,"train"))->second;const Metric&b=metrics.find(Key(choices[j].id,"validation"))->second;
    TGraph*line=new TGraph(2);line->SetPoint(0,a.retention,a.purity);line->SetPoint(1,b.retention,b.purity);line->SetLineColor(Color(j));line->SetLineStyle(3);line->Draw("L SAME");
    for(int q=0;q<2;++q){const Metric&m=q?b:a;TGraphErrors*g=new TGraphErrors(1);g->SetName(Form("Choice%u_%s",(unsigned)j,q?"Validation":"Training"));g->SetPoint(0,m.retention,m.purity);g->SetPointError(0,0,m.purityError);
      g->SetMarkerColor(Color(j));g->SetLineColor(Color(j));g->SetMarkerStyle(q?20:24);g->SetMarkerSize(1.5);g->Draw("P SAME");legend->AddEntry(g,(choices[j].label+(q?" / validation":" / training")).c_str(),"pe");out.cd();g->Write();canvas->cd();}
  }legend->Draw();TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.028);t.DrawLatex(.12,.95,c.title.c_str());t.SetTextSize(.022);
  t.DrawLatex(.12,.725,"Open: training; filled: held-out validation. Dashed links join the same fixed selection, not a fit.");
  t.SetTextSize(.018);t.DrawLatex(.12,.035,"Vertical bars: sideband Poisson approximation. Retention errors/correlations are not shown; no significance claim from point separation.");
  Save(canvas,stem,out);
}
static void Flag(Fit&r,const char*s){if(!r.flags.empty())r.flags+=';';r.flags+=s;}
static double QuadProduct(const double*a,const TMatrixDSym&v,const double*b){double s=0;for(int i=0;i<6;++i)for(int j=0;j<6;++j)s+=a[i]*v(i,j)*b[j];return s;}
static Fit FitOne(TH1D*h,const Range&r,const Settings&c,int index,TFile&out){
  Fit a;int first=Edge(h,r.lo),last=Edge(h,r.hi)-1;double total=h->Integral(first,last);a.ndf=last-first+1-6;
  if(total<c.minEntries||total<=0||last-first+1<=6){Flag(a,"insufficient_fit_entries");return a;}
  const double center=(r.lo+r.hi)/2,scale=(r.hi-r.lo)/2;
  std::string formula=Form("[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*((x-%.17g)/%.17g)+[5]*((x-%.17g)/%.17g)^2",center,scale,center,scale);
  a.model=new TF1(Form("PurityMassFit_%d",index),formula.c_str(),r.lo,r.hi);
  a.model->SetParNames("Height","Mean","Sigma","B0","B1Scaled","B2Scaled");
  double side=0;int sideBins=0,peak=Edge(h,c.windowLo);
  for(int b=first;b<=last;++b){double x=h->GetBinCenter(b);if(x<c.windowLo||x>=c.windowHi){side+=h->GetBinContent(b);++sideBins;}
    if(x>=c.meanLo&&x<c.meanHi&&h->GetBinContent(b)>h->GetBinContent(peak))peak=b;}
  double background=sideBins?side/sideBins:0,mean=std::max(c.meanLo,std::min(c.meanHi,h->GetBinCenter(peak)));
  a.model->SetParameters(std::max(1.,h->GetBinContent(peak)-background),mean,std::sqrt(c.sigmaLo*c.sigmaHi),background,0.,0.);
  a.model->SetParLimits(0,0.,total*1000);a.model->SetParLimits(1,c.meanLo,c.meanHi);a.model->SetParLimits(2,c.sigmaLo,c.sigmaHi);
  // I integrates the model over each bin and divides by its width. Therefore
  // model components are counts/bin and their mass integrals divide by bw.
  TFitResultPtr result=h->Fit(a.model,"RSQLI0N");a.status=int(result);
  if(!result.Get()){Flag(a,"missing_fit_result");return a;}
  a.covariance=result->CovMatrixStatus();out.cd();result->Write(Form("PurityFitResult_%d",index));a.model->Write();
  if(a.status!=0)Flag(a,"nonzero_fit_status");if(a.covariance!=3)Flag(a,"covariance_not_accurate");
  double p[6];a.model->GetParameters(p);for(int k=0;k<6;++k)if(!TMath::Finite(p[k])){Flag(a,"nonfinite_parameters");return a;}
  a.mean=p[1];a.sigma=p[2];if(p[0]<=0||p[2]<=0){Flag(a,"nonpositive_signal");return a;}
  if(std::min(p[1]-c.meanLo,c.meanHi-p[1])<=1.e-5*(c.meanHi-c.meanLo))Flag(a,"mean_at_bound");
  if(std::min(p[2]-c.sigmaLo,c.sigmaHi-p[2])<=1.e-5*(c.sigmaHi-c.sigmaLo))Flag(a,"sigma_at_bound");
  const double xlo=(r.lo-center)/scale,xhi=(r.hi-center)/scale;
  double minB=std::min(p[3]+p[4]*xlo+p[5]*xlo*xlo,p[3]+p[4]*xhi+p[5]*xhi*xhi);
  if(p[5]>0){double vertex=-p[4]/(2*p[5]);if(vertex>xlo&&vertex<xhi)minB=std::min(minB,p[3]+p[4]*vertex+p[5]*vertex*vertex);}
  if(!TMath::Finite(minB)||minB<0)Flag(a,"negative_background_in_fit_range");
  a.background=new TF1(Form("PurityMassBackground_%d",index),Form("[0]+[1]*((x-%.17g)/%.17g)+[2]*((x-%.17g)/%.17g)^2",center,scale,center,scale),r.lo,r.hi);
  a.background->SetParameters(p[3],p[4],p[5]);a.background->Write();
  const double zlo=(c.windowLo-p[1])/p[2],zhi=(c.windowHi-p[1])/p[2];
  const double elo=std::exp(-.5*zlo*zlo),ehi=std::exp(-.5*zhi*zhi);
  const double base=p[2]*std::sqrt(TMath::Pi()/2)*(TMath::Erf(zhi/std::sqrt(2.))-TMath::Erf(zlo/std::sqrt(2.)))/c.bw;
  const double dl=c.windowLo-center,dh=c.windowHi-center;
  double gs[6]={base,p[0]*(elo-ehi)/c.bw,0,0,0,0};
  double gb[6]={0,0,0,(c.windowHi-c.windowLo)/c.bw,(dh*dh-dl*dl)/(2*scale*c.bw),(dh*dh*dh-dl*dl*dl)/(3*scale*scale*c.bw)};
  a.s=p[0]*base;a.b=gb[3]*p[3]+gb[4]*p[4]+gb[5]*p[5];gs[2]=a.s/p[2]+p[0]*(zlo*elo-zhi*ehi)/c.bw;
  if(!TMath::Finite(a.s)||!TMath::Finite(a.b)||a.s<=0||a.b<=0)Flag(a,"invalid_window_components");
  for(int b=first;b<=last;++b){double expected=a.model->Integral(h->GetBinLowEdge(b),h->GetBinLowEdge(b)+c.bw)/c.bw;double observed=h->GetBinContent(b);
    if(!TMath::Finite(expected)||expected<=0){Flag(a,"nonpositive_expected_bin");break;}
    a.deviance+=2*(expected-observed+(observed>0?observed*std::log(observed/expected):0));}
  a.devianceP=TMath::Prob(std::max(0.,a.deviance),a.ndf);
  TMatrixDSym covariance=result->GetCovarianceMatrix();if(covariance.GetNrows()!=6){Flag(a,"wrong_covariance_dimension");return a;}
  if(!a.valid())return a;
  double gp[6];for(int i=0;i<6;++i)gp[i]=(a.b*gs[i]-a.s*gb[i])/((a.s+a.b)*(a.s+a.b));
  double vs=QuadProduct(gs,covariance,gs),vb=QuadProduct(gb,covariance,gb),vp=QuadProduct(gp,covariance,gp);
  if(!TMath::Finite(vs)||!TMath::Finite(vb)||!TMath::Finite(vp)||vs<0||vb<0||vp<0){Flag(a,"invalid_propagated_covariance");return a;}
  a.sError=std::sqrt(vs);a.bError=std::sqrt(vb);a.purity=a.s/(a.s+a.b);a.purityError=std::sqrt(vp);return a;
}
static void Fits(const std::vector<Choice>&choices,const std::map<std::string,TH1D*>&hist,const Settings&c,const std::string&stem,TFile&out){
  std::ofstream csv((stem+"_fit_metrics.csv").c_str());Require(bool(csv),"cannot create fit CSV");csv<<std::setprecision(17);
  csv<<"choice_id,split,fitLow,fitHigh,accepted,status,covariance_quality,flags,mean,sigma,poisson_deviance,ndf,deviance_asymptotic_p,S,SError,B,BError,purity,purityError,signalRetention\n";
  const char*splits[]={"full","train","validation"};int index=0;
  for(size_t r=0;r<c.fits.size();++r)for(int q=0;q<3;++q){std::vector<Fit>fits;
    TCanvas*canvas=new TCanvas(Form("PurityFitRange%u_%s",(unsigned)r,splits[q]),"Independent Gaussian + quadratic Poisson fits",1800,choices.size()>3?1100:650);
    int cols=choices.size()>3?3:int(choices.size());int rows=int((choices.size()+cols-1)/cols);canvas->Divide(cols,rows);
    for(size_t j=0;j<choices.size();++j){TH1D*h=hist.find(Key(choices[j].id,splits[q]))->second;Fit f=FitOne(h,c.fits[r],c,index++,out);fits.push_back(f);
      csv<<CsvQuote(choices[j].id)<<','<<splits[q]<<','<<c.fits[r].lo<<','<<c.fits[r].hi<<','<<int(f.valid())<<','<<f.status<<','<<f.covariance<<','<<CsvQuote(f.valid()?"ok":f.flags)<<',';
      if(f.model)csv<<f.mean<<','<<f.sigma<<','<<f.deviance<<','<<f.ndf<<','<<f.devianceP;else csv<<",,,,";csv<<',';
      if(f.valid()){csv<<f.s<<','<<f.sError<<','<<f.b<<','<<f.bError<<','<<f.purity<<','<<f.purityError<<',';if(fits[0].valid())csv<<f.s/fits[0].s;}else csv<<",,,,,,";csv<<'\n';
      canvas->cd(j+1);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.13);gPad->SetTopMargin(.18);gPad->SetTicks(1,1);
      TH1D*draw=(TH1D*)h->Clone(Form("FitDraw_%u_%d_%u",(unsigned)r,q,(unsigned)j));draw->SetDirectory(0);HistStyle(draw,j);draw->GetXaxis()->SetRange(Edge(draw,c.fits[r].lo),Edge(draw,c.fits[r].hi)-1);Axis(draw,Form("Candidates / %.3g MeV/c^{2}",1000*c.bw));draw->SetMinimum(0);draw->SetMaximum(draw->GetMaximum()*1.15);draw->Draw("E");
      if(f.model){f.model->SetLineColor(kRed+1);f.model->SetLineWidth(2);f.model->Draw("SAME");}if(f.background){f.background->SetLineColor(kBlue+1);f.background->SetLineStyle(2);f.background->Draw("SAME");}
      TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.036);text.DrawLatex(.14,.94,Form("%s / %s",choices[j].label.c_str(),splits[q]));
      text.SetTextSize(.033);text.SetTextColor(f.valid()?kBlack:kRed+1);text.DrawLatex(.14,.875,f.valid()?Form("Valid; P = %.3f #pm %.3f (fit-model cross-check)",f.purity,f.purityError):Form("INVALID: status %d, covariance %d",f.status,f.covariance));
    }
    Save(canvas,Form("%s_fit_range%02u_%s",stem.c_str(),(unsigned)r,splits[q]),out);
  }csv.close();Require(bool(csv),"failed fit CSV write");
}
} // namespace LambdaPurityQA

Int_t drawLambdaPurityStudy(const char*scanDirectory,const char*outputPrefix){
  using namespace LambdaPurityQA;
  try{
    Require(scanDirectory&&*scanDirectory&&outputPrefix&&*outputPrefix,"missing arguments");
    const std::string dir(scanDirectory),stem(outputPrefix);Settings c=ReadManifest(Join(dir,"plot_manifest.tsv"));
    std::vector<Row> selectionRows=ReadCsv(Join(dir,"chosen_selections.csv"));std::vector<Choice> choices;std::map<std::string,size_t> ids;
    for(size_t i=0;i<selectionRows.size();++i){Choice a;a.id=Get(selectionRows[i],"choice_id");a.label=Get(selectionRows[i],"label");a.role=Get(selectionRows[i],"role");Require(ids.insert(std::make_pair(a.id,choices.size())).second,"duplicate choice id");choices.push_back(a);}
    Require(!choices.empty()&&choices[0].role=="baseline","baseline must be first");
    const int bins=int(std::floor((c.hi-c.lo)/c.bw+.5));Require(bins>0&&std::fabs(c.lo+bins*c.bw-c.hi)<1.e-9,"nonintegral bin count");
    std::map<std::string,TH1D*> hist;std::map<std::string,std::vector<bool> > seen;std::map<std::string,double> norm;
    const char*splits[]={"full","train","validation"};
    for(size_t j=0;j<choices.size();++j)for(int q=0;q<3;++q){std::string key=Key(choices[j].id,splits[q]);TH1D*h=new TH1D(Form("Counts_%u_%s",(unsigned)j,splits[q]),"",bins,c.lo,c.hi);h->SetDirectory(0);h->Sumw2();hist[key]=h;seen[key]=std::vector<bool>(bins,false);}
    std::vector<Row> rows=ReadCsv(Join(dir,"histograms.csv"));
    for(size_t j=0;j<rows.size();++j){const Row&r=rows[j];std::string key=Key(Get(r,"choice_id"),Get(r,"split"));Require(hist.count(key),"unknown histogram choice/split");TH1D*h=hist[key];
      double bval=N(r,"bin");int b=int(bval);Require(bval==b&&b>=1&&b<=bins&&!seen[key][b-1],"invalid/duplicate bin");seen[key][b-1]=true;
      Require(std::fabs(N(r,"low")-h->GetBinLowEdge(b))<1.e-9&&std::fabs(N(r,"high")-h->GetBinLowEdge(b+1))<1.e-9,"bin-edge mismatch");
      double count=N(r,"count"),error=N(r,"error"),normalization=N(r,"normalization_count");Require(count>=0&&std::fabs(count-std::floor(count+.5))<1.e-8&&error>=0&&std::fabs(error*error-count)<1.e-7*std::max(1.,count),"not unweighted Poisson counts");
      if(norm.count(key))Require(norm[key]==normalization,"inconsistent normalization_count");else norm[key]=normalization;
      h->SetBinContent(b,count);h->SetBinError(b,error);
    }
    for(std::map<std::string,TH1D*>::iterator it=hist.begin();it!=hist.end();++it){for(int b=0;b<bins;++b)Require(seen[it->first][b],"missing histogram bin");Require(it->second->Integral(1,bins)==norm[it->first],"normalization_count mismatch");}
    for(size_t j=0;j<choices.size();++j)for(int b=1;b<=bins;++b)Require(hist[Key(choices[j].id,"full")]->GetBinContent(b)==hist[Key(choices[j].id,"train")]->GetBinContent(b)+hist[Key(choices[j].id,"validation")]->GetBinContent(b),"full histogram != train + validation");
    rows=ReadCsv(Join(dir,"chosen_metrics.csv"));std::map<std::string,Metric> metrics;
    for(size_t j=0;j<rows.size();++j){const Row&r=rows[j];if(Get(r,"sideband")!="default")continue;Metric m;m.n=N(r,"Nwin");m.s=N(r,"S");m.b=N(r,"B");m.purity=N(r,"purity");m.purityError=N(r,"purityError");m.retention=N(r,"signalRetention");std::string key=Key(Get(r,"choice_id"),Get(r,"split"));Require(hist.count(key)&&metrics.insert(std::make_pair(key,m)).second,"invalid/duplicate metric");
      Require(m.n>0&&m.s>0&&m.b>=0&&m.purity>=0&&m.purity<=1&&m.purityError>=0&&m.retention>0,"invalid selected metric");
      Require(std::fabs(m.s+m.b-m.n)<1.e-7&&std::fabs(m.purity-m.s/m.n)<1.e-10,"metric arithmetic mismatch");TH1D*h=hist[key];Require(h->Integral(Edge(h,c.windowLo),Edge(h,c.windowHi)-1)==m.n,"window count does not match metrics");}
    Require(metrics.size()==hist.size(),"missing default metrics for selected split");
    int overlays=std::max(1,int((choices.size()-1+2)/3));std::vector<std::string>products;products.push_back(stem+".root");products.push_back(stem+"_fit_metrics.csv");products.push_back(stem+"_audit.txt");
    const char*ext[]={".pdf",".png"};for(int e=0;e<2;++e){products.push_back(stem+"_retention"+ext[e]);for(int k=0;k<overlays;++k)products.push_back(Form("%s_overlay_%02d%s",stem.c_str(),k,ext[e]));
      for(size_t r=0;r<c.fits.size();++r)for(int q=0;q<3;++q)products.push_back(Form("%s_fit_range%02u_%s%s",stem.c_str(),(unsigned)r,splits[q],ext[e]));}
    for(size_t k=0;k<products.size();++k)Require(gSystem->AccessPathName(products[k].c_str()),"refusing to overwrite "+products[k]);
    gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetTextFont(42);gStyle->SetGridColor(TColor::GetColor("#e0e3e7"));gStyle->SetGridStyle(3);
    TFile out((stem+".root").c_str(),"CREATE");Require(!out.IsZombie(),"cannot create ROOT output");for(std::map<std::string,TH1D*>::iterator it=hist.begin();it!=hist.end();++it)it->second->Write();
    for(int k=0;k<overlays;++k){std::vector<size_t>which;which.push_back(0);for(size_t j=1+3*k;j<choices.size()&&j<size_t(4+3*k);++j)which.push_back(j);Overlay(choices,which,hist,metrics,c,Form("%s_overlay_%02d",stem.c_str(),k),k,out);}
    Frontier(choices,metrics,c,stem+"_retention",out);Fits(choices,hist,c,stem,out);
    std::ostringstream audit;audit<<"Saved CSV-only QA; source="<<dir<<"\n"<<c.rawManifest;
    audit<<"PASS: all CSV bins complete, integer unweighted counts, Poisson errors, full=train+validation, complete-range normalization and fixed-window metrics agree.\n"
      <<"Selections and ranking are imported, never recomputed here. Full-sample overlays are descriptive. Held-out metrics shown separately from training.\n"
      <<"Fits use RSQLI0N, Gaussian plus scaled quadratic, integrated bin expectation; model mass integral/bin width gives yield. Fits are cross-checks, not selection scores.\n"
      <<"Accepted fit requires status0, covariance3, parameters inside explicit bounds, nonnegative quadratic over whole fit interval, positive S/B components and valid propagated full covariance. Invalid physical estimates are blank in CSV.\n"
      <<"Deviance NDF uses fitted bins minus six free parameters. Its chi-square p-value is asymptotic, unreliable for sparse bins, and a model diagnostic rather than an extra selection gate. Numerical acceptance does not prove model adequacy.\n"
      <<"Uncertainties are conditional on fit model; model dependence and shared-candidate correlations are not included. No true reconstruction efficiency is measured.\n";
    out.cd();TObjString(audit.str().c_str()).Write("QAProvenance");out.Flush();Require(!out.TestBit(TFile::kWriteError),"ROOT write failure");out.Close();
    std::ofstream log((stem+"_audit.txt").c_str());log<<audit.str();log.close();Require(bool(log),"audit write failure");
    std::cout<<"[drawLambdaPurityStudy] PASS: "<<choices.size()<<" choices; "<<overlays<<" overlays; independent fit validity in CSV.\n";return 0;
  }catch(const std::exception&e){std::cerr<<"[drawLambdaPurityStudy] ERROR: "<<e.what()<<std::endl;return 1;}
}

#endif // compiled implementation; CINT only needs the public entry point

