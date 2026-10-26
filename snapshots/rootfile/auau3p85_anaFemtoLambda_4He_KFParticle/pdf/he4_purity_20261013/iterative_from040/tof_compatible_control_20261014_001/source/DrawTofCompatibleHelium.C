// READ-only physical dE/dx controls. D TOF window also accepts He4 m/q.
#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#ifndef __CINT__
namespace TofCompatibleHelium {
typedef std::map<std::string,double> Settings;
struct Reference {int index,accepted,quality;double lo,hi,r,legacy[4],old[3];};
double Need(const Settings& c,const char* key) {Settings::const_iterator i=c.find(key);if(i==c.end()||!TMath::Finite(i->second))throw std::runtime_error(std::string("Missing setting: ")+key);return i->second;}
void Require(bool x,const char* what) {if(!x)throw std::runtime_error(what);}
Settings ReadSettings(const char* path) {std::ifstream f(path);Require(bool(f),"Cannot read settings TSV");Settings c;std::string k;double x;while(f>>k>>x)Require(c.insert(std::make_pair(k,x)).second,"Duplicate setting");Require(f.eof(),"Malformed setting TSV");return c;}
std::vector<Reference> ReadReferences(const char* path) {
  std::ifstream f(path);Require(bool(f),"Cannot read references TSV");std::vector<Reference> out;Reference r;
  while(f>>r.index>>r.lo>>r.hi>>r.r>>r.accepted>>r.quality>>r.legacy[0]>>r.legacy[1]>>r.legacy[2]>>r.legacy[3]>>r.old[0]>>r.old[1]>>r.old[2])out.push_back(r);
  Require(f.eof()&&!out.empty(),"Malformed/empty reference TSV");return out;
}
bool Close(double a,double b,double tol) {return std::fabs(a-b)<=tol*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
void Axis(TAxis* a,TAxis* b,double tol) {Require(a->GetNbins()==b->GetNbins(),"TOF histogram axis dimensions differ");for(int i=1;i<=a->GetNbins()+1;++i)Require(Close(a->GetBinLowEdge(i),b->GetBinLowEdge(i),tol),"TOF histogram axis boundaries differ");}
void Subsets(TH2* all,TH2* tag3,TH2* tagD,const Settings& c,const TString& output) {
  Axis(all->GetXaxis(),tag3->GetXaxis(),Need(c,"axisRelativeTolerance"));Axis(all->GetYaxis(),tag3->GetYaxis(),Need(c,"axisRelativeTolerance"));Axis(all->GetXaxis(),tagD->GetXaxis(),Need(c,"axisRelativeTolerance"));Axis(all->GetYaxis(),tagD->GetYaxis(),Need(c,"axisRelativeTolerance"));
  double counts[3]={0,0,0},variances[3]={0,0,0},minRemaining=1.e300;unsigned cells=0;
  for(int x=0;x<=all->GetNbinsX()+1;++x)for(int y=0;y<=all->GetNbinsY()+1;++y) {
    TH2* h[]={all,tag3,tagD};double n[3];
    for(int j=0;j<3;++j) {n[j]=h[j]->GetBinContent(x,y);const double v=std::pow(h[j]->GetBinError(x,y),2);Require(TMath::Finite(n[j])&&n[j]>=0&&Close(n[j],v,Need(c,"varianceRelativeTolerance")),"Histogram is not finite unweighted unit counts");counts[j]+=n[j];variances[j]+=v;}
    Require(n[1]<=n[0]&&n[2]<=n[0],"TOF tag is not an ALL cell subset");Require(n[1]+n[2]<=n[0],"Disjoint TOF tag sum exceeds ALL");minRemaining=std::min(minRemaining,n[0]-n[1]-n[2]);++cells;
  }
  std::ofstream f(output.Data());Require(bool(f),"Cannot write subset audit");f<<std::setprecision(17)<<"{\"cellsIncludingFlows\":"<<cells<<",\"axesEqual\":true,\"tag3SubsetOfAll\":true,\"tagDSubsetOfAll\":true,\"disjointTagSumSubsetOfAll\":true,\"minimumRemainingCount\":"<<minRemaining<<",\"allCount\":"<<counts[0]<<",\"tag3Count\":"<<counts[1]<<",\"tagDCount\":"<<counts[2]<<",\"allVariance\":"<<variances[0]<<",\"tag3Variance\":"<<variances[1]<<",\"tagDVariance\":"<<variances[2]<<",\"sumw2Sizes\":["<<all->GetSumw2N()<<','<<tag3->GetSumw2N()<<','<<tagD->GetSumw2N()<<"],\"scope\":\"Mass2 windows are disjoint, TAG3/TAGD are correlated ALL subsets. TAGD is D/He4-compatible, NOT pure He4. No independent ALL+TAG likelihood.\"}\n";
}
double Interval(TH1D* h,double lo,double hi) {double sum=0;for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=lo&&h->GetBinCenter(b)<=hi)sum+=h->GetBinContent(b);return sum;}
struct Moment {double mode,smoothedHeight,count,mean,rms,lo,hi;bool enough;Moment():mode(-999),smoothedHeight(0),count(0),mean(-999),rms(-999),lo(0),hi(0),enough(false){}};
Moment Diagnose(TH1D* h,const Settings& c) {
  Moment m;const int half=int(Need(c,"modeSmoothingHalfBins"));
  for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=Need(c,"nuclearDedxMin")&&h->GetBinCenter(b)<=Need(c,"nuclearDedxMax")) {
    double count=0;int bins=0;for(int j=std::max(1,b-half);j<=std::min(h->GetNbinsX(),b+half);++j){count+=h->GetBinContent(j);++bins;}count/=bins;
    if(count>m.smoothedHeight){m.smoothedHeight=count;m.mode=h->GetBinCenter(b);}
  }
  if(m.mode<-900)return m;m.lo=std::max(Need(c,"nuclearDedxMin"),m.mode-Need(c,"momentHalfWindowKeV"));m.hi=std::min(Need(c,"nuclearDedxMax"),m.mode+Need(c,"momentHalfWindowKeV"));double first=0,second=0;
  for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=m.lo&&h->GetBinCenter(b)<=m.hi){double n=h->GetBinContent(b),x=h->GetBinCenter(b);m.count+=n;first+=n*x;second+=n*x*x;}
  if(m.count>0){m.mean=first/m.count;m.rms=std::sqrt(std::max(0.,second/m.count-m.mean*m.mean));}m.enough=m.count>=Need(c,"minimumMomentCounts");return m;
}
void ReferenceLines(const Reference& r,double lo,double hi,double ymin,double ymax,TLegend* legend) {
  const int colors[]={kGray+2,kOrange+7,kBlue+1,kRed+1};const char* names[]={"D","T","^{3}He","^{4}He"};
  for(int s=0;s<4;++s)if(r.legacy[s]>=lo&&r.legacy[s]<=hi) {TLine* line=new TLine(r.legacy[s],ymin,r.legacy[s],ymax);line->SetLineColor(colors[s]);line->SetLineStyle(2);line->Draw();legend->AddEntry(line,(std::string("Legacy ")+names[s]).c_str(),"l");}
  for(int s=0;s<3;++s)if(r.old[s]>=lo&&r.old[s]<=hi) {TLine* line=new TLine(r.old[s],ymin,r.old[s],ymax);line->SetLineColor(colors[s+1]);line->SetLineStyle(1);line->SetLineWidth(2);line->Draw();legend->AddEntry(line,(std::string("from040 ")+names[s+1]).c_str(),"l");}
}
void Draw(TH1D* observed,const Reference& r,const Settings& c,int sample,bool nuclear,bool log,const Moment& m) {
  gPad->SetLeftMargin(Need(c,"marginLeft"));gPad->SetRightMargin(Need(c,"marginRight"));gPad->SetBottomMargin(Need(c,"marginBottom"));gPad->SetTopMargin(Need(c,"marginTop"));gPad->SetLogy(log);
  TH1D* h=static_cast<TH1D*>(observed->Clone(Form("display_%d_%d_%d_%d",r.index,sample,int(nuclear),int(log))));h->SetDirectory(0);
  const double lo=std::max(h->GetXaxis()->GetXmin(),Need(c,nuclear?"nuclearDedxMin":"fullDedxMin")),hi=std::min(h->GetXaxis()->GetXmax(),Need(c,nuclear?"nuclearDedxMax":"fullDedxMax"));h->GetXaxis()->SetRangeUser(lo,hi);h->GetXaxis()->SetRange(std::max(1,h->GetXaxis()->GetFirst()),std::min(h->GetNbinsX(),h->GetXaxis()->GetLast()));double max=1;
  for(int b=h->GetXaxis()->GetFirst();b<=h->GetXaxis()->GetLast();++b)max=std::max(max,h->GetBinContent(b)+h->GetBinError(b));const double ymin=log?Need(c,"logMinimum"):0.,ymax=max*Need(c,log?"logMaximumFactor":"linearMaximumFactor");h->SetMinimum(ymin);h->SetMaximum(ymax);
  h->SetTitle(Form(";TPC dE/dx (keV/cm);Tracks / %.2f keV/cm",h->GetBinWidth(1)));h->GetXaxis()->SetTitleSize(Need(c,"axisTitleSize"));h->GetYaxis()->SetTitleSize(Need(c,"axisTitleSize"));h->GetXaxis()->SetLabelSize(Need(c,"axisLabelSize"));h->GetYaxis()->SetLabelSize(Need(c,"axisLabelSize"));h->GetXaxis()->SetTitleOffset(Need(c,"xTitleOffset"));h->GetYaxis()->SetTitleOffset(Need(c,"yTitleOffset"));h->SetMarkerStyle(20);h->SetMarkerSize(.25);h->Draw("E1");
  TLegend* legend=new TLegend(.18,.74,.96,.91);legend->SetNColumns(3);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextFont(42);legend->SetTextSize(Need(c,"legendTextSize"));ReferenceLines(r,lo,hi,ymin,ymax,legend);legend->Draw();
  TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(Need(c,"headerTextSize"));const char* labels[]={"ALL pre-PID","TOF ^{3}He window","TOF D / ^{4}He-compatible"};t.DrawLatex(.18,.965,Form("R=%.3f | %s | %s",r.r,labels[sample],log?"log":"linear"));
  t.SetTextSize(Need(c,"captionTextSize"));if(sample) t.DrawLatex(.18,.09,Form("%.5f < (m/q)^{2} < %.5f; pre-TPC subset, NOT pure isotope",Need(c,sample==1?"tag3Mass2Min":"tagDMass2Min"),Need(c,sample==1?"tag3Mass2Max":"tagDMass2Max")));
  else t.DrawLatex(.18,.09,"Unselected ALL; same basic quality as TOF controls");
  t.DrawLatex(.18,.055,Form("from040 centers are references (%s); no fit / no PID changes",r.accepted?"old accepted":"old unaccepted"));
  t.DrawLatex(.18,.023,Form("Nuclear window mode %.2f, N %.0f%s; window mean/RMS are NOT Gaussian fit",m.mode,m.count,m.enough?"":" [low count]"));
}
void Map(TCanvas& canvas,TH2* input,const std::vector<Reference>& refs,const Settings& c,int sample,const TString& pdf,const TString& png) {
  canvas.Clear();gPad->SetLeftMargin(.14);gPad->SetRightMargin(.16);gPad->SetBottomMargin(.24);gPad->SetTopMargin(.25);gPad->SetLogz();TH2* h=static_cast<TH2*>(input->Clone(Form("map_%d",sample)));h->SetDirectory(0);h->RebinY(int(Need(c,"dedxRebinFactor")));
  const double tol=Need(c,"axisRelativeTolerance");h->GetXaxis()->SetRange(h->GetXaxis()->FindFixBin(Need(c,"rigidityMin")+tol),h->GetXaxis()->FindFixBin(Need(c,"rigidityMax")-tol));h->GetYaxis()->SetRangeUser(Need(c,"nuclearDedxMin"),Need(c,"nuclearDedxMax"));h->SetTitle(";R = |p/q| (GeV/c);TPC dE/dx (keV/cm)");h->GetXaxis()->SetTitleOffset(1.3);h->GetYaxis()->SetTitleOffset(1.6);h->GetXaxis()->SetTitleSize(.037);h->GetYaxis()->SetTitleSize(.037);h->GetXaxis()->SetLabelSize(.030);h->GetYaxis()->SetLabelSize(.030);h->SetMinimum(1);h->Draw("COLZ");
  TLegend* legend=new TLegend(.14,.79,.94,.925);legend->SetNColumns(3);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.022);legend->SetTextFont(42);const int colors[]={kOrange+7,kBlue+1,kRed+1};const char* names[]={"T","^{3}He","^{4}He"};
  for(int s=0;s<3;++s)for(int old=0;old<2;++old){TGraph* graph=new TGraph();for(unsigned i=0;i<refs.size();++i)graph->SetPoint(i,refs[i].r,old?refs[i].old[s]:refs[i].legacy[s+1]);graph->SetLineColor(colors[s]);graph->SetLineStyle(old?1:2);graph->SetLineWidth(2);graph->Draw("L same");legend->AddEntry(graph,(std::string(old?"from040 ":"Legacy ")+names[s]).c_str(),"l");}legend->Draw();TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.028);const char* titles[]={"ALL pre-PID","TOF 3He window (TPC-independent)","TOF D / 4He-compatible window (NOT pure He4)"};text.DrawLatex(.14,.965,titles[sample]);text.SetTextSize(.019);text.DrawLatex(.14,.07,"Line centers: legacy and old from040 Gaussian references, NOT isotope identification");text.DrawLatex(.14,.037,Form("TOF windows: 3He %.5f-%.5f; D/He4-compatible %.5f-%.5f in raw (m/q)^{2}",Need(c,"tag3Mass2Min"),Need(c,"tag3Mass2Max"),Need(c,"tagDMass2Min"),Need(c,"tagDMass2Max")));canvas.Print(pdf);canvas.Print(png);delete h;
}
}
#endif

void DrawTofCompatibleHelium(const char* input,const char* output,const char* settings,const char* references,const char* allKey,const char* tag3Key,const char* tagDKey) {
#ifndef __CINT__
  using namespace TofCompatibleHelium;
  try {
    const Settings c=ReadSettings(settings);const std::vector<Reference> refs=ReadReferences(references);const TString folder(output);TFile* file=TFile::Open(input,"READ");Require(file&&!file->IsZombie(),"Cannot read original ROOT");
    TH2* source[]={dynamic_cast<TH2*>(file->Get(allKey)),dynamic_cast<TH2*>(file->Get(tag3Key)),dynamic_cast<TH2*>(file->Get(tagDKey))};Require(source[0]&&source[1]&&source[2],"Required ALL, TOF3 or TOFD histogram is missing");Subsets(source[0],source[1],source[2],c,folder+"/subset_audit.json");
    const TString full=folder+"/tof_controls_full_dedx.pdf",nuclear=folder+"/tof_controls_nuclear_dedx.pdf",maps=folder+"/tof_controls_nuclear_colz.pdf";Require(gSystem->AccessPathName(full)&&gSystem->AccessPathName(nuclear)&&gSystem->AccessPathName(maps),"Refusing to overwrite TOF control PDF");
    gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetPalette(1);TCanvas map("tof_control_maps","Independent TOF compatible controls",1500,1100);map.Print(maps+"[");const char* mapNames[]={"all_nuclear_colz.png","tag3_nuclear_colz.png","tagD_nuclear_colz.png"};for(int s=0;s<3;++s)Map(map,source[s],refs,c,s,maps,folder+"/"+mapNames[s]);map.Print(maps+"]");
    TCanvas detail("tof_control_slices","Physical dEdx controls",int(Need(c,"canvasWidth")),int(Need(c,"canvasHeight")));detail.Print(full+"[");detail.Print(nuclear+"[");TFile product(folder+"/control_projections.root","NEW");Require(!product.IsZombie(),"Cannot create new projection ROOT");
    std::ofstream rows((folder+"/control_moments.tsv").Data());rows<<std::setprecision(17)<<"rigidity\tsample\tregularCount\tflowCount\tnuclearCount\tnuclearSmoothedMode\tmodeWindowLow\tmodeWindowHigh\tmodeWindowCount\tmodeWindowMean\tmodeWindowRms\tminimumCountPassed\toldReferenceHe4Mean\tlegacyHe4Mean\n";
    for(unsigned i=0;i<refs.size();++i) {
      const Reference& r=refs[i];const double tol=Need(c,"axisRelativeTolerance");const int first=source[0]->GetXaxis()->FindFixBin(r.lo+tol),last=source[0]->GetXaxis()->FindFixBin(r.hi-tol);Require(first==last,"Requested control must use one native R bin");TH1D* h[3];Moment moments[3];
      for(int sample=0;sample<3;++sample) {
        h[sample]=source[sample]->ProjectionY(Form("control_native_%d_%03d",sample,r.index),first,last,"e");h[sample]->SetDirectory(0);product.cd();h[sample]->Write();const double count=h[sample]->Integral(0,h[sample]->GetNbinsX()+1);h[sample]->Rebin(int(Need(c,"dedxRebinFactor")));Require(h[sample]->Integral(0,h[sample]->GetNbinsX()+1)==count,"Control rebin changed counts");h[sample]->SetName(Form("control_display_%d_%03d",sample,r.index));h[sample]->Write();moments[sample]=Diagnose(h[sample],c);const Moment& m=moments[sample];
        rows<<r.r<<'\t'<<sample<<'\t'<<h[sample]->Integral(1,h[sample]->GetNbinsX())<<'\t'<<h[sample]->GetBinContent(0)+h[sample]->GetBinContent(h[sample]->GetNbinsX()+1)<<'\t'<<Interval(h[sample],Need(c,"nuclearDedxMin"),Need(c,"nuclearDedxMax"))<<'\t'<<m.mode<<'\t'<<m.lo<<'\t'<<m.hi<<'\t'<<m.count<<'\t'<<m.mean<<'\t'<<m.rms<<'\t'<<int(m.enough)<<'\t'<<r.old[2]<<'\t'<<r.legacy[3]<<'\n';
      }
      for(int focus=0;focus<2;++focus){detail.Clear();detail.Divide(3,2);for(int sample=0;sample<3;++sample){detail.cd(sample+1);Draw(h[sample],r,c,sample,focus!=0,true,moments[sample]);detail.cd(sample+4);Draw(h[sample],r,c,sample,focus!=0,false,moments[sample]);}detail.Print(focus?nuclear:full);detail.Print(folder+Form("/%s_R_%03d.png",focus?"nuclear":"full",int(r.r*1000+.5)));}
      for(int s=0;s<3;++s)delete h[s];
    }detail.Print(full+"]");detail.Print(nuclear+"]");rows.close();product.Close();file->Close();delete file;Require(rows.good(),"Control moments write failed");std::cout<<"TOF_COMPATIBLE_HELIUM_COMPLETE Rbins="<<refs.size()<<"; no fit, no PID edits"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"TOF_COMPATIBLE_HELIUM_ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
