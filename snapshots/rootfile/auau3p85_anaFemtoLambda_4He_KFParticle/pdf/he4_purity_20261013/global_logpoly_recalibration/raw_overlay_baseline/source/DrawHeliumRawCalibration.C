// Histogram-only diagnostics: input READ, no event loop, no PID mutation/refit.
#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "TPad.h"
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
namespace HeliumRawOverlay {
struct Point { double r,mean,sigma;int custom; };
struct Curve { int mode,species;std::vector<Point> points; };
struct Knot { int species,quality;double r,mean,sigma; };
struct Moment { double count,mode,mean,sigma,meanError,rmsError,low,high;Moment():count(0),mode(-999),mean(-999),sigma(-999),meanError(0),rmsError(0),low(0),high(0){} };
typedef std::map<std::string,double> Settings;
double Need(const Settings& c,const char* key) {
  Settings::const_iterator i=c.find(key);if(i==c.end()||!TMath::Finite(i->second))throw std::runtime_error(std::string("Missing setting: ")+key);return i->second;
}
Settings ReadSettings(const char* path) {
  std::ifstream f(path);if(!f)throw std::runtime_error("Unreadable explicit settings");Settings c;std::string key;double value;
  while(f>>key>>value)if(!c.insert(std::make_pair(key,value)).second)throw std::runtime_error("Duplicate setting");
  if(!f.eof())throw std::runtime_error("Malformed settings");return c;
}
std::vector<Curve> ReadCurves(const char* path) {
  std::ifstream f(path);if(!f)throw std::runtime_error("Unreadable curve table");std::vector<Curve> c;int mode,species;Point p;
  while(f>>mode>>species>>p.r>>p.mean>>p.sigma>>p.custom) {
    if(mode<0||mode>2||species<0||species>3||!(p.r>0&&p.mean>0&&p.sigma>0))throw std::runtime_error("Invalid curve point");
    if(c.empty()||c.back().mode!=mode||c.back().species!=species){Curve item;item.mode=mode;item.species=species;c.push_back(item);}
    if(!c.back().points.empty()&&p.r<=c.back().points.back().r)throw std::runtime_error("Unordered curve");c.back().points.push_back(p);
  }
  if(!f.eof()||c.size()<6)throw std::runtime_error("Incomplete curves");return c;
}
std::vector<Knot> ReadKnots(const char* path) {
  std::ifstream f(path);if(!f)throw std::runtime_error("Unreadable knot table");std::vector<Knot> out;Knot k;
  while(f>>k.species>>k.r>>k.mean>>k.sigma>>k.quality)out.push_back(k);
  if(!f.eof()||out.empty())throw std::runtime_error("Invalid knot table");return out;
}
Point At(const std::vector<Curve>& curves,int mode,int species,double r) {
  for(unsigned c=0;c<curves.size();++c)if(curves[c].mode==mode&&curves[c].species==species) {
    const std::vector<Point>& p=curves[c].points;
    for(unsigned i=0;i<p.size();++i)if(p[i].r>=r) {
      if(i==0||std::fabs(p[i].r-r)<1.e-11)return p[i];Point out=p[i];const double f=(r-p[i-1].r)/(p[i].r-p[i-1].r);out.r=r;
      out.mean=p[i-1].mean+f*(p[i].mean-p[i-1].mean);out.sigma=p[i-1].sigma+f*(p[i].sigma-p[i-1].sigma);return out;
    }return p.back();
  }throw std::runtime_error("Requested curve absent");
}
int Color(int mode,int s) { const int old[]={kGray+2,kOrange+7,kBlue+1,kRed+1};return !mode?old[s]:mode==1?(s==2?kCyan+2:kMagenta+2):(s==2?kGreen+2:kViolet+2); }
std::string Label(int mode,int s) { const char* names[]={"D","T","^{3}He","^{4}He"};return std::string(!mode?"Legacy ":mode==1?"Current knots ":"Global candidate ")+names[s]; }
bool SameAxis(const TAxis& a,const TAxis& b) { if(a.GetNbins()!=b.GetNbins())return false;for(int i=1;i<=a.GetNbins()+1;++i)if(a.GetBinLowEdge(i)!=b.GetBinLowEdge(i))return false;return true; }
void AuditSubset(TH2* all,TH2* tag,const char* path) {
  std::ofstream out(path);out<<std::setprecision(17);if(!tag){out<<"{\"tofAvailable\":false}\n";return;}
  if(!SameAxis(*all->GetXaxis(),*tag->GetXaxis())||!SameAxis(*all->GetYaxis(),*tag->GetYaxis()))throw std::runtime_error("ALL/TAG axes differ");
  double na=0,nt=0,ra=0,rt=0,va=0,vt=0,minCount=1.e300,minVar=1.e300;unsigned cells=0;
  for(int x=0;x<=all->GetNbinsX()+1;++x)for(int y=0;y<=all->GetNbinsY()+1;++y) {
    const double a=all->GetBinContent(x,y),t=tag->GetBinContent(x,y),ea=std::pow(all->GetBinError(x,y),2),et=std::pow(tag->GetBinError(x,y),2);
    if(!TMath::Finite(a)||!TMath::Finite(t)||a<0||t<0||t>a)throw std::runtime_error("TAG is not an exact nonnegative ALL cell subset");
    if(et>ea+1.e-10*std::max(1.,ea))throw std::runtime_error("TAG variance exceeds ALL");
    na+=a;nt+=t;va+=ea;vt+=et;minCount=std::min(minCount,a-t);minVar=std::min(minVar,ea-et);++cells;
    if(x>0&&x<=all->GetNbinsX()&&y>0&&y<=all->GetNbinsY()){ra+=a;rt+=t;}
  }
  out<<"{\"tofAvailable\":true,\"axesExactlyEqual\":true,\"cellsIncludingFlows\":"<<cells<<",\"allClass\":\""<<all->ClassName()<<"\",\"tagClass\":\""<<tag->ClassName()<<"\",\"allEntries\":"<<all->GetEntries()<<",\"tagEntries\":"<<tag->GetEntries()<<",\"allStoredCount\":"<<na<<",\"tagStoredCount\":"<<nt<<",\"notTagStoredCount\":"<<na-nt<<",\"allRegularStoredCount\":"<<ra<<",\"tagRegularStoredCount\":"<<rt<<",\"allFlowStoredCount\":"<<na-ra<<",\"tagFlowStoredCount\":"<<nt-rt<<",\"allVarianceSum\":"<<va<<",\"tagVarianceSum\":"<<vt<<",\"notTagVarianceSumAssumingUnitWeightSubset\":"<<va-vt<<",\"allSumw2Size\":"<<all->GetSumw2N()<<",\"tagSumw2Size\":"<<tag->GetSumw2N()<<",\"minimumCountDifference\":"<<minCount<<",\"minimumVarianceDifference\":"<<minVar<<",\"scope\":\"TAG subset of ALL. No independent ALL+TAG likelihood. Unit-weight NOT_TAG variance is Var(ALL)-Var(TAG), not sum. Float32 accumulation can make stored integral differ from entries.\"}\n";
}
void DrawMap(TCanvas& canvas,TH2* input,const std::vector<Curve>& curves,const std::vector<Knot>& knots,
             const Settings& cfg,double lo,double hi,bool tag,const TString& pdf,const TString& png) {
  canvas.Clear();canvas.cd();gPad->SetLeftMargin(.11);gPad->SetRightMargin(.15);gPad->SetBottomMargin(.13);gPad->SetTopMargin(.20);gPad->SetLogz(Need(cfg,"logZ")!=0);
  TH2* h=static_cast<TH2*>(input->Clone(Form("map_%d_%g",int(tag),lo)));h->SetDirectory(0);
  const int rebin=int(Need(cfg,"displayDedxRebinFactor"));if(h->GetNbinsY()%rebin)throw std::runtime_error("Display rebin does not divide axis");h->RebinY(rebin);
  h->GetXaxis()->SetRangeUser(lo,hi);h->GetYaxis()->SetRangeUser(Need(cfg,"dedxDisplayMin"),Need(cfg,"dedxDisplayMax"));
  h->SetTitle(";R = |p/q| (GeV/c), NOT physical He momentum;TPC dE/dx (keV/cm)");h->SetMinimum(1.);h->GetXaxis()->SetTitleSize(.042);h->GetYaxis()->SetTitleSize(.042);h->Draw("COLZ");
  TLegend* legend=new TLegend(.11,.83,.88,.975);legend->SetNColumns(3);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.024);
  for(unsigned c=0;c<curves.size();++c)for(int band=-1;band<=1;++band) {
    const Curve& curve=curves[c];if(band&&curve.species<2)continue;TGraph* graph=new TGraph();
    for(unsigned i=0;i<curve.points.size();++i)graph->SetPoint(i,curve.points[i].r,curve.points[i].mean+band*Need(cfg,"curveNSigmaBand")*curve.points[i].sigma);
    graph->SetLineColor(Color(curve.mode,curve.species));graph->SetLineStyle(band?3:(curve.mode?1:2));graph->SetLineWidth(band?1:3);graph->Draw("L same");
    if(!band)legend->AddEntry(graph,Label(curve.mode,curve.species).c_str(),"l");
  }
  for(int s=2;s<4;++s)for(int quality=0;quality<2;++quality) {
    TGraph* points=new TGraph();int n=0;for(unsigned i=0;i<knots.size();++i)if(knots[i].species==s&&knots[i].quality==quality)points->SetPoint(n++,knots[i].r,knots[i].mean);
    points->SetMarkerColor(Color(1,s));points->SetMarkerStyle(quality?20:24);points->SetMarkerSize(quality?.8:.45);points->Draw("P same");
    if(s==3)legend->AddEntry(points,quality?"Filled: quality free-fit point":"Open: reference free-fit point","p");
  }
  legend->Draw();TLatex text;text.SetNDC();text.SetTextSize(.032);text.DrawLatex(.11,.985,tag?"TOF-tagged ^{3}He QA (pre-TPC; not isotope truth)":"ALL pre-PID tracks: raw dE/dx and calibration centers");
  text.SetTextSize(.021);text.DrawLatex(.11,.065,Form("R %.2f-%.2f; displayed dE/dx bin %.2f keV/cm; dotted: mean #pm %.1f#sigma",lo,hi,h->GetYaxis()->GetBinWidth(1),Need(cfg,"curveNSigmaBand")));
  text.DrawLatex(.11,.035,tag?Form("Existing tag: %.5f < (m/q)^{2} < %.5f (GeV/c^{2})^{2}; subset of ALL, no double counting",Need(cfg,"tofMass2Min"),Need(cfg,"tofMass2Max")):"Charge-inclusive, centrality-integrated; native dE/dx storage ends at 120 keV/cm");
  canvas.Print(pdf);canvas.Print(png);delete h;
}
Moment TaggedMoment(TH1D* h,const Point& legacy,const Settings& cfg) {
  Moment out;if(!h)return out;int best=-1;double maximum=-1.;
  const double low=legacy.mean-Need(cfg,"tofModeSearchNSigma")*legacy.sigma,high=legacy.mean+Need(cfg,"tofModeSearchNSigma")*legacy.sigma;
  for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=low&&h->GetBinCenter(b)<=high&&h->GetBinContent(b)>maximum){maximum=h->GetBinContent(b);best=b;}
  if(best<0||maximum<=0)return out;out.mode=h->GetBinCenter(best);out.low=out.mode-Need(cfg,"tofMomentNSigma")*legacy.sigma;out.high=out.mode+Need(cfg,"tofMomentNSigma")*legacy.sigma;double sum=0,square=0;
  for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=out.low&&h->GetBinCenter(b)<=out.high) {
    const double n=h->GetBinContent(b),x=h->GetBinCenter(b);out.count+=n;sum+=n*x;square+=n*x*x;
  }
  if(out.count>0) {
    out.mean=sum/out.count;out.sigma=std::sqrt(std::max(0.,square/out.count-out.mean*out.mean));double meanVariance=0,rmsVariance=0;
    for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=out.low&&h->GetBinCenter(b)<=out.high) {
      const double delta=h->GetBinCenter(b)-out.mean,variance=std::pow(h->GetBinError(b),2);
      meanVariance+=variance*delta*delta;rmsVariance+=variance*std::pow(delta*delta-out.sigma*out.sigma,2);
    }
    out.meanError=std::sqrt(meanVariance)/out.count;if(out.sigma>0)out.rmsError=std::sqrt(rmsVariance)/(2*out.sigma*out.count);
  }return out;
}
double IntervalCount(TH1D* h,double low,double high) {
  if(!h)return 0.;double sum=0;for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=low&&h->GetBinCenter(b)<=high)sum+=h->GetBinContent(b);return sum;
}
void DrawProjection(TH1D* original,TH1D* scaledTag,const std::vector<Curve>& curves,const Settings& cfg,
                    double r,bool tag,bool logarithmic,const Moment& moment,double tagScale=0.) {
  gPad->SetLeftMargin(.14);gPad->SetRightMargin(.04);gPad->SetBottomMargin(.18);gPad->SetTopMargin(.22);gPad->SetLogy(logarithmic);
  TH1D* h=static_cast<TH1D*>(original->Clone(Form("proj_%d_%d_%g_%g",int(tag),int(logarithmic),r,tagScale)));h->SetDirectory(0);
  const Point old4=At(curves,0,3,r);const double requestedLow=old4.mean+Need(cfg,"projectionDisplayNSigmaMin")*old4.sigma,requestedHigh=old4.mean+Need(cfg,"projectionDisplayNSigmaMax")*old4.sigma;
  const double lo=std::max(requestedLow,h->GetXaxis()->GetXmin()),hi=std::min(requestedHigh,h->GetXaxis()->GetXmax());h->GetXaxis()->SetRangeUser(lo,hi);double maximum=1.;
  for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=lo&&h->GetBinCenter(b)<=hi)maximum=std::max(maximum,h->GetBinContent(b)+h->GetBinError(b));
  const double lower=logarithmic?Need(cfg,"projectionLogMinimum"):0.,upper=maximum*Need(cfg,logarithmic?"projectionLogMaximumFactor":"projectionLinearMaximumFactor");
  h->SetMinimum(lower);h->SetMaximum(upper);h->SetTitle(Form(";TPC dE/dx (keV/cm);Stored tracks / %.2f keV/cm",h->GetBinWidth(1)));h->SetMarkerStyle(20);h->SetMarkerSize(.35);h->Draw("E1");
  TLegend* legend=new TLegend(.14,.82,.96,.955);legend->SetNColumns(3);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.026);
  if(scaledTag&&tagScale>0){TH1D* overlay=static_cast<TH1D*>(scaledTag->Clone(Form("scaled_tag_%g_%d",r,int(logarithmic))));overlay->SetDirectory(0);overlay->Scale(tagScale);overlay->SetLineColor(kGreen+2);overlay->SetLineWidth(3);overlay->SetMarkerSize(0);overlay->Draw("HIST same");legend->AddEntry(overlay,"TOF tag scaled in core","l");}
  for(unsigned c=0;c<curves.size();++c) {
    const Curve& curve=curves[c];if(curve.species<1)continue;const Point p=At(curves,curve.mode,curve.species,r);
    for(int sign=-1;sign<=1;++sign) {
      if(sign&&curve.species<2)continue;const double x=p.mean+sign*Need(cfg,"curveNSigmaBand")*p.sigma;if(x<lo||x>hi)continue;
      TLine* line=new TLine(x,lower,x,upper);line->SetLineColor(Color(curve.mode,curve.species));line->SetLineStyle(sign?3:(curve.mode?1:2));line->SetLineWidth(sign?1:2);line->Draw();if(!sign)legend->AddEntry(line,Label(curve.mode,curve.species).c_str(),"l");
    }
  }
  if(tag&&moment.count>0&&moment.mode>=lo&&moment.mode<=hi){TLine* line=new TLine(moment.mode,lower,moment.mode,upper);line->SetLineStyle(9);line->SetLineWidth(2);line->Draw();legend->AddEntry(line,"TOF local mode (no fit)","l");}
  legend->Draw();TLatex text;text.SetNDC();text.SetTextSize(.034);text.DrawLatex(.14,.975,Form("R=%.3f: %s / %s",r,tag?"TOF ^{3}He tag":"ALL pre-PID",logarithmic?"log":"linear"));
  text.SetTextSize(.021);text.DrawLatex(.14,.11,Form("Display: legacy ^{4}He %.0f to %+.0f#sigma -> %.2f-%.2f keV/cm%s",Need(cfg,"projectionDisplayNSigmaMin"),Need(cfg,"projectionDisplayNSigmaMax"),lo,hi,(requestedHigh>h->GetXaxis()->GetXmax()||requestedLow<h->GetXaxis()->GetXmin())?" (axis clipped)":""));
  if(scaledTag&&tagScale>0){text.DrawLatex(.14,.07,Form("TAG scale = %.4g, ALL/TAG core integral; normalization recorded in TSV",tagScale));text.DrawLatex(.14,.035,"Not a purity/background coefficient; TOF efficiency may bias tag shape");}
  else {text.DrawLatex(.14,.07,tag&&moment.count>0?Form("Window mean/RMS %.3f / %.3f; mode %.2f; NOT Gaussian fit",moment.mean,moment.sigma,moment.mode):"Centers at R-bin midpoint; dotted bounds are widths, not fit errors");text.DrawLatex(.14,.035,"Unselected data; TOF tag is a subset, not a second independent sample");}
}
}
#endif
void DrawHeliumRawCalibration(const char* input,const char* output,const char* settings,const char* curves,const char* knots,const char* rawKey,const char* tofKey) {
#ifndef __CINT__
  using namespace HeliumRawOverlay;
  try {
    const Settings cfg=ReadSettings(settings);const std::vector<Curve> c=ReadCurves(curves);const std::vector<Knot> k=ReadKnots(knots);
    const TString folder(output),maps=folder+"/raw_dedx_calibration_overlay.pdf",projections=folder+"/low_rigidity_raw_projections.pdf",shapes=folder+"/tof3_vs_all_shape_overlay.pdf";
    if(!gSystem->AccessPathName(maps)||!gSystem->AccessPathName(projections)||!gSystem->AccessPathName(shapes))throw std::runtime_error("Refusing to overwrite diagnostic PDFs");
    TFile* source=TFile::Open(input,"READ");if(!source||source->IsZombie())throw std::runtime_error("Unreadable raw ROOT");
    TH2* all=dynamic_cast<TH2*>(source->Get(rawKey));TH2* tag=dynamic_cast<TH2*>(source->Get(tofKey));
    if(!all||(Need(cfg,"requireTofHistogram")!=0&&!tag))throw std::runtime_error("Required histogram missing");AuditSubset(all,tag,(folder+"/subset_audit.json").Data());
    gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetPalette(1);
    TCanvas canvas("helium_raw_map","Helium raw diagnostics",1500,1100);canvas.Print(maps+"[");
    DrawMap(canvas,all,c,k,cfg,Need(cfg,"wideRigidityMin"),Need(cfg,"wideRigidityMax"),false,maps,folder+"/raw_wide.png");
    DrawMap(canvas,all,c,k,cfg,Need(cfg,"lowRigidityMin"),Need(cfg,"lowRigidityMax"),false,maps,folder+"/raw_low.png");
    if(tag){DrawMap(canvas,tag,c,k,cfg,Need(cfg,"wideRigidityMin"),Need(cfg,"wideRigidityMax"),true,maps,folder+"/tof3_wide.png");DrawMap(canvas,tag,c,k,cfg,Need(cfg,"lowRigidityMin"),Need(cfg,"lowRigidityMax"),true,maps,folder+"/tof3_low.png");}canvas.Print(maps+"]");
    std::ofstream rows((folder+"/projection_audit.tsv").Data());rows<<std::setprecision(17)
      <<"rigidity\tlo\thi\tall_regular_count\tall_flow_count\ttag_regular_count\ttag_flow_count\told_he3_mean\told_he3_sigma\tcurrent_he3_mean\tcurrent_he3_sigma\told_he4_mean\told_he4_sigma\tcurrent_he4_mean\tcurrent_he4_sigma\ttag_mode\ttag_window_count\ttag_window_mean\ttag_window_rms\ttag_mean_stat_error_proxy\ttag_rms_stat_error_proxy\ttag_window_low\ttag_window_high\told_he4_mean_bin_drift\tshape_core_low\tshape_core_high\tall_core_count\ttag_core_count\ttag_scale_to_all_core\tshape_high_low\tshape_high_high\tall_high_count\ttag_high_count\ttag_scaled_high_count\tcurrent_custom\n";
    TFile product((folder+"/raw_projections.root").Data(),"CREATE");if(product.IsZombie())throw std::runtime_error("Cannot create projections ROOT");
    TCanvas detail("helium_raw_slices","Raw and TOF projections",1600,1250);detail.Print(projections+"[");
    TCanvas shape("helium_shape_slices","TAG versus ALL shape",1700,850);shape.Print(shapes+"[");unsigned count=0,allSlices=0;
    for(int bin=1;bin<=all->GetNbinsX();++bin) {
      const double lo=all->GetXaxis()->GetBinLowEdge(bin),hi=all->GetXaxis()->GetBinUpEdge(bin),r=all->GetXaxis()->GetBinCenter(bin);
      if(lo<Need(cfg,"wideRigidityMin")-1.e-10||hi>Need(cfg,"wideRigidityMax")+1.e-10)continue;
      TH1D* raw=all->ProjectionY(Form("raw_native_R_%03d",int(1000*r+.5)),bin,bin,"e");raw->SetDirectory(0);
      TH1D* tagged=tag?tag->ProjectionY(Form("tof3_native_R_%03d",int(1000*r+.5)),bin,bin,"e"):0;if(tagged)tagged->SetDirectory(0);
      const double n=raw->Integral(1,raw->GetNbinsX()),flow=raw->Integral(0,raw->GetNbinsX()+1)-n,nt=tagged?tagged->Integral(1,tagged->GetNbinsX()):0.,ft=tagged?tagged->Integral(0,tagged->GetNbinsX()+1)-nt:0.;
      product.cd();raw->Write();if(tagged)tagged->Write();const int factor=int(Need(cfg,"projectionDedxRebinFactor"));if(raw->GetNbinsX()%factor)throw std::runtime_error("Projection rebin does not divide axis");
      raw->Rebin(factor);raw->SetName(Form("raw_display_R_%03d",int(1000*r+.5)));if(tagged){tagged->Rebin(factor);tagged->SetName(Form("tof3_display_R_%03d",int(1000*r+.5)));}
      if(raw->Integral(0,raw->GetNbinsX()+1)!=n+flow||(tagged&&tagged->Integral(0,tagged->GetNbinsX()+1)!=nt+ft))throw std::runtime_error("Projection rebin changed counts");raw->Write();if(tagged)tagged->Write();
      const Point old3=At(c,0,2,r),current3=At(c,1,2,r),old4=At(c,0,3,r),current4=At(c,1,3,r);const Moment m=TaggedMoment(tagged,old3,cfg);
      const double coreLo=m.mode-Need(cfg,"tofShapeNormalizeNSigma")*old3.sigma,coreHi=m.mode+Need(cfg,"tofShapeNormalizeNSigma")*old3.sigma;
      const double highLo=m.mode+Need(cfg,"tofShapeHighMinNSigma")*old3.sigma,highHi=m.mode+Need(cfg,"tofShapeHighMaxNSigma")*old3.sigma;
      const double allCore=IntervalCount(raw,coreLo,coreHi),tagCore=IntervalCount(tagged,coreLo,coreHi),scale=tagCore>0?allCore/tagCore:0,allHigh=IntervalCount(raw,highLo,highHi),tagHigh=IntervalCount(tagged,highLo,highHi);
      rows<<r<<'\t'<<lo<<'\t'<<hi<<'\t'<<n<<'\t'<<flow<<'\t'<<nt<<'\t'<<ft<<'\t'<<old3.mean<<'\t'<<old3.sigma<<'\t'<<current3.mean<<'\t'<<current3.sigma<<'\t'<<old4.mean<<'\t'<<old4.sigma<<'\t'<<current4.mean<<'\t'<<current4.sigma<<'\t'<<m.mode<<'\t'<<m.count<<'\t'<<m.mean<<'\t'<<m.sigma<<'\t'<<m.meanError<<'\t'<<m.rmsError<<'\t'<<m.low<<'\t'<<m.high<<'\t'<<At(c,0,3,hi).mean-At(c,0,3,lo).mean<<'\t'<<coreLo<<'\t'<<coreHi<<'\t'<<allCore<<'\t'<<tagCore<<'\t'<<scale<<'\t'<<highLo<<'\t'<<highHi<<'\t'<<allHigh<<'\t'<<tagHigh<<'\t'<<scale*tagHigh<<'\t'<<current4.custom<<'\n';++allSlices;
      if(lo>=Need(cfg,"projectionRigidityMin")-1.e-10&&hi<=Need(cfg,"projectionRigidityMax")+1.e-10) {
        detail.Clear();detail.Divide(2,2);detail.cd(1);DrawProjection(raw,0,c,cfg,r,false,true,m);detail.cd(2);DrawProjection(raw,0,c,cfg,r,false,false,m);
        if(tagged){detail.cd(3);DrawProjection(tagged,0,c,cfg,r,true,true,m);detail.cd(4);DrawProjection(tagged,0,c,cfg,r,true,false,m);}detail.Print(projections);detail.Print(folder+Form("/raw_tof3_R_%03d.png",int(1000*r+.5)));
        if(tagged){shape.Clear();shape.Divide(2,1);shape.cd(1);DrawProjection(raw,tagged,c,cfg,r,false,true,m,scale);shape.cd(2);DrawProjection(raw,tagged,c,cfg,r,false,false,m,scale);shape.Print(shapes);shape.Print(folder+Form("/all_tag_shape_R_%03d.png",int(1000*r+.5)));}++count;
      }delete raw;delete tagged;
    }
    detail.Print(projections+"]");shape.Print(shapes+"]");product.Close();source->Close();delete source;rows.close();if(!rows.good()||count==0)throw std::runtime_error("Empty/unwritable projection audit");
    std::cout<<"HELIUM_RAW_OVERLAY_COMPLETE maps="<<(tag?4:2)<<" projections="<<count<<" momentSlices="<<allSlices<<"; READ input, no fits or PID edits"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"HELIUM_RAW_OVERLAY_ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
