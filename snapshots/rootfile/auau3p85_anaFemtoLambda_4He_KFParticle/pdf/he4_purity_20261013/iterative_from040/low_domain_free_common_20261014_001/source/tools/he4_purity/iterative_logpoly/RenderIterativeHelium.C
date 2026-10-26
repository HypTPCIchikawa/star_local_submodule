// Saved histogram/model renderer only: READ phase ROOT, no fitting or PID writes.
#include "TFile.h"
#include "TH1D.h"
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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef __CINT__
namespace IterativeHeliumViews {
typedef std::map<std::string,std::string> Fields;
const char* components[] = {"total", "smooth", "t", "he3", "he4"};
const int colors[] = {kBlack, kMagenta+1, kOrange+7, kCyan+2, kGreen+2};
const int styles[] = {1, 2, 4, 7, 1};

std::string RequiredString(const Fields& fields,const std::string& key) {
  Fields::const_iterator i=fields.find(key);
  if(i==fields.end()||i->second.empty())throw std::runtime_error("Missing explicit field: "+key);
  return i->second;
}
double Number(const Fields& fields,const std::string& key) {
  const std::string value=RequiredString(fields,key);std::istringstream in(value);double x=0.;in>>x;
  in>>std::ws;if(!in.eof()||!TMath::Finite(x))throw std::runtime_error("Invalid numeric field: "+key);return x;
}
int Integer(const Fields& fields,const std::string& key) {
  const double x=Number(fields,key);if(std::fabs(x-std::floor(x+.5))>1.e-9)throw std::runtime_error("Noninteger field: "+key);return int(x);
}
bool Boolean(const Fields& fields,const std::string& key) {
  const std::string value=RequiredString(fields,key);
  if(value=="1"||value=="true"||value=="True")return true;
  if(value=="0"||value=="false"||value=="False")return false;
  throw std::runtime_error("Invalid explicit boolean: "+key);
}
void Require(bool condition,const std::string& message) {if(!condition)throw std::runtime_error(message);}
bool Near(double a,double b,double tolerance) {return std::fabs(a-b)<=tolerance*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
std::string Quote(const std::string& input) {
  std::ostringstream out;out<<'"';for(unsigned i=0;i<input.size();++i) {
    if(input[i]=='"'||input[i]=='\\')out<<'\\';
    if(input[i]=='\n')out<<"\\n";else if(input[i]=='\r')out<<"\\r";else out<<input[i];
  }out<<'"';return out.str();
}
std::vector<std::string> Split(const std::string& line) {
  std::vector<std::string> output;unsigned start=0;
  for(unsigned i=0;i<=line.size();++i)if(i==line.size()||line[i]=='\t') {output.push_back(line.substr(start,i-start));start=i+1;}
  return output;
}
Fields ReadSettings(const char* path) {
  std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read explicit renderer settings TSV");
  Fields out;std::string line;
  while(std::getline(in,line)) {
    if(!line.empty()&&line[line.size()-1]=='\r')line.erase(line.size()-1);
    if(line.empty())continue;const std::vector<std::string> columns=Split(line);
    Require(columns.size()==2&&out.insert(std::make_pair(columns[0],columns[1])).second,"Malformed/duplicate renderer settings TSV");
  }Require(in.eof(),"Renderer settings read failed");return out;
}
std::vector<Fields> ReadManifest(const char* path) {
  std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read explicit renderer manifest TSV");
  std::string line;Require(bool(std::getline(in,line)),"Empty manifest TSV");
  if(!line.empty()&&line[line.size()-1]=='\r')line.erase(line.size()-1);
  const std::vector<std::string> header=Split(line);std::set<std::string> unique(header.begin(),header.end());
  Require(unique.size()==header.size(),"Duplicate manifest column names");std::vector<Fields> rows;
  while(std::getline(in,line)) {
    if(!line.empty()&&line[line.size()-1]=='\r')line.erase(line.size()-1);
    if(line.empty())continue;const std::vector<std::string> values=Split(line);
    Require(values.size()==header.size(),"Manifest field count differs from header");Fields row;
    for(unsigned i=0;i<header.size();++i)row.insert(std::make_pair(header[i],values[i]));rows.push_back(row);
  }Require(in.eof()&&!rows.empty(),"Empty/unreadable manifest rows");return rows;
}
std::string HistogramName(const Fields& settings,const std::string& key,int index) {
  const std::string pattern=RequiredString(settings,key);const size_t position=pattern.find("%03d");
  Require(position!=std::string::npos&&pattern.find('%',position+1)==std::string::npos&&pattern.find('%')==position,"Histogram pattern must have exactly one %03d");
  return pattern.substr(0,position)+Form("%03d",index)+pattern.substr(position+4);
}
TH1D* Histogram(TFile* input,const std::string& name) {
  TH1D* h=dynamic_cast<TH1D*>(input->Get(name.c_str()));if(!h)throw std::runtime_error("Missing saved TH1D: "+name);return h;
}
void SameAxis(const TH1D* a,const TH1D* b,double tolerance) {
  Require(a->GetNbinsX()==b->GetNbinsX(),"Saved histogram bin counts differ");
  for(int i=1;i<=a->GetNbinsX()+1;++i)Require(Near(a->GetXaxis()->GetBinLowEdge(i),b->GetXaxis()->GetBinLowEdge(i),tolerance),"Saved raw model axes differ");
}

struct Slice {
  Fields fields;int index,firstFit,lastFit;double lo,hi,r,anchorMean,anchorSigma,fitLow,fitHigh;
  TH1D *native,*raw,*nsigma,*mask,*prediction[5];
  int firstDisplay,lastDisplay;double maximumClosureRelativeError;bool diagnostic;
  Slice():index(-1),firstFit(0),lastFit(0),lo(0),hi(0),r(0),anchorMean(0),anchorSigma(0),fitLow(0),fitHigh(0),native(0),raw(0),nsigma(0),mask(0),firstDisplay(0),lastDisplay(0),maximumClosureRelativeError(0),diagnostic(false) {for(int i=0;i<5;++i)prediction[i]=0;}
};
void ValidateSettings(const Fields& c) {
  Require(Integer(c,"schemaVersion")==1,"Unsupported renderer settings schema");
  Boolean(c,"allowSparseSlices");
  Require(Integer(c,"panelsPerPage")==4,"Renderer uses exactly four slices per linear page");
  Require(Integer(c,"expectedSlices")>0&&Integer(c,"dedxRebinFactor")>0,"Invalid explicit slice/rebin dimensions");
  Require(Number(c,"displayNSigmaMax")>Number(c,"displayNSigmaMin"),"Reversed display interval");
  Require(Number(c,"diagnosticRigidityMax")>Number(c,"diagnosticRigidityMin"),"Reversed diagnostic R interval");
  for(int d=0;d<3;++d) {
    const char* keys[]={"axisRelativeTolerance","componentRelativeTolerance","predictionNonnegativeTolerance"};Require(Number(c,keys[d])>0,"Invalid numeric audit tolerance");
  }
  Require(Number(c,"linearMaximumFactor")>1&&Number(c,"logMaximumFactor")>1&&Number(c,"logMinimum")>0,"Invalid y display factors");
  for(int i=0;i<4;++i) {const char* sides[]={"marginLeft","marginRight","marginBottom","marginTop"};Require(Number(c,sides[i])>0&&Number(c,sides[i])<.5,"Invalid pad margins");}
  Require(Number(c,"marginBottom")+Number(c,"marginTop")<.8,"Plotting area too small");
  Require(Number(c,"captionThirdY")>0&&Number(c,"captionThirdY")<Number(c,"captionSecondY")&&Number(c,"captionSecondY")<Number(c,"captionFirstY")&&Number(c,"captionFirstY")<Number(c,"marginBottom"),"Invalid caption order");
  Require(Number(c,"legendYMin")<Number(c,"legendYMax")&&Number(c,"legendYMax")<Number(c,"headerY")&&Number(c,"headerY")<1,"Invalid header/legend order");
}
Slice LoadSlice(TFile* input,const Fields& fields,const Fields& c) {
  Slice s;s.fields=fields;s.index=Integer(fields,"index");s.lo=Number(fields,"lo");s.hi=Number(fields,"hi");s.r=Number(fields,"r");
  s.anchorMean=Number(fields,"anchorMean");s.anchorSigma=Number(fields,"anchorSigma");s.fitLow=Number(fields,"fitLoDedx");s.fitHigh=Number(fields,"fitHiDedx");s.firstFit=Integer(fields,"firstFitBin");s.lastFit=Integer(fields,"lastFitBin");
  Require(s.index>=0&&s.hi>s.lo&&Near(s.r,(s.lo+s.hi)/2,Number(c,"axisRelativeTolerance"))&&s.anchorSigma>0&&s.fitHigh>s.fitLow,"Invalid slice coordinates");
  const std::string mode=RequiredString(fields,"stageMode");Require(mode=="local_free"||mode=="regularized"||mode=="shared_fixed","Unknown stageMode");
  const std::string model=RequiredString(fields,"model");Require(model=="linear"||model=="exponential","Unknown background model");
  const char* integerFields[]={"valid","quality","status","covarianceStatus","atBoundary","ndf","rawLocalMaximum","oldAccepted"};for(unsigned i=0;i<sizeof(integerFields)/sizeof(integerFields[0]);++i)Integer(fields,integerFields[i]);
  Require(Number(fields,"chi2")>=0||Number(fields,"chi2")==-1.,"Invalid chi2 (use explicit -1 when unrecorded)");
  const char* floatFields[]={"deviance","penalty"};for(unsigned i=0;i<sizeof(floatFields)/sizeof(floatFields[0]);++i)Require(Number(fields,floatFields[i])>=0,"Negative data statistic/penalty");
  const char* triStateFields[]={"valid","quality","atBoundary","rawLocalMaximum"};for(unsigned i=0;i<sizeof(triStateFields)/sizeof(triStateFields[0]);++i)Require(Integer(fields,triStateFields[i])>=-1&&Integer(fields,triStateFields[i])<=1,"Tri-state metadata must be -1 (unrecorded), 0 or 1");
  const char* states[]={"tShapeState","he3ShapeState","he4ShapeState"};for(int i=0;i<3;++i) {const std::string state=RequiredString(fields,states[i]);Require(state=="free"||state=="regularized"||state=="shared_fixed"||state=="fixed"||state=="inactive"||state=="unidentified","Unknown component shape state");}
  s.native=Histogram(input,HistogramName(c,"nativeHistogramPattern",s.index));s.raw=Histogram(input,HistogramName(c,"rawHistogramPattern",s.index));s.nsigma=Histogram(input,HistogramName(c,"nsigmaHistogramPattern",s.index));s.mask=Histogram(input,HistogramName(c,"fitMaskPattern",s.index));
  for(int i=0;i<5;++i)s.prediction[i]=Histogram(input,HistogramName(c,std::string("prediction_")+components[i]+"Pattern",s.index));
  SameAxis(s.raw,s.mask,Number(c,"axisRelativeTolerance"));for(int i=0;i<5;++i)SameAxis(s.raw,s.prediction[i],Number(c,"axisRelativeTolerance"));
  Require(s.raw->GetNbinsX()==s.nsigma->GetNbinsX(),"Observed raw/nSigma bin counts differ");
  Require(s.firstFit>=1&&s.lastFit>=s.firstFit&&s.lastFit<=s.raw->GetNbinsX(),"Invalid used fit bin span");
  const int factor=Integer(c,"dedxRebinFactor");Require(s.native->GetNbinsX()==s.raw->GetNbinsX()*factor,"Native/fit rebin dimensions differ");
  for(int b=1;b<=s.raw->GetNbinsX()+1;++b) {
    Require(Near(s.raw->GetXaxis()->GetBinLowEdge(b),s.native->GetXaxis()->GetBinLowEdge(1+(b-1)*factor),Number(c,"axisRelativeTolerance")),"Native/fit rebin edges differ");
    Require(Near(s.nsigma->GetXaxis()->GetBinLowEdge(b),(s.raw->GetXaxis()->GetBinLowEdge(b)-s.anchorMean)/s.anchorSigma,Number(c,"axisRelativeTolerance")),"nSigma does not use the fixed saved anchor");
  }
  for(int b=0;b<=s.raw->GetNbinsX()+1;++b) {
    const double raw=s.raw->GetBinContent(b),error=s.raw->GetBinError(b);
    Require(TMath::Finite(raw)&&TMath::Finite(error)&&raw>=0&&error>=0,"Invalid observed count/error");
    Require(raw==s.nsigma->GetBinContent(b)&&error==s.nsigma->GetBinError(b),"raw/nSigma observed count/error/flow differs");
    double original=0,variance=0;
    if(b==0||b==s.raw->GetNbinsX()+1) {const int i=b==0?0:s.native->GetNbinsX()+1;original=s.native->GetBinContent(i);variance=std::pow(s.native->GetBinError(i),2);}
    else for(int i=1+(b-1)*factor;i<=b*factor;++i) {original+=s.native->GetBinContent(i);variance+=std::pow(s.native->GetBinError(i),2);}
    // These are saved unweighted integer counts, so no relative count tolerance
    // may silently absorb a one-track change in a very populated bin.
    Require(raw==original&&Near(error*error,variance,Number(c,"axisRelativeTolerance")),"Native rebin count/variance/flow differs");
    const double mask=s.mask->GetBinContent(b);Require(mask==double(b>=s.firstFit&&b<=s.lastFit),"Saved fit mask differs from declared contiguous used bins");
    double sum=0.;for(int component=1;component<5;++component)sum+=s.prediction[component]->GetBinContent(b);
    const double total=s.prediction[0]->GetBinContent(b),relative=std::fabs(sum-total)/std::max(1.,std::max(std::fabs(total),std::fabs(sum)));
    s.maximumClosureRelativeError=std::max(s.maximumClosureRelativeError,relative);Require(relative<=Number(c,"componentRelativeTolerance"),"Saved component bin-integral sum does not close");
    for(int component=0;component<5;++component) {const double value=s.prediction[component]->GetBinContent(b);Require(TMath::Finite(value)&&value>=-Number(c,"predictionNonnegativeTolerance"),"Nonfinite/negative saved bin expectation");}
  }
  const double low=s.anchorMean+s.anchorSigma*Number(c,"displayNSigmaMin"),high=s.anchorMean+s.anchorSigma*Number(c,"displayNSigmaMax");
  Require(high>s.raw->GetXaxis()->GetXmin()&&low<s.raw->GetXaxis()->GetXmax(),"Display outside stored data axis");
  TAxis range(*s.raw->GetXaxis());range.SetRangeUser(std::max(low,range.GetXmin()),std::min(high,range.GetXmax()));
  s.firstDisplay=std::max(1,range.GetFirst());s.lastDisplay=std::min(s.raw->GetNbinsX(),range.GetLast());
  Require(s.lastDisplay>=s.firstDisplay,"Empty normal-bin display interval");
  s.diagnostic=s.lo>=Number(c,"diagnosticRigidityMin")-Number(c,"axisRelativeTolerance")&&s.hi<=Number(c,"diagnosticRigidityMax")+Number(c,"axisRelativeTolerance");
  return s;
}

std::string ShapeName(const Slice& s,int component) {
  if(component==0)return "Total bin-integrated expectation";
  if(component==1)return "Positive smooth background";
  const char* labels[]={"", "", "T", "^{3}He", "^{4}He"};
  const char* keys[]={"", "", "tShapeState", "he3ShapeState", "he4ShapeState"};
  return std::string(labels[component])+" ("+RequiredString(s.fields,keys[component])+")";
}
void Draw(const Slice& s,const Fields& c,bool physical,bool logarithmic,const char* phaseLabel) {
  gPad->SetLeftMargin(Number(c,"marginLeft"));gPad->SetRightMargin(Number(c,"marginRight"));gPad->SetBottomMargin(Number(c,"marginBottom"));gPad->SetTopMargin(Number(c,"marginTop"));gPad->SetLogy(logarithmic);
  TH1D* source=physical?s.raw:s.nsigma;
  TH1D* data=static_cast<TH1D*>(source->Clone(Form("view_%d_%d_%d",s.index,int(physical),int(logarithmic))));data->SetDirectory(0);
  // Both axes select exactly the same normal-bin indices, never a rounded overflow.
  data->GetXaxis()->SetRange(s.firstDisplay,s.lastDisplay);double maximum=1.;
  for(int b=s.firstDisplay;b<=s.lastDisplay;++b)maximum=std::max(maximum,data->GetBinContent(b)+data->GetBinError(b));
  const double ymin=logarithmic?Number(c,"logMinimum"):0.,ymax=maximum*Number(c,logarithmic?"logMaximumFactor":"linearMaximumFactor");
  data->SetMinimum(ymin);data->SetMaximum(ymax);data->SetTitle(Form(";%s;Tracks / %.2f keV/cm fit bin",physical?"TPC dE/dx (keV/cm)":"Legacy n#sigma_{^{4}He} (fixed anchor)",s.raw->GetBinWidth(1)));
  data->SetMarkerStyle(20);data->SetMarkerSize(Number(c,"markerSize"));data->SetLineColor(kBlack);
  data->GetXaxis()->SetTitleSize(Number(c,"axisTitleSize"));data->GetYaxis()->SetTitleSize(Number(c,"axisTitleSize"));data->GetXaxis()->SetLabelSize(Number(c,"axisLabelSize"));data->GetYaxis()->SetLabelSize(Number(c,"axisLabelSize"));
  data->GetXaxis()->SetTitleOffset(Number(c,"xTitleOffset"));data->GetYaxis()->SetTitleOffset(Number(c,"yTitleOffset"));
  data->GetXaxis()->SetTitleFont(Integer(c,"axisFont"));data->GetYaxis()->SetTitleFont(Integer(c,"axisFont"));data->GetXaxis()->SetLabelFont(Integer(c,"axisFont"));data->GetYaxis()->SetLabelFont(Integer(c,"axisFont"));data->Draw("E1");
  TLegend* legend=new TLegend(Number(c,"legendXMin"),Number(c,"legendYMin"),Number(c,"legendXMax"),Number(c,"legendYMax"));legend->SetNColumns(Integer(c,"legendColumns"));legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextFont(Integer(c,"axisFont"));legend->SetTextSize(Number(c,"legendTextSize"));legend->AddEntry(data,"Saved pre-PID tracks","lep");
  for(int component=0;component<5;++component) {
    TGraph* graph=new TGraph();int point=0;
    for(int b=std::max(s.firstFit,s.firstDisplay);b<=std::min(s.lastFit,s.lastDisplay);++b)graph->SetPoint(point++,source->GetBinCenter(b),s.prediction[component]->GetBinContent(b));
    graph->SetLineColor(colors[component]);graph->SetLineStyle(styles[component]);graph->SetLineWidth(component==0?3:2);if(point>0)graph->Draw("L same");legend->AddEntry(graph,ShapeName(s,component).c_str(),"l");
  }
  const double actualLow=s.raw->GetXaxis()->GetBinLowEdge(s.firstFit),actualHigh=s.raw->GetXaxis()->GetBinUpEdge(s.lastFit);
  const double displayLow=source->GetXaxis()->GetBinLowEdge(s.firstDisplay),displayHigh=source->GetXaxis()->GetBinUpEdge(s.lastDisplay);
  for(int side=0;side<2;++side) {
    const double raw=side?actualHigh:actualLow,x=physical?raw:(raw-s.anchorMean)/s.anchorSigma;
    if(x>=displayLow&&x<=displayHigh) {TLine* edge=new TLine(x,ymin,x,ymax);edge->SetLineColor(kGray+2);edge->SetLineStyle(3);edge->SetLineWidth(2);edge->Draw();}
  }
  legend->Draw();TLatex text;text.SetNDC();text.SetTextFont(Integer(c,"axisFont"));text.SetTextSize(Number(c,"headerTextSize"));
  text.DrawLatex(Number(c,"marginLeft"),Number(c,"headerY"),Form("R %.3f-%.3f | %s | %s",s.lo,s.hi,physical?"raw dE/dx":"legacy n#sigma",logarithmic?"log":"linear"));
  text.SetTextSize(Number(c,"statusTextSize"));const int valid=Integer(s.fields,"valid"),quality=Integer(s.fields,"quality"),boundary=Integer(s.fields,"atBoundary");
  text.SetTextColor(valid==1&&quality==1&&boundary==0?kBlack:kRed+1);
  text.DrawLatex(Number(c,"marginLeft"),Number(c,"statusFirstY"),Form("%s: %s | %s%s%s",phaseLabel,RequiredString(s.fields,"stageMode").c_str(),valid==1?"fit finite":valid==0?"FIT FAILED":"VALIDITY UNTESTED",quality==1?"":quality==0?" / POOR OR REFERENCE":" / QUALITY UNTESTED",boundary==1?" / BOUNDARY":boundary==0?"":" / BOUNDARY UNTESTED"));
  text.SetTextColor(kBlack);const int ndf=Integer(s.fields,"ndf");
  const std::string ratio=ndf>0&&Number(s.fields,"chi2")>=0?std::string(Form("%.2f",Number(s.fields,"chi2")/ndf)):std::string("unavailable");
  text.DrawLatex(Number(c,"marginLeft"),Number(c,"statusSecondY"),Form("Raw #chi^{2}/ndf %s; deviance %.3g; penalty %.3g; status %d, cov %d",ratio.c_str(),Number(s.fields,"deviance"),Number(s.fields,"penalty"),Integer(s.fields,"status"),Integer(s.fields,"covarianceStatus")));
  text.SetTextSize(Number(c,"captionTextSize"));
  text.DrawLatex(Number(c,"marginLeft"),Number(c,"captionFirstY"),Form("Actual fit bins: %.2f-%.2f keV/cm; requested support: %.2f-%.2f",actualLow,actualHigh,s.fitLow,s.fitHigh));
  text.DrawLatex(Number(c,"marginLeft"),Number(c,"captionSecondY"),"Same observed bin counts/errors; curves: saved bin expectations, NO density/Jacobian rescaling");
  const int rawPeak=Integer(s.fields,"rawLocalMaximum");text.SetTextColor(rawPeak==1?kBlack:kRed+1);
  text.DrawLatex(Number(c,"marginLeft"),Number(c,"captionThirdY"),Form("%s; seed %s; gray lines: fit bins, NOT PID cuts",rawPeak==1?"Local-peak heuristic present; isotope labels remain model-dependent":rawPeak==0?"NO SEPARATE RAW PEAK: labels are model-dependent":"Local-peak diagnostic unrecorded; labels are model-dependent",Integer(s.fields,"oldAccepted")?"accepted":"unaccepted / reference"));text.SetTextColor(kBlack);
}
void WriteAudit(const std::vector<Slice>& slices,const Fields& c,const char* output,const char* input,const char* phase) {
  std::ofstream f((std::string(output)+"/render_audit.json").c_str());Require(f.good(),"Cannot create render audit");f<<std::setprecision(17);
  unsigned diagnostics=0;double maxClosure=0;
  for(unsigned i=0;i<slices.size();++i){diagnostics+=slices[i].diagnostic?1:0;maxClosure=std::max(maxClosure,slices[i].maximumClosureRelativeError);}
  f<<"{\"schemaVersion\":1,\"phase\":\"complete\",\"inputPhaseRoot\":"<<Quote(input)<<",\"phaseLabel\":"<<Quote(phase)<<",\"noRefit\":true,\"productionChanged\":false,\"slices\":"<<slices.size()<<",\"nsigmaLinearPages\":"<<(slices.size()+3)/4<<",\"physicalLinearPages\":"<<(slices.size()+3)/4<<",\"lowRigidityPages\":"<<diagnostics<<",\"nSigmaAnchor\":\"original legacy He4 at saved R center\",\"physicalYConvention\":\"counts per saved physical fit bin; NO Jacobian\",\"nativeRebinCountsVarianceFlowsPreserved\":true,\"rawNSigmaCountsErrorsFlowsEqual\":true,\"sameVisibleBinIndices\":true,\"normalBinsOnly\":true,\"componentClosureChecked\":true,\"maximumComponentClosureRelativeError\":"<<maxClosure<<",\"slicesAudit\":[";
  for(unsigned i=0;i<slices.size();++i) {
    const Slice& s=slices[i];if(i)f<<',';
    f<<"{\"index\":"<<s.index<<",\"rigidity\":"<<s.r<<",\"diagnostic\":"<<(s.diagnostic?"true":"false")<<",\"visibleBins\":["<<s.firstDisplay<<','<<s.lastDisplay<<"],\"visiblePhysicalRange\":["<<s.raw->GetXaxis()->GetBinLowEdge(s.firstDisplay)<<','<<s.raw->GetXaxis()->GetBinUpEdge(s.lastDisplay)<<"],\"visibleNSigmaRange\":["<<s.nsigma->GetXaxis()->GetBinLowEdge(s.firstDisplay)<<','<<s.nsigma->GetXaxis()->GetBinUpEdge(s.lastDisplay)<<"],\"visibleCounts\":"<<s.raw->Integral(s.firstDisplay,s.lastDisplay)<<",\"allNormalCount\":"<<s.raw->Integral(1,s.raw->GetNbinsX())<<",\"underflow\":"<<s.raw->GetBinContent(0)<<",\"overflow\":"<<s.raw->GetBinContent(s.raw->GetNbinsX()+1)<<",\"usedFitBins\":["<<s.firstFit<<','<<s.lastFit<<"],\"actualFitPhysicalRange\":["<<s.raw->GetXaxis()->GetBinLowEdge(s.firstFit)<<','<<s.raw->GetXaxis()->GetBinUpEdge(s.lastFit)<<"],\"valid\":"<<Integer(s.fields,"valid")<<",\"quality\":"<<Integer(s.fields,"quality")<<",\"atBoundary\":"<<Integer(s.fields,"atBoundary")<<",\"rawLocalMaximum\":"<<Integer(s.fields,"rawLocalMaximum")<<",\"maximumClosureRelativeError\":"<<s.maximumClosureRelativeError<<'}';
  }f<<"]}\n";f.close();Require(f.good(),"Render audit write failed");
}
}
#endif

void RenderIterativeHelium(const char* inputPhaseRoot,const char* manifestTsv,const char* settingsTsv,const char* outputDirectory,const char* phaseLabel) {
#ifndef __CINT__
  using namespace IterativeHeliumViews;
  try {
    const Fields c=ReadSettings(settingsTsv);ValidateSettings(c);const std::vector<Fields> rows=ReadManifest(manifestTsv);
    Require(rows.size()==unsigned(Integer(c,"expectedSlices")),"Manifest slice count differs from explicit expectedSlices");
    const TString out(outputDirectory),nsigma=out+"/fits_nsigma_linear.pdf",physical=out+"/fits_raw_dedx_linear.pdf",low=out+"/low_rigidity_components.pdf";
    const TString products[]={nsigma,physical,low,out+"/render_audit.json"};
    for(unsigned i=0;i<sizeof(products)/sizeof(products[0]);++i)Require(gSystem->AccessPathName(products[i]),"Refusing to overwrite an existing renderer product");
    Require(!gSystem->AccessPathName(outputDirectory),"Renderer output directory must be explicitly created by its launcher");
    TFile* input=TFile::Open(inputPhaseRoot,"READ");Require(input&&!input->IsZombie(),"Cannot read saved iterative-fit ROOT");std::vector<Slice> slices;
    for(unsigned i=0;i<rows.size();++i) {
      Slice s=LoadSlice(input,rows[i],c);
      if(!Boolean(c,"allowSparseSlices"))Require(s.index==int(i),"Complete-manifest indices must start at zero and be contiguous");
      if(i) {
        Require(s.index>slices.back().index,"Manifest indices must be unique and strictly increasing");
        if(Boolean(c,"allowSparseSlices"))Require(s.lo>=slices.back().hi-Number(c,"axisRelativeTolerance"),"Sparse R intervals overlap or are not ascending");
        else Require(Near(s.lo,slices.back().hi,Number(c,"axisRelativeTolerance")),"Manifest R intervals are not contiguous");
      }slices.push_back(s);
    }
    gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetPaperSize(Number(c,"paperWidthCm"),Number(c,"paperHeightCm"));
    TCanvas pages("iterative_helium_linear","Saved iterative helium fits",Integer(c,"linearCanvasWidth"),Integer(c,"linearCanvasHeight"));
    for(int axis=0;axis<2;++axis) {
      const TString book=axis?physical:nsigma;pages.Print(book+"[");
      for(unsigned start=0;start<slices.size();start+=4) {
        pages.Clear();pages.Divide(2,2);for(unsigned j=0;j<4&&start+j<slices.size();++j){pages.cd(j+1);Draw(slices[start+j],c,axis!=0,false,phaseLabel);}pages.Print(book);
        if(start==0)pages.Print(out+(axis?"/fits_raw_dedx_linear_first.png":"/fits_nsigma_linear_first.png"));
      }pages.Print(book+"]");
    }
    TCanvas diagnostics("iterative_helium_low","Same-bin nSigma/raw log/linear",Integer(c,"diagnosticCanvasWidth"),Integer(c,"diagnosticCanvasHeight"));diagnostics.Print(low+"[");
    for(unsigned i=0;i<slices.size();++i)if(slices[i].diagnostic) {
      diagnostics.Clear();diagnostics.Divide(2,2);diagnostics.cd(1);Draw(slices[i],c,false,true,phaseLabel);diagnostics.cd(2);Draw(slices[i],c,false,false,phaseLabel);diagnostics.cd(3);Draw(slices[i],c,true,true,phaseLabel);diagnostics.cd(4);Draw(slices[i],c,true,false,phaseLabel);
      diagnostics.Print(low);diagnostics.Print(out+Form("/low_components_R_%03d.png",int(slices[i].r*1000+.5)));
    }diagnostics.Print(low+"]");WriteAudit(slices,c,outputDirectory,inputPhaseRoot,phaseLabel);input->Close();delete input;
    std::cout<<"ITERATIVE_HELIUM_RENDER_COMPLETE slices="<<slices.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"ITERATIVE_HELIUM_RENDER_ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
