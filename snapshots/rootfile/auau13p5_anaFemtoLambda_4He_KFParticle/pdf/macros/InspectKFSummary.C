// ROOT 5 / ACLiC: inspect an existing merged Lambda--4He ROOT, never rewrite it.
// Inclusive purity is a proxy, NOT a measurement of pair- or k*-dependent purity.
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TParameter.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace KFSummaryInspect {
std::string Quote(const std::string &s) {
  std::ostringstream o; o << '"';
  for (std::string::const_iterator i=s.begin(); i!=s.end(); ++i) {
    if (*i=='"' || *i=='\\') o << '\\' << *i;
    else if (*i=='\n') o << "\\n";
    else if (*i=='\r') o << "\\r";
    else if (*i=='\t') o << "\\t";
    else o << *i;
  }
  o << '"'; return o.str();
}
void Num(std::ostream &o, double x) {
  if (TMath::Finite(x)) o << std::setprecision(17) << x;
  else o << "null";
}
struct Window {
  int first, last;
  double lower, upper, sum, error;
  Window():first(0),last(0),lower(0),upper(0),sum(0),error(0) {}
};
// Use centers to select complete bins; all adjacent windows share the same cut.
// No fractional-bin assumptions and no overlap/double counting.
int FirstCenterAtLeast(TH1 *h, double edge) {
  int b=1;
  while (b<=h->GetNbinsX() && h->GetBinCenter(b)<edge) ++b;
  return b;
}
Window MakeWindow(TH1 *h, double low, double high) {
  Window w;
  w.first=FirstCenterAtLeast(h,low);
  w.last=FirstCenterAtLeast(h,high)-1;
  if (w.first<1 || w.last>h->GetNbinsX() || w.last<w.first) return w;
  w.lower=h->GetXaxis()->GetBinLowEdge(w.first);
  w.upper=h->GetXaxis()->GetBinUpEdge(w.last);
  w.sum=h->IntegralAndError(w.first,w.last,w.error);
  return w;
}
void WriteWindow(std::ostream &o,const Window &w) {
  o << "{\"firstBin\":" << w.first << ",\"lastBin\":" << w.last
    << ",\"actualLower\":"; Num(o,w.lower);
  o << ",\"actualUpper\":"; Num(o,w.upper);
  o << ",\"count\":"; Num(o,w.sum);
  o << ",\"error\":"; Num(o,w.error); o << '}';
}
TH1 *GetHistogram(TFile *f,const char *name) {
  TObject *o=f->Get(name);
  if (!o && std::string(name)=="hLambda_InvMass") o=f->Get("lambda/hLambda_InvMass");
  return dynamic_cast<TH1*>(o);
}
void AddNames(std::vector<std::string> &names) {
  const char *base[]={"hLambda_InvMass","hN","hVz","hRefMult","hRawMult",
    "hCentrality","hCentralityRaw","hCentrality16","hRefMultCorr","hRefMultWeight",
    "hLambda_PtVsYLab_signal","hLambdaProton_PtVsYLab_signal",
    "hLambdaPion_PtVsYLab_signal","hNucleus_PtVsYLab_he4",
    "true/hDedxP","true/hDedxP_4He","true/hM2P"};
  for (unsigned i=0;i<sizeof(base)/sizeof(base[0]);++i) names.push_back(base[i]);
  for (int c=0;c<9;++c) {
    names.push_back(Form("hLambda_InvMass_CentBin%d",c));
    for(int mix=0;mix<2;++mix) {
      const char *d=mix?"mix":"true"; const char *m=mix?"Mixed_":"";
      names.push_back(Form("%s/hKstar_%s4He_CentBin%d",d,m,c));
      names.push_back(Form("%s/hKstar_%s4He_SBNeg_CentBin%d",d,m,c));
      names.push_back(Form("%s/hKstar_%s4He_SBPos_CentBin%d",d,m,c));
      names.push_back(Form("%s/hKstarMass_%s4He_CentBin%d",d,m,c));
      names.push_back(Form("%s/hKstarMt_%s4He_CentBin%d",d,m,c));
    }
  }
}
void Fail(const std::string &s) { std::cerr << "InspectKFSummary: " << s << std::endl; gSystem->Exit(2); }
}

void InspectKFSummary(const char *input_root, const char *output_dir,
                     double mean=1.11596, double sigma=.00296,
                     double nsigma=3., double outer=4.,
                     double fit_min=1.08, double fit_max=1.15) {
  using namespace KFSummaryInspect;
  if (!input_root || !*input_root || !output_dir || !*output_dir)
    { Fail("input_root and output_dir are required"); return; }
  if (!(sigma>0 && nsigma>0 && outer>1 && fit_min<fit_max))
    { Fail("invalid explicit mass-window / fit arguments"); return; }
  const TString out(output_dir);
  const TString jsonPath=out+"/inspect.json";
  const TString pdfPath=out+"/Lambda_InvMass.pdf";
  const TString pngPath=out+"/Lambda_InvMass.png";
  if (!gSystem->AccessPathName(jsonPath) || !gSystem->AccessPathName(pdfPath) || !gSystem->AccessPathName(pngPath))
    { Fail("refusing to overwrite existing inspect.json / Lambda_InvMass outputs"); return; }
  if (gSystem->AccessPathName(out) && gSystem->mkdir(out,kTRUE)!=0)
    { Fail("cannot create output directory"); return; }
  TFile *f=TFile::Open(input_root,"READ");
  if (!f || f->IsZombie()) { Fail("cannot open input ROOT"); return; }
  TH1 *h=GetHistogram(f,"hLambda_InvMass");
  if (!h || h->GetDimension()!=1 || h->Integral()<=0)
    { Fail("missing or empty one-dimensional hLambda_InvMass"); return; }
  const double low=mean-nsigma*sigma, high=mean+nsigma*sigma;
  const double sbLow=mean-outer*nsigma*sigma, sbHigh=mean+outer*nsigma*sigma;
  if (sbLow<h->GetXaxis()->GetXmin() || sbHigh>h->GetXaxis()->GetXmax())
    { Fail("requested sidebands exceed stored mass axis"); return; }
  const Window signal=MakeWindow(h,low,high);
  const Window left=MakeWindow(h,sbLow,low), right=MakeWindow(h,high,sbHigh);
  const double sigWidth=signal.upper-signal.lower;
  const double sbWidth=left.upper-left.lower+right.upper-right.lower;
  if (!(sigWidth>0 && sbWidth>0 && signal.sum>0))
    { Fail("empty signal or sideband window"); return; }
  const double alpha=sigWidth/sbWidth;
  const double background=alpha*(left.sum+right.sum);
  const double backgroundError=alpha*std::sqrt(left.error*left.error+right.error*right.error);
  const double purity=1.-background/signal.sum;
  const double purityError=std::sqrt(std::pow(backgroundError/signal.sum,2)+
                         std::pow(background*signal.error/(signal.sum*signal.sum),2));
  const bool purityUsable=TMath::Finite(purity) && purity>0 && purity<1;

  // Compare with gaus+pol2, but retain the transparent sideband estimate as the
  // primary inclusive proxy. No fit-derived change to the production mass gate.
  const double bgPerBin=(left.sum+right.sum)/(left.last-left.first+right.last-right.first+2);
  TF1 fit("kfSummaryInclusiveMassFit","gaus(0)+pol2(3)",fit_min,fit_max);
  fit.SetParameters(TMath::Max(1.,h->GetMaximum()-bgPerBin),mean,sigma,bgPerBin,0,0);
  fit.SetParLimits(0,0.,10.*h->GetMaximum());
  fit.SetParLimits(1,mean-.008,mean+.008);
  fit.SetParLimits(2,.0005,.010);
  TFitResultPtr result=h->Fit(&fit,"SQNR0","",fit_min,fit_max);
  int fitStatus=int(result), covStatus=result.Get()?result->CovMatrixStatus():-1;
  bool fitValid=result.Get() && result->IsValid() && fitStatus==0 && covStatus>=2;
  double fitBackground=0, fitBackgroundError=0, fitPurity=0;
  TF1 bg("kfSummaryInclusiveBackground","pol2",fit_min,fit_max);
  bg.SetParameters(fit.GetParameter(3),fit.GetParameter(4),fit.GetParameter(5));
  if (fitValid) {
    TMatrixDSym covariance=result->GetCovarianceMatrix();
    double covBg[9];
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) covBg[3*i+j]=covariance(i+3,j+3);
    fitBackground=bg.Integral(signal.lower,signal.upper)/h->GetBinWidth(signal.first);
    fitBackgroundError=bg.IntegralError(signal.lower,signal.upper,bg.GetParameters(),covBg)/h->GetBinWidth(signal.first);
    fitPurity=1.-fitBackground/signal.sum;
    fitValid=TMath::Finite(fitBackground) && TMath::Finite(fitBackgroundError) &&
             fitBackground>=0 && fitPurity>0 && fitPurity<1;
  }

  std::ofstream o(jsonPath.Data());
  if(!o) { Fail("cannot open JSON output"); return; }
  o << "{\n\"schemaVersion\":1,\n\"purity_for_cfsub\":";
  if(purityUsable) Num(o,purity); else o<<"null";
  o << ",\n\"inputRoot\":" << Quote(input_root)
    << ",\n\"histogram\":\"hLambda_InvMass\",\n\"massWindow\":{\"mean\":"; Num(o,mean);
  o << ",\"sigma\":"; Num(o,sigma); o << ",\"nSigma\":"; Num(o,nsigma);
  o << ",\"sidebandOuterFactor\":"; Num(o,outer);
  o << ",\"nominalSignalLower\":"; Num(o,low); o << ",\"nominalSignalUpper\":"; Num(o,high);
  o << ",\"nominalSidebandLower\":"; Num(o,sbLow); o << ",\"nominalSidebandUpper\":"; Num(o,sbHigh);
  o << ",\"wholeBinPolicy\":\"center in [lower,upper), no overlap; actual bin edges below\",\"signal\":";
  WriteWindow(o,signal); o << ",\"left\":"; WriteWindow(o,left); o << ",\"right\":"; WriteWindow(o,right); o << "},\n";
  o << "\"inclusivePurityProxy\":{\"method\":\"width_scaled_two_sidebands\",\"usableForCFsub\":" << (purityUsable?"true":"false")
    << ",\"purity\":"; Num(o,purity); o << ",\"purityError\":"; Num(o,purityError);
  o << ",\"background\":"; Num(o,background); o << ",\"backgroundError\":"; Num(o,backgroundError);
  o << ",\"signalEstimate\":"; Num(o,signal.sum-background);
  o << ",\"signalEstimateError\":"; Num(o,std::sqrt(signal.error*signal.error+backgroundError*backgroundError));
  o << ",\"signalToBackground\":"; if(background>0) Num(o,(signal.sum-background)/background); else o<<"null";
  o << ",\"actualWidthScale\":"; Num(o,alpha);
  o << ",\"nominalWidthScale\":"; Num(o,1./(outer-1.));
  o << ",\"caveat\":\"Inclusive candidate proxy; not pair/kstar/centrality-dependent purity. Flat-background interpolation over sidebands. Binned window approximates, not exactly reproduces, the unbinned Maker mass gate. Statistical errors exclude model systematics.\"},\n";
  o << "\"gausPol2Comparison\":{\"status\":" << fitStatus << ",\"covarianceStatus\":" << covStatus
    << ",\"usable\":false,\"numericallyValid\":" << (fitValid?"true":"false") << ",\"fitMin\":"; Num(o,fit_min);
  o << ",\"fitMax\":"; Num(o,fit_max); o << ",\"mean\":"; Num(o,fit.GetParameter(1));
  o << ",\"sigma\":"; Num(o,fit.GetParameter(2)); o << ",\"chi2\":"; Num(o,fit.GetChisquare());
  o << ",\"ndf\":" << fit.GetNDF() << ",\"chi2Probability\":";
  if(fit.GetNDF()>0) Num(o,TMath::Prob(fit.GetChisquare(),fit.GetNDF())); else o<<"null";
  o << ",\"chi2PerNdf\":";
  if(fit.GetNDF()>0) Num(o,fit.GetChisquare()/fit.GetNDF()); else o<<"null";
  o << ",\"goodnessThresholdApplied\":false,\"backgroundInBinnedSignal\":";
  if(fitValid) Num(o,fitBackground); else o<<"null";
  o << ",\"backgroundError\":"; if(fitValid) Num(o,fitBackgroundError); else o<<"null";
  o << ",\"purity\":"; if(fitValid) Num(o,fitPurity); else o<<"null";
  o << ",\"primaryEstimator\":false,\"caveat\":\"Comparison only, never used for CF correction. Numerical convergence does not establish a good model fit. Inspect chi2/ndf and chi2Probability; no hidden goodness threshold applied.\"},\n\"eventCounters\":{";
  const char *counters[]={"inputEvents","reconstructedEvents","selectedLambdaCandidates","mixingVzClampedEvents"};
  for(unsigned i=0;i<4;++i) {
    if(i) o<<','; o<<Quote(counters[i])<<':';
    TParameter<Long64_t> *p=dynamic_cast<TParameter<Long64_t>*>(f->Get(counters[i]));
    if(p) o<<p->GetVal(); else o<<"null";
  }
  o << "},\n\"histograms\":[\n";
  std::vector<std::string> names; AddNames(names);
  for(unsigned i=0;i<names.size();++i) {
    if(i) o<<",\n";
    TH1 *v=GetHistogram(f,names[i].c_str());
    o << "{\"name\":"<<Quote(names[i])<<",\"present\":"<<(v?"true":"false");
    if(v) {
      o<<",\"class\":"<<Quote(v->ClassName())<<",\"dimension\":"<<v->GetDimension()<<",\"entries\":"; Num(o,v->GetEntries());
      o<<",\"sumWeights\":"; Num(o,v->GetSumOfWeights());
      o<<",\"xBins\":"<<v->GetNbinsX()<<",\"xMin\":"; Num(o,v->GetXaxis()->GetXmin());
      o<<",\"xMax\":"; Num(o,v->GetXaxis()->GetXmax());
      if(v->GetDimension()==1) {
        o<<",\"binContents\":[";
        for(int b=1;b<=v->GetNbinsX();++b) { if(b>1)o<<','; Num(o,v->GetBinContent(b)); }
        o<<']';
      } else { o<<",\"yBins\":"<<v->GetNbinsY()<<",\"yMin\":"; Num(o,v->GetYaxis()->GetXmin()); o<<",\"yMax\":"; Num(o,v->GetYaxis()->GetXmax()); }
    }
    o << '}';
  }
  o << "\n]}\n"; o.close();

  gStyle->SetOptStat(0); gStyle->SetOptFit(0);
  TCanvas canvas("kfSummaryMassCanvas","KFParticle Lambda invariant mass",1100,900);
  TPad plotPad("kfSummaryMassPlot","mass",0,.30,1,1);
  TPad notePad("kfSummaryMassNotes","notes",0,0,1,.30);
  plotPad.SetLeftMargin(.12); plotPad.SetBottomMargin(.13);
  plotPad.Draw(); notePad.Draw(); plotPad.cd();
  TH1 *draw=dynamic_cast<TH1*>(h->Clone("kfSummaryMassDraw")); draw->SetDirectory(0);
  draw->SetTitle("KFParticle #Lambda: inclusive selected candidates");
  draw->GetXaxis()->SetTitle("M_{p#pi^{-}} [GeV/c^{2}]");
  draw->GetYaxis()->SetTitle(Form("Candidates / %.1f MeV/c^{2}",1000.*h->GetBinWidth(1)));
  draw->GetXaxis()->SetRangeUser(sbLow-.006,sbHigh+.006);
  draw->SetLineColor(kBlack); draw->Draw("E1");
  fit.SetLineColor(kOrange+7); fit.SetLineWidth(2);
  bg.SetLineColor(kMagenta+1); bg.SetLineStyle(3);
  if(fitValid) { fit.Draw("same"); bg.Draw("same"); }
  TLine lines[4]; const double edges[]={sbLow,low,high,sbHigh};
  for(int i=0;i<4;++i) { lines[i].SetX1(edges[i]); lines[i].SetX2(edges[i]); lines[i].SetY1(0); lines[i].SetY2(h->GetMaximum()*.80); lines[i].SetLineColor((i==0 || i==3)?kBlue:kRed); lines[i].SetLineStyle(2); lines[i].Draw(); }
  notePad.cd();
  TPaveText note(.08,.03,.97,.97,"NDC"); note.SetFillColor(0); note.SetBorderSize(0); note.SetTextAlign(12); note.SetTextSize(.092);
  note.AddText(Form("Nominal signal: %.5f -- %.5f",low,high));
  note.AddText(Form("Binned signal: %.5f -- %.5f",signal.lower,signal.upper));
  note.AddText(Form("Width-scaled SB purity: %.4f #pm %.4f",purity,purityError));
  if(fitValid) note.AddText(Form("Gaus+pol2 comparison: P=%.4f; #chi^{2}/ndf=%.1f/%d; Prob=%.2g",fitPurity,fit.GetChisquare(),fit.GetNDF(),fit.GetNDF()>0?TMath::Prob(fit.GetChisquare(),fit.GetNDF()):0.));
  else note.AddText(Form("Gaus+pol2 unavailable (status %d, cov %d)",fitStatus,covStatus));
  note.AddText("Fit convergence #neq good fit; no fit-based correction is applied");
  note.AddText("Inclusive proxy; no pair-dependent correction claim"); note.Draw();
  canvas.Print(pdfPath); canvas.Print(pngPath);
  std::cout << "KFSUMMARY_INSPECT " << jsonPath << " purity=" << purity
            << " usableForCFsub=" << purityUsable << std::endl;
  delete draw; f->Close(); delete f;
}
