// Read-only extraction for low/full/high Lambda-mass CF comparisons.
// ROOT 5 / ACLiC. All analysis choices arrive in the explicit settings file.
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
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

namespace MassHalfExtract {

void Require(bool yes, const std::string &message) {
  if (!yes) throw std::runtime_error(message);
}
bool Finite(double x) { return TMath::Finite(x); }
bool Close(double a, double b) {
  return std::fabs(a-b) <= 1.e-12 * std::max(1., std::max(std::fabs(a),std::fabs(b)));
}
std::string Quote(const std::string &s) {
  std::ostringstream o; o << '"';
  for (std::string::const_iterator i=s.begin();i!=s.end();++i) {
    if (*i=='"' || *i=='\\') o << '\\' << *i;
    else if (*i=='\n') o << "\\n";
    else if (*i=='\r') o << "\\r";
    else if (*i=='\t') o << "\\t";
    else o << *i;
  }
  o << '"'; return o.str();
}
void Num(std::ostream &o, double x) {
  if (Finite(x)) o << std::setprecision(17) << x; else o << "null";
}
void Array(std::ostream &o, const std::vector<double> &v) {
  o << '[';
  for (unsigned i=0;i<v.size();++i) { if(i) o << ','; Num(o,v[i]); }
  o << ']';
}
void Pair(std::ostream &o, const double *v, bool available=true) {
  o << '[';
  for (int i=0;i<2;++i) { if(i) o << ','; if(available) Num(o,v[i]); else o << "null"; }
  o << ']';
}
std::string Bool(bool x) { return x ? "true" : "false"; }

struct Group { std::string id; int lo,hi; };
struct Settings {
  double mean,sigma,nsigma,outer,normLo,normHi,fitLo,fitHi;
  double meanLo,meanHi,sigmaLo,sigmaHi,sigmaInitial,minEntries;
  int rebin,purityRebin,halfWidth,minCov;
  std::vector<Group> groups;
  Settings():mean(0),sigma(0),nsigma(0),outer(0),normLo(0),normHi(0),
    fitLo(0),fitHi(0),meanLo(0),meanHi(0),sigmaLo(0),sigmaHi(0),
    sigmaInitial(0),minEntries(0),rebin(0),purityRebin(0),halfWidth(0),minCov(0) {}
};
Settings ReadSettings(const char *path) {
  std::ifstream in(path); Require(in.good(),std::string("cannot read settings: ")+path);
  Settings c; std::set<std::string> seen, groupIds; std::string line; int lineNo=0;
  while(std::getline(in,line)) {
    ++lineNo; const std::string::size_type comment=line.find('#');
    if(comment!=std::string::npos) line.resize(comment);
    std::istringstream s(line); std::string key; if(!(s>>key)) continue;
    if(key!="group") Require(seen.insert(key).second,"duplicate settings key: "+key);
    if(key=="gate") s>>c.mean>>c.sigma>>c.nsigma>>c.outer;
    else if(key=="rebin") s>>c.rebin;
    else if(key=="purityRebin") s>>c.purityRebin;
    else if(key=="slidingHalfWidth") s>>c.halfWidth;
    else if(key=="normalization") s>>c.normLo>>c.normHi;
    else if(key=="fitRange") s>>c.fitLo>>c.fitHi;
    else if(key=="fitMeanLimits") s>>c.meanLo>>c.meanHi;
    else if(key=="fitSigmaLimits") s>>c.sigmaLo>>c.sigmaHi;
    else if(key=="fitSigmaInitial") s>>c.sigmaInitial;
    else if(key=="minimumEntries") s>>c.minEntries;
    else if(key=="minimumCovarianceStatus") s>>c.minCov;
    else if(key=="group") {
      Group g; g.lo=g.hi=-1; s>>g.id>>g.lo>>g.hi;
      Require(!s.fail(),Form("malformed group on settings line %d",lineNo));
      Require(!g.id.empty() && g.id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")==std::string::npos,
              "group id must use letters, digits, underscore or dash");
      Require(groupIds.insert(g.id).second,"duplicate group id: "+g.id);
      Require(g.lo>=0 && g.hi<=8 && g.lo<=g.hi,"invalid cent9 group: "+g.id);
      c.groups.push_back(g);
    } else throw std::runtime_error("unknown settings key: "+key);
    Require(!s.fail(),Form("malformed settings line %d",lineNo));
    std::string extra; Require(!(s>>extra),Form("extra tokens on settings line %d",lineNo));
  }
  const char *required[]={"gate","rebin","purityRebin","slidingHalfWidth","normalization","fitRange",
    "fitMeanLimits","fitSigmaLimits","fitSigmaInitial","minimumEntries","minimumCovarianceStatus"};
  for(unsigned i=0;i<sizeof(required)/sizeof(required[0]);++i)
    Require(seen.count(required[i])!=0,std::string("missing settings key: ")+required[i]);
  Require(c.mean>0 && c.sigma>0 && c.nsigma>0 && c.outer>1,"invalid gate settings");
  Require(c.rebin>0 && c.purityRebin>0 && c.halfWidth>=0 && c.minCov>=0 && c.minCov<=3 && c.minEntries>=0,
          "invalid integer/statistical settings");
  Require(c.normLo>=0 && c.normHi>c.normLo && c.fitHi>c.fitLo,"invalid fit/normalization range");
  Require(c.meanLo<c.mean && c.mean<c.meanHi && c.sigmaLo>0 && c.sigmaLo<c.sigmaHi &&
          c.sigmaInitial>=c.sigmaLo && c.sigmaInitial<=c.sigmaHi,"invalid fit shape settings");
  Require(!c.groups.empty(),"no explicit centrality groups");
  return c;
}
void WriteSettings(std::ostream &o,const Settings &c) {
  o << "{\"mean\":"; Num(o,c.mean); o << ",\"sigma\":"; Num(o,c.sigma);
  o << ",\"nSigma\":"; Num(o,c.nsigma); o << ",\"sidebandOuterFactor\":"; Num(o,c.outer);
  o << ",\"rebin\":"<<c.rebin<<",\"purityRebin\":"<<c.purityRebin
    <<",\"slidingHalfWidth\":"<<c.halfWidth<<",\"normalization\":[";
  Num(o,c.normLo);o<<',';Num(o,c.normHi);o<<"],\"fitRange\":[";Num(o,c.fitLo);o<<',';Num(o,c.fitHi);
  o<<"],\"fitMeanLimits\":[";Num(o,c.meanLo);o<<',';Num(o,c.meanHi);
  o<<"],\"fitSigmaLimits\":[";Num(o,c.sigmaLo);o<<',';Num(o,c.sigmaHi);
  o<<"],\"fitSigmaInitial\":";Num(o,c.sigmaInitial);o<<",\"minimumEntries\":";Num(o,c.minEntries);
  o<<",\"minimumCovarianceStatus\":"<<c.minCov<<",\"groups\":[";
  for(unsigned i=0;i<c.groups.size();++i) { if(i)o<<',';
    o<<"{\"id\":"<<Quote(c.groups[i].id)<<",\"centMin\":"<<c.groups[i].lo<<",\"centMax\":"<<c.groups[i].hi<<'}'; }
  o<<"]}";
}

TH1 *Histogram(TFile *f,const std::string &name,int dimension) {
  TH1 *h=dynamic_cast<TH1*>(f->Get(name.c_str()));
  Require(h && h->GetDimension()==dimension,"missing/wrong histogram: "+name); return h;
}
void SameAxis(const TAxis *a,const TAxis *b,const std::string &context) {
  Require(a->GetNbins()==b->GetNbins(),"different bin counts: "+context);
  for(int i=1;i<=a->GetNbins()+1;++i)
    Require(Close(a->GetBinLowEdge(i),b->GetBinLowEdge(i)),"different axis edges: "+context);
}
void UniformAxis(const TAxis *a,const std::string &context) {
  for(int b=2;b<=a->GetNbins();++b)
    Require(Close(a->GetBinWidth(1),a->GetBinWidth(b)),"nonuniform axis unsupported: "+context);
}
double Variance(const TH1 *h,int bin) {
  const double n=h->GetBinContent(bin);
  Require(Finite(n) && n>=0,std::string("negative/nonfinite count in ")+h->GetName());
  if(h->GetSumw2N()==0) return n;
  const double e=h->GetBinError(bin); Require(Finite(e),"nonfinite bin error"); return e*e;
}
std::string VariancePolicy(const TH1 *h) {
  return h->GetSumw2N()==0 ? "unit_fill_no_Sumw2_error2_equals_content" : "stored_Sumw2";
}
// Double-precision clone avoids TH2F Add rounding when summing centralities.
TH2D *DoubleCopy(TH2 *h,const std::string &name) {
  UniformAxis(h->GetXaxis(),h->GetName());UniformAxis(h->GetYaxis(),h->GetName());
  TH2D *d=new TH2D(name.c_str(),h->GetTitle(),h->GetNbinsX(),h->GetXaxis()->GetXmin(),
    h->GetXaxis()->GetXmax(),h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  d->SetDirectory(0); d->Sumw2();
  for(int x=0;x<=h->GetNbinsX()+1;++x) for(int y=0;y<=h->GetNbinsY()+1;++y) {
    const int bin=h->GetBin(x,y); d->SetBinContent(x,y,h->GetBinContent(bin));
    d->SetBinError(x,y,std::sqrt(Variance(h,bin)));
  }
  d->SetEntries(h->GetEntries()); return d;
}
struct Window { int first,last; double lower,upper; };
Window MakeWindow(const TAxis *a,double lo,double hi) {
  Require(lo>=a->GetXmin() && hi<=a->GetXmax() && hi>lo,"mass interval outside histogram axis");
  int first=1;while(first<=a->GetNbins() && a->GetBinCenter(first)<lo)++first;
  int after=first;while(after<=a->GetNbins() && a->GetBinCenter(after)<hi)++after;
  Require(first<=a->GetNbins() && after>first,"empty selected whole-bin mass interval");
  Window w;w.first=first;w.last=after-1;w.lower=a->GetBinLowEdge(first);w.upper=a->GetBinUpEdge(w.last);return w;
}
struct Partition { Window low,high,full,left,right; };
Partition MakePartition(const TAxis *axis,const Settings &c) {
  UniformAxis(axis,"mass");const double d=c.sigma*c.nsigma;Partition p;
  p.low=MakeWindow(axis,c.mean-d,c.mean);p.high=MakeWindow(axis,c.mean,c.mean+d);
  p.full=MakeWindow(axis,c.mean-d,c.mean+d);
  p.left=MakeWindow(axis,c.mean-c.outer*d,c.mean-d);p.right=MakeWindow(axis,c.mean+d,c.mean+c.outer*d);
  Require(p.low.last+1==p.high.first && p.low.first==p.full.first && p.high.last==p.full.last,
          "inconsistent low/high/full mass partition");return p;
}
void WriteWindow(std::ostream &o,const Window &w) {
  o<<"{\"firstBin\":"<<w.first<<",\"lastBin\":"<<w.last<<",\"lower\":";Num(o,w.lower);
  o<<",\"upper\":";Num(o,w.upper);o<<'}';
}
void WritePartition(std::ostream &o,const Partition &p,const Settings &c) {
  const double d=c.sigma*c.nsigma;o<<"{\"policy\":\"whole bins whose centers are in [lower,upper); no split-bin interpolation\",\"nominal\":{\"mean\":";
  Num(o,c.mean);o<<",\"low\":";Num(o,c.mean-d);o<<",\"high\":";Num(o,c.mean+d);
  o<<",\"sideLow\":";Num(o,c.mean-c.outer*d);o<<",\"sideHigh\":";Num(o,c.mean+c.outer*d);
  o<<"},\"binned\":{\"low\":";WriteWindow(o,p.low);o<<",\"high\":";WriteWindow(o,p.high);
  o<<",\"full\":";WriteWindow(o,p.full);o<<",\"leftSB\":";WriteWindow(o,p.left);
  o<<",\"rightSB\":";WriteWindow(o,p.right);o<<"}}";
}
void Counts(const TH1 *h,const Window &w,double &count,double &error2) {
  count=error2=0;
  for(int b=w.first;b<=w.last;++b) { count+=h->GetBinContent(b);error2+=Variance(h,b); }
}
void WriteInclusive(std::ostream &o,TH1 *h,const Partition &p) {
  o<<"{\"histogram\":"<<Quote(h->GetName())<<",\"variancePolicy\":"<<Quote(VariancePolicy(h));
  const Window windows[]={p.low,p.high,p.full,p.left,p.right};
  const char *names[]={"low","high","full","leftSB","rightSB"};
  for(int i=0;i<5;++i) { double n,v;Counts(h,windows[i],n,v);o<<','<<Quote(names[i])<<":{\"count\":";
    Num(o,n);o<<",\"error2\":";Num(o,v);o<<'}'; }
  o<<'}';
}

struct Spectrum { std::vector<double> counts,error2; std::string policy; };
Spectrum RawSpectrum(TH2 *h,const Group &g,int rebin) {
  Require(h->GetNbinsY()==9 && Close(h->GetYaxis()->GetXmin(),-.5) &&
    Close(h->GetYaxis()->GetXmax(),8.5),"raw pair histogram must have canonical cent9 Y axis");
  Require(h->GetNbinsX()%rebin==0,"rebin must divide raw kstar bin count");
  TH2D *copy=DoubleCopy(h,"mh_raw_copy");
  TH1D *pr=copy->ProjectionX("mh_raw_projection",g.lo+1,g.hi+1,"e");pr->SetDirectory(0);pr->Rebin(rebin);
  Spectrum a;a.policy=VariancePolicy(h);
  for(int b=1;b<=pr->GetNbinsX();++b) {
    a.counts.push_back(pr->GetBinContent(b));a.error2.push_back(std::pow(pr->GetBinError(b),2));
  }
  delete pr;delete copy;return a;
}
void WriteSpectrum(std::ostream &o,const Spectrum &a) {
  o<<"{\"counts\":";Array(o,a.counts);o<<",\"error2\":";Array(o,a.error2);
  o<<",\"variancePolicy\":"<<Quote(a.policy)<<'}';
}
TH2D *GroupedMass(TFile *f,const Group &g,bool mixed,int rebin,const TAxis *massAxis,const TAxis *kstarAxis) {
  TH2D *sum=0;
  for(int c=g.lo;c<=g.hi;++c) {
    const std::string name=Form("%s/hKstarMass_%s4He_CentBin%d",mixed?"mix":"true",mixed?"Mixed_":"",c);
    TH2 *h=static_cast<TH2*>(Histogram(f,name,2));
    SameAxis(massAxis,h->GetYaxis(),name);SameAxis(kstarAxis,h->GetXaxis(),name);
    TH2D *one=DoubleCopy(h,Form("mh_mass_copy_%d_%d",mixed,c));
    if(!sum) { sum=one;sum->SetName(Form("mh_mass_sum_%s_%d",g.id.c_str(),mixed)); }
    else { sum->Add(one);delete one; }
  }
  Require(sum && sum->GetNbinsX()%rebin==0,"invalid grouped mass histogram/rebin");
  sum->RebinX(rebin);return sum;
}

struct FitDiagnostics {
  int status,covStatus,ndf;bool isValid,boundary,finiteCov,numericallyValid;
  double chi2,probability,parameters[6],covariance[6][6];std::string reason;
  FitDiagnostics():status(-1),covStatus(-1),ndf(0),isValid(false),boundary(false),finiteCov(false),
    numericallyValid(false),chi2(0),probability(0),reason("not_attempted") {
    for(int i=0;i<6;++i) {parameters[i]=0;for(int j=0;j<6;++j)covariance[i][j]=0;}
  }
};
bool Boundary(const TF1 &f) {
  for(int i=0;i<f.GetNpar();++i) {
    double lo=0,hi=0;f.GetParLimits(i,lo,hi);
    if(hi>lo && (f.GetParameter(i)<=lo+(hi-lo)*1.e-4 || f.GetParameter(i)>=hi-(hi-lo)*1.e-4))return true;
  }
  return false;
}
FitDiagnostics FitMass(TH1D *h,const Settings &c) {
  FitDiagnostics d;
  // Unit-fill projections: the counted in-range observations, not propagated
  // histogram entry metadata, control the explicit minimumEntries requirement.
  if(h->Integral()<c.minEntries) {d.reason="insufficient_entries";return d;}
  TF1 fit("mh_total_fit","gaus(0)+pol2(3)",c.fitLo,c.fitHi);
  fit.SetParameters(h->GetMaximum(),c.mean,c.sigmaInitial,0.,0.,0.);
  fit.SetParLimits(1,c.meanLo,c.meanHi);fit.SetParLimits(2,c.sigmaLo,c.sigmaHi);
  TFitResultPtr result=h->Fit(&fit,"SQNR0","",c.fitLo,c.fitHi);
  d.status=int(result);d.covStatus=result.Get()?result->CovMatrixStatus():-1;
  d.isValid=result.Get() && result->IsValid();d.boundary=Boundary(fit);
  d.chi2=fit.GetChisquare();d.ndf=fit.GetNDF();
  d.probability=d.ndf>0 && Finite(d.chi2)?TMath::Prob(d.chi2,d.ndf):0;
  bool finite=true;for(int i=0;i<6;++i) {d.parameters[i]=fit.GetParameter(i);finite=finite && Finite(d.parameters[i]);}
  if(result.Get()) {
    TMatrixDSym cov=result->GetCovarianceMatrix();d.finiteCov=cov.GetNrows()==6;
    for(int i=0;i<6 && cov.GetNrows()==6;++i) for(int j=0;j<6;++j) {
      d.covariance[i][j]=cov(i,j);d.finiteCov=d.finiteCov && Finite(cov(i,j));
    }
  }
  d.numericallyValid=d.status==0 && d.isValid && d.covStatus>=c.minCov && finite && d.finiteCov;
  d.reason=d.numericallyValid?(d.boundary?"numerically_valid_boundary_warning":"numerically_valid"):
    "invalid_fit_or_covariance";return d;
}
void WriteFit(std::ostream &o,const FitDiagnostics &d,double observations) {
  o<<"{\"status\":"<<d.status<<",\"covarianceStatus\":"<<d.covStatus<<",\"isValid\":"<<Bool(d.isValid)
   <<",\"numericallyValid\":"<<Bool(d.numericallyValid)<<",\"boundary\":"<<Bool(d.boundary)
   <<",\"reason\":"<<Quote(d.reason)<<",\"chi2\":";
  if(d.status!=-1)Num(o,d.chi2);else o<<"null";
  o<<",\"ndf\":"<<d.ndf<<",\"probability\":";if(d.status!=-1)Num(o,d.probability);else o<<"null";
  o<<",\"observationsInStoredMassRange\":";Num(o,observations);o<<",\"parameters\":[";
  for(int i=0;i<6;++i) {if(i)o<<',';if(d.status!=-1)Num(o,d.parameters[i]);else o<<"null";}
  o<<"],\"covariance\":[";
  for(int i=0;i<6;++i) {if(i)o<<',';o<<'[';for(int j=0;j<6;++j) {if(j)o<<',';
    if(d.finiteCov)Num(o,d.covariance[i][j]);else o<<"null";}o<<']';}o<<"]}";
}
void WriteFitBin(std::ostream &o,TH2D *h,int bx,int halfWidth,int originalRebin,
                 bool isPurityGrid,const Settings &c,const Partition &p) {
  const int first=std::max(1,bx-halfWidth),last=std::min(h->GetNbinsX(),bx+halfWidth);
  TH1D *center=h->ProjectionY("mh_central_projection",bx,bx,"e");center->SetDirectory(0);
  TH1D *window=h->ProjectionY("mh_window_projection",first,last,"e");window->SetDirectory(0);
  double centralN[2],centralV[2],windowN[2],windowV[2];const Window regions[]={p.low,p.high};
  for(int r=0;r<2;++r) {Counts(center,regions[r],centralN[r],centralV[r]);Counts(window,regions[r],windowN[r],windowV[r]);}
  FitDiagnostics fit=FitMass(window,c);double background[2]={0,0},covB[2][2]={{0,0},{0,0}};
  long double gradient[2][3];
  bool bgAvailable=fit.numericallyValid;
  const double binWidth=h->GetYaxis()->GetBinWidth(p.full.first);
  // Stable integrals of 1, x, x^2 over a narrow mass interval. The raw pol2
  // coefficients/covariances are strongly correlated, so double contractions
  // lose significant digits. Do not evaluate the mathematically identical
  // LH and HL contractions in different floating-point summation orders.
  for(int r=0;r<2;++r) {
    const long double lo=regions[r].lower,hi=regions[r].upper;
    gradient[r][0]=(hi-lo)/static_cast<long double>(binWidth);
    gradient[r][1]=gradient[r][0]*(hi+lo)/2.L;
    gradient[r][2]=gradient[r][0]*(hi*hi+hi*lo+lo*lo)/3.L;
  }
  if(bgAvailable) {
    for(int r=0;r<2;++r) {
      long double value=0;
      for(int a=0;a<3;++a)value+=gradient[r][a]*static_cast<long double>(fit.parameters[a+3]);
      background[r]=static_cast<double>(value);
    }
    for(int r=0;r<2;++r)for(int s=r;s<2;++s) {
      long double value=0;
      for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        value+=gradient[r][a]*static_cast<long double>(fit.covariance[a+3][b+3])*gradient[s][b];
      covB[r][s]=covB[s][r]=static_cast<double>(value);
    }
    for(int r=0;r<2;++r) for(int s=0;s<2;++s) bgAvailable=bgAvailable && Finite(covB[r][s]);
    bgAvailable=bgAvailable && Finite(background[0]) && Finite(background[1]) && covB[0][0]>=0 && covB[1][1]>=0;
  }
  const double fullB=background[0]+background[1],fullBV=covB[0][0]+2*covB[0][1]+covB[1][1];
  bool valid[2];std::string regionReason[2];
  for(int r=0;r<2;++r) {
    valid[r]=bgAvailable && windowN[r]>0 && background[r]>=0 && background[r]<=windowN[r];
    regionReason[r]=!fit.numericallyValid?fit.reason:!bgAvailable?"nonfinite_or_negative_background_covariance":
      windowN[r]<=0?"empty_window_region":background[r]<0 || background[r]>windowN[r]?"nonphysical_background":
      fit.boundary?"valid_boundary_warning":"valid";
  }
  const bool fullValid=bgAvailable && fullBV>=0 && windowN[0]+windowN[1]>0 && fullB>=0 && fullB<=windowN[0]+windowN[1];
  o<<"{\"bin\":"<<bx<<",\"x\":";Num(o,h->GetXaxis()->GetBinCenter(bx));
  o<<",\"grid\":"<<Quote(isPurityGrid?"method5_coarse_purity":"method3_cf")
   <<",\"sourceRebin\":"<<originalRebin;
  if(isPurityGrid)o<<",\"purityBin\":"<<bx;
  o<<",\"centralOriginalFirstBin\":"<<(bx-1)*originalRebin+1
   <<",\"centralOriginalLastBin\":"<<bx*originalRebin
   <<",\"windowOriginalFirstBin\":"<<(first-1)*originalRebin+1
   <<",\"windowOriginalLastBin\":"<<last*originalRebin;
  o<<",\"windowFirstBin\":"<<first<<",\"windowLastBin\":"<<last<<",\"centralNtot\":";Pair(o,centralN);
  o<<",\"centralError2\":";Pair(o,centralV);o<<",\"windowNtot\":";Pair(o,windowN);o<<",\"windowError2\":";Pair(o,windowV);
  o<<",\"background\":";Pair(o,background,bgAvailable);o<<",\"backgroundCovariance\":[";
  Pair(o,covB[0],bgAvailable);o<<',';Pair(o,covB[1],bgAvailable);o<<"],\"valid\":["<<Bool(valid[0])<<','<<Bool(valid[1])
    <<"],\"regionReason\":["<<Quote(regionReason[0])<<','<<Quote(regionReason[1])<<"],\"full\":{\"centralNtot\":";
  Num(o,centralN[0]+centralN[1]);o<<",\"centralError2\":";Num(o,centralV[0]+centralV[1]);
  o<<",\"windowNtot\":";Num(o,windowN[0]+windowN[1]);o<<",\"windowError2\":";Num(o,windowV[0]+windowV[1]);
  o<<",\"background\":";if(bgAvailable)Num(o,fullB);else o<<"null";
  o<<",\"backgroundVariance\":";if(bgAvailable)Num(o,fullBV);else o<<"null";
  o<<",\"valid\":"<<Bool(fullValid)<<"},\"fit\":";WriteFit(o,fit,window->Integral());o<<'}';
  delete center;delete window;
}

// Original (5 MeV/c for these files) disjoint cells are the random-variable
// basis for Method 5. Apply each coarse-bin purity to its original cells, then
// aggregate cells into the independently configured CF bins downstream.
void WriteMethod5Base(std::ostream &o,TH2D *h,const Partition &p) {
  o<<'[';const Window regions[]={p.low,p.high};
  for(int bx=1;bx<=h->GetNbinsX();++bx) {
    if(bx>1)o<<',';double n[2]={0,0},v[2]={0,0};
    for(int r=0;r<2;++r)for(int by=regions[r].first;by<=regions[r].last;++by) {
      const int cell=h->GetBin(bx,by);n[r]+=h->GetBinContent(cell);v[r]+=Variance(h,cell);
    }
    o<<"{\"bin\":"<<bx<<",\"x\":";Num(o,h->GetXaxis()->GetBinCenter(bx));
    o<<",\"counts\":";Pair(o,n);o<<",\"error2\":";Pair(o,v);o<<'}';
  }
  o<<']';
}

std::vector<double> RebinEdges(const TAxis *axis,int rebin) {
  Require(axis->GetNbins()%rebin==0,"configured rebin must divide original kstar bin count");
  std::vector<double> edges;
  for(int b=1;b<=axis->GetNbins()+1;b+=rebin)edges.push_back(axis->GetBinLowEdge(b));
  return edges;
}

void Run(const char *inputRoot,const char *outputDir,const char *settingsFile) {
  Require(inputRoot && *inputRoot && outputDir && *outputDir && settingsFile && *settingsFile,
          "inputRoot, outputDir and settingsFile are all required");
  const Settings c=ReadSettings(settingsFile);const std::string out(outputDir);
  const std::string target=out+"/mass_half_extract.json",temporary=target+".tmp";
  Require(gSystem->AccessPathName(target.c_str()) && gSystem->AccessPathName(temporary.c_str()),
          "refusing to overwrite existing mass_half_extract.json or .tmp");
  if(gSystem->AccessPathName(out.c_str()))Require(gSystem->mkdir(out.c_str(),kTRUE)==0,"cannot create output directory");
  TFile *f=TFile::Open(inputRoot,"READ");Require(f && !f->IsZombie(),"cannot open input ROOT read-only");
  TH1 *inclusive=Histogram(f,"hLambda_InvMass",1);const Partition p=MakePartition(inclusive->GetXaxis(),c);
  Require(c.fitLo>=inclusive->GetXaxis()->GetXmin() && c.fitHi<=inclusive->GetXaxis()->GetXmax(),"fit interval outside stored mass axis");
  const char *modes[]={"SE","ME"};const char *regions[]={"signal","signalLow","signalHigh","leftSB","rightSB"};
  TH2 *raw[2][5];
  for(int m=0;m<2;++m)for(int r=0;r<5;++r) {
    raw[m][r]=static_cast<TH2*>(Histogram(f,Form("hKstar%sVsCent_lambda_he4_%s",modes[m],regions[r]),2));
    SameAxis(raw[0][0]->GetXaxis(),raw[m][r]->GetXaxis(),raw[m][r]->GetName());
  }
  Require(c.normLo>=raw[0][0]->GetXaxis()->GetXmin() && c.normHi<=raw[0][0]->GetXaxis()->GetXmax(),"normalization interval outside stored kstar axis");
  std::ofstream o(temporary.c_str());Require(o.good(),"cannot open output JSON");
  o<<"{\n\"schemaVersion\":2,\"inputRoot\":"<<Quote(inputRoot)<<",\"settingsFile\":"<<Quote(settingsFile)<<",\"settings\":";
  WriteSettings(o,c);o<<",\n\"massPartition\":";WritePartition(o,p,c);o<<",\n\"inclusive\":";WriteInclusive(o,inclusive,p);
  o<<",\n\"caveats\":[\"Raw signalLow/signalHigh/full use exact Maker candidate gates; mass-fit projections use the documented whole-bin partition. These are not identical gates.\","
   "\"Low/high background estimates come from the SAME fit and their covariance is retained. Full variance includes the off-diagonal terms.\","
   "\"Polynomial integral gradients and covariance contractions use long double; one symmetric off-diagonal contraction is mirrored. No covariance eigenvalue clipping or fit-threshold relaxation is applied.\","
   "\"Covariance of observed counts with fitted background and sliding-window cross-kstar covariance are not evaluated here.\","
   "\"Method5 rows use the independent coarse purity grid, not the CF grid. Each coarse window is fitted once; apply that purity to original cells, then sum cells into CF bins.\","
   "\"windowFirstBin/windowLastBin index the row's own grid. windowOriginalFirstBin/windowOriginalLastBin index the original kstar cells, inclusive and one-based.\","
   "\"Numerical fit acceptance is not a goodness-of-fit criterion. Fit probability is diagnostic, not an implicit selection cut.\","
   "\"Raw closure deviations are recorded unchanged; full signal is never replaced by low+high. Unit-fill Poisson variances are assumed only when original Sumw2 is absent.\"],\n\"groups\":[\n";
  for(unsigned gi=0;gi<c.groups.size();++gi) {
    if(gi)o<<",\n";const Group &g=c.groups[gi];std::cout<<"Mass-half extraction group "<<g.id<<std::endl;
    Spectrum spectra[2][5];for(int m=0;m<2;++m)for(int r=0;r<5;++r)spectra[m][r]=RawSpectrum(raw[m][r],g,c.rebin);
    const std::vector<double> edges=RebinEdges(raw[0][0]->GetXaxis(),c.rebin);
    const std::vector<double> purityEdges=RebinEdges(raw[0][0]->GetXaxis(),c.purityRebin);
    const std::vector<double> originalEdges=RebinEdges(raw[0][0]->GetXaxis(),1);
    o<<"{\"id\":"<<Quote(g.id)<<",\"centMin\":"<<g.lo<<",\"centMax\":"<<g.hi<<",\"kstarEdges\":";Array(o,edges);
    o<<",\"method5PurityEdges\":";Array(o,purityEdges);o<<",\"originalKstarEdges\":";Array(o,originalEdges);o<<",\"raw\":{";
    for(int m=0;m<2;++m) {if(m)o<<',';o<<Quote(modes[m])<<":{";
      for(int r=0;r<5;++r) {if(r)o<<',';o<<Quote(regions[r])<<':';WriteSpectrum(o,spectra[m][r]);}o<<'}';}
    o<<"},\"rawClosure\":{";
    for(int m=0;m<2;++m) {if(m)o<<',';std::vector<double> delta;double maxAbs=0,totalDiff=0,fullTotal=0;
      for(unsigned b=0;b<spectra[m][0].counts.size();++b) {const double d=spectra[m][0].counts[b]-spectra[m][1].counts[b]-spectra[m][2].counts[b];
        delta.push_back(d);maxAbs=std::max(maxAbs,std::fabs(d));totalDiff+=d;fullTotal+=spectra[m][0].counts[b];}
      o<<Quote(modes[m])<<":{\"difference\":";Array(o,delta);o<<",\"maxAbsDifference\":";Num(o,maxAbs);
      o<<",\"totalDifference\":";Num(o,totalDiff);o<<",\"fullCount\":";Num(o,fullTotal);o<<'}';}
    o<<"},\"method5Base\":{";
    TH2D *original[2],*fine[2],*coarse[2];
    for(int m=0;m<2;++m) {
      original[m]=GroupedMass(f,g,m==1,1,inclusive->GetXaxis(),raw[0][0]->GetXaxis());
      if(m)o<<',';o<<Quote(modes[m])<<':';WriteMethod5Base(o,original[m],p);
      fine[m]=static_cast<TH2D*>(original[m]->Clone(Form("mh_fine_%s_%d",g.id.c_str(),m)));
      coarse[m]=static_cast<TH2D*>(original[m]->Clone(Form("mh_coarse_%s_%d",g.id.c_str(),m)));
      fine[m]->SetDirectory(0);coarse[m]->SetDirectory(0);fine[m]->RebinX(c.rebin);coarse[m]->RebinX(c.purityRebin);
    }
    o<<"},\"methods\":{";
    for(int method=0;method<2;++method) {if(method)o<<',';o<<Quote(method?"method5":"method3")<<":{";
      for(int m=0;m<2;++m) {if(m)o<<',';o<<Quote(modes[m])<<":[";
        TH2D *h=method?coarse[m]:fine[m];
        for(int b=1;b<=h->GetNbinsX();++b) {if(b>1)o<<',';
          WriteFitBin(o,h,b,method?c.halfWidth:0,method?c.purityRebin:c.rebin,method!=0,c,p);}o<<']';}o<<'}';}
    for(int m=0;m<2;++m) {delete original[m];delete fine[m];delete coarse[m];}
    o<<"}}";o.flush();Require(o.good(),"JSON write failed");
  }
  o<<"\n]}\n";o.close();Require(!o.fail(),"cannot close output JSON");f->Close();delete f;
  Require(gSystem->Rename(temporary.c_str(),target.c_str())==0,"cannot finalize output JSON");
  std::cout<<"EXTRACT_MASS_HALF_COMPLETE "<<target<<std::endl;
}
} // namespace MassHalfExtract

void ExtractMassHalfCF(const char *inputRoot,const char *outputDir,const char *settingsFile) {
  try { MassHalfExtract::Run(inputRoot,outputDir,settingsFile); }
  catch(const std::exception &e) {
    std::cerr<<"ExtractMassHalfCF: "<<e.what()<<std::endl;gSystem->Exit(2);
  }
}
