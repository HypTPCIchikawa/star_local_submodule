// Joint bin-integrated Poisson study. He3/He4 have 14 COMMON coefficients;
// no independently fitted He3/He4 centers or widths exist in this program.
void GlobalHeliumFit(const char* input,const char* output,const char* settings,const char* seeds);
#ifndef __CINT__
#include "../../../StMaker/common/NuclearIdDeDxVsMom.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TString.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
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
namespace GlobalHelium {
typedef std::vector<double> Vec;
typedef std::vector<Vec> Matrix;
typedef std::map<std::string,std::string> Settings;
double Value(const Settings& cfg,const std::string& key) {Settings::const_iterator i=cfg.find(key);if(i==cfg.end())throw std::runtime_error("missing setting "+key);std::istringstream s(i->second);double v;std::string extra;s>>v;if(!s||!TMath::Finite(v)||(s>>extra))throw std::runtime_error("bad setting "+key);return v;}
std::string Text(const Settings& cfg,const std::string& key){Settings::const_iterator i=cfg.find(key);if(i==cfg.end())throw std::runtime_error("missing text setting "+key);return i->second;}
Settings ReadSettings(const char* path){std::ifstream f(path);if(!f)throw std::runtime_error("missing settings");Settings c;std::string line;while(std::getline(f,line)){size_t i=line.find('\t');if(i==std::string::npos||i==0||!c.insert(std::make_pair(line.substr(0,i),line.substr(i+1))).second)throw std::runtime_error("bad/duplicate settings TSV");}return c;}
std::string Quote(const std::string& s){std::ostringstream o;o<<'"';for(unsigned j=0;j<s.size();++j){if(s[j]=='"'||s[j]=='\\')o<<'\\';o<<s[j];}o<<'"';return o.str();}
double Poly(const double* c,int n,double x){double y=0;for(int i=n-1;i>=0;--i)y=y*x+c[i];return y;}
double Extreme(const double* c,int degree,double lo,double hi,bool minimum){Vec x;x.push_back(lo);x.push_back(hi);if(degree==2&&std::fabs(c[2])>1.e-14)x.push_back(-c[1]/(2*c[2]));if(degree==3){double a=3*c[3],b=2*c[2],d=c[1];if(std::fabs(a)<1.e-14){if(std::fabs(b)>1.e-14)x.push_back(-d/b);}else if(b*b-4*a*d>=0){double q=std::sqrt(b*b-4*a*d);x.push_back((-b+q)/(2*a));x.push_back((-b-q)/(2*a));}}double v=Poly(c,degree+1,lo);for(unsigned j=0;j<x.size();++j)if(x[j]>=lo&&x[j]<=hi){double p=Poly(c,degree+1,x[j]);v=minimum?std::min(v,p):std::max(v,p);}return v;}
double Mean(const Vec& g,int species,double r){return Poly(&g[(species-2)*7],4,std::log10(r));}
double Sigma(const Vec& g,int species,double r){return Poly(&g[(species-2)*7+4],3,std::log10(r));}
double Gaussian(double lo,double hi,double mean,double sigma){if(!(sigma>0))return 0;const double a=(lo-mean)/(std::sqrt(2.)*sigma),b=(hi-mean)/(std::sqrt(2.)*sigma);if(a>=0)return .5*(std::erfc(a)-std::erfc(b));if(b<=0)return .5*(std::erfc(-b)-std::erfc(-a));return .5*(std::erf(b)-std::erf(a));}
Vec Legacy(){Vec g;for(int s=2;s<4;++s){const NuclearIdDeDxVsMom::LogPolyParams& m=NuclearIdDeDxVsMom::Mean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s));const NuclearIdDeDxVsMom::LogPolyParams& w=NuclearIdDeDxVsMom::Sigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s));g.push_back(m.c0);g.push_back(m.c1);g.push_back(m.c2);g.push_back(m.c3);g.push_back(w.c0);g.push_back(w.c1);g.push_back(w.c2);}return g;}
struct Bin {double lo,hi,n;};
struct Slice {double lo,hi,r,count,fitLo,fitHi;int classId,index;std::vector<Bin> bins;Vec nodes,weights,shape,amplitudes;int localStatus,amplitudeIterations;double localDeviance,amplitudeKktResidual,localEdm;bool amplitudeConverged,localMinimizerConverged,localPointAccepted;Slice():count(0),classId(0),index(0),localStatus(-1),amplitudeIterations(0),localDeviance(0),amplitudeKktResidual(1.e20),localEdm(1.e20),amplitudeConverged(false),localMinimizerConverged(false),localPointAccepted(false){shape.assign(5,0);shape[1]=shape[3]=1;}};
struct Profile {Vec amplitude;double deviance,kktResidual;int iterations;bool converged;Profile():deviance(0),kktResidual(1.e20),iterations(0),converged(false){}};
bool Solve(Matrix h,Vec b,Vec& x){int n=b.size();x.assign(n,0);for(int i=0;i<n;++i){int pivot=i;for(int j=i+1;j<n;++j)if(std::fabs(h[j][i])>std::fabs(h[pivot][i]))pivot=j;if(std::fabs(h[pivot][i])<1.e-16)return false;std::swap(h[i],h[pivot]);std::swap(b[i],b[pivot]);double d=h[i][i];for(int k=i;k<n;++k)h[i][k]/=d;b[i]/=d;for(int j=0;j<n;++j)if(j!=i){double a=h[j][i];for(int k=i;k<n;++k)h[j][k]-=a*h[i][k];b[j]-=a*b[i];}}x=b;return true;}
double Deviance(const Slice& s,const Matrix& t,const Vec& a){double d=0;for(unsigned b=0;b<t.size();++b){double v=0;for(unsigned k=0;k<a.size();++k)v+=a[k]*t[b][k];double expectation=std::max(1.e-20,s.count*v),n=s.bins[b].n;if(n>0){double q=(expectation-n)/n;d+=2*n*(std::fabs(q)<.5?q-std::log1p(q):expectation/n-1+std::log(n/expectation));}else d+=2*expectation;}return std::max(0.,d);}
const double coefficientScales[]={1.,3.,5.,10.,.1,.3,.5};
struct InputAudit {unsigned long long cells;double maximumIntegerResidual,maximumVarianceResidual;InputAudit():cells(0),maximumIntegerResidual(0),maximumVarianceResidual(0){}};
void AuditCount(double n,double variance,const Settings& cfg,InputAudit& audit){
  if(!TMath::Finite(n)||!TMath::Finite(variance)||n<0||variance<0)throw std::runtime_error("nonfinite/negative count or variance");
  double integerResidual=std::fabs(n-std::floor(n+.5)),varianceResidual=std::fabs(variance-n)/std::max(1.,n);
  ++audit.cells;audit.maximumIntegerResidual=std::max(audit.maximumIntegerResidual,integerResidual);audit.maximumVarianceResidual=std::max(audit.maximumVarianceResidual,varianceResidual);
  if(integerResidual>Value(cfg,"countIntegerTolerance"))throw std::runtime_error("Poisson study requires integer counts");
  if(varianceResidual>Value(cfg,"countVarianceRelativeTolerance"))throw std::runtime_error("Poisson study requires variance=count (unit weights)");
}
void AuditAxes(const TAxis& all,const TAxis& tag,const Settings& cfg){
  if(all.GetNbins()!=tag.GetNbins())throw std::runtime_error("ALL/TAG axis number of bins differs");
  for(int i=1;i<=all.GetNbins()+1;++i)if(std::fabs(all.GetBinLowEdge(i)-tag.GetBinLowEdge(i))>Value(cfg,"axisTolerance"))throw std::runtime_error("ALL/TAG bin edges differ");
}
void AuditSourceCells(const TH2& all,const TH2* tag,int first,int last,const Settings& cfg,InputAudit& audit){
  for(int x=first;x<=last;++x)for(int y=0;y<=all.GetNbinsY()+1;++y){
    double n=all.GetBinContent(x,y),v=std::pow(all.GetBinError(x,y),2);AuditCount(n,v,cfg,audit);
    if(tag){double nt=tag->GetBinContent(x,y),vt=std::pow(tag->GetBinError(x,y),2);AuditCount(nt,vt,cfg,audit);if(nt>n+Value(cfg,"countIntegerTolerance"))throw std::runtime_error("TOF TAG not a native-cell subset of ALL");AuditCount(std::max(0.,n-nt),v-vt,cfg,audit);}
  }
}
void ProtectOutput(const TString& out,const char* settings,const char* seeds){
  const char* names[]={"global_fit.json","input_slices.root","completed.json"};
  for(unsigned j=0;j<sizeof(names)/sizeof(names[0]);++j)if(!gSystem->AccessPathName(out+"/"+names[j]))throw std::runtime_error("existing generated output will not be replaced");
  if(!gSystem->AccessPathName(out)){if(gSystem->AccessPathName(out+"/context.json")||TString(settings)!=out+"/settings.tsv"||TString(seeds)!=out+"/seeds.tsv")throw std::runtime_error("existing output is not a prepared launcher study");}
  else if(gSystem->mkdir(out,kTRUE))throw std::runtime_error("cannot create new study output");
}
class Engine {
 public:
  Settings cfg;std::vector<Slice> slices;bool exponential;int quadrature;Vec legacy;
  Engine(const Settings& c):cfg(c),exponential(Text(c,"backgroundModel")=="exponential"),quadrature(int(Value(c,"rQuadratureOrder"))),legacy(Legacy()){}
  bool Feasible(const Vec& g)const{
    if(g.size()!=14)return false;for(unsigned j=0;j<g.size();++j)if(!TMath::Finite(g[j]))return false;
    double low=std::log10(Value(cfg,"rigidityMin")),high=std::log10(Value(cfg,"rigidityMax"));
    for(int s=0;s<2;++s){const double* m=&g[s*7];const double* w=&g[s*7+4];if(Extreme(m,3,low,high,true)<=0||Extreme(w,2,low,high,true)<Value(cfg,"minimumSigma"))return false;if(Value(cfg,"requireMonotonicMean")){double d[]={m[1],2*m[2],3*m[3]};if(Extreme(d,2,low,high,false)>0)return false;}double shiftLow[4],shiftHigh[4],width[3];for(int j=0;j<4;++j){double sig=j<3?legacy[s*7+4+j]:0;shiftLow[j]=m[j]-legacy[s*7+j]-Value(cfg,"meanShiftMinSigma")*sig;shiftHigh[j]=m[j]-legacy[s*7+j]-Value(cfg,"meanShiftMaxSigma")*sig;}for(int j=0;j<3;++j)width[j]=w[j]-Value(cfg,"maximumSigmaScale")*legacy[s*7+4+j];if(Extreme(shiftLow,3,low,high,true)<0||Extreme(shiftHigh,3,low,high,false)>0||Extreme(width,2,low,high,false)>0)return false;}
    for(int j=0;j<14;++j)if(std::fabs(g[j]-legacy[j])>Value(cfg,"coefficientBoundSteps")*coefficientScales[j%7])return false;
    double order[4];for(int j=0;j<4;++j)order[j]=g[7+j]-g[j];return Extreme(order,3,low,high,true)>=Value(cfg,"minimumMeanGap");
  }
  bool ShapeFeasible(const Slice& s,const Vec& g,const Vec& shape)const{
    if(shape.size()!=5)return false;for(unsigned j=0;j<shape.size();++j)if(!TMath::Finite(shape[j]))return false;
    for(int sp=0;sp<2;++sp)if(std::fabs(shape[2*sp])>Value(cfg,"neighborMeanShiftMax")||shape[2*sp+1]<Value(cfg,"neighborWidthMin")||shape[2*sp+1]>Value(cfg,"neighborWidthMax"))return false;
    if(std::fabs(shape[4])>Value(cfg,"exponentialSlopeMax"))return false;
    double m[2][4];for(int sp=0;sp<2;++sp){const NuclearIdDeDxVsMom::LogPolyParams& a=NuclearIdDeDxVsMom::Mean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(sp));const NuclearIdDeDxVsMom::LogPolyParams& w=NuclearIdDeDxVsMom::Sigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(sp));const double ma[]={a.c0,a.c1,a.c2,a.c3},wa[]={w.c0,w.c1,w.c2,0};for(int j=0;j<4;++j)m[sp][j]=ma[j]+shape[2*sp]*wa[j];}double orderDT[4],orderT3[4];for(int j=0;j<4;++j){orderDT[j]=m[1][j]-m[0][j];orderT3[j]=g[j]-m[1][j];}double low=std::log10(s.lo),high=std::log10(s.hi);return Extreme(orderDT,3,low,high,true)>=Value(cfg,"minimumMeanGap")&&Extreme(orderT3,3,low,high,true)>=Value(cfg,"minimumMeanGap");
  }
  void Nodes(Slice& s)const{const double x3[]={-.7745966692414834,0,.7745966692414834},w3[]={5./18,4./9,5./18};const double x5[]={-.9061798459386640,-.5384693101056831,0,.5384693101056831,.9061798459386640},w5[]={.1184634425280945,.2393143352496832,.2844444444444444,.2393143352496832,.1184634425280945};s.nodes.clear();s.weights.clear();if(quadrature==1){s.nodes.push_back(s.r);s.weights.push_back(1);}else if(quadrature==3||quadrature==5){for(int j=0;j<quadrature;++j){s.nodes.push_back(s.r+(s.hi-s.lo)/2*(quadrature==3?x3[j]:x5[j]));s.weights.push_back(quadrature==3?w3[j]:w5[j]);}}else throw std::runtime_error("quadrature must be 1,3,5");}
  Matrix Templates(const Slice& s,const Vec& g,const Vec& shape)const{
    const int n=exponential?5:6;Matrix t(s.bins.size(),Vec(n,0));double width=s.fitHi-s.fitLo;
    for(unsigned b=0;b<s.bins.size();++b){const Bin& bin=s.bins[b];for(int species=0;species<4;++species){for(unsigned q=0;q<s.nodes.size();++q){double r=s.nodes[q],m,w;if(species<2){NuclearIdDeDxVsMom::SpeciesIndex sp=static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(species);w=NuclearIdDeDxVsMom::GetSigma(sp,r);m=NuclearIdDeDxVsMom::GetMean(sp,r)+shape[2*species]*w;w*=shape[2*species+1];}else{m=Mean(g,species,r);w=Sigma(g,species,r);}t[b][species]+=s.weights[q]*Gaussian(bin.lo,bin.hi,m,w);}}
      if(exponential){double k=shape[4];if(std::fabs(k)<1.e-8)t[b][4]=(bin.hi-bin.lo)/width;else t[b][4]=std::exp(k*(bin.lo-s.fitLo)/width)*std::expm1(k*(bin.hi-bin.lo)/width)/std::expm1(k);}
      else{double center=(bin.lo+bin.hi)/2;t[b][4]=2*(bin.hi-bin.lo)*(s.fitHi-center)/(width*width);t[b][5]=2*(bin.hi-bin.lo)*(center-s.fitLo)/(width*width);}}
    return t;
  }
  Profile Amplitudes(const Slice& s,const Matrix& t,const Vec& initial)const{
    const int n=t[0].size();Profile out;out.amplitude=initial.size()==unsigned(n)?initial:Vec(n,1./n);Vec& a=out.amplitude;double cap=Value(cfg,"maximumAmplitudeFraction");for(int j=0;j<n;++j)a[j]=std::min(cap,std::max(0.,a[j]));double previous=Deviance(s,t,a);
    for(int iteration=0;iteration<int(Value(cfg,"amplitudeMaxIterations"));++iteration){Vec gradient(n,0);Matrix hessian(n,Vec(n,0));for(unsigned b=0;b<t.size();++b){double v=0;for(int j=0;j<n;++j)v+=a[j]*t[b][j];v=std::max(1.e-20,v);double observed=s.bins[b].n/s.count;for(int j=0;j<n;++j){gradient[j]+=t[b][j]*(1-observed/v);for(int k=0;k<n;++k)hessian[j][k]+=observed*t[b][j]*t[b][k]/(v*v);}}
      std::vector<int> active;double violation=0;for(int j=0;j<n;++j){bool blocked=(a[j]<1.e-12&&gradient[j]>0)||(a[j]>cap-1.e-12&&gradient[j]<0);if(!blocked){active.push_back(j);violation=std::max(violation,std::fabs(gradient[j]));}}
      out.iterations=iteration+1;out.kktResidual=violation;if(violation<Value(cfg,"amplitudeGradientTolerance")){out.converged=true;break;}if(active.empty()){out.converged=true;break;}
      Matrix h(active.size(),Vec(active.size()));Vec rhs(active.size());for(unsigned j=0;j<active.size();++j){rhs[j]=-gradient[active[j]];for(unsigned k=0;k<active.size();++k)h[j][k]=hessian[active[j]][active[k]];h[j][j]+=1.e-10;}
      Vec direction;bool solved=Solve(h,rhs,direction),accepted=false;
      for(int strategy=0;strategy<2&&!accepted;++strategy){
        if(strategy==1||!solved){direction.assign(active.size(),0);for(unsigned j=0;j<active.size();++j)direction[j]=-gradient[active[j]]/std::max(1.e-8,hessian[active[j]][active[j]]);}
        double step=1.;for(unsigned j=0;j<active.size();++j){double d=direction[j],v=a[active[j]];if(d<0)step=std::min(step,-v/d);if(d>0)step=std::min(step,(cap-v)/d);}step=std::max(0.,step);
        for(int line=0;line<30&&step>1.e-16;++line){Vec next=a;for(unsigned j=0;j<active.size();++j)next[active[j]]=std::min(cap,std::max(0.,a[active[j]]+step*direction[j]));double dev=Deviance(s,t,next);if(dev<previous){a=next;previous=dev;accepted=true;break;}step*=.5;}
      }
      if(!accepted)break;
    }
    // Recompute the projected gradient at the actual returned point, including
    // iteration-limit / failed-line-search exits. Never imply a profile optimum.
    Vec gradient(n,0);for(unsigned b=0;b<t.size();++b){double v=0;for(int j=0;j<n;++j)v+=a[j]*t[b][j];v=std::max(1.e-20,v);for(int j=0;j<n;++j)gradient[j]+=t[b][j]*(1-s.bins[b].n/(s.count*v));}
    out.kktResidual=0;for(int j=0;j<n;++j){bool blocked=(a[j]<1.e-12&&gradient[j]>0)||(a[j]>cap-1.e-12&&gradient[j]<0);if(!blocked)out.kktResidual=std::max(out.kktResidual,std::fabs(gradient[j]));}out.converged=out.kktResidual<Value(cfg,"amplitudeGradientTolerance");out.deviance=Deviance(s,t,a);return out;
  }
  double ProfileSlice(const Slice& s,const Vec& g,const Vec& shape)const{if(!ShapeFeasible(s,g,shape))return 1.e20;Profile p=Amplitudes(s,Templates(s,g,shape),s.amplitudes);return p.converged?p.deviance:1.e20;}
  double Objective(const Vec& g)const{if(!Feasible(g))return 1.e20;double loss=0;for(unsigned j=0;j<slices.size();++j){double d=ProfileSlice(slices[j],g,slices[j].shape);if(d>=1.e19)return d;loss+=d;}return loss;}
};
struct GlobalObjective{const Engine* engine;double operator()(const double* p)const{return engine->Objective(Vec(p,p+14));}};
struct LocalObjective{const Engine* engine;const Slice* slice;const Vec* global;std::vector<int> free;double operator()(const double* p)const{Vec shape=slice->shape;for(unsigned j=0;j<free.size();++j)shape[free[j]]=p[j];return engine->ProfileSlice(*slice,*global,shape);}};
void RefreshAmplitude(const Engine& e,Slice& s,const Vec& g){Profile p=e.Amplitudes(s,e.Templates(s,g,s.shape),s.amplitudes);s.amplitudes=p.amplitude;s.amplitudeIterations=p.iterations;s.localDeviance=p.deviance;s.amplitudeKktResidual=p.kktResidual;s.amplitudeConverged=p.converged;}
void FitLocal(Engine& e,Slice& s,const Vec& g){
  LocalObjective objective;objective.engine=&e;objective.slice=&s;objective.global=&g;
  if(Value(e.cfg,"floatDeuteronShape")){objective.free.push_back(0);objective.free.push_back(1);}
  if(Value(e.cfg,"floatTritonShape")){objective.free.push_back(2);objective.free.push_back(3);}
  if(e.exponential)objective.free.push_back(4);
  if(!objective.free.empty()){
    ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    if(!m)throw std::runtime_error("cannot create local Minuit2");
    ROOT::Math::Functor f(objective,objective.free.size());m->SetFunction(f);
    m->SetMaxFunctionCalls(int(Value(e.cfg,"localMaxCalls")));m->SetTolerance(Value(e.cfg,"localTolerance"));m->SetPrintLevel(-1);
    for(unsigned j=0;j<objective.free.size();++j){int index=objective.free[j];bool width=index==1||index==3;double low=width?Value(e.cfg,"neighborWidthMin"):-Value(e.cfg,index==4?"exponentialSlopeMax":"neighborMeanShiftMax");double high=width?Value(e.cfg,"neighborWidthMax"):Value(e.cfg,index==4?"exponentialSlopeMax":"neighborMeanShiftMax");m->SetLimitedVariable(j,Form("shape_%d",index),s.shape[index],.05,low,high);}
    Vec packed;for(unsigned j=0;j<objective.free.size();++j)packed.push_back(s.shape[objective.free[j]]);
    double before=objective(&packed[0]);bool ok=m->Minimize();double proposed=objective(m->X());
    s.localPointAccepted=proposed<1.e19&&proposed<=before;
    if(s.localPointAccepted)for(unsigned j=0;j<objective.free.size();++j)s.shape[objective.free[j]]=m->X()[j];
    s.localStatus=m->Status();s.localEdm=m->Edm();s.localMinimizerConverged=ok&&s.localStatus==0&&s.localPointAccepted;delete m;
  }else{s.localStatus=0;s.localEdm=0;s.localMinimizerConverged=true;s.localPointAccepted=true;}
  RefreshAmplitude(e,s,g);
}
bool AllAmplitudeConverged(const Engine& e){for(unsigned j=0;j<e.slices.size();++j)if(!e.slices[j].amplitudeConverged)return false;return true;}
bool AllLocalConverged(const Engine& e){for(unsigned j=0;j<e.slices.size();++j)if(!e.slices[j].localMinimizerConverged)return false;return true;}
double CoefficientDistance(const Vec& a,const Vec& b){double d=0;for(unsigned j=0;j<a.size();++j)d=std::max(d,std::fabs(a[j]-b[j])/coefficientScales[j%7]);return d;}
struct Stationarity {double residual,count;int central,oneSided,blocked;bool evaluable;Stationarity():residual(0),count(0),central(0),oneSided(0),blocked(0),evaluable(true){}};
// Axis-direction diagnostic only. At nonlinear constraint intersections this
// does not prove stationarity in every feasible tangent direction.
Stationarity GlobalStationarity(const Engine& e,const Vec& g){
  Stationarity s;double v=e.Objective(g),step=Value(e.cfg,"stationarityStep");
  for(unsigned j=0;j<e.slices.size();++j)s.count+=e.slices[j].count;
  if(v>=1.e19||!(step>0)||!(s.count>0)){s.evaluable=false;s.residual=1.e20;return s;}
  for(unsigned j=0;j<g.size();++j){Vec up=g,down=g;up[j]+=step*coefficientScales[j%7];down[j]-=step*coefficientScales[j%7];double vp=e.Objective(up),vm=e.Objective(down);bool plus=vp<1.e19,minus=vm<1.e19;
    if(plus&&minus){s.residual=std::max(s.residual,std::fabs(vp-vm)/(2*step));++s.central;}
    else if(plus){s.residual=std::max(s.residual,std::max(0.,v-vp)/step);++s.oneSided;}
    else if(minus){s.residual=std::max(s.residual,std::max(0.,v-vm)/step);++s.oneSided;}
    else{++s.blocked;s.evaluable=false;}
  }return s;
}
void WriteVec(std::ostream& o,const Vec& v){o<<'[';for(unsigned j=0;j<v.size();++j){if(j)o<<',';o<<v[j];}o<<']';}
}
#endif

void GlobalHeliumFit(const char* input,const char* output,const char* settings,const char* seeds) {
#ifndef __CINT__
  using namespace GlobalHelium;
  try {
    const Settings cfg=ReadSettings(settings);const TString out(output);
    ProtectOutput(out,settings,seeds);
    if(gSystem->Load("libMinuit2")<0)throw std::runtime_error("cannot load Minuit2");
    if(Text(cfg,"sampleMode")!="all"&&Text(cfg,"sampleMode")!="tag_not_tag")throw std::runtime_error("unsupported sampleMode");
    if(Text(cfg,"backgroundModel")!="linear"&&Text(cfg,"backgroundModel")!="exponential")throw std::runtime_error("unsupported backgroundModel");
    Value(cfg,"outerCoefficientTolerance");Value(cfg,"stationarityRelativeTolerance");Value(cfg,"stationarityStep");
    TFile source(input,"READ");if(source.IsZombie())throw std::runtime_error("cannot open source ROOT");
    TH2* all=dynamic_cast<TH2*>(source.Get(Text(cfg,"histogram").c_str()));
    bool disjoint=Text(cfg,"sampleMode")=="tag_not_tag";
    TH2* tag=disjoint?dynamic_cast<TH2*>(source.Get(Text(cfg,"tagHistogram").c_str())):0;
    if(!all)throw std::runtime_error("missing prePID TH2");
    if(disjoint&&!tag)throw std::runtime_error("missing TOF tag histogram");
    if(tag){AuditAxes(*all->GetXaxis(),*tag->GetXaxis(),cfg);AuditAxes(*all->GetYaxis(),*tag->GetYaxis(),cfg);}
    const double rmin=Value(cfg,"rigidityMin"),rmax=Value(cfg,"rigidityMax"),step=Value(cfg,"rigidityStep"),tol=Value(cfg,"axisTolerance");
    if(!(rmin>0&&rmax>rmin&&step>0))throw std::runtime_error("invalid R domain");
    int total=int((rmax-rmin)/step+.5);
    if(std::fabs(rmin+total*step-rmax)>tol)throw std::runtime_error("noninteger R slice count");
    std::vector<int> indices;std::set<int> used;std::istringstream list(Text(cfg,"sliceIndices"));std::string token;
    while(std::getline(list,token,',')){std::istringstream item(token);int index;std::string extra;if(!(item>>index)||(item>>extra)||index<0||index>=total||!used.insert(index).second)throw std::runtime_error("invalid/duplicate slice index");indices.push_back(index);}
    if(indices.empty())for(int j=0;j<total;++j)indices.push_back(j);
    Engine base(cfg);InputAudit inputAudit;TFile saved(out+"/input_slices.root","CREATE");
    if(saved.IsZombie())throw std::runtime_error("cannot CREATE protected slice ROOT");
    for(unsigned a=0;a<indices.size();++a){
      int index=indices[a];double lo=rmin+index*step,hi=lo+step,r=(lo+hi)/2;
      int first=all->GetXaxis()->FindFixBin(lo+tol),last=all->GetXaxis()->FindFixBin(hi-tol);
      if(first<1||last>all->GetNbinsX()||last<first||std::fabs(all->GetXaxis()->GetBinLowEdge(first)-lo)>tol||std::fabs(all->GetXaxis()->GetBinUpEdge(last)-hi)>tol)throw std::runtime_error("R slice not native-bin aligned");
      AuditSourceCells(*all,tag,first,last,cfg,inputAudit);
      TH1D* raw=all->ProjectionY(Form("all_native_%03d",index),first,last,"e");raw->SetDirectory(0);saved.cd();raw->Write();
      TH1D* tagged=0;if(disjoint){tagged=tag->ProjectionY(Form("tag_native_%03d",index),first,last,"e");tagged->SetDirectory(0);tagged->Write();}
      const int classes=disjoint?2:1;
      for(int classIndex=0;classIndex<classes;++classIndex){
        TH1D* h=static_cast<TH1D*>((disjoint&&classIndex==0?tagged:raw)->Clone(Form("fitraw_%03d_class%d",index,classIndex)));h->SetDirectory(0);
        if(disjoint&&classIndex==1)for(int b=0;b<=h->GetNbinsX()+1;++b){
          double n=raw->GetBinContent(b)-tagged->GetBinContent(b),variance=std::pow(raw->GetBinError(b),2)-std::pow(tagged->GetBinError(b),2);
          AuditCount(n,variance,cfg,inputAudit);h->SetBinContent(b,n);h->SetBinError(b,std::sqrt(std::max(0.,variance)));
        }
        double factorValue=Value(cfg,"dedxRebinFactor");int factor=int(factorValue);
        if(factor<1||factor!=factorValue||h->GetNbinsX()%factor)throw std::runtime_error("dE/dx rebin not an exact partition");
        h->Rebin(factor);h->Write();Slice s;s.lo=lo;s.hi=hi;s.r=r;s.index=index;s.classId=disjoint?classIndex+1:0;base.Nodes(s);
        double oldMean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),oldSigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);
        if(!(oldMean>0&&oldSigma>0))throw std::runtime_error("R outside legacy fit-window validity");
        double requestLow=std::max(h->GetXaxis()->GetXmin(),oldMean+Value(cfg,"fitNSigmaMin")*oldSigma),requestHigh=std::min(h->GetXaxis()->GetXmax(),oldMean+Value(cfg,"fitNSigmaMax")*oldSigma);
        int bfirst=std::max(1,h->FindBin(requestLow+tol)),blast=std::min(h->GetNbinsX(),h->FindBin(requestHigh-tol));
        if(blast<bfirst)throw std::runtime_error("empty fit window");
        s.fitLo=h->GetBinLowEdge(bfirst);s.fitHi=h->GetBinLowEdge(blast)+h->GetBinWidth(blast);
        for(int b=bfirst;b<=blast;++b){Bin v;v.lo=h->GetBinLowEdge(b);v.hi=v.lo+h->GetBinWidth(b);v.n=h->GetBinContent(b);AuditCount(v.n,std::pow(h->GetBinError(b),2),cfg,inputAudit);s.count+=v.n;s.bins.push_back(v);}
        if(s.count>=Value(cfg,"minimumSliceCount"))base.slices.push_back(s);delete h;
      }delete raw;delete tagged;
    }
    saved.Close();if(base.slices.empty())throw std::runtime_error("no usable slices");
    std::ifstream seedInput(seeds);if(!seedInput)throw std::runtime_error("missing explicit shared seeds");
    std::ofstream json((out+"/global_fit.json").Data());if(!json)throw std::runtime_error("cannot create fit JSON");
    json<<std::setprecision(17)<<"{\"schemaVersion\":1,\"globalParameterCount\":14,\"input\":"<<Quote(input)
      <<",\"sampleMode\":"<<Quote(Text(cfg,"sampleMode"))<<",\"background\":"<<Quote(Text(cfg,"backgroundModel"))<<",\"rQuadratureOrder\":"<<base.quadrature
      <<",\"rWithinBinPolicy\":\"uniform GL-average; source has no sub-bin track momenta\""
      <<",\"inputAuditScope\":\"native selected regular X bins, every Y including flow; unused X flow excluded\""
      <<",\"fitWindowPolicy\":\"fixed legacy-relative request clipped to native regular Y axis and expanded to full bins; outside-axis data unavailable\""
      <<",\"amplitudeDefinition\":\"Gaussian full-line yield / observed fit-window count, not fraction within window or sum-to-one; smooth background uses window yield / window count\""
      <<",\"stationarityScope\":\"axis-direction projected finite differences in scaled coefficients; nonlinear-boundary tangent stationarity not proven\""
      <<",\"covarianceScope\":\"global14 conditional on retained final T/d and nonlinear background shapes; amplitudes profiled; NOT full nuisance covariance\""
      <<",\"inputAudit\":{\"cells\":"<<inputAudit.cells<<",\"maximumIntegerResidual\":"<<inputAudit.maximumIntegerResidual<<",\"maximumVarianceRelativeResidual\":"<<inputAudit.maximumVarianceResidual<<"},\"runs\":[";
    std::string seedName;int runIndex=0;std::set<std::string> seedNames;
    while(seedInput>>seedName){
      if(!seedNames.insert(seedName).second)throw std::runtime_error("duplicate seed name");
      Vec g(14);for(int j=0;j<14;++j)if(!(seedInput>>g[j]))throw std::runtime_error("invalid seed");
      if(!base.Feasible(g))throw std::runtime_error("infeasible shared seed "+seedName);
      Engine engine=base;Matrix covariance(14,Vec(14,0));Vec covariancePoint,trace;std::vector<std::string> outerRecords;
      int globalStatus=-1,covStatus=-1;bool outerConverged=false,globalAccepted=false,globalMinimizerConverged=false,covariancePointMatches=false;
      double previous=1.e20,globalEdm=1.e20,coefficientDrift=1.e20;
      for(int outer=0;outer<int(Value(cfg,"maximumOuterIterations"));++outer){
        for(unsigned j=0;j<engine.slices.size();++j)FitLocal(engine,engine.slices[j],g);
        Vec priorCoefficients=g;double before=engine.Objective(g);
        GlobalObjective objective;objective.engine=&engine;
        ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");if(!m)throw std::runtime_error("cannot create global Minuit2");
        ROOT::Math::Functor f(objective,14);m->SetFunction(f);m->SetPrintLevel(-1);m->SetMaxFunctionCalls(int(Value(cfg,"globalMaxCalls")));m->SetTolerance(Value(cfg,"globalTolerance"));
        for(int j=0;j<14;++j){double span=Value(cfg,"coefficientBoundSteps")*coefficientScales[j%7];m->SetLimitedVariable(j,Form("coefficient_%d",j),g[j],coefficientScales[j%7],engine.legacy[j]-span,engine.legacy[j]+span);}
        bool ok=m->Minimize();Vec next(m->X(),m->X()+14);double proposed=engine.Objective(next);
        globalAccepted=engine.Feasible(next)&&proposed<1.e19&&proposed<=before;
        if(globalAccepted)g=next;
        globalStatus=m->Status();globalEdm=m->Edm();globalMinimizerConverged=ok&&globalStatus==0&&globalAccepted;
        covStatus=m->CovMatrixStatus();covariancePointMatches=globalAccepted&&covStatus>0;covariancePoint=next;
        if(covariancePointMatches)for(int j=0;j<14;++j)for(int k=0;k<14;++k){covariance[j][k]=m->CovMatrix(j,k);if(!TMath::Finite(covariance[j][k]))covariancePointMatches=false;}
        delete m;
        for(unsigned j=0;j<engine.slices.size();++j)RefreshAmplitude(engine,engine.slices[j],g);
        double dev=engine.Objective(g);coefficientDrift=CoefficientDistance(g,priorCoefficients);
        bool stableDeviance=previous<1.e19&&std::fabs(previous-dev)<=Value(cfg,"outerRelativeTolerance")*std::max(1.,std::fabs(dev));
        bool stableCoefficients=coefficientDrift<=Value(cfg,"outerCoefficientTolerance");
        outerConverged=stableDeviance&&stableCoefficients&&AllAmplitudeConverged(engine)&&AllLocalConverged(engine)&&globalMinimizerConverged;
        trace.push_back(dev);std::ostringstream record;record<<std::setprecision(17)
          <<"{\"iteration\":"<<outer<<",\"deviance\":"<<dev<<",\"globalStatus\":"<<globalStatus<<",\"globalPointAccepted\":"<<(globalAccepted?"true":"false")
          <<",\"globalMinimizerConverged\":"<<(globalMinimizerConverged?"true":"false")<<",\"globalEdm\":"<<globalEdm
          <<",\"allAmplitudeConverged\":"<<(AllAmplitudeConverged(engine)?"true":"false")
          <<",\"allLocalMinimizerConvergedBeforeGlobalStep\":"<<(AllLocalConverged(engine)?"true":"false")
          <<",\"scaledCoefficientDrift\":"<<coefficientDrift<<",\"stableDeviance\":"<<(stableDeviance?"true":"false")<<",\"stableCoefficients\":"<<(stableCoefficients?"true":"false")<<"}";
        outerRecords.push_back(record.str());
        std::cout<<"GLOBAL seed="<<seedName<<" outer="<<outer<<" deviance="<<dev<<" status="<<globalStatus<<" AMP="<<AllAmplitudeConverged(engine)<<" local="<<AllLocalConverged(engine)<<" step="<<coefficientDrift<<" He4mean(.41)="<<Mean(g,3,.41)<<" He4sigma(.41)="<<Sigma(g,3,.41)<<std::endl;
        if(outerConverged)break;previous=dev;
      }
      // Final local stationarity probe does not mutate retained local shapes:
      // covariance remains evaluated at exactly the exported global+shape point.
      Engine localProbe=engine;double largestLocalImprovement=0;
      for(unsigned j=0;j<localProbe.slices.size();++j){double before=engine.slices[j].localDeviance;FitLocal(localProbe,localProbe.slices[j],g);largestLocalImprovement=std::max(largestLocalImprovement,std::max(0.,before-localProbe.slices[j].localDeviance)/std::max(1.,engine.slices[j].count));}
      bool localStationarity=AllLocalConverged(localProbe)&&largestLocalImprovement<=Value(cfg,"stationarityRelativeTolerance");
      Stationarity gradient=GlobalStationarity(engine,g);double normalizedGradient=gradient.residual/std::max(1.,gradient.count);
      bool globalStationarity=gradient.evaluable&&normalizedGradient<=Value(cfg,"stationarityRelativeTolerance");
      bool converged=outerConverged&&AllAmplitudeConverged(engine)&&localStationarity&&globalStationarity;
      if(runIndex++)json<<',';
      json<<"{\"seed\":"<<Quote(seedName)<<",\"converged\":"<<(converged?"true":"false")<<",\"outerConverged\":"<<(outerConverged?"true":"false")
        <<",\"allAmplitudeConverged\":"<<(AllAmplitudeConverged(engine)?"true":"false")
        <<",\"allLocalMinimizerConvergedBeforeGlobalStep\":"<<(AllLocalConverged(engine)?"true":"false")
        <<",\"finalLocalProbeMinimizerConverged\":"<<(AllLocalConverged(localProbe)?"true":"false")
        <<",\"finalLocalStationarity\":"<<(localStationarity?"true":"false")<<",\"largestFinalLocalImprovementPerCount\":"<<largestLocalImprovement
        <<",\"finalGlobalStationarity\":"<<(globalStationarity?"true":"false")<<",\"finalGlobalStationarityEvaluable\":"<<(gradient.evaluable?"true":"false")
        <<",\"finalGlobalStationarityAbsoluteResidual\":"<<gradient.residual<<",\"finalGlobalStationarityResidualPerCount\":"<<normalizedGradient
        <<",\"stationarityNormalizationCount\":"<<gradient.count<<",\"stationarityCentralDirections\":"<<gradient.central<<",\"stationarityOneSidedDirections\":"<<gradient.oneSided<<",\"stationarityBlockedDirections\":"<<gradient.blocked
        <<",\"globalStatus\":"<<globalStatus<<",\"globalPointAccepted\":"<<(globalAccepted?"true":"false")<<",\"globalMinimizerConverged\":"<<(globalMinimizerConverged?"true":"false")
        <<",\"globalEdm\":"<<globalEdm<<",\"conditionalCovarianceStatus\":"<<covStatus<<",\"conditionalCovariancePointMatches\":"<<(covariancePointMatches?"true":"false")
        <<",\"conditionalCovarianceAccurate\":"<<((covariancePointMatches&&covStatus==3&&globalMinimizerConverged)?"true":"false")
        <<",\"deviance\":"<<engine.Objective(g)<<",\"coefficients\":";WriteVec(json,g);json<<",\"trace\":";WriteVec(json,trace);
      json<<",\"outerIterations\":[";for(unsigned j=0;j<outerRecords.size();++j){if(j)json<<',';json<<outerRecords[j];}
      json<<"],\"conditionalCovarianceEvaluationPoint\":";if(covariancePointMatches)WriteVec(json,covariancePoint);else json<<"null";
      json<<",\"conditionalCovariance\":";if(covariancePointMatches){json<<'[';for(int j=0;j<14;++j){if(j)json<<',';WriteVec(json,covariance[j]);}json<<']';}else json<<"null";
      json<<",\"slices\":[";
      for(unsigned j=0;j<engine.slices.size();++j){
        const Slice& s=engine.slices[j];if(j)json<<',';Matrix templates=engine.Templates(s,g,s.shape);Vec acceptance(s.amplitudes.size(),0),windowYield(s.amplitudes.size(),0);
        for(unsigned b=0;b<templates.size();++b)for(unsigned k=0;k<acceptance.size();++k)acceptance[k]+=templates[b][k];
        for(unsigned k=0;k<acceptance.size();++k)windowYield[k]=s.count*s.amplitudes[k]*acceptance[k];
        json<<"{\"index\":"<<s.index<<",\"classId\":"<<s.classId<<",\"lo\":"<<s.lo<<",\"hi\":"<<s.hi<<",\"r\":"<<s.r<<",\"fitLo\":"<<s.fitLo<<",\"fitHi\":"<<s.fitHi<<",\"count\":"<<s.count<<",\"deviance\":"<<s.localDeviance
          <<",\"localStatus\":"<<s.localStatus<<",\"localEdm\":"<<s.localEdm<<",\"localPointAccepted\":"<<(s.localPointAccepted?"true":"false")
          <<",\"localMinimizerConvergedBeforeGlobalStep\":"<<(s.localMinimizerConverged?"true":"false")<<",\"amplitudeIterations\":"<<s.amplitudeIterations
          <<",\"amplitudeKktResidual\":"<<s.amplitudeKktResidual<<",\"amplitudeConverged\":"<<(s.amplitudeConverged?"true":"false")
          <<",\"finalLocalProbeStatus\":"<<localProbe.slices[j].localStatus<<",\"finalLocalProbeConverged\":"<<(localProbe.slices[j].localMinimizerConverged?"true":"false")
          <<",\"He3meanDerived\":"<<Mean(g,2,s.r)<<",\"He3sigmaDerived\":"<<Sigma(g,2,s.r)<<",\"He4meanDerived\":"<<Mean(g,3,s.r)<<",\"He4sigmaDerived\":"<<Sigma(g,3,s.r)
          <<",\"amplitudesFraction\":";WriteVec(json,s.amplitudes);json<<",\"componentWindowAcceptance\":";WriteVec(json,acceptance);json<<",\"componentExpectedWindowYield\":";WriteVec(json,windowYield);
        json<<",\"componentNumericallyInvisible\":[";for(unsigned k=0;k<acceptance.size();++k){if(k)json<<',';json<<(acceptance[k]<1.e-15?"true":"false");}json<<"],\"neighborBackgroundShape\":";WriteVec(json,s.shape);json<<'}';
      }json<<"]}";
    }
    if(!runIndex)throw std::runtime_error("no shared seeds supplied");
    json<<"]}\n";json.close();std::cout<<"GLOBAL_HELIUM_COMPLETE runs="<<runIndex<<" slices="<<base.slices.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"GlobalHeliumFit: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
