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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace GlobalHelium {
typedef std::vector<double> Vec;
typedef std::vector<Vec> Matrix;
typedef std::map<std::string,std::string> Settings;
double Value(const Settings& cfg,const std::string& key) {Settings::const_iterator i=cfg.find(key);if(i==cfg.end())throw std::runtime_error("missing setting "+key);std::istringstream s(i->second);double v;s>>v;if(!s||!TMath::Finite(v))throw std::runtime_error("bad setting "+key);return v;}
std::string Text(const Settings& cfg,const std::string& key){Settings::const_iterator i=cfg.find(key);if(i==cfg.end())throw std::runtime_error("missing text setting "+key);return i->second;}
Settings ReadSettings(const char* path){std::ifstream f(path);if(!f)throw std::runtime_error("missing settings");Settings c;std::string line;while(std::getline(f,line)){size_t i=line.find('\t');if(i==std::string::npos)throw std::runtime_error("bad settings TSV");c[line.substr(0,i)]=line.substr(i+1);}return c;}
std::string Quote(const std::string& s){std::ostringstream o;o<<'"';for(unsigned j=0;j<s.size();++j){if(s[j]=='"'||s[j]=='\\')o<<'\\';o<<s[j];}o<<'"';return o.str();}
double Poly(const double* c,int n,double x){double y=0;for(int i=n-1;i>=0;--i)y=y*x+c[i];return y;}
double Extreme(const double* c,int degree,double lo,double hi,bool minimum){Vec x;x.push_back(lo);x.push_back(hi);if(degree==2&&std::fabs(c[2])>1.e-14)x.push_back(-c[1]/(2*c[2]));if(degree==3){double a=3*c[3],b=2*c[2],d=c[1];if(std::fabs(a)<1.e-14){if(std::fabs(b)>1.e-14)x.push_back(-d/b);}else if(b*b-4*a*d>=0){double q=std::sqrt(b*b-4*a*d);x.push_back((-b+q)/(2*a));x.push_back((-b-q)/(2*a));}}double v=Poly(c,degree+1,lo);for(unsigned j=0;j<x.size();++j)if(x[j]>=lo&&x[j]<=hi){double p=Poly(c,degree+1,x[j]);v=minimum?std::min(v,p):std::max(v,p);}return v;}
double Mean(const Vec& g,int species,double r){return Poly(&g[(species-2)*7],4,std::log10(r));}
double Sigma(const Vec& g,int species,double r){return Poly(&g[(species-2)*7+4],3,std::log10(r));}
double Gaussian(double lo,double hi,double mean,double sigma){if(!(sigma>0))return 0;return .5*(std::erf((hi-mean)/(std::sqrt(2.)*sigma))-std::erf((lo-mean)/(std::sqrt(2.)*sigma)));}
Vec Legacy(){Vec g;for(int s=2;s<4;++s){const NuclearIdDeDxVsMom::LogPolyParams& m=NuclearIdDeDxVsMom::Mean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s));const NuclearIdDeDxVsMom::LogPolyParams& w=NuclearIdDeDxVsMom::Sigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s));g.push_back(m.c0);g.push_back(m.c1);g.push_back(m.c2);g.push_back(m.c3);g.push_back(w.c0);g.push_back(w.c1);g.push_back(w.c2);}return g;}
struct Bin {double lo,hi,n;};
struct Slice {double lo,hi,r,count,fitLo,fitHi;int classId,index;std::vector<Bin> bins;Vec nodes,weights,shape,amplitudes;int localStatus,amplitudeIterations;double localDeviance;Slice():count(0),classId(0),index(0),localStatus(-1),amplitudeIterations(0),localDeviance(0){shape.assign(5,0);shape[1]=shape[3]=1;}};
struct Profile {Vec amplitude;double deviance;int iterations;bool converged;Profile():deviance(0),iterations(0),converged(false){}};
bool Solve(Matrix h,Vec b,Vec& x){int n=b.size();x.assign(n,0);for(int i=0;i<n;++i){int pivot=i;for(int j=i+1;j<n;++j)if(std::fabs(h[j][i])>std::fabs(h[pivot][i]))pivot=j;if(std::fabs(h[pivot][i])<1.e-16)return false;std::swap(h[i],h[pivot]);std::swap(b[i],b[pivot]);double d=h[i][i];for(int k=i;k<n;++k)h[i][k]/=d;b[i]/=d;for(int j=0;j<n;++j)if(j!=i){double a=h[j][i];for(int k=i;k<n;++k)h[j][k]-=a*h[i][k];b[j]-=a*b[i];}}x=b;return true;}
double Deviance(const Slice& s,const Matrix& t,const Vec& a){double d=0;for(unsigned b=0;b<t.size();++b){double v=0;for(unsigned k=0;k<a.size();++k)v+=a[k]*t[b][k];double expectation=std::max(1.e-20,s.count*v),n=s.bins[b].n;d+=2*(expectation-n+(n>0?n*std::log(n/expectation):0));}return d;}
class Engine {
 public:
  Settings cfg;std::vector<Slice> slices;bool exponential;int quadrature;Vec legacy;
  Engine(const Settings& c):cfg(c),exponential(Text(c,"backgroundModel")=="exponential"),quadrature(int(Value(c,"rQuadratureOrder"))),legacy(Legacy()){}
  bool Feasible(const Vec& g)const{
    double low=std::log10(Value(cfg,"rigidityMin")),high=std::log10(Value(cfg,"rigidityMax"));
    for(int s=0;s<2;++s){const double* m=&g[s*7];const double* w=&g[s*7+4];if(Extreme(m,3,low,high,true)<=0||Extreme(w,2,low,high,true)<Value(cfg,"minimumSigma"))return false;if(Value(cfg,"requireMonotonicMean")){double d[]={m[1],2*m[2],3*m[3]};if(Extreme(d,2,low,high,false)>0)return false;}double shiftLow[4],shiftHigh[4],width[3];for(int j=0;j<4;++j){double sig=j<3?legacy[s*7+4+j]:0;shiftLow[j]=m[j]-legacy[s*7+j]-Value(cfg,"meanShiftMinSigma")*sig;shiftHigh[j]=m[j]-legacy[s*7+j]-Value(cfg,"meanShiftMaxSigma")*sig;}for(int j=0;j<3;++j)width[j]=w[j]-Value(cfg,"maximumSigmaScale")*legacy[s*7+4+j];if(Extreme(shiftLow,3,low,high,true)<0||Extreme(shiftHigh,3,low,high,false)>0||Extreme(width,2,low,high,false)>0)return false;}
    double order[4];for(int j=0;j<4;++j)order[j]=g[7+j]-g[j];return Extreme(order,3,low,high,true)>=Value(cfg,"minimumMeanGap");
  }
  bool ShapeFeasible(const Slice& s,const Vec& g,const Vec& shape)const{for(unsigned j=0;j<s.nodes.size();++j){double r=s.nodes[j],md=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kDeuteron,r)+shape[0]*NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kDeuteron,r),mt=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kTriton,r)+shape[2]*NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kTriton,r);if(mt-md<Value(cfg,"minimumMeanGap")||Mean(g,2,r)-mt<Value(cfg,"minimumMeanGap"))return false;}return true;}
  void Nodes(Slice& s)const{const double x3[]={-.7745966692414834,0,.7745966692414834},w3[]={5./18,4./9,5./18};const double x5[]={-.9061798459386640,-.5384693101056831,0,.5384693101056831,.9061798459386640},w5[]={.1184634425280945,.2393143352496832,.2844444444444444,.2393143352496832,.1184634425280945};if(quadrature==1){s.nodes.push_back(s.r);s.weights.push_back(1);}else if(quadrature==3||quadrature==5){for(int j=0;j<quadrature;++j){s.nodes.push_back(s.r+(s.hi-s.lo)/2*(quadrature==3?x3[j]:x5[j]));s.weights.push_back(quadrature==3?w3[j]:w5[j]);}}else throw std::runtime_error("quadrature must be 1,3,5");}
  Matrix Templates(const Slice& s,const Vec& g,const Vec& shape)const{
    const int n=exponential?5:6;Matrix t(s.bins.size(),Vec(n,0));double width=s.fitHi-s.fitLo;
    for(unsigned b=0;b<s.bins.size();++b){const Bin& bin=s.bins[b];for(int species=0;species<4;++species){for(unsigned q=0;q<s.nodes.size();++q){double r=s.nodes[q],m,w;if(species<2){NuclearIdDeDxVsMom::SpeciesIndex sp=static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(species);w=NuclearIdDeDxVsMom::GetSigma(sp,r);m=NuclearIdDeDxVsMom::GetMean(sp,r)+shape[2*species]*w;w*=shape[2*species+1];}else{m=Mean(g,species,r);w=Sigma(g,species,r);}t[b][species]+=s.weights[q]*Gaussian(bin.lo,bin.hi,m,w);}}
      if(exponential){double k=shape[4];if(std::fabs(k)<1.e-8)t[b][4]=(bin.hi-bin.lo)/width;else t[b][4]=(std::exp(k*(bin.hi-s.fitLo)/width)-std::exp(k*(bin.lo-s.fitLo)/width))/(std::exp(k)-1.);}
      else{double center=(bin.lo+bin.hi)/2;t[b][4]=2*(bin.hi-bin.lo)*(s.fitHi-center)/(width*width);t[b][5]=2*(bin.hi-bin.lo)*(center-s.fitLo)/(width*width);}}
    return t;
  }
  Profile Amplitudes(const Slice& s,const Matrix& t,const Vec& initial)const{
    const int n=t[0].size();Profile out;out.amplitude=initial.size()==unsigned(n)?initial:Vec(n,1./n);Vec& a=out.amplitude;double cap=Value(cfg,"maximumAmplitudeFraction");for(int j=0;j<n;++j)a[j]=std::min(cap,std::max(1.e-8,a[j]));double previous=Deviance(s,t,a);
    for(int iteration=0;iteration<int(Value(cfg,"amplitudeMaxIterations"));++iteration){Vec gradient(n,0);Matrix hessian(n,Vec(n,0));for(unsigned b=0;b<t.size();++b){double v=0;for(int j=0;j<n;++j)v+=a[j]*t[b][j];v=std::max(1.e-20,v);double observed=s.bins[b].n/s.count;for(int j=0;j<n;++j){gradient[j]+=t[b][j]*(1-observed/v);for(int k=0;k<n;++k)hessian[j][k]+=observed*t[b][j]*t[b][k]/(v*v);}}
      std::vector<int> active;double violation=0;for(int j=0;j<n;++j){bool blocked=(a[j]<1.e-9&&gradient[j]>0)||(a[j]>cap-1.e-9&&gradient[j]<0);if(!blocked){active.push_back(j);violation=std::max(violation,std::fabs(gradient[j]));}}
      out.iterations=iteration+1;if(violation<Value(cfg,"amplitudeGradientTolerance")){out.converged=true;break;}if(active.empty()){out.converged=true;break;}
      Matrix h(active.size(),Vec(active.size()));Vec rhs(active.size());for(unsigned j=0;j<active.size();++j){rhs[j]=-gradient[active[j]];for(unsigned k=0;k<active.size();++k)h[j][k]=hessian[active[j]][active[k]];h[j][j]+=1.e-10;}
      Vec direction;if(!Solve(h,rhs,direction))break;double step=1.;for(unsigned j=0;j<active.size();++j){double d=direction[j],v=a[active[j]];if(d<0)step=std::min(step,-v/d);if(d>0)step=std::min(step,(cap-v)/d);}step=std::max(0.,step);bool accepted=false;
      for(int line=0;line<25;++line){Vec next=a;for(unsigned j=0;j<active.size();++j)next[active[j]]=std::min(cap,std::max(0.,a[active[j]]+step*direction[j]));double dev=Deviance(s,t,next);if(dev<previous){a=next;previous=dev;accepted=true;break;}step*=.5;}
      if(!accepted)break;
    }
    out.deviance=Deviance(s,t,a);return out;
  }
  double ProfileSlice(const Slice& s,const Vec& g,const Vec& shape)const{if(!ShapeFeasible(s,g,shape))return 1.e20;return Amplitudes(s,Templates(s,g,shape),s.amplitudes).deviance;}
  double Objective(const Vec& g)const{if(!Feasible(g))return 1.e20;double loss=0;for(unsigned j=0;j<slices.size();++j){double d=ProfileSlice(slices[j],g,slices[j].shape);if(d>=1.e19)return d;loss+=d;}return loss;}
};
struct GlobalObjective{const Engine* engine;double operator()(const double* p)const{return engine->Objective(Vec(p,p+14));}};
struct LocalObjective{const Engine* engine;const Slice* slice;const Vec* global;std::vector<int> free;double operator()(const double* p)const{Vec shape=slice->shape;for(unsigned j=0;j<free.size();++j)shape[free[j]]=p[j];return engine->ProfileSlice(*slice,*global,shape);}};
void FitLocal(Engine& e,Slice& s,const Vec& g){LocalObjective objective;objective.engine=&e;objective.slice=&s;objective.global=&g;if(Value(e.cfg,"floatDeuteronShape")){objective.free.push_back(0);objective.free.push_back(1);}if(Value(e.cfg,"floatTritonShape")){objective.free.push_back(2);objective.free.push_back(3);}if(e.exponential)objective.free.push_back(4);
  if(!objective.free.empty()){ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");ROOT::Math::Functor f(objective,objective.free.size());m->SetFunction(f);m->SetMaxFunctionCalls(int(Value(e.cfg,"localMaxCalls")));m->SetTolerance(Value(e.cfg,"localTolerance"));m->SetPrintLevel(-1);for(unsigned j=0;j<objective.free.size();++j){int index=objective.free[j];bool width=index==1||index==3;double low=width?Value(e.cfg,"neighborWidthMin"):-Value(e.cfg,index==4?"exponentialSlopeMax":"neighborMeanShiftMax");double high=width?Value(e.cfg,"neighborWidthMax"):Value(e.cfg,index==4?"exponentialSlopeMax":"neighborMeanShiftMax");m->SetLimitedVariable(j,Form("shape_%d",index),s.shape[index],.05,low,high);}double before=objective(&s.shape[0]); // actual baseline calculated below with correctly packed variables
    Vec packed;for(unsigned j=0;j<objective.free.size();++j)packed.push_back(s.shape[objective.free[j]]);before=objective(&packed[0]);m->Minimize();s.localStatus=m->Status();if(m->MinValue()<=before)for(unsigned j=0;j<objective.free.size();++j)s.shape[objective.free[j]]=m->X()[j];delete m;
  }else s.localStatus=0;
  Profile p=e.Amplitudes(s,e.Templates(s,g,s.shape),s.amplitudes);s.amplitudes=p.amplitude;s.amplitudeIterations=p.iterations;s.localDeviance=p.deviance;
}
void WriteVec(std::ostream& o,const Vec& v){o<<'[';for(unsigned j=0;j<v.size();++j){if(j)o<<',';o<<v[j];}o<<']';}
}
#endif

void GlobalHeliumFit(const char* input,const char* output,const char* settings,const char* seeds) {
#ifndef __CINT__
  using namespace GlobalHelium;
  try {
    const Settings cfg=ReadSettings(settings);const TString out(output);if(!gSystem->AccessPathName(out+"/global_fit.json"))throw std::runtime_error("output already exists");gSystem->mkdir(out,kTRUE);gSystem->Load("libMinuit2");
    TFile source(input,"READ");if(source.IsZombie())throw std::runtime_error("cannot open source ROOT");TH2* all=dynamic_cast<TH2*>(source.Get(Text(cfg,"histogram").c_str()));TH2* tag=dynamic_cast<TH2*>(source.Get(Text(cfg,"tagHistogram").c_str()));if(!all)throw std::runtime_error("missing prePID TH2");bool disjoint=Text(cfg,"sampleMode")=="tag_not_tag";if(disjoint&&!tag)throw std::runtime_error("missing TOF tag histogram");
    std::vector<int> indices;std::istringstream list(Text(cfg,"sliceIndices"));std::string token;while(std::getline(list,token,','))indices.push_back(std::atoi(token.c_str()));const double rmin=Value(cfg,"rigidityMin"),rmax=Value(cfg,"rigidityMax"),step=Value(cfg,"rigidityStep");int total=int((rmax-rmin)/step+.5);if(indices.empty())for(int j=0;j<total;++j)indices.push_back(j);Engine base(cfg);TFile saved(out+"/input_slices.root","RECREATE");
    for(unsigned a=0;a<indices.size();++a){int index=indices[a];if(index<0||index>=total)throw std::runtime_error("invalid slice index");double lo=rmin+index*step,hi=lo+step,r=(lo+hi)/2;int first=all->GetXaxis()->FindFixBin(lo+1.e-9),last=all->GetXaxis()->FindFixBin(hi-1.e-9);TH1D* raw=all->ProjectionY(Form("all_native_%03d",index),first,last,"e");raw->SetDirectory(0);saved.cd();raw->Write();TH1D* tagged=0;if(disjoint){tagged=tag->ProjectionY(Form("tag_native_%03d",index),first,last,"e");tagged->SetDirectory(0);tagged->Write();for(int b=0;b<=raw->GetNbinsX()+1;++b)if(tagged->GetBinContent(b)>raw->GetBinContent(b)+1.e-6)throw std::runtime_error("TOF tag not a subset of ALL");}
      const int classes=disjoint?2:1;for(int classIndex=0;classIndex<classes;++classIndex){TH1D* h=static_cast<TH1D*>((disjoint&&classIndex==0?tagged:raw)->Clone(Form("fitraw_%03d_class%d",index,classIndex)));h->SetDirectory(0);if(disjoint&&classIndex==1){h->Add(tagged,-1);for(int b=0;b<=h->GetNbinsX()+1;++b)h->SetBinError(b,std::sqrt(std::max(0.,h->GetBinContent(b))));}h->Rebin(int(Value(cfg,"dedxRebinFactor")));h->Write();Slice s;s.lo=lo;s.hi=hi;s.r=r;s.index=index;s.classId=disjoint?classIndex+1:0;base.Nodes(s);double oldMean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),oldSigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);double requestLow=std::max(0.,oldMean+Value(cfg,"fitNSigmaMin")*oldSigma),requestHigh=std::min(120.,oldMean+Value(cfg,"fitNSigmaMax")*oldSigma);int bfirst=h->FindBin(requestLow+1.e-9),blast=h->FindBin(requestHigh-1.e-9);s.fitLo=h->GetBinLowEdge(bfirst);s.fitHi=h->GetBinLowEdge(blast)+h->GetBinWidth(blast);for(int b=bfirst;b<=blast;++b){Bin v;v.lo=h->GetBinLowEdge(b);v.hi=v.lo+h->GetBinWidth(b);v.n=h->GetBinContent(b);if(v.n<0)throw std::runtime_error("negative count");s.count+=v.n;s.bins.push_back(v);}if(s.count>=Value(cfg,"minimumSliceCount"))base.slices.push_back(s);delete h;}delete raw;delete tagged;}
    saved.Close();if(base.slices.empty())throw std::runtime_error("no slices");std::ifstream seedInput(seeds);if(!seedInput)throw std::runtime_error("missing explicit shared seeds");std::ofstream json((out+"/global_fit.json").Data());json<<std::setprecision(17)<<"{\"schemaVersion\":1,\"globalParameterCount\":14,\"input\":"<<Quote(input)<<",\"sampleMode\":"<<Quote(Text(cfg,"sampleMode"))<<",\"background\":"<<Quote(Text(cfg,"backgroundModel"))<<",\"rQuadratureOrder\":"<<base.quadrature<<",\"rWithinBinPolicy\":\"uniform GL-average; source has no sub-bin track momenta\",\"covarianceScope\":\"global14 conditional on T/d and nonlinear background shapes; amplitudes profiled; NOT full nuisance covariance\",\"runs\":[";
    std::string seedName;int runIndex=0;while(seedInput>>seedName){Vec g(14);for(int j=0;j<14;++j)if(!(seedInput>>g[j]))throw std::runtime_error("invalid seed");if(!base.Feasible(g))throw std::runtime_error("infeasible shared seed "+seedName);Engine engine=base;Matrix covariance(14,Vec(14,0));std::vector<double> trace;int globalStatus=-1,covStatus=-1;bool converged=false;double previous=1.e20;
      for(int outer=0;outer<int(Value(cfg,"maximumOuterIterations"));++outer){for(unsigned j=0;j<engine.slices.size();++j)FitLocal(engine,engine.slices[j],g);double before=engine.Objective(g);GlobalObjective objective;objective.engine=&engine;ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");ROOT::Math::Functor f(objective,14);m->SetFunction(f);m->SetPrintLevel(-1);m->SetMaxFunctionCalls(int(Value(cfg,"globalMaxCalls")));m->SetTolerance(Value(cfg,"globalTolerance"));const double scales[]={1.,3.,5.,10.,.1,.3,.5};for(int j=0;j<14;++j){double span=Value(cfg,"coefficientBoundSteps")*scales[j%7];m->SetLimitedVariable(j,Form("coefficient_%d",j),g[j],scales[j%7],engine.legacy[j]-span,engine.legacy[j]+span);}m->Minimize();Vec next(m->X(),m->X()+14);if(engine.Feasible(next)&&m->MinValue()<=before)g=next;globalStatus=m->Status();covStatus=m->CovMatrixStatus();for(int j=0;j<14;++j)for(int k=0;k<14;++k)covariance[j][k]=m->CovMatrix(j,k);delete m;for(unsigned j=0;j<engine.slices.size();++j){Profile p=engine.Amplitudes(engine.slices[j],engine.Templates(engine.slices[j],g,engine.slices[j].shape),engine.slices[j].amplitudes);engine.slices[j].amplitudes=p.amplitude;engine.slices[j].localDeviance=p.deviance;}
        double dev=engine.Objective(g);trace.push_back(dev);std::cout<<"GLOBAL seed="<<seedName<<" outer="<<outer<<" deviance="<<dev<<" status="<<globalStatus<<" He4mean(.41)="<<Mean(g,3,.41)<<" He4sigma(.41)="<<Sigma(g,3,.41)<<std::endl;if(std::fabs(previous-dev)<=Value(cfg,"outerRelativeTolerance")*std::max(1.,std::fabs(dev))){converged=true;break;}previous=dev;
      }
      if(runIndex++)json<<',';json<<"{\"seed\":"<<Quote(seedName)<<",\"converged\":"<<(converged?"true":"false")<<",\"globalStatus\":"<<globalStatus<<",\"conditionalCovarianceStatus\":"<<covStatus<<",\"deviance\":"<<engine.Objective(g)<<",\"coefficients\":";WriteVec(json,g);json<<",\"trace\":";WriteVec(json,trace);json<<",\"conditionalCovariance\":[";for(int j=0;j<14;++j){if(j)json<<',';WriteVec(json,covariance[j]);}json<<"],\"slices\":[";
      for(unsigned j=0;j<engine.slices.size();++j){const Slice& s=engine.slices[j];if(j)json<<',';json<<"{\"index\":"<<s.index<<",\"classId\":"<<s.classId<<",\"lo\":"<<s.lo<<",\"hi\":"<<s.hi<<",\"r\":"<<s.r<<",\"fitLo\":"<<s.fitLo<<",\"fitHi\":"<<s.fitHi<<",\"count\":"<<s.count<<",\"deviance\":"<<s.localDeviance<<",\"localStatus\":"<<s.localStatus<<",\"He3meanDerived\":"<<Mean(g,2,s.r)<<",\"He3sigmaDerived\":"<<Sigma(g,2,s.r)<<",\"He4meanDerived\":"<<Mean(g,3,s.r)<<",\"He4sigmaDerived\":"<<Sigma(g,3,s.r)<<",\"amplitudesFraction\":";WriteVec(json,s.amplitudes);json<<",\"neighborBackgroundShape\":";WriteVec(json,s.shape);json<<'}';}json<<"]}";
    }
    json<<"]}\n";json.close();std::cout<<"GLOBAL_HELIUM_COMPLETE runs="<<runIndex<<" slices="<<base.slices.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"GlobalHeliumFit: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
