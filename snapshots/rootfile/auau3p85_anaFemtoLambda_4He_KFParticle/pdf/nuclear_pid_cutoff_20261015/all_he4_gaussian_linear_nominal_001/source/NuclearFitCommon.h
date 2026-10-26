#ifndef STAR_NUCLEAR_FIT_COMMON_H
#define STAR_NUCLEAR_FIT_COMMON_H
// Study-only bin-integral Poisson fits. No legacy bounds, yield caps or PID writes.
#include "TH1D.h"
#include "TMath.h"
#include "TSystem.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace NuclearFit {
typedef std::vector<double> Vec;
typedef std::vector<Vec> Matrix;
typedef std::map<std::string,std::string> Settings;
const double invalid=1.e100;
inline void Require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
inline bool Finite(double value){return TMath::Finite(value);}
inline double Number(const Settings& c,const std::string& key){
  Settings::const_iterator i=c.find(key);Require(i!=c.end(),"Missing numeric setting: "+key);
  std::istringstream in(i->second);double value=0;Require(bool(in>>value),"Bad setting: "+key);in>>std::ws;
  Require(in.eof()&&Finite(value),"Nonfinite/non-numeric setting: "+key);return value;
}
inline std::string Text(const Settings& c,const std::string& key){Settings::const_iterator i=c.find(key);Require(i!=c.end()&&!i->second.empty(),"Missing text setting: "+key);return i->second;}
inline int Integer(const Settings& c,const std::string& key){double n=Number(c,key);Require(n==std::floor(n)&&n>=0&&n<1.e9,"Noninteger setting: "+key);return int(n);}
inline Settings ReadSettings(const char* path){std::ifstream in(path);Require(bool(in),"Missing explicit settings TSV");Settings out;std::string line;
  while(std::getline(in,line)){size_t p=line.find('\t');Require(p!=std::string::npos&&line.find('\t',p+1)==std::string::npos,"Malformed settings row");
    Require(out.insert(std::make_pair(line.substr(0,p),line.substr(p+1))).second,"Duplicate setting key");}
  Require(in.eof(),"Settings read failure");return out;}
inline void ValidateSettings(const Settings& c){
  const char* positive[]={"amplitudeGradientTolerance","amplitudeHessianRidge","shapeTolerance","minimumSigma","maximumSigma","minimumGaussianAcceptance","meanStep","sigmaLogStep","backgroundSlopeStep","exponentialSlopeLimit","boundaryTolerance"};
  for(unsigned i=0;i<sizeof(positive)/sizeof(positive[0]);++i)Require(Number(c,positive[i])>0,"Nonpositive mathematical setting");
  const char* counts[]={"amplitudeMaxIterations","amplitudeLineSearchIterations","shapeMaxCalls"};
  for(unsigned i=0;i<sizeof(counts)/sizeof(counts[0]);++i)Require(Integer(c,counts[i])>0,"Zero iteration/call limit");
  Require(Number(c,"meanAxisPadding")>=0&&Number(c,"minimumSigma")<Number(c,"maximumSigma"),"Reversed numeric support");
}
inline std::string Quote(const std::string& s){std::ostringstream out;out<<'"';for(unsigned i=0;i<s.size();++i){if(s[i]=='"'||s[i]=='\\')out<<'\\';if(s[i]=='\n')out<<"\\n";else out<<s[i];}out<<'"';return out.str();}
inline void WriteVec(std::ostream& out,const Vec& v){out<<'[';for(unsigned i=0;i<v.size();++i){if(i)out<<',';out<<v[i];}out<<']';}
inline void WriteMatrix(std::ostream& out,const Matrix& m){out<<'[';for(unsigned i=0;i<m.size();++i){if(i)out<<',';WriteVec(out,m[i]);}out<<']';}
struct Bin {double lo,hi,count;Bin(double l=0,double h=0,double n=0):lo(l),hi(h),count(n){}};
struct Node {double weight,mean,sigma;Node(double w=1,double m=0,double s=1):weight(w),mean(m),sigma(s){}};
struct Component {
  std::string name;double mean,sigma,startWindowYield;bool active,freeMean,freeSigma;std::vector<Node> nodes;
  Component():mean(0),sigma(1),startWindowYield(0),active(true),freeMean(true),freeSigma(true){}
};
struct Request {
  std::vector<Bin> bins;std::vector<Component> components;std::string backgroundModel;
  double backgroundSlope,backgroundWindowYield,nativeAxisMin,nativeAxisMax;bool freeBackgroundSlope;
  Request():backgroundModel("linear"),backgroundSlope(0),backgroundWindowYield(0),nativeAxisMin(0),nativeAxisMax(120),freeBackgroundSlope(true){}
};
inline double GaussianIntegral(double lo,double hi,double mean,double sigma){
  if(!(sigma>0)||!Finite(mean)||!Finite(sigma)||!(hi>lo))return 0;
  double a=(lo-mean)/(std::sqrt(2.)*sigma),b=(hi-mean)/(std::sqrt(2.)*sigma),p;
  if(a>=0)p=.5*(std::erfc(a)-std::erfc(b));else if(b<=0)p=.5*(std::erfc(-b)-std::erfc(-a));else p=.5*(std::erf(b)-std::erf(a));
  return std::max(0.,p);
}
inline double BackgroundIntegral(double lo,double hi,double fitLo,double fitHi,const std::string& model,double slope,int basis){
  lo=std::max(lo,fitLo);hi=std::min(hi,fitHi);if(!(hi>lo&&fitHi>fitLo))return 0;
  double width=fitHi-fitLo,a=(lo-fitLo)/width,b=(hi-fitLo)/width;
  if(model=="linear")return basis==0?2*(b-a)-(b*b-a*a):(b*b-a*a);
  Require(model=="exponential"&&basis==0,"Unknown background basis");
  return std::fabs(slope)<1.e-8?(b-a):std::exp(slope*a)*std::expm1(slope*(b-a))/std::expm1(slope);
}
inline double KernelIntegral(const Component& c,double lo,double hi){
  if(!c.active)return 0;if(c.nodes.empty())return GaussianIntegral(lo,hi,c.mean,c.sigma);
  double result=0;for(unsigned i=0;i<c.nodes.size();++i)result+=c.nodes[i].weight*GaussianIntegral(lo,hi,c.nodes[i].mean,c.nodes[i].sigma);return result;
}
inline double TotalCount(const std::vector<Bin>& bins){double n=0;for(unsigned i=0;i<bins.size();++i)n+=bins[i].count;return n;}
inline bool Solve(Matrix a,Vec b,Vec& result,double threshold=1.e-14){
  const unsigned n=b.size();if(a.size()!=n)return false;
  for(unsigned i=0;i<n;++i){if(a[i].size()!=n)return false;unsigned pivot=i;
    for(unsigned j=i+1;j<n;++j)if(std::fabs(a[j][i])>std::fabs(a[pivot][i]))pivot=j;
    if(!(std::fabs(a[pivot][i])>threshold)||!Finite(a[pivot][i]))return false;
    std::swap(a[i],a[pivot]);std::swap(b[i],b[pivot]);double scale=a[i][i];for(unsigned k=i;k<n;++k)a[i][k]/=scale;b[i]/=scale;
    for(unsigned j=0;j<n;++j)if(j!=i){double x=a[j][i];for(unsigned k=i;k<n;++k)a[j][k]-=x*a[i][k];b[j]-=x*b[i];}}
  result=b;for(unsigned i=0;i<n;++i)if(!Finite(result[i]))return false;return true;
}
inline double LogRemainder(double x){if(std::fabs(x)>1.e-4)return x-std::log1p(x);double p=x*x,result=0;for(int k=2;k<=10;++k){result+=(k%2? -1.:1.)*p/k;p*=x;}return result;}
inline double Deviance(const std::vector<Bin>& bins,const Matrix& templates,const Vec& fractions){
  const double count=TotalCount(bins);double result=0;
  for(unsigned b=0;b<bins.size();++b){double mu=0;for(unsigned j=0;j<fractions.size();++j)mu+=count*fractions[j]*templates[b][j];
    if(!Finite(mu)||mu<0)return invalid;if(bins[b].count==0){result+=2*mu;continue;}if(!(mu>0))return invalid;
    double x=(mu-bins[b].count)/bins[b].count;if(x<=-1)return invalid;result+=2*bins[b].count*LogRemainder(x);}
  return Finite(result)?result:invalid;
}
inline Matrix BuildTemplates(const Request& request,const Settings& cfg,Vec* acceptance=0){
  if(request.bins.empty())return Matrix();const double lo=request.bins.front().lo,hi=request.bins.back().hi;
  const unsigned n=request.components.size(),bg=request.backgroundModel=="linear"?2:1;Matrix t(request.bins.size(),Vec(n+bg,0));Vec accepted(n,0);
  for(unsigned j=0;j<n;++j)if(request.components[j].active){accepted[j]=KernelIntegral(request.components[j],lo,hi);
    if(!(accepted[j]>Number(cfg,"minimumGaussianAcceptance")))return Matrix();
    for(unsigned b=0;b<t.size();++b)t[b][j]=KernelIntegral(request.components[j],request.bins[b].lo,request.bins[b].hi)/accepted[j];}
  for(unsigned b=0;b<t.size();++b)for(unsigned j=0;j<bg;++j)t[b][n+j]=BackgroundIntegral(request.bins[b].lo,request.bins[b].hi,lo,hi,request.backgroundModel,request.backgroundSlope,j);
  for(unsigned j=0;j<n+bg;++j){double sum=0;for(unsigned b=0;b<t.size();++b){if(!Finite(t[b][j])||t[b][j]<0)return Matrix();sum+=t[b][j];}
    if(sum>0)for(unsigned b=0;b<t.size();++b)t[b][j]/=sum;}
  if(acceptance)*acceptance=accepted;return t;
}
struct ProfileResult {Vec fractions,gradient;Matrix hessian;double dataDeviance,kkt;int iterations;bool converged;
  ProfileResult():dataDeviance(invalid),kkt(invalid),iterations(0),converged(false){} };
inline void Derivatives(const std::vector<Bin>& bins,const Matrix& t,const Vec& a,Vec& gradient,Matrix& h){
  const unsigned n=a.size();const double count=TotalCount(bins);gradient.assign(n,0);h.assign(n,Vec(n,0));
  for(unsigned b=0;b<t.size();++b){double v=0;for(unsigned j=0;j<n;++j)v+=a[j]*t[b][j];v=std::max(1.e-100,v);double observed=bins[b].count/count;
    for(unsigned j=0;j<n;++j){gradient[j]+=t[b][j]*(1-observed/v);for(unsigned k=0;k<n;++k)h[j][k]+=observed*t[b][j]*t[b][k]/(v*v);}}
}
inline ProfileResult Profile(const std::vector<Bin>& bins,const Matrix& t,const Vec& initial,const Settings& cfg){
  ProfileResult p;if(t.empty()||!(TotalCount(bins)>0))return p;const unsigned n=t[0].size();p.fractions=initial.size()==n?initial:Vec(n,1./n);
  Vec sums(n,0);for(unsigned b=0;b<t.size();++b)for(unsigned j=0;j<n;++j)sums[j]+=t[b][j];
  for(unsigned j=0;j<n;++j)p.fractions[j]=sums[j]>0?std::max(0.,p.fractions[j]):0;
  double initialSum=0;unsigned nonzeroColumns=0;for(unsigned j=0;j<n;++j){initialSum+=p.fractions[j];nonzeroColumns+=sums[j]>0;}
  if(initialSum==0&&nonzeroColumns)for(unsigned j=0;j<n;++j)p.fractions[j]=sums[j]>0?1./nonzeroColumns:0;
  double previous=Deviance(bins,t,p.fractions);
  for(int iteration=0;iteration<Integer(cfg,"amplitudeMaxIterations");++iteration){
    Derivatives(bins,t,p.fractions,p.gradient,p.hessian);std::vector<unsigned> active;p.kkt=0;
    for(unsigned j=0;j<n;++j)if(sums[j]>0&&!(p.fractions[j]<1.e-12&&p.gradient[j]>0)){active.push_back(j);p.kkt=std::max(p.kkt,std::fabs(p.gradient[j]));}
    p.iterations=iteration+1;if(p.kkt<Number(cfg,"amplitudeGradientTolerance"))break;
    Matrix h(active.size(),Vec(active.size()));Vec rhs(active.size());
    for(unsigned j=0;j<active.size();++j){rhs[j]=-p.gradient[active[j]];for(unsigned k=0;k<active.size();++k)h[j][k]=p.hessian[active[j]][active[k]];h[j][j]+=Number(cfg,"amplitudeHessianRidge");}
    Vec direction;bool solved=Solve(h,rhs,direction),accepted=false;
    for(int strategy=0;strategy<2&&!accepted;++strategy){if(strategy||!solved){direction.resize(active.size());for(unsigned j=0;j<active.size();++j)direction[j]=-p.gradient[active[j]]/std::max(1.e-12,p.hessian[active[j]][active[j]]);}
      double step=1;for(unsigned j=0;j<active.size();++j)if(direction[j]<0)step=std::min(step,-p.fractions[active[j]]/direction[j]);
      for(int line=0;line<Integer(cfg,"amplitudeLineSearchIterations")&&step>1.e-16;++line){Vec next=p.fractions;
        for(unsigned j=0;j<active.size();++j)next[active[j]]=std::max(0.,next[active[j]]+step*direction[j]);double d=Deviance(bins,t,next);
        if(d<previous){p.fractions=next;previous=d;accepted=true;break;}step*=.5;}}
    if(!accepted)break;
  }
  Derivatives(bins,t,p.fractions,p.gradient,p.hessian);p.kkt=0;
  for(unsigned j=0;j<n;++j)if(sums[j]>0&&!(p.fractions[j]<1.e-12&&p.gradient[j]>0))p.kkt=std::max(p.kkt,std::fabs(p.gradient[j]));
  p.dataDeviance=Deviance(bins,t,p.fractions);p.converged=Finite(p.dataDeviance)&&p.dataDeviance<invalid&&p.kkt<Number(cfg,"amplitudeGradientTolerance");return p;
}
struct Result {
  std::vector<Component> components;std::vector<std::string> columnNames,shapeNames;std::vector<bool> yieldCovarianceColumnActive;Vec columnWindowYields,windowAcceptance,effectiveCount,shapePoint,meanError,sigmaError,conditionalDropDeltaD;
  Matrix totalTemplates,componentCounts,fixedShapeYieldCovariance,shapeCovariance;Vec totalCounts,smoothCounts;
  double dataDeviance,pearson,edm,amplitudeKkt,backgroundSlope;int status,covarianceStatus,amplitudeIterations,conditionalNdf;
  bool finite,amplitudeConverged,pointAccepted,atBoundary,yieldCovarianceValid,shapeCovarianceValid;
  Result():dataDeviance(invalid),pearson(0),edm(invalid),amplitudeKkt(invalid),backgroundSlope(0),status(-1),covarianceStatus(-1),amplitudeIterations(0),conditionalNdf(0),finite(false),amplitudeConverged(false),pointAccepted(false),atBoundary(false),yieldCovarianceValid(false),shapeCovarianceValid(false){}
};
class Objective {
 public:
  Request original;Settings cfg;Vec initial;std::vector<int> species,kind;std::vector<std::string> names;
  Objective(const Request& request,const Settings& settings):original(request),cfg(settings){
    const double count=TotalCount(request.bins);for(unsigned j=0;j<request.components.size();++j){const Component& c=request.components[j];initial.push_back(c.active?c.startWindowYield/count:0);
      if(c.active&&c.freeMean){species.push_back(j);kind.push_back(0);names.push_back(c.name+"_mean");}
      if(c.active&&c.freeSigma){species.push_back(j);kind.push_back(1);names.push_back(c.name+"_logSigma");}}
    initial.push_back(request.backgroundWindowYield/count/(request.backgroundModel=="linear"?2:1));if(request.backgroundModel=="linear")initial.push_back(initial.back());
    if(request.backgroundModel=="exponential"&&request.freeBackgroundSlope){species.push_back(-1);kind.push_back(2);names.push_back("backgroundSlope");}
  }
  Vec Point()const{Vec p;for(unsigned j=0;j<kind.size();++j)p.push_back(kind[j]==2?original.backgroundSlope:(kind[j]?std::log(original.components[species[j]].sigma):original.components[species[j]].mean));return p;}
  Request Model(const double* p)const{Request r=original;
    for(unsigned j=0;j<kind.size();++j){if(kind[j]==2)r.backgroundSlope=p[j];else if(kind[j])r.components[species[j]].sigma=std::exp(p[j]);else r.components[species[j]].mean=p[j];}
    for(unsigned j=0;j<r.components.size();++j){const Component& old=original.components[j];for(unsigned q=0;q<r.components[j].nodes.size();++q){r.components[j].nodes[q].mean=old.nodes[q].mean+r.components[j].mean-old.mean;r.components[j].nodes[q].sigma=old.nodes[q].sigma*r.components[j].sigma/old.sigma;}}
    return r;
  }
  bool Feasible(const Request& r)const{
    if(!Finite(r.backgroundSlope)||std::fabs(r.backgroundSlope)>Number(cfg,"exponentialSlopeLimit"))return false;
    for(unsigned j=0;j<r.components.size();++j)if(r.components[j].active){const Component& c=r.components[j];
      if(!Finite(c.mean)||!Finite(c.sigma)||c.mean<r.nativeAxisMin-Number(cfg,"meanAxisPadding")||c.mean>r.nativeAxisMax+Number(cfg,"meanAxisPadding")||c.sigma<Number(cfg,"minimumSigma")||c.sigma>Number(cfg,"maximumSigma"))return false;
      for(unsigned q=0;q<c.nodes.size();++q)if(!Finite(c.nodes[q].mean)||!Finite(c.nodes[q].sigma)||c.nodes[q].sigma<Number(cfg,"minimumSigma")||c.nodes[q].sigma>Number(cfg,"maximumSigma")||c.nodes[q].mean<r.nativeAxisMin-Number(cfg,"meanAxisPadding")||c.nodes[q].mean>r.nativeAxisMax+Number(cfg,"meanAxisPadding")||c.nodes[q].weight<0||!Finite(c.nodes[q].weight))return false;}
    return true;
  }
  double operator()(const double* p)const{Request r=Model(p);if(!Feasible(r))return invalid;Matrix t=BuildTemplates(r,cfg);ProfileResult fit=Profile(r.bins,t,initial,cfg);return fit.converged?fit.dataDeviance:invalid;}
};
inline Result Fit(const Request& request,const Settings& cfg){
  ValidateSettings(cfg);Require(!request.bins.empty(),"No fit bins supplied");Require(request.backgroundModel=="linear"||request.backgroundModel=="exponential","Unknown explicit background model");
  Require(Finite(request.backgroundWindowYield)&&request.backgroundWindowYield>=0&&request.nativeAxisMax>request.nativeAxisMin,"Bad request normalization/axis");
  for(unsigned b=0;b<request.bins.size();++b){const Bin& x=request.bins[b];Require(Finite(x.lo)&&Finite(x.hi)&&Finite(x.count)&&x.hi>x.lo&&x.count>=0,"Invalid finite nonnegative raw bin");if(b)Require(x.lo==request.bins[b-1].hi,"Fit bins are not adjacent whole native bins");}
  for(unsigned j=0;j<request.components.size();++j){const Component& c=request.components[j];Require(!c.name.empty()&&Finite(c.startWindowYield)&&c.startWindowYield>=0,"Invalid component seed");double weight=0;for(unsigned q=0;q<c.nodes.size();++q)weight+=c.nodes[q].weight;
    if(!c.nodes.empty())Require(std::fabs(weight-1)<1.e-9,"R kernel weights must sum to one");}
  Result result;if(!(TotalCount(request.bins)>0)){result.status=-2;return result;}
  Objective objective(request,cfg);Require(objective.Feasible(request),"Seed is outside explicit mathematical support; no silent clipping");
  Vec point=objective.Point(),chosen=point;double start=objective(point.empty()?0:&point[0]);Require(start<invalid,"Seed AMP/template profile failed");
  ROOT::Math::Minimizer* minimizer=0;
  if(!point.empty()){
    Require(gSystem->Load("libMinuit2")>=0,"Cannot load ROOT Minuit2");minimizer=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");Require(minimizer!=0,"Cannot create Minuit2");
    minimizer->SetMaxFunctionCalls(Integer(cfg,"shapeMaxCalls"));minimizer->SetTolerance(Number(cfg,"shapeTolerance"));minimizer->SetPrintLevel(-1);minimizer->SetErrorDef(1.);
    ROOT::Math::Functor f(objective,point.size());minimizer->SetFunction(f);
    for(unsigned j=0;j<point.size();++j){int k=objective.kind[j];double low=k==0?request.nativeAxisMin-Number(cfg,"meanAxisPadding"):(k==1?std::log(Number(cfg,"minimumSigma")):-Number(cfg,"exponentialSlopeLimit"));
      double high=k==0?request.nativeAxisMax+Number(cfg,"meanAxisPadding"):(k==1?std::log(Number(cfg,"maximumSigma")):Number(cfg,"exponentialSlopeLimit"));
      double step=Number(cfg,k==0?"meanStep":(k==1?"sigmaLogStep":"backgroundSlopeStep"));minimizer->SetLimitedVariable(j,objective.names[j],point[j],step,low,high);}
    minimizer->Minimize();Vec candidate(minimizer->X(),minimizer->X()+point.size());double next=objective(&candidate[0]);
    result.pointAccepted=next<invalid&&next<=start;chosen=result.pointAccepted?candidate:point;result.status=minimizer->Status();result.edm=minimizer->Edm();
    result.covarianceStatus=minimizer->CovMatrixStatus();
    if(result.pointAccepted){result.shapeCovariance.assign(point.size(),Vec(point.size()));result.shapeCovarianceValid=result.covarianceStatus==3;
      for(unsigned j=0;j<point.size();++j)for(unsigned k=0;k<point.size();++k){double v=minimizer->CovMatrix(j,k);result.shapeCovariance[j][k]=v;result.shapeCovarianceValid=result.shapeCovarianceValid&&Finite(v);}}
  }else{result.status=0;result.edm=0;result.pointAccepted=true;result.covarianceStatus=0;}
  Request fitted=objective.Model(chosen.empty()?0:&chosen[0]);result.components=fitted.components;result.shapePoint=chosen;result.shapeNames=objective.names;result.backgroundSlope=fitted.backgroundSlope;
  Matrix t=BuildTemplates(fitted,cfg,&result.windowAcceptance);ProfileResult profile=Profile(fitted.bins,t,objective.initial,cfg);
  result.dataDeviance=profile.dataDeviance;result.amplitudeKkt=profile.kkt;result.amplitudeConverged=profile.converged;result.amplitudeIterations=profile.iterations;result.totalTemplates=t;
  const double count=TotalCount(fitted.bins);const unsigned n=fitted.components.size();result.meanError.assign(n,-1);result.sigmaError.assign(n,-1);
  for(unsigned j=0;j<n;++j)result.columnNames.push_back(fitted.components[j].name);
  if(fitted.backgroundModel=="linear"){result.columnNames.push_back("BG_left");result.columnNames.push_back("BG_right");}else result.columnNames.push_back("BG_exp");
  for(unsigned j=0;j<profile.fractions.size();++j)result.columnWindowYields.push_back(count*profile.fractions[j]);
  if(result.shapeCovarianceValid)for(unsigned j=0;j<chosen.size();++j){int sp=objective.species[j];double variance=result.shapeCovariance[j][j];if(variance<0){result.shapeCovarianceValid=false;continue;}
    if(sp>=0){if(objective.kind[j]==0)result.meanError[sp]=std::sqrt(variance);else result.sigmaError[sp]=result.components[sp].sigma*std::sqrt(variance);}}
  result.componentCounts.assign(fitted.bins.size(),Vec(n,0));result.totalCounts.assign(fitted.bins.size(),0);result.smoothCounts.assign(fitted.bins.size(),0);result.effectiveCount.assign(n,0);
  for(unsigned b=0;b<t.size();++b){for(unsigned j=0;j<profile.fractions.size();++j){double v=result.columnWindowYields[j]*t[b][j];result.totalCounts[b]+=v;if(j<n)result.componentCounts[b][j]=v;else result.smoothCounts[b]+=v;}
    double mu=result.totalCounts[b];result.finite=Finite(mu)&&mu>=0;if(!result.finite)break;if(mu>0)result.pearson+=std::pow(fitted.bins[b].count-mu,2)/mu;
    for(unsigned j=0;j<n;++j)result.effectiveCount[j]+=result.componentCounts[b][j]*result.componentCounts[b][j]/std::max(1.e-100,mu);}
  result.finite=result.finite&&Finite(result.dataDeviance)&&result.dataDeviance<invalid;
  // Yield covariance is conditional on final shapes, for interior active yields only.
  std::vector<unsigned> active;for(unsigned j=0;j<profile.fractions.size();++j)if(profile.fractions[j]>1.e-12)active.push_back(j);
  result.yieldCovarianceColumnActive.assign(profile.fractions.size(),false);for(unsigned j=0;j<active.size();++j)result.yieldCovarianceColumnActive[active[j]]=true;
  Matrix h(active.size(),Vec(active.size()));for(unsigned j=0;j<active.size();++j)for(unsigned k=0;k<active.size();++k)h[j][k]=profile.hessian[active[j]][active[k]];
  result.fixedShapeYieldCovariance.assign(profile.fractions.size(),Vec(profile.fractions.size(),0));result.yieldCovarianceValid=!active.empty();
  for(unsigned k=0;k<active.size();++k){Vec unit(active.size(),0),inverse;unit[k]=1;if(!Solve(h,unit,inverse)){result.yieldCovarianceValid=false;break;}
    for(unsigned j=0;j<active.size();++j)result.fixedShapeYieldCovariance[active[j]][active[k]]=count*inverse[j];}
  if(result.yieldCovarianceValid)for(unsigned j=0;j<active.size();++j)result.yieldCovarianceValid=result.yieldCovarianceValid&&Finite(result.fixedShapeYieldCovariance[active[j]][active[j]])&&result.fixedShapeYieldCovariance[active[j]][active[j]]>0;
  if(!result.yieldCovarianceValid)result.fixedShapeYieldCovariance.clear();
  for(unsigned j=0;j<chosen.size();++j){int k=objective.kind[j];double low=k==0?request.nativeAxisMin-Number(cfg,"meanAxisPadding"):(k==1?std::log(Number(cfg,"minimumSigma")):-Number(cfg,"exponentialSlopeLimit"));
    double high=k==0?request.nativeAxisMax+Number(cfg,"meanAxisPadding"):(k==1?std::log(Number(cfg,"maximumSigma")):Number(cfg,"exponentialSlopeLimit"));
    result.atBoundary=result.atBoundary||std::min(chosen[j]-low,high-chosen[j])<=Number(cfg,"boundaryTolerance")*std::max(1.,high-low);}
  result.conditionalDropDeltaD.assign(n,-1);for(unsigned j=0;j<n;++j)if(fitted.components[j].active){Request dropped=fitted;dropped.components[j].active=false;Matrix dt=BuildTemplates(dropped,cfg);Vec a=profile.fractions;a[j]=0;ProfileResult d=Profile(dropped.bins,dt,a,cfg);if(d.converged)result.conditionalDropDeltaD[j]=std::max(0.,d.dataDeviance-result.dataDeviance);}
  result.conditionalNdf=int(fitted.bins.size())-int(chosen.size())-int(active.size());delete minimizer;return result;
}
} // namespace NuclearFit
#endif
