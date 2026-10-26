#ifndef STAR_RESTART_ALL_FIXED_RESPONSE_FIT_H
#define STAR_RESTART_ALL_FIXED_RESPONSE_FIT_H
// Study only. Errors are conditional on fixed curves and the R-bin mixture.
#include "SplitNormalFit.h"
#include "Float16BinMapping.h"
namespace RestartAllFixed {
using NuclearFit::Bin; using NuclearFit::Component; using NuclearFit::Matrix;
using NuclearFit::Request; using NuclearFit::Settings; using NuclearFit::Vec;
struct Shape { Component component; Vec ratios; bool supportProvided;
  Shape():supportProvided(false){} };
struct Interval { double lo,hi;bool includeLo,includeHi;
  Interval(double a=0.,double b=0.,bool l=false,bool h=false):lo(a),hi(b),includeLo(l),includeHi(h){} };
typedef std::vector<Interval> Intervals;
struct Gate {
  std::vector<Intervals> nominal,effective; Vec fractions,gateYields,gradient;
  Matrix countCovariance; double signal,background,total,purity,error;
  bool completeSupport,insideFit,covarianceValid;
  Gate():signal(0),background(0),total(0),purity(-1),error(-1),completeSupport(false),insideFit(false),covarianceValid(false){}
};
struct Result {
  NuclearFit::Result fit;Gate winner,nsigma;double informationCondition,targetRelativeError;
  bool purityAvailable;std::vector<std::string> reasons;
  Result():informationCondition(NuclearFit::invalid),targetRelativeError(NuclearFit::invalid),purityAvailable(false){}
};
inline void ValidateSettings(const Settings& cfg) {
  NuclearFit::ValidateSettings(cfg);
  const char* keys[]={"gateNSigma","maximumPurityPearsonPerNdf","maximumYieldInformationCondition",
    "minimumPurityTargetEffectiveCount","minimumPurityGateModelCount","minimumPurityTargetDropDeltaD","maximumPurityTargetRelativeError"};
  for(unsigned k=0;k<sizeof(keys)/sizeof(keys[0]);++k)NuclearFit::Require(NuclearFit::Number(cfg,keys[k])>0.,"Invalid explicit purity-quality setting");
  NuclearFit::Require(NuclearFit::Text(cfg,"gatePolicy")=="winner","Require production winner gate");
  const std::string m=NuclearFit::Text(cfg,"storageMapping");
  NuclearFit::Require(m=="continuous"||m=="root_default12","Unknown storage mapping");
}
inline double Kernel(const Shape& s,double lo,double hi) {
  if(!s.component.active)return 0.;
  NuclearFit::Require(s.supportProvided&&!s.component.nodes.empty()&&s.component.nodes.size()==s.ratios.size(),"Missing fixed response nodes");
  double v=0.;for(unsigned q=0;q<s.ratios.size();++q){const NuclearFit::Node& n=s.component.nodes[q];
    v+=n.weight*RestartSplitNormal::Integral(lo,hi,n.mean,n.sigma,s.ratios[q]);}return v;
}
inline Matrix Templates(const Request& r,const std::vector<Shape>& shapes,const Settings& cfg,Vec* acceptance=0) {
  NuclearFit::Require(r.components.size()==4&&shapes.size()==4,"Need four ordered species");
  if(r.bins.empty())return Matrix();const double lo=r.bins.front().lo,hi=r.bins.back().hi;
  const unsigned bg=r.backgroundModel=="linear"?2:1;Matrix t(r.bins.size(),Vec(4+bg,0.));Vec accepted(4,0.);
  for(unsigned j=0;j<4;++j)if(r.components[j].active){Shape s=shapes[j];s.component=r.components[j];accepted[j]=Kernel(s,lo,hi);
    if(!(accepted[j]>NuclearFit::Number(cfg,"minimumGaussianAcceptance")))return Matrix();
    for(unsigned b=0;b<t.size();++b)t[b][j]=Kernel(s,r.bins[b].lo,r.bins[b].hi)/accepted[j];}
  for(unsigned b=0;b<t.size();++b)for(unsigned j=0;j<bg;++j)t[b][4+j]=NuclearFit::BackgroundIntegral(r.bins[b].lo,r.bins[b].hi,lo,hi,r.backgroundModel,r.backgroundSlope,j);
  for(unsigned j=0;j<4+bg;++j){double sum=0.;for(unsigned b=0;b<t.size();++b){if(!NuclearFit::Finite(t[b][j])||t[b][j]<0.)return Matrix();sum+=t[b][j];}
    if(sum>0.){if(j<4)accepted[j]*=sum;for(unsigned b=0;b<t.size();++b)t[b][j]/=sum;}}
  if(acceptance)*acceptance=accepted;return t;
}
inline bool CommonSupport(const std::vector<Shape>& s,unsigned target) {
  if(s.size()!=4||target>=4||!s[target].supportProvided)return false;
  const unsigned nq=s[target].component.nodes.size();if(!nq)return false;
  for(unsigned j=0;j<4;++j){if(!s[j].supportProvided||s[j].component.nodes.size()!=nq||s[j].ratios.size()!=nq)return false;
    for(unsigned q=0;q<nq;++q)if(std::fabs(s[j].component.nodes[q].weight-s[target].component.nodes[q].weight)>1.e-12)return false;}return true;
}
inline bool Selected(double x,const Vec& mean,const Vec& sigma,unsigned target,double gate,bool winner) {
  if(!(std::fabs((x-mean[target])/sigma[target])<gate))return false;if(!winner)return true;
  unsigned chosen=0;double best=std::fabs((x-mean[0])/sigma[0]);
  for(unsigned j=1;j<mean.size();++j){double p=std::fabs((x-mean[j])/sigma[j]);if(p<best){best=p;chosen=j;}}return chosen==target;
}
inline Intervals NodeGate(const Vec& mean,const Vec& sigma,unsigned target,double gate,bool winner) {
  NuclearFit::Require(mean.size()==sigma.size()&&target<mean.size(),"Invalid gate curves");
  for(unsigned j=0;j<sigma.size();++j)NuclearFit::Require(sigma[j]>0.&&NuclearFit::Finite(sigma[j])&&NuclearFit::Finite(mean[j]),"Invalid gate mean/width");
  double lo=std::max(0.,mean[target]-gate*sigma[target]),hi=mean[target]+gate*sigma[target];if(!(hi>lo))return Intervals();
  Vec roots;roots.push_back(lo);roots.push_back(hi);
  if(winner)for(unsigned j=0;j<mean.size();++j)if(j!=target){double den=sigma[j]-sigma[target];if(den!=0.){
    double x=(mean[target]*sigma[j]-mean[j]*sigma[target])/den;if(x>lo&&x<hi)roots.push_back(x);}
    double x=(mean[target]*sigma[j]+mean[j]*sigma[target])/(sigma[j]+sigma[target]);if(x>lo&&x<hi)roots.push_back(x);}
  std::sort(roots.begin(),roots.end());roots.erase(std::unique(roots.begin(),roots.end()),roots.end());
  std::vector<bool> point(roots.size(),false),piece(roots.size()-1,false);
  for(unsigned i=0;i<roots.size();++i)point[i]=Selected(roots[i],mean,sigma,target,gate,winner);
  for(unsigned i=0;i<piece.size();++i)piece[i]=Selected(.5*(roots[i]+roots[i+1]),mean,sigma,target,gate,winner);
  Intervals result;for(unsigned i=0;i<piece.size();++i)if(piece[i])result.push_back(Interval(roots[i],roots[i+1],point[i]&&!(i&&piece[i-1]),point[i+1]));
  for(unsigned i=0;i<roots.size();++i)if(point[i]&&!(i&&piece[i-1])&&!(i<piece.size()&&piece[i]))result.push_back(Interval(roots[i],roots[i],true,true));return result;
}
inline Interval EffectiveInterval(const Interval& x,const std::string& mapping) {
  if(mapping=="continuous")return x;
  const double l=x.includeLo?x.lo:std::nextafter(x.lo,std::numeric_limits<double>::infinity());
  const double h=x.includeHi?std::nextafter(x.hi,std::numeric_limits<double>::infinity()):x.hi;
  return Interval(RestartFloat16::ContinuumEdge(l),RestartFloat16::ContinuumEdge(h),true,false);
}
inline double InfinityNorm(const Matrix& a) {
  double v=0.;for(unsigned j=0;j<a.size();++j){double row=0.;for(unsigned k=0;k<a[j].size();++k)row+=std::fabs(a[j][k]);v=std::max(v,row);}return v;
}
inline void Covariance(const NuclearFit::ProfileResult& p,const Matrix& t,double count,NuclearFit::Result& fit,double& condition) {
  const unsigned n=p.fractions.size();std::vector<unsigned> used;fit.yieldCovarianceColumnActive.assign(n,false);
  // Keep nonzero columns even at zero fitted yield: do not manufacture P=1+-0.
  for(unsigned j=0;j<n;++j){double sum=0.;for(unsigned b=0;b<t.size();++b)sum+=t[b][j];if(sum>0.){used.push_back(j);fit.yieldCovarianceColumnActive[j]=true;}}
  fit.yieldCovarianceValid=false;fit.fixedShapeYieldCovariance.clear();if(!p.converged||used.empty())return;
  const unsigned m=used.size();Matrix h(m,Vec(m)),scaled(m,Vec(m)),inverse(m,Vec(m)),invScaled(m,Vec(m));
  for(unsigned j=0;j<m;++j)if(!(p.hessian[used[j]][used[j]]>0.))return;
  for(unsigned j=0;j<m;++j)for(unsigned k=0;k<m;++k){h[j][k]=p.hessian[used[j]][used[k]];scaled[j][k]=h[j][k]/std::sqrt(p.hessian[used[j]][used[j]]*p.hessian[used[k]][used[k]]);}
  for(unsigned k=0;k<m;++k){Vec rhs(m,0.),col;rhs[k]=1.;if(!NuclearFit::Solve(h,rhs,col))return;for(unsigned j=0;j<m;++j)inverse[j][k]=count*col[j];
    if(!NuclearFit::Solve(scaled,rhs,col))return;for(unsigned j=0;j<m;++j)invScaled[j][k]=col[j];}
  if(!RestartResponse::PositiveDefinite(inverse))return;condition=InfinityNorm(scaled)*InfinityNorm(invScaled);
  fit.fixedShapeYieldCovariance.assign(n,Vec(n,0.));for(unsigned j=0;j<m;++j)for(unsigned k=0;k<m;++k)fit.fixedShapeYieldCovariance[used[j]][used[k]]=inverse[j][k];
  fit.yieldCovarianceValid=NuclearFit::Finite(condition);
}
inline Gate EvaluateGate(const Request& r,const std::vector<Shape>& shapes,const NuclearFit::Result& fit,unsigned target,const Settings& cfg,bool winner) {
  Gate g;g.completeSupport=CommonSupport(shapes,target);if((winner&&!g.completeSupport)||!shapes[target].supportProvided)return g;
  const unsigned nq=shapes[target].component.nodes.size(),n=fit.columnWindowYields.size();g.fractions.assign(n,0.);g.insideFit=true;
  const double fitLo=r.bins.front().lo,fitHi=r.bins.back().hi;
  for(unsigned q=0;q<nq;++q){Vec mean(4,0.),sigma(4,1.);for(unsigned j=0;j<4;++j)if(shapes[j].supportProvided&&shapes[j].component.nodes.size()==nq){mean[j]=shapes[j].component.nodes[q].mean;sigma[j]=shapes[j].component.nodes[q].sigma;}
    const Intervals nominal=NodeGate(mean,sigma,target,NuclearFit::Number(cfg,"gateNSigma"),winner);Intervals effective;
    for(unsigned k=0;k<nominal.size();++k){const Interval x=EffectiveInterval(nominal[k],NuclearFit::Text(cfg,"storageMapping"));effective.push_back(x);
      if(x.lo<fitLo||x.hi>fitHi)g.insideFit=false;if(!(x.hi>x.lo))continue;
      for(unsigned j=0;j<4;++j)if(r.components[j].active){const NuclearFit::Node& node=shapes[j].component.nodes[q];g.fractions[j]+=node.weight*RestartSplitNormal::Integral(x.lo,x.hi,node.mean,node.sigma,shapes[j].ratios[q])/fit.windowAcceptance[j];}
      for(unsigned j=4;j<n;++j)g.fractions[j]+=shapes[target].component.nodes[q].weight*NuclearFit::BackgroundIntegral(x.lo,x.hi,fitLo,fitHi,r.backgroundModel,r.backgroundSlope,j-4);}
    g.nominal.push_back(nominal);g.effective.push_back(effective);}
  if(!g.insideFit)return g;g.gateYields.resize(n);
  for(unsigned j=0;j<n;++j){g.gateYields[j]=fit.columnWindowYields[j]*g.fractions[j];if(j==target)g.signal+=g.gateYields[j];else g.background+=g.gateYields[j];}
  g.total=g.signal+g.background;if(!(g.total>0.))return g;g.purity=g.signal/g.total;g.gradient.resize(n);
  for(unsigned j=0;j<n;++j)g.gradient[j]=g.fractions[j]*((j==target?1.:0.)-g.purity)/g.total;
  if(!fit.yieldCovarianceValid)return g;double variance=0.;g.countCovariance.assign(2,Vec(2,0.));
  for(unsigned j=0;j<n;++j)for(unsigned k=0;k<n;++k){const double c=fit.fixedShapeYieldCovariance[j][k];variance+=g.gradient[j]*c*g.gradient[k];g.countCovariance[j==target?0:1][k==target?0:1]+=g.fractions[j]*c*g.fractions[k];}
  if(variance>=-1.e-14&&NuclearFit::Finite(variance)){g.error=std::sqrt(std::max(0.,variance));g.covarianceValid=true;}return g;
}
inline Result Fit(const Request& r,const std::vector<Shape>& shapes,unsigned target,const Settings& cfg) {
  ValidateSettings(cfg);NuclearFit::Require(target<4&&shapes.size()==4&&r.components.size()==4&&!r.bins.empty(),"Invalid fit dimensions");
  NuclearFit::Require(shapes[target].supportProvided&&r.components[target].active,"Target must be supported and active");
  NuclearFit::Require(!r.freeBackgroundSlope&&(r.backgroundModel=="linear"||r.backgroundModel=="exponential"),"Background shape must be fixed");
  const unsigned nq=shapes[target].component.nodes.size();
  for(unsigned j=0;j<4;++j){NuclearFit::Require(!r.components[j].freeMean&&!r.components[j].freeSigma,"Nuclear centers/widths must not float");
    if(r.components[j].active){NuclearFit::Require(shapes[j].supportProvided&&shapes[j].component.nodes.size()==nq&&shapes[j].ratios.size()==nq,"Active response node dimensions differ");
      for(unsigned q=0;q<nq;++q)NuclearFit::Require(std::fabs(shapes[j].component.nodes[q].weight-shapes[target].component.nodes[q].weight)<1.e-12,"Active quadrature weights differ");
      NuclearFit::Require(r.components[j].mean>=r.bins.front().lo&&r.components[j].mean<=r.bins.back().hi,"Active nuisance peak center must be inside local fit");}}
  Result result;NuclearFit::Result& f=result.fit;f.components=r.components;f.backgroundSlope=r.backgroundSlope;f.meanError.assign(4,-1.);f.sigmaError.assign(4,-1.);
  const double count=NuclearFit::TotalCount(r.bins);if(!(count>0.)){f.status=-2;result.reasons.push_back("empty_roi");return result;}
  Matrix t=Templates(r,shapes,cfg,&f.windowAcceptance);NuclearFit::Require(!t.empty(),"Fixed templates failed");
  Vec initial;for(unsigned j=0;j<4;++j)initial.push_back(r.components[j].active?r.components[j].startWindowYield/count:0.);
  const unsigned bg=r.backgroundModel=="linear"?2:1;for(unsigned j=0;j<bg;++j)initial.push_back(r.backgroundWindowYield/count/bg);
  const NuclearFit::ProfileResult p=NuclearFit::Profile(r.bins,t,initial,cfg);
  f.status=p.converged?0:-3;f.pointAccepted=p.converged;f.amplitudeConverged=p.converged;f.amplitudeIterations=p.iterations;f.amplitudeKkt=p.kkt;f.dataDeviance=p.dataDeviance;
  f.totalTemplates=t;f.totalCounts.assign(t.size(),0.);f.smoothCounts.assign(t.size(),0.);f.componentCounts.assign(t.size(),Vec(4,0.));f.effectiveCount.assign(4,0.);
  for(unsigned j=0;j<p.fractions.size();++j){f.columnWindowYields.push_back(count*p.fractions[j]);f.columnNames.push_back(j<4?r.components[j].name:(bg==2?(j==4?"BG_left":"BG_right"):"BG_exp"));
    if(p.fractions[j]<=NuclearFit::Number(cfg,"boundaryTolerance")&&(j>=4||r.components[j].active))f.atBoundary=true;}
  f.finite=p.converged&&NuclearFit::Finite(f.dataDeviance);
  for(unsigned b=0;b<t.size();++b){for(unsigned j=0;j<p.fractions.size();++j){double v=f.columnWindowYields[j]*t[b][j];f.totalCounts[b]+=v;if(j<4)f.componentCounts[b][j]=v;else f.smoothCounts[b]+=v;}
    double mu=f.totalCounts[b],obs=r.bins[b].count;if(!(mu>0.)){if(obs>0.)f.finite=false;continue;}f.pearson+=(obs-mu)*(obs-mu)/mu;for(unsigned j=0;j<4;++j)f.effectiveCount[j]+=f.componentCounts[b][j]*f.componentCounts[b][j]/mu;}
  Covariance(p,t,count,f,result.informationCondition);f.covarianceStatus=f.yieldCovarianceValid?3:0;
  unsigned active=0;for(unsigned j=0;j<p.fractions.size();++j)if(p.fractions[j]>1.e-12)++active;f.conditionalNdf=int(r.bins.size())-int(active);
  f.conditionalDropDeltaD.assign(4,-1.);for(unsigned j=0;j<4;++j)if(r.components[j].active){Request reduced=r;reduced.components[j].active=false;Matrix rt=Templates(reduced,shapes,cfg);Vec start=p.fractions;start[j]=0.;
    NuclearFit::ProfileResult drop=NuclearFit::Profile(r.bins,rt,start,cfg);if(drop.converged)f.conditionalDropDeltaD[j]=std::max(0.,drop.dataDeviance-f.dataDeviance);}
  result.nsigma=EvaluateGate(r,shapes,f,target,cfg,false);result.winner=EvaluateGate(r,shapes,f,target,cfg,true);
  if(f.yieldCovarianceValid&&f.columnWindowYields[target]>0.)result.targetRelativeError=std::sqrt(std::max(0.,f.fixedShapeYieldCovariance[target][target]))/f.columnWindowYields[target];
  if(!f.finite)result.reasons.push_back("nonconverged_or_nonfinite_fit");
  if(!result.winner.completeSupport)result.reasons.push_back("missing_common_curve_for_winner_gate");
  if(!result.winner.insideFit)result.reasons.push_back("winner_gate_outside_fit_roi");
  if(!f.yieldCovarianceValid||!result.winner.covarianceValid)result.reasons.push_back("invalid_conditional_covariance");
  if(result.informationCondition>NuclearFit::Number(cfg,"maximumYieldInformationCondition"))result.reasons.push_back("ill_conditioned_yield_information");
  if(f.conditionalNdf<=0||f.pearson/f.conditionalNdf>NuclearFit::Number(cfg,"maximumPurityPearsonPerNdf"))result.reasons.push_back("poor_local_fit_quality");
  if(f.effectiveCount[target]<NuclearFit::Number(cfg,"minimumPurityTargetEffectiveCount"))result.reasons.push_back("insufficient_target_effective_count");
  if(result.winner.total<NuclearFit::Number(cfg,"minimumPurityGateModelCount"))result.reasons.push_back("insufficient_fitted_gate_count");
  if(f.conditionalDropDeltaD[target]<NuclearFit::Number(cfg,"minimumPurityTargetDropDeltaD"))result.reasons.push_back("target_not_identifiable_against_null");
  if(result.targetRelativeError>NuclearFit::Number(cfg,"maximumPurityTargetRelativeError"))result.reasons.push_back("large_conditional_target_yield_error");
  if(!(f.columnWindowYields[target]>0.))result.reasons.push_back("target_yield_at_boundary");
  bool bgZero=true;for(unsigned j=0;j<f.columnWindowYields.size();++j)if(j!=target&&f.columnWindowYields[j]>count*NuclearFit::Number(cfg,"boundaryTolerance"))bgZero=false;
  if(bgZero)result.reasons.push_back("all_background_at_nonnegative_boundary");result.purityAvailable=result.reasons.empty();return result;
}
}
#endif
