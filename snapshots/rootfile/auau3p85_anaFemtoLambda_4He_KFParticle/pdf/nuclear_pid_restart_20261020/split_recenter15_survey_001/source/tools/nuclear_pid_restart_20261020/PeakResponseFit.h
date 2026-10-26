#ifndef STAR_RESTART_PEAK_RESPONSE_FIT_H
#define STAR_RESTART_PEAK_RESPONSE_FIT_H
// Study-only response shapes. No PID coefficients, interpolation or cutoffs.
#include "../nuclear_pid_tof/NuclearFitCommon.h"
#include <limits>
#include <memory>
#include <set>

namespace RestartResponse {
using NuclearFit::Component;
using NuclearFit::Matrix;
using NuclearFit::Request;
using NuclearFit::Settings;
using NuclearFit::Vec;

inline bool IsTail(const std::string& model) {
  NuclearFit::Require(model=="gaussian" || model=="gaussian_left_exp" ||
                      model=="gaussian_right_exp","Unknown explicit response model: "+model);
  return model!="gaussian";
}
inline void CheckAlpha(double alpha) {
  NuclearFit::Require(NuclearFit::Finite(alpha) && alpha>0.,"Response alpha must be finite and positive");
}
inline double Normalization(double alpha) {
  CheckAlpha(alpha);
  const double value=std::sqrt(TMath::Pi()/2.)*std::erfc(-alpha/std::sqrt(2.))+
                     std::exp(-.5*alpha*alpha)/alpha;
  NuclearFit::Require(NuclearFit::Finite(value) && value>0.,"Response normalization is not finite/positive");
  return value;
}
inline double StandardDensity(double u,const std::string& model,double alpha) {
  const bool tail=IsTail(model);
  NuclearFit::Require(u==u,"NaN response coordinate");
  if(!tail)return std::exp(-.5*u*u)/std::sqrt(2.*TMath::Pi());
  CheckAlpha(alpha);
  if(model=="gaussian_right_exp")u=-u;
  const double value=u < -alpha ? std::exp(-.5*alpha*alpha+alpha*(u+alpha)) : std::exp(-.5*u*u);
  return value/Normalization(alpha);
}
inline double StandardDensityDerivative(double u,const std::string& model,double alpha) {
  const bool tail=IsTail(model);
  NuclearFit::Require(NuclearFit::Finite(u),"Derivative requires a finite response coordinate");
  if(!tail)return -u*StandardDensity(u,model,alpha);
  CheckAlpha(alpha);
  const double direction=model=="gaussian_right_exp" ? -1. : 1.;
  const double v=direction*u;
  return direction*(v < -alpha ? alpha : -v)*StandardDensity(u,model,alpha);
}
inline double LeftIntegral(double lo,double hi,double alpha) {
  CheckAlpha(alpha);
  if(!(hi>lo))return 0.;
  const double join=-alpha,tailScale=std::exp(-.5*alpha*alpha)/alpha;
  double value=0.;
  if(lo<join) {
    const double upper=std::min(hi,join);
    // Difference of exponentials evaluated at the larger endpoint.  This also
    // handles lo=-infinity without subtracting nearly equal CDF values.
    value+=tailScale*std::exp(alpha*(upper+alpha))*(-std::expm1(alpha*(lo-upper)));
  }
  if(hi>join) {
    const double lower=std::max(lo,join);
    value+=std::sqrt(2.*TMath::Pi())*NuclearFit::GaussianIntegral(lower,hi,0.,1.);
  }
  return value/Normalization(alpha);
}
inline double StandardIntegral(double lo,double hi,const std::string& model,double alpha) {
  const bool tail=IsTail(model);
  NuclearFit::Require(lo==lo && hi==hi,"NaN response integration endpoint");
  if(!(hi>lo))return 0.;
  if(!tail)return NuclearFit::GaussianIntegral(lo,hi,0.,1.);
  return model=="gaussian_left_exp" ? LeftIntegral(lo,hi,alpha) : LeftIntegral(-hi,-lo,alpha);
}
inline double Integral(double lo,double hi,double mean,double sigma,
                       const std::string& model,double alpha) {
  NuclearFit::Require(NuclearFit::Finite(mean) && NuclearFit::Finite(sigma) && sigma>0.,
                      "Response mean/core width must be finite with positive width");
  return StandardIntegral((lo-mean)/sigma,(hi-mean)/sigma,model,alpha);
}
inline double Cdf(double x,double mean,double sigma,const std::string& model,double alpha) {
  return Integral(-std::numeric_limits<double>::infinity(),x,mean,sigma,model,alpha);
}
inline double KernelIntegral(const Component& component,double lo,double hi,
                             const std::string& model,double alpha) {
  if(!component.active)return 0.;
  if(component.nodes.empty())return Integral(lo,hi,component.mean,component.sigma,model,alpha);
  double value=0.;
  for(unsigned i=0;i<component.nodes.size();++i)
    value+=component.nodes[i].weight*Integral(lo,hi,component.nodes[i].mean,
                                           component.nodes[i].sigma,model,alpha);
  return value;
}
inline Matrix BuildTemplates(const Request& request,const Settings& cfg,
                              const std::string& model,double alpha,Vec* acceptance=0) {
  // Exact baseline route: no separately rounded Gaussian implementation.
  if(!IsTail(model))return NuclearFit::BuildTemplates(request,cfg,acceptance);
  CheckAlpha(alpha);
  if(request.bins.empty())return Matrix();
  const double lo=request.bins.front().lo,hi=request.bins.back().hi;
  const unsigned n=request.components.size(),bg=request.backgroundModel=="linear" ? 2 : 1;
  Matrix result(request.bins.size(),Vec(n+bg,0.));Vec accepted(n,0.);
  for(unsigned j=0;j<n;++j)if(request.components[j].active) {
    accepted[j]=KernelIntegral(request.components[j],lo,hi,model,alpha);
    if(!(accepted[j]>NuclearFit::Number(cfg,"minimumGaussianAcceptance")))return Matrix();
    for(unsigned b=0;b<result.size();++b)
      result[b][j]=KernelIntegral(request.components[j],request.bins[b].lo,
                                 request.bins[b].hi,model,alpha)/accepted[j];
  }
  for(unsigned b=0;b<result.size();++b)for(unsigned j=0;j<bg;++j)
    result[b][n+j]=NuclearFit::BackgroundIntegral(request.bins[b].lo,request.bins[b].hi,
                     lo,hi,request.backgroundModel,request.backgroundSlope,j);
  for(unsigned j=0;j<n+bg;++j) {
    double sum=0.;
    for(unsigned b=0;b<result.size();++b) {
      if(!NuclearFit::Finite(result[b][j]) || result[b][j]<0.)return Matrix();
      sum+=result[b][j];
    }
    if(sum>0.)for(unsigned b=0;b<result.size();++b)result[b][j]/=sum;
  }
  if(acceptance)*acceptance=accepted;
  return result;
}

struct Result {
  NuclearFit::Result fit;
  std::string model;
  double alpha,alphaError;
  int freeParameters;
  // The legacy Gaussian delegate does not expose the separate Minimize/Hesse
  // phases.  Keep availability false there instead of inventing observations.
  bool optimizerDiagnosticsAvailable,minimizeReturned,minimizeConverged;
  int migradStatus,hesseStatus;
  double migradEdm;
  bool hesseAttempted,hesseSucceeded;
  Result():alpha(0.),alphaError(-1.),freeParameters(0),
    optimizerDiagnosticsAvailable(false),minimizeReturned(false),minimizeConverged(false),
    migradStatus(-1),hesseStatus(-1),migradEdm(NuclearFit::invalid),
    hesseAttempted(false),hesseSucceeded(false){}
};

inline bool MigradConverged(bool returned,int status,double edm) {
  return returned && status==0 && NuclearFit::Finite(edm) && edm>=0.;
}
inline int PhaseStatus(const Result& response) {
  if(!response.optimizerDiagnosticsAvailable)return response.fit.status;
  // A successful Hesse can never erase a failed/nonconverged MIGRAD phase.
  if(response.migradStatus!=0)return response.migradStatus;
  if(!response.minimizeConverged)return -10; // false return or invalid MIGRAD EDM despite status 0
  if(!response.fit.pointAccepted)return -12; // no accepted profile point
  if(!response.hesseAttempted || !response.hesseSucceeded)
    return response.hesseStatus>0 ? response.hesseStatus : -11;
  return response.hesseStatus; // 0 only after both independent phases passed
}

inline void ValidateSettings(const Settings& cfg) {
  NuclearFit::ValidateSettings(cfg);
  const double lo=NuclearFit::Number(cfg,"responseAlphaMin");
  const double hi=NuclearFit::Number(cfg,"responseAlphaMax");
  NuclearFit::Require(lo>0. && hi>lo && NuclearFit::Number(cfg,"responseAlphaLogStep")>0.,
                      "Invalid explicit response-alpha bounds/log step");
  Normalization(lo);Normalization(hi);
}
inline unsigned ValidateRequest(const Request& request) {
  NuclearFit::Require(request.components.size()==4,"Response study requires all four component records");
  NuclearFit::Require(request.backgroundModel=="linear" || request.backgroundModel=="exponential",
                      "Unknown explicit response background model");
  NuclearFit::Require(!request.bins.empty(),"No whole native response bins");
  NuclearFit::Require(NuclearFit::Finite(request.backgroundWindowYield) &&
                      request.backgroundWindowYield>=0. &&
                      NuclearFit::Finite(request.nativeAxisMin) && NuclearFit::Finite(request.nativeAxisMax) &&
                      request.nativeAxisMax>request.nativeAxisMin,"Invalid response normalization/native axis");
  for(unsigned b=0;b<request.bins.size();++b) {
    const NuclearFit::Bin& bin=request.bins[b];
    NuclearFit::Require(NuclearFit::Finite(bin.lo) && NuclearFit::Finite(bin.hi) &&
                        NuclearFit::Finite(bin.count) && bin.hi>bin.lo && bin.count>=0.,
                        "Invalid whole native response bin");
    if(b)NuclearFit::Require(bin.lo==request.bins[b-1].hi,"Response bins must be adjacent native bins");
  }
  unsigned active=0,target=0;std::set<std::string> componentNames;
  for(unsigned j=0;j<request.components.size();++j) {
    const Component& c=request.components[j];
    NuclearFit::Require(!c.name.empty() && NuclearFit::Finite(c.mean) && NuclearFit::Finite(c.sigma) &&
                        NuclearFit::Finite(c.startWindowYield) && c.startWindowYield>=0.,
                        "Invalid preserved response component/seed");
    NuclearFit::Require(componentNames.insert(c.name).second,"Duplicate response component name");
    if(c.active) {
      ++active;target=j;
      NuclearFit::Require(c.freeMean && c.freeSigma && c.sigma>0.,
                          "The single active target must have free mean and core width");
    } else NuclearFit::Require(!c.freeMean && !c.freeSigma,
                               "Inactive foreign species cannot have free peak parameters");
    double sum=0.;
    for(unsigned q=0;q<c.nodes.size();++q) {
      const NuclearFit::Node& node=c.nodes[q];
      NuclearFit::Require(NuclearFit::Finite(node.weight) && node.weight>=0. &&
                          NuclearFit::Finite(node.mean) && NuclearFit::Finite(node.sigma),
                          "Invalid native-R kernel node");
      if(c.active)NuclearFit::Require(node.sigma>0.,"Active native-R node has nonpositive core width");
      sum+=node.weight;
    }
    if(!c.nodes.empty())NuclearFit::Require(std::fabs(sum-1.)<1.e-9,"Native-R node weights must sum to one");
  }
  NuclearFit::Require(active==1,"Exactly one target may be active/free; no foreign-peak fitting");
  NuclearFit::Require(componentNames.count("d") && componentNames.count("t") &&
                      componentNames.count("he3") && componentNames.count("he4"),
                      "Response component records must be d/t/he3/he4");
  return target;
}

class Objective {
 public:
  NuclearFit::Objective base;
  std::string model;
  double alphaStart;
  bool tail;
  Objective(const Request& request,const Settings& cfg,const std::string& family,double start):
    base(request,cfg),model(family),alphaStart(start),tail(IsTail(family)){}
  Vec Point()const {
    Vec point=base.Point();if(tail)point.push_back(std::log(alphaStart));return point;
  }
  double Alpha(const double* point)const {return tail ? std::exp(point[base.kind.size()]) : 0.;}
  Request Model(const double* point)const {
    Request request=base.Model(point);
    // Base::Model rescales nodes for every component.  Restore inactive records
    // exactly, including unsupported zero-width placeholders, rather than let
    // an irrelevant 0/0 transformation alter a foreign species.
    for(unsigned j=0;j<request.components.size();++j)
      if(!base.original.components[j].active)request.components[j]=base.original.components[j];
    return request;
  }
  bool Feasible(const Request& request,double alpha)const {
    return base.Feasible(request) && (!tail ||
      (NuclearFit::Finite(alpha) && alpha>=NuclearFit::Number(base.cfg,"responseAlphaMin") &&
       alpha<=NuclearFit::Number(base.cfg,"responseAlphaMax")));
  }
  double operator()(const double* point)const {
    const Request request=Model(point);const double alpha=Alpha(point);
    if(!Feasible(request,alpha))return NuclearFit::invalid;
    const Matrix templates=BuildTemplates(request,base.cfg,model,alpha);
    const NuclearFit::ProfileResult profile=NuclearFit::Profile(request.bins,templates,base.initial,base.cfg);
    return profile.converged ? profile.dataDeviance : NuclearFit::invalid;
  }
};

inline bool PositiveDefinite(const Matrix& covariance) {
  if(covariance.empty())return false;
  const unsigned n=covariance.size();Matrix lower(n,Vec(n,0.));
  for(unsigned i=0;i<n;++i) {
    if(covariance[i].size()!=n)return false;
    for(unsigned j=0;j<=i;++j) {
      if(!NuclearFit::Finite(covariance[i][j]) || !NuclearFit::Finite(covariance[j][i]))return false;
      double value=.5*(covariance[i][j]+covariance[j][i]);
      for(unsigned k=0;k<j;++k)value-=lower[i][k]*lower[j][k];
      if(i==j) {if(!(value>0.) || !NuclearFit::Finite(value))return false;lower[i][j]=std::sqrt(value);}
      else {lower[i][j]=value/lower[j][j];if(!NuclearFit::Finite(lower[i][j]))return false;}
    }
  }
  return true;
}
inline void ParameterLimits(const Objective& objective,unsigned index,
                             double& low,double& high,double& step,std::string& name) {
  const Settings& cfg=objective.base.cfg;
  if(index==objective.base.kind.size()) {
    low=std::log(NuclearFit::Number(cfg,"responseAlphaMin"));
    high=std::log(NuclearFit::Number(cfg,"responseAlphaMax"));
    step=NuclearFit::Number(cfg,"responseAlphaLogStep");name="response_logAlpha";return;
  }
  const int kind=objective.base.kind[index];name=objective.base.names[index];
  low=kind==0 ? objective.base.original.nativeAxisMin-NuclearFit::Number(cfg,"meanAxisPadding") :
       (kind==1 ? std::log(NuclearFit::Number(cfg,"minimumSigma")) : -NuclearFit::Number(cfg,"exponentialSlopeLimit"));
  high=kind==0 ? objective.base.original.nativeAxisMax+NuclearFit::Number(cfg,"meanAxisPadding") :
        (kind==1 ? std::log(NuclearFit::Number(cfg,"maximumSigma")) : NuclearFit::Number(cfg,"exponentialSlopeLimit"));
  step=NuclearFit::Number(cfg,kind==0 ? "meanStep" : (kind==1 ? "sigmaLogStep" : "backgroundSlopeStep"));
}

inline Result Fit(const Request& request,const Settings& cfg,
                   const std::string& model,double alphaStart) {
  ValidateSettings(cfg);ValidateRequest(request);
  const bool tail=IsTail(model);
  Result result;result.model=model;
  if(!tail) {
    result.fit=NuclearFit::Fit(request,cfg);result.alpha=0.;result.alphaError=0.;
    for(unsigned j=0;j<result.fit.components.size();++j)
      if(!request.components[j].active)result.fit.components[j]=request.components[j];
    result.freeParameters=result.fit.shapePoint.size();return result;
  }
  CheckAlpha(alphaStart);
  NuclearFit::Require(alphaStart>=NuclearFit::Number(cfg,"responseAlphaMin") &&
                      alphaStart<=NuclearFit::Number(cfg,"responseAlphaMax"),
                      "Explicit response alpha start is outside its bounds");
  if(!(NuclearFit::TotalCount(request.bins)>0.)) {
    result.fit.status=-2;result.alpha=alphaStart;return result;
  }
  Objective objective(request,cfg,model,alphaStart);
  NuclearFit::Require(objective.Feasible(request,alphaStart),"Response seed is outside mathematical support");
  Vec point=objective.Point(),chosen=point;
  const double start=objective(&point[0]);
  NuclearFit::Require(start<NuclearFit::invalid,"Response seed template/amplitude profile failed");
  NuclearFit::Require(gSystem->Load("libMinuit2")>=0,"Cannot load response ROOT Minuit2");
  std::unique_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
  NuclearFit::Require(minimizer.get()!=0,"Cannot create response Minuit2");
  minimizer->SetMaxFunctionCalls(NuclearFit::Integer(cfg,"shapeMaxCalls"));
  minimizer->SetTolerance(NuclearFit::Number(cfg,"shapeTolerance"));
  minimizer->SetPrintLevel(-1);minimizer->SetErrorDef(1.);
  ROOT::Math::Functor functor(objective,point.size());minimizer->SetFunction(functor);
  for(unsigned j=0;j<point.size();++j) {
    double low,high,step;std::string name;ParameterLimits(objective,j,low,high,step,name);
    minimizer->SetLimitedVariable(j,name,point[j],step,low,high);
  }
  result.optimizerDiagnosticsAvailable=true;
  result.minimizeReturned=minimizer->Minimize();
  result.migradStatus=minimizer->Status();result.migradEdm=minimizer->Edm();
  result.minimizeConverged=MigradConverged(result.minimizeReturned,result.migradStatus,result.migradEdm);
  Vec candidate(minimizer->X(),minimizer->X()+point.size());
  const double next=objective(&candidate[0]);
  NuclearFit::Result& fit=result.fit;
  fit.pointAccepted=next<NuclearFit::invalid && next<=start;
  chosen=fit.pointAccepted ? candidate : point;
  if(result.minimizeConverged && fit.pointAccepted) {
    result.hesseAttempted=true;result.hesseSucceeded=minimizer->Hesse();
    result.hesseStatus=minimizer->Status();
  }
  fit.status=PhaseStatus(result);fit.edm=result.migradEdm;fit.covarianceStatus=minimizer->CovMatrixStatus();
  bool covarianceAtChosen=result.minimizeConverged && fit.pointAccepted;
  for(unsigned j=0;j<chosen.size();++j)covarianceAtChosen=covarianceAtChosen && chosen[j]==minimizer->X()[j];
  if(covarianceAtChosen && result.hesseAttempted && result.hesseSucceeded &&
     result.hesseStatus==0 && fit.covarianceStatus==3) {
    fit.shapeCovariance.assign(chosen.size(),Vec(chosen.size(),0.));
    for(unsigned j=0;j<chosen.size();++j)for(unsigned k=0;k<chosen.size();++k)
      fit.shapeCovariance[j][k]=minimizer->CovMatrix(j,k);
    fit.shapeCovarianceValid=PositiveDefinite(fit.shapeCovariance);
  }
  if(!fit.shapeCovarianceValid)fit.shapeCovariance.clear();
  const Request fitted=objective.Model(&chosen[0]);result.alpha=objective.Alpha(&chosen[0]);
  result.freeParameters=chosen.size();fit.components=fitted.components;fit.shapePoint=chosen;
  fit.shapeNames=objective.base.names;fit.shapeNames.push_back("response_logAlpha");
  fit.backgroundSlope=fitted.backgroundSlope;
  const Matrix templates=BuildTemplates(fitted,cfg,model,result.alpha,&fit.windowAcceptance);
  const NuclearFit::ProfileResult profile=NuclearFit::Profile(fitted.bins,templates,objective.base.initial,cfg);
  fit.dataDeviance=profile.dataDeviance;fit.amplitudeKkt=profile.kkt;
  fit.amplitudeConverged=profile.converged;fit.amplitudeIterations=profile.iterations;fit.totalTemplates=templates;
  const double count=NuclearFit::TotalCount(fitted.bins);const unsigned n=fitted.components.size();
  fit.meanError.assign(n,-1.);fit.sigmaError.assign(n,-1.);
  for(unsigned j=0;j<n;++j)fit.columnNames.push_back(fitted.components[j].name);
  if(fitted.backgroundModel=="linear"){fit.columnNames.push_back("BG_left");fit.columnNames.push_back("BG_right");}
  else fit.columnNames.push_back("BG_exp");
  for(unsigned j=0;j<profile.fractions.size();++j)fit.columnWindowYields.push_back(count*profile.fractions[j]);
  if(fit.shapeCovarianceValid) {
    for(unsigned j=0;j<objective.base.kind.size();++j) {
      const int species=objective.base.species[j];
      if(species>=0) {
        const double error=std::sqrt(fit.shapeCovariance[j][j]);
        if(objective.base.kind[j]==0)fit.meanError[species]=error;
        else if(objective.base.kind[j]==1)fit.sigmaError[species]=fit.components[species].sigma*error;
      }
    }
    result.alphaError=result.alpha*std::sqrt(fit.shapeCovariance.back().back());
    bool finiteErrors=NuclearFit::Finite(result.alphaError) && result.alphaError>0.;
    for(unsigned j=0;j<n;++j)if(fitted.components[j].active)
      finiteErrors=finiteErrors && NuclearFit::Finite(fit.meanError[j]) && fit.meanError[j]>0. &&
                   NuclearFit::Finite(fit.sigmaError[j]) && fit.sigmaError[j]>0.;
    if(!finiteErrors) {
      fit.shapeCovarianceValid=false;fit.shapeCovariance.clear();
      fit.meanError.assign(n,-1.);fit.sigmaError.assign(n,-1.);result.alphaError=-1.;
    }
  }
  fit.componentCounts.assign(fitted.bins.size(),Vec(n,0.));fit.totalCounts.assign(fitted.bins.size(),0.);
  fit.smoothCounts.assign(fitted.bins.size(),0.);fit.effectiveCount.assign(n,0.);
  fit.finite=profile.converged && !templates.empty();fit.pearson=0.;
  for(unsigned b=0;b<templates.size();++b) {
    for(unsigned j=0;j<profile.fractions.size();++j) {
      const double value=fit.columnWindowYields[j]*templates[b][j];
      fit.totalCounts[b]+=value;if(j<n)fit.componentCounts[b][j]=value;else fit.smoothCounts[b]+=value;
    }
    const double mu=fit.totalCounts[b];
    if(!NuclearFit::Finite(mu) || mu<0. || (mu==0. && fitted.bins[b].count>0.))fit.finite=false;
    if(mu>0.) {
      fit.pearson+=std::pow(fitted.bins[b].count-mu,2)/mu;
      for(unsigned j=0;j<n;++j)fit.effectiveCount[j]+=fit.componentCounts[b][j]*fit.componentCounts[b][j]/mu;
    }
  }
  fit.finite=fit.finite && NuclearFit::Finite(fit.dataDeviance) &&
             fit.dataDeviance<NuclearFit::invalid && NuclearFit::Finite(fit.pearson);
  // This is only the conditional fixed-shape amplitude covariance.  It is not
  // an uncomputed joint shape/yield covariance, and inactive columns stay masked.
  std::vector<unsigned> active;
  for(unsigned j=0;j<profile.fractions.size();++j)if(profile.fractions[j]>1.e-12)active.push_back(j);
  fit.yieldCovarianceColumnActive.assign(profile.fractions.size(),false);
  for(unsigned j=0;j<active.size();++j)fit.yieldCovarianceColumnActive[active[j]]=true;
  Matrix hessian(active.size(),Vec(active.size(),0.));
  for(unsigned j=0;j<active.size();++j)for(unsigned k=0;k<active.size();++k)
    hessian[j][k]=profile.hessian[active[j]][active[k]];
  fit.fixedShapeYieldCovariance.assign(profile.fractions.size(),Vec(profile.fractions.size(),0.));
  fit.yieldCovarianceValid=profile.converged && !active.empty();
  for(unsigned k=0;k<active.size() && fit.yieldCovarianceValid;++k) {
    Vec unit(active.size(),0.),inverse;unit[k]=1.;
    if(!NuclearFit::Solve(hessian,unit,inverse)){fit.yieldCovarianceValid=false;break;}
    for(unsigned j=0;j<active.size();++j)
      fit.fixedShapeYieldCovariance[active[j]][active[k]]=count*inverse[j];
  }
  Matrix activeCovariance(active.size(),Vec(active.size(),0.));
  if(fit.yieldCovarianceValid) {
    for(unsigned j=0;j<active.size();++j)for(unsigned k=0;k<active.size();++k)
      activeCovariance[j][k]=fit.fixedShapeYieldCovariance[active[j]][active[k]];
    fit.yieldCovarianceValid=PositiveDefinite(activeCovariance);
  }
  if(!fit.yieldCovarianceValid)fit.fixedShapeYieldCovariance.clear();
  for(unsigned j=0;j<chosen.size();++j) {
    double low,high,step;std::string name;ParameterLimits(objective,j,low,high,step,name);
    fit.atBoundary=fit.atBoundary || std::min(chosen[j]-low,high-chosen[j])<=
      NuclearFit::Number(cfg,"boundaryTolerance")*std::max(1.,high-low);
  }
  fit.conditionalDropDeltaD.assign(n,-1.);
  for(unsigned j=0;j<n;++j)if(fitted.components[j].active) {
    Request dropped=fitted;dropped.components[j].active=false;
    const Matrix droppedTemplates=BuildTemplates(dropped,cfg,model,result.alpha);
    Vec initial=profile.fractions;if(j<initial.size())initial[j]=0.;
    const NuclearFit::ProfileResult reduced=NuclearFit::Profile(dropped.bins,droppedTemplates,initial,cfg);
    if(reduced.converged)fit.conditionalDropDeltaD[j]=std::max(0.,reduced.dataDeviance-fit.dataDeviance);
  }
  fit.conditionalNdf=int(fitted.bins.size())-int(chosen.size())-int(active.size());
  return result;
}
} // namespace RestartResponse
#endif
