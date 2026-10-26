// Synthetic mathematical tests, not analysis cuts or physics acceptance criteria.
int TestPeakResponseFit(const char*,const char*);
#ifndef __CINT__
#include "PeakResponseFit.h"
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
unsigned ResponseChecks=0;
void ResponseCheck(bool ok,const std::string& label) {
  ++ResponseChecks;if(!ok)throw std::runtime_error(label);
}
void ResponseNear(double actual,double expected,double tolerance,const std::string& label) {
  ResponseCheck(NuclearFit::Finite(actual) && NuclearFit::Finite(expected) &&
                std::fabs(actual-expected)<=tolerance*std::max(1.,std::fabs(expected)),label);
}
NuclearFit::Request ResponseToy(const std::string& family,double alpha) {
  NuclearFit::Request request;request.nativeAxisMin=0.;request.nativeAxisMax=120.;
  request.backgroundModel="linear";request.backgroundWindowYield=4000.;
  const char* names[]={"d","t","he3","he4"};
  for(unsigned j=0;j<4;++j) {
    NuclearFit::Component component;component.name=names[j];component.mean=10.*(j+1);
    component.sigma=2.;component.active=j==3;component.freeMean=j==3;component.freeSigma=j==3;
    component.startWindowYield=j==3 ? 20000. : 0.;request.components.push_back(component);
  }
  NuclearFit::Component& target=request.components[3];
  target.nodes.push_back(NuclearFit::Node(5./18.,39.5,1.9));
  target.nodes.push_back(NuclearFit::Node(8./18.,40.,2.));
  target.nodes.push_back(NuclearFit::Node(5./18.,40.5,2.1));
  // Construct adjacent endpoints once; no floating-point gaps between native bins.
  double lo=18.;
  for(unsigned b=0;b<220;++b) {
    const double hi=18.+.2*(b+1);request.bins.push_back(NuclearFit::Bin(lo,hi,0.));lo=hi;
  }
  // The same exact bin integrals generate and fit each Asimov response.
  for(unsigned b=0;b<request.bins.size();++b) {
    const NuclearFit::Bin& bin=request.bins[b];
    const double acceptance=RestartResponse::KernelIntegral(target,18.,62.,family,alpha);
    const double probability=RestartResponse::KernelIntegral(target,bin.lo,bin.hi,family,alpha)/acceptance;
    const double bg=NuclearFit::BackgroundIntegral(bin.lo,bin.hi,18.,62.,"linear",0.,0)+
                    NuclearFit::BackgroundIntegral(bin.lo,bin.hi,18.,62.,"linear",0.,1);
    request.bins[b].count=20000.*probability+2000.*bg;
  }
  return request;
}
void ResponseSeedShift(NuclearFit::Request& request) {
  NuclearFit::Component& target=request.components[3];
  const double oldMean=target.mean,oldSigma=target.sigma;
  target.mean=39.7;target.sigma=2.2;target.startWindowYield=18000.;request.backgroundWindowYield=4500.;
  for(unsigned q=0;q<target.nodes.size();++q) {
    target.nodes[q].mean+=target.mean-oldMean;target.nodes[q].sigma*=target.sigma/oldSigma;
  }
}
void ResponseIntegralTests() {
  const double infinity=std::numeric_limits<double>::infinity();
  const char* families[]={"gaussian","gaussian_left_exp","gaussian_right_exp"};
  const double alphas[]={.2,.5,1.,2.,5.,30.};
  for(unsigned f=0;f<3;++f)for(unsigned a=0;a<6;++a) {
    const std::string family=families[f];const double alpha=alphas[a];
    ResponseNear(RestartResponse::StandardIntegral(-infinity,infinity,family,alpha),1.,2.e-14,"full-line normalization");
    ResponseNear(RestartResponse::Cdf(-infinity,40.,2.,family,alpha),0.,1.e-14,"CDF lower infinity");
    ResponseNear(RestartResponse::Cdf(infinity,40.,2.,family,alpha),1.,2.e-14,"CDF upper infinity");
    for(int i=-40;i<40;++i) {
      const double lo=.25*i,mid=lo+.11,hi=lo+.25;
      const double whole=RestartResponse::StandardIntegral(lo,hi,family,alpha);
      ResponseCheck(whole>=0. && whole<=1.,"nonnegative bin probability");
      ResponseNear(RestartResponse::StandardIntegral(lo,mid,family,alpha)+
                   RestartResponse::StandardIntegral(mid,hi,family,alpha),whole,3.e-15,"bin additivity");
      if(f>0)ResponseNear(whole,RestartResponse::StandardIntegral(-hi,-lo,
                          f==1 ? "gaussian_right_exp" : "gaussian_left_exp",alpha),3.e-15,"left/right mirror");
    }
    if(f>0) {
      const double join=f==1 ? -alpha : alpha,epsilon=1.e-8;
      ResponseNear(RestartResponse::StandardDensity(join-epsilon,family,alpha),
                   RestartResponse::StandardDensity(join+epsilon,family,alpha),1.e-7,"join value continuity");
      ResponseNear(RestartResponse::StandardDensityDerivative(join-epsilon,family,alpha),
                   RestartResponse::StandardDensityDerivative(join+epsilon,family,alpha),1.e-7,"join derivative continuity");
      const double peak=RestartResponse::StandardDensity(0.,family,alpha);
      ResponseCheck(peak>=RestartResponse::StandardDensity(-.1,family,alpha) &&
                    peak>=RestartResponse::StandardDensity(.1,family,alpha),"mu is the response mode");
    }
  }
  for(int b=-40;b<40;++b) {
    const double lo=.25*b,hi=lo+.25;
    ResponseNear(RestartResponse::StandardIntegral(lo,hi,"gaussian_left_exp",30.),
                 NuclearFit::GaussianIntegral(lo,hi,0.,1.),1.e-14,"large alpha Gaussian limit");
    // A very narrow far-tail bin exercises expm1 rather than CDF subtraction.
    const double probability=RestartResponse::StandardIntegral(-30.,-30.+1.e-9,"gaussian_left_exp",1.);
    ResponseCheck(probability>0. && NuclearFit::Finite(probability),"stable small far-tail bin");
  }
}
void ResponseFitTests(const NuclearFit::Settings& cfg) {
  NuclearFit::Request gaussian=ResponseToy("gaussian",2.);ResponseSeedShift(gaussian);
  const NuclearFit::Matrix originalTemplates=NuclearFit::BuildTemplates(gaussian,cfg);
  const NuclearFit::Matrix responseTemplates=RestartResponse::BuildTemplates(gaussian,cfg,"gaussian",2.);
  ResponseCheck(originalTemplates==responseTemplates,"Gaussian templates are exactly the original implementation");
  const NuclearFit::Result original=NuclearFit::Fit(gaussian,cfg);
  const RestartResponse::Result delegated=RestartResponse::Fit(gaussian,cfg,"gaussian",2.);
  ResponseNear(delegated.fit.pearson,original.pearson,1.e-12,"free Gaussian Pearson parity");
  ResponseNear(delegated.fit.dataDeviance,original.dataDeviance,1.e-12,"free Gaussian deviance parity");
  ResponseCheck(delegated.fit.totalTemplates==original.totalTemplates,"fitted Gaussian template parity");
  ResponseCheck(delegated.fit.conditionalNdf==original.conditionalNdf,"Gaussian NDF parity");
  ResponseCheck(!delegated.optimizerDiagnosticsAvailable && delegated.migradStatus==-1 &&
                !delegated.hesseAttempted,"Gaussian delegate does not invent separate optimizer-phase observations");
  const char* families[]={"gaussian_left_exp","gaussian_right_exp"};
  for(unsigned f=0;f<2;++f) {
    const double truthAlpha=1.6;
    ResponseCheck(NuclearFit::Number(cfg,"responseAlphaMin")<truthAlpha &&
                  truthAlpha<NuclearFit::Number(cfg,"responseAlphaMax"),"synthetic truth is inside configured alpha support");
    NuclearFit::Request request=ResponseToy(families[f],truthAlpha);ResponseSeedShift(request);
    const RestartResponse::Result fit=RestartResponse::Fit(request,cfg,families[f],2.);
    ResponseCheck(fit.optimizerDiagnosticsAvailable && fit.minimizeReturned &&
                  fit.minimizeConverged && fit.migradStatus==0 &&
                  NuclearFit::Finite(fit.migradEdm) && fit.migradEdm>=0.,"actual successful MIGRAD diagnostics retained");
    ResponseCheck(fit.hesseAttempted && fit.hesseSucceeded && fit.hesseStatus==0 &&
                  fit.fit.status==0,"successful Hesse is recorded separately from MIGRAD");
    ResponseCheck(fit.fit.finite && fit.fit.amplitudeConverged && fit.fit.pointAccepted,"Asimov profile solved");
    ResponseNear(fit.fit.components[3].mean,40.,2.e-3,"Asimov mode recovery");
    ResponseNear(fit.fit.components[3].sigma,2.,2.e-2,"Asimov core-width recovery");
    ResponseNear(fit.alpha,truthAlpha,3.e-2,"Asimov tail-alpha recovery");
    ResponseCheck(fit.fit.dataDeviance<1.,"Asimov deviance near zero");
    ResponseCheck(fit.fit.shapeCovarianceValid && fit.alphaError>0.,"actual profile-Hesse alpha uncertainty");
    ResponseCheck(fit.fit.yieldCovarianceValid,"actual conditional amplitude covariance");
    ResponseCheck(fit.freeParameters==3,"tail adds exactly one shape coordinate to linear-BG Gaussian");
    unsigned activeYields=0;
    for(unsigned j=0;j<fit.fit.yieldCovarianceColumnActive.size();++j)
      activeYields+=fit.fit.yieldCovarianceColumnActive[j];
    ResponseCheck(fit.fit.conditionalNdf==int(request.bins.size())-fit.freeParameters-int(activeYields),"tail NDF counts actual shape and yield dimensions");
    for(unsigned j=0;j<3;++j) {
      ResponseCheck(fit.fit.components[j].mean==request.components[j].mean &&
                    fit.fit.components[j].sigma==request.components[j].sigma &&
                    !fit.fit.components[j].active && !fit.fit.components[j].freeMean &&
                    !fit.fit.components[j].freeSigma,"foreign original shapes are unchanged and inactive");
    }
    for(unsigned b=0;b<fit.fit.totalCounts.size();++b) {
      double sum=fit.fit.smoothCounts[b];
      for(unsigned j=0;j<4;++j)sum+=fit.fit.componentCounts[b][j];
      ResponseNear(sum,fit.fit.totalCounts[b],2.e-14,"saved prediction is sum of actual components and BG");
    }
  }
  NuclearFit::Request invalid=gaussian;invalid.components[2].freeMean=true;
  bool rejected=false;try{RestartResponse::Fit(invalid,cfg,"gaussian",2.);}catch(const std::exception&){rejected=true;}
  ResponseCheck(rejected,"inactive foreign peak cannot become free");
  invalid=gaussian;invalid.components[2].active=true;invalid.components[2].freeMean=true;invalid.components[2].freeSigma=true;
  rejected=false;try{RestartResponse::Fit(invalid,cfg,"gaussian_left_exp",2.);}catch(const std::exception&){rejected=true;}
  ResponseCheck(rejected,"a second active peak is rejected");
  // Deliberately inadequate computational budget, not a physics selection.
  NuclearFit::Settings limited=cfg;limited["shapeMaxCalls"]="1";
  NuclearFit::Request difficult=ResponseToy("gaussian_left_exp",1.6);ResponseSeedShift(difficult);
  const RestartResponse::Result unfinished=RestartResponse::Fit(difficult,limited,"gaussian_left_exp",2.);
  ResponseCheck(unfinished.optimizerDiagnosticsAvailable && !unfinished.minimizeConverged &&
                unfinished.fit.status!=0,"call-limited MIGRAD cannot become numerical success");
  ResponseCheck(!unfinished.hesseAttempted && !unfinished.hesseSucceeded &&
                !unfinished.fit.shapeCovarianceValid,"Hesse/covariance skipped for unconverged MIGRAD");
}
void ResponsePhaseTests() {
  ResponseCheck(RestartResponse::MigradConverged(true,0,0.),"finite zero-EDM MIGRAD may converge");
  ResponseCheck(!RestartResponse::MigradConverged(false,0,0.),"false Minimize return is not convergence");
  ResponseCheck(!RestartResponse::MigradConverged(true,4,0.),"nonzero MIGRAD status is not convergence");
  ResponseCheck(!RestartResponse::MigradConverged(true,0,-1.),"negative MIGRAD EDM is invalid");
  ResponseCheck(!RestartResponse::MigradConverged(true,0,std::numeric_limits<double>::quiet_NaN()),
                "NaN MIGRAD EDM is invalid");
  RestartResponse::Result response;response.optimizerDiagnosticsAvailable=true;
  response.migradStatus=4;response.minimizeConverged=false;
  response.fit.pointAccepted=true;response.hesseAttempted=true;response.hesseSucceeded=true;response.hesseStatus=0;
  ResponseCheck(RestartResponse::PhaseStatus(response)==4,"Hesse status zero cannot overwrite MIGRAD failure");
  response.migradStatus=0;
  ResponseCheck(RestartResponse::PhaseStatus(response)==-10,"status zero without converged Minimize stays nonzero");
  response.minimizeConverged=true;response.hesseSucceeded=false;
  ResponseCheck(RestartResponse::PhaseStatus(response)==-11,"false Hesse return with status zero stays nonzero");
  response.hesseStatus=2;
  ResponseCheck(RestartResponse::PhaseStatus(response)==2,"observed nonzero Hesse status is retained");
  response.hesseStatus=0;response.hesseSucceeded=true;
  ResponseCheck(RestartResponse::PhaseStatus(response)==0,"both successful phases permit status zero");
  response.fit.pointAccepted=false;
  ResponseCheck(RestartResponse::PhaseStatus(response)==-12,"rejected candidate cannot receive status zero");
}
}

int TestPeakResponseFit(const char* settingsPath,const char* reportPath) {
  try {
    ResponseChecks=0;
    ResponseCheck(settingsPath && settingsPath[0] && reportPath && reportPath[0],"explicit test settings/report required");
    ResponseCheck(gSystem->AccessPathName(reportPath),"test report must be new");
    const NuclearFit::Settings cfg=NuclearFit::ReadSettings(settingsPath);
    RestartResponse::ValidateSettings(cfg);ResponseIntegralTests();ResponsePhaseTests();ResponseFitTests(cfg);
    std::ofstream out(reportPath);ResponseCheck(out.good(),"create test report");
    out<<"{\n  \"status\": \"passed\",\n  \"checks\": "<<(ResponseChecks+1)
       <<",\n  \"syntheticOnly\": true,\n  \"productionEnabled\": false,\n  \"coefficientsUpdated\": false\n}\n";
    out.close();ResponseCheck(out.good(),"write test report");
    std::cout<<"PEAK_RESPONSE_TEST_PASS checks="<<ResponseChecks<<std::endl;return 0;
  }catch(const std::exception& error) {
    std::cerr<<"PEAK_RESPONSE_TEST_FAIL "<<error.what()<<std::endl;return 1;
  }
}
#endif
