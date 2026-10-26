// Numerical fixtures, not analysis parameters.
void SelfTestLocalNuclearCutoff();
#ifndef __CINT__
#include "FitLocalNuclearCutoff.C"
namespace LocalCutoffTest {
using namespace NuclearLocalCutoff;
int checks=0;
void Check(bool ok,const char* message){++checks;Require(ok,message);}
void Near(double value,double expected,double tolerance,const char* message){Check(Finite(value)&&std::fabs(value-expected)<=tolerance*std::max(1.,std::fabs(expected)),message);}
Settings Config(){Settings c;const char* keys[]={"amplitudeGradientTolerance","amplitudeHessianRidge","shapeTolerance","minimumSigma","maximumSigma","minimumGaussianAcceptance","minimumResponseAcceptance","minimumMode","meanStep","sigmaLogStep","backgroundSlopeStep","exponentialSlopeLimit","boundaryTolerance","amplitudeMaxIterations","amplitudeLineSearchIterations","shapeMaxCalls","meanAxisPadding"};
  const char* values[]={"1e-9","1e-10","1e-6","0.05","120","1e-12","1e-12","0.001","0.05","0.02","0.05","8","1e-6","1000","80","6000","120"};for(unsigned i=0;i<sizeof(keys)/sizeof(keys[0]);++i)c[keys[i]]=values[i];return c;}
void Gates(){std::vector<Choice> c(4);for(int s=0;s<4;++s){c[s].active=true;c[s].mean=5+4*s;c[s].sigma=1+.3*s;}
  for(int target=0;target<4;++target){Regions r=SelectionRegions(c,target,2,0,25);for(int b=0;b<25000;++b){double x=(b+.5)/1000.;bool in=false;for(unsigned k=0;k<r.size();++k)in=in||(x>r[k].lo&&x<r[k].hi);Check(in==(Select(c,x,2)==target),"Analytic best-species partition disagrees with direct selector");}}
  c[1]=c[0];Check(Select(c,c[0].mean,2)==0,"Tie order changed");c[1].active=false;c[2].active=false;c[3].active=false;Check(Select(c,c[0].mean+2*c[0].sigma,2)==-1,"Strict gate boundary changed");
  Regions r;r.push_back(Interval(1.02,1.17));Near(Overlap(r,1,1.1),.8,1e-12,"Fractional raw bin overlap wrong");Near(Overlap(r,1.1,1.2),.7,1e-12,"Fractional raw bin overlap wrong");
  Near(GaussianEngine::Integral(-100,100,10,2),1,1e-12,"Gaussian normalization");Near(LognormalEngine::Integral(0,1e6,10,2),1,1e-12,"Lognormal normalization");Check(std::fabs(GaussianEngine::Integral(6,14,10,2)-LognormalEngine::Integral(6,14,10,2))>.001,"Lognormal integral incorrectly replaced with Gaussian shortcut");
}
template<class E> void Closure(const std::string& background){
  const double means[]={8,16,27,40},sigmas[]={.8,1.1,1.4,2},yields[]={3000,2000,1200,800};const double lo=0,hi=60,step=.1;typename E::Request request;request.nativeAxisMin=lo;request.nativeAxisMax=hi;request.backgroundModel=background;request.backgroundSlope=background=="linear"?0:.7;request.freeBackgroundSlope=background!="linear";
  Slice slice;slice.index=0;slice.rBin=1;slice.lo=.8;slice.hi=.81;slice.fitLo=lo;slice.fitHi=hi;
  const double weights[]={5./18.,8./18.,5./18.};
  for(int q=0;q<3;++q)for(int s=0;s<4;++s){Row r;r.index=0;r.rBin=1;r.lo=slice.lo;r.hi=slice.hi;r.node=q;r.species=s;r.weight=weights[q];r.r=.801+q*.004;r.mean=means[s];r.sigma=sigmas[s];r.supported=true;r.fitLo=lo;r.fitHi=hi;r.gate.mean=means[s]+(q-1)*.07;r.gate.sigma=sigmas[s];r.gate.active=true;slice.rows.push_back(r);}
  for(int s=0;s<4;++s){typename E::Component c;c.name=Species[s];c.mean=means[s];c.sigma=sigmas[s];c.freeMean=false;c.freeSigma=false;c.startWindowYield=1500;for(int q=0;q<3;++q)c.nodes.push_back(typename E::Node(weights[q],means[s],sigmas[s]));request.components.push_back(c);}
  request.backgroundWindowYield=1000;TH1D raw("cutoff_test_raw","",600,lo,hi);raw.SetDirectory(0);raw.Sumw2();
  for(int b=1;b<=600;++b){double a=raw.GetBinLowEdge(b),z=raw.GetXaxis()->GetBinUpEdge(b),count=0;
    for(int s=0;s<4;++s)count+=yields[s]*E::Integral(a,z,means[s],sigmas[s])/E::Integral(lo,hi,means[s],sigmas[s]);
    if(background=="linear")count+=500*NuclearFit::BackgroundIntegral(a,z,lo,hi,background,0,0)+700*NuclearFit::BackgroundIntegral(a,z,lo,hi,background,0,1);else count+=1200*NuclearFit::BackgroundIntegral(a,z,lo,hi,background,.7,0);
    raw.SetBinContent(b,count);raw.SetBinError(b,std::sqrt(count));request.bins.push_back(typename E::Bin(a,z,count));
  }
  if(background=="exponential")request.backgroundSlope=0;
  typename E::Result fit=E::Fit(request,Config());Check(fit.finite&&fit.amplitudeConverged&&fit.pointAccepted&&fit.status==0&&!fit.atBoundary,"Exact mixture fit did not converge");
  for(int s=0;s<4;++s)Near(fit.columnWindowYields[s],yields[s],2e-5,"Known mixture yield not recovered");
  for(int target=0;target<4;++target){Estimate estimate=Integrate<E>(fit,request,slice,raw,1,target,2,1e-9,1e-9);Check(estimate.available&&estimate.errorAvailable&&estimate.covered&&estimate.axisCovered&&estimate.identitySupport&&!estimate.boundaryTarget,"Exact mixture purity unavailable");
    double expectedSignal=0,expectedTotal=0;for(int q=0;q<3;++q){std::vector<Choice> choices(4);for(int s=0;s<4;++s)choices[s]=slice.rows[4*q+s].gate;Regions regions=SelectionRegions(choices,target,2,lo,hi);
      for(unsigned k=0;k<regions.size();++k){double a=regions[k].lo,z=regions[k].hi;for(int s=0;s<4;++s){double n=weights[q]*yields[s]*E::Integral(a,z,means[s],sigmas[s])/E::Integral(lo,hi,means[s],sigmas[s]);expectedTotal+=n;if(s==target)expectedSignal+=n;}
        if(background=="linear")expectedTotal+=weights[q]*(500*NuclearFit::BackgroundIntegral(a,z,lo,hi,background,0,0)+700*NuclearFit::BackgroundIntegral(a,z,lo,hi,background,0,1));else expectedTotal+=weights[q]*1200*NuclearFit::BackgroundIntegral(a,z,lo,hi,background,.7,0);
      }}
    Near(estimate.signal,expectedSignal,3e-5,"Selected signal closure");Near(estimate.total,expectedTotal,3e-5,"Selected total closure");Near(estimate.purity,expectedSignal/expectedTotal,3e-5,"Selected purity closure");Check(estimate.rawLower<=estimate.rawCount+1e-8&&estimate.rawCount<=estimate.rawUpper+1e-8,"Raw count bin-boundary bounds fail");Check(estimate.localResidual<1e-4,"Exact closure local residual nonzero");
  }
}
}
#endif
void SelfTestLocalNuclearCutoff(){
#ifndef __CINT__
  try{LocalCutoffTest::Gates();LocalCutoffTest::Closure<NuclearLocalCutoff::GaussianEngine>("linear");LocalCutoffTest::Closure<NuclearLocalCutoff::GaussianEngine>("exponential");LocalCutoffTest::Closure<NuclearLocalCutoff::LognormalEngine>("linear");LocalCutoffTest::Closure<NuclearLocalCutoff::LognormalEngine>("exponential");bool ratio=false;LocalCutoffTest::Near(NuclearLocalCutoff::OmittedBound(1000,1e-14,.5,ratio),2e-11,1e-12,"Omission bound ratio");LocalCutoffTest::Check(ratio,"Omission ratio availability");LocalCutoffTest::Near(NuclearLocalCutoff::OmittedBound(1000,0,0,ratio),1000,1e-12,"Underflow must not imply zero omitted background");LocalCutoffTest::Check(!ratio,"Omission underflow availability");std::cout<<"NUCLEAR_LOCAL_CUTOFF_SELFTEST_PASS checks="<<LocalCutoffTest::checks<<std::endl;}catch(const std::exception& e){std::cerr<<"Local cutoff selftest: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
