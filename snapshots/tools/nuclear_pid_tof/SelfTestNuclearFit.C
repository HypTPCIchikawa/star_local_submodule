void SelfTestNuclearFit(const char*,const char*);
#ifndef __CINT__
#include "NuclearFitCommon.h"
namespace NuclearFitTest {
using namespace NuclearFit;
unsigned checks=0;
void Check(bool ok,const char* message){++checks;Require(ok,message);}
bool Near(double a,double b,double tolerance){return std::fabs(a-b)<=tolerance*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
Request Toy(){Request r;r.nativeAxisMin=0;r.nativeAxisMax=100;r.backgroundModel="linear";r.freeBackgroundSlope=false;r.backgroundWindowYield=3000;
  Component c;c.name="toy";c.mean=52;c.sigma=4;c.startWindowYield=40000;r.components.push_back(c);
  for(int b=0;b<1000;++b){double lo=.1*b,hi=.1*(b+1),n=40000*GaussianIntegral(lo,hi,55,3)+3000*.001;r.bins.push_back(Bin(lo,hi,n));}return r;}
}
#endif
void SelfTestNuclearFit(const char* settings,const char* output){
#ifndef __CINT__
  using namespace NuclearFitTest;
  try{
    Require(gSystem->AccessPathName(output),"Existing numeric result will not be replaced");Settings cfg=ReadSettings(settings);ValidateSettings(cfg);
    Check(Near(GaussianIntegral(-100,100,0,1),1,1.e-15),"Gaussian normalization");
    Check(GaussianIntegral(30,31,0,1)>0,"Stable erfc tail");
    Check(GaussianIntegral(0,1,0,0)==0,"Nonpositive sigma guard");
    const char* models[]={"linear","exponential"};
    for(int m=0;m<2;++m)for(int k=-5;k<=5;++k)for(int basis=0;basis<(m?1:2);++basis){double sum=0;
      for(int b=0;b<100;++b)sum+=BackgroundIntegral(.1*b,.1*(b+1),0,10,models[m],k,basis);
      Check(Near(sum,1,1.e-12),"Positive BG integral normalization");Check(BackgroundIntegral(3,4,0,10,models[m],k,basis)>0,"Positive BG basis");}
    std::vector<Bin> two;two.push_back(Bin(0,1,5));two.push_back(Bin(1,2,9));Matrix identity(2,Vec(2,0));identity[0][0]=identity[1][1]=1;
    ProfileResult p=Profile(two,identity,Vec(2,0),cfg);Check(p.converged,"All-zero starting amplitudes converge");Check(Near(p.fractions[0]*14,5,1.e-6)&&Near(p.fractions[1]*14,9,1.e-6),"Independent Poisson yields");
    Check(Near(p.hessian[0][0],14./5,1.e-6)&&Near(p.hessian[1][1],14./9,1.e-6),"Normalized amplitude Hessian units");
    std::vector<Bin> zero;zero.push_back(Bin(0,1,0));zero.push_back(Bin(1,2,10));ProfileResult boundary=Profile(zero,identity,Vec(2,.5),cfg);
    Check(boundary.converged,"Zero yield AMP boundary KKT");Check(boundary.fractions[0]<1.e-12,"Zero yield is not forced positive");
    Request request=Toy();Matrix t=BuildTemplates(request,cfg);for(unsigned j=0;j<t[0].size();++j){double sum=0;for(unsigned b=0;b<t.size();++b)sum+=t[b][j];Check(Near(sum,1,1.e-12),"Fit-window template normalization");}
    Result fit=Fit(request,cfg);Check(fit.finite&&fit.amplitudeConverged&&fit.pointAccepted&&fit.status==0,"Asimov free shape fit");
    Check(std::fabs(fit.components[0].mean-55)<.02&&std::fabs(fit.components[0].sigma-3)<.02,"Free mean/sigma recovered");
    Check(fit.dataDeviance<.05,"Asimov deviance");Check(fit.yieldCovarianceValid,"Conditional AMP covariance");
    Check(fit.conditionalDropDeltaD[0]>100,"Conditional component-removal diagnostic works");
    Check(Near(fit.shapePoint[0],fit.components[0].mean,1.e-12)&&Near(std::exp(fit.shapePoint[1]),fit.components[0].sigma,1.e-12),"Export covariance and component point match");
    Check(fit.shapeCovarianceValid&&fit.meanError[0]>0&&fit.sigmaError[0]>0,"Profiled shape covariance is available");
    for(unsigned b=0;b<fit.totalCounts.size();++b)Check(Near(fit.totalCounts[b],fit.componentCounts[b][0]+fit.smoothCounts[b],1.e-13),"Prediction component count closure");
    Request fixed=request;fixed.components[0].mean=55;fixed.components[0].sigma=3;fixed.components[0].freeMean=false;fixed.components[0].freeSigma=false;Result f=Fit(fixed,cfg);
    Check(f.status==0&&f.shapePoint.empty()&&f.amplitudeConverged,"Fixed shapes need no Minuit");Check(f.yieldCovarianceColumnActive.size()==3,"BG-left/right covariance column map");
    Request averaged=fixed;averaged.components[0].nodes.push_back(Node(.5,55,3));averaged.components[0].nodes.push_back(Node(.5,55,3));Matrix avg=BuildTemplates(averaged,cfg);
    for(unsigned b=0;b<avg.size();++b)Check(Near(avg[b][0],f.totalTemplates[b][0],1.e-13),"Constant GL kernel equals centre template");
    Request empty=fixed;for(unsigned b=0;b<empty.bins.size();++b)empty.bins[b].count=0;Result none=Fit(empty,cfg);Check(none.status==-2&&!none.finite,"Empty window remains unsupported");
    bool rejected=false;try{Request bad=request;bad.components[0].sigma=Number(cfg,"minimumSigma")/2;Fit(bad,cfg);}catch(const std::exception&){rejected=true;}Check(rejected,"Invalid width seed is not silently clipped");
    std::ofstream out(output);Require(bool(out),"Cannot write new numeric selftest");out<<"{\"passed\":true,\"checks\":"<<checks<<",\"scope\":\"Numerical templates/Poisson AMP/MIGRAD/conditional covariance only; NOT isotope identity or physics validation\"}\n";
    std::cout<<"NUCLEAR_FIT_SELFTEST_PASS checks="<<checks<<std::endl;
  }catch(const std::exception& e){std::cerr<<"SelfTestNuclearFit: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
