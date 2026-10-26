void SelfTestCoreTail(const char*,const char*);
#ifndef __CINT__
#include "NuclearCoreTailFitCommon.h"
namespace CoreTailTest {
using namespace NuclearCoreTailFit;
unsigned long checks=0;
void Check(bool ok,const char* message){++checks;Require(ok,message);}
bool Near(double a,double b,double eps=1.e-12){return std::fabs(a-b)<=eps*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
Component Core(){Component c;c.name="test_isotope";c.mean=30;c.sigma=2;c.startWindowYield=1000;return c;}
}
#endif
void SelfTestCoreTail(const char* settingsPath,const char* outputPath) {
#ifndef __CINT__
  using namespace CoreTailTest;
  try {
    Require(gSystem->AccessPathName(outputPath),"Existing core+tail test result will not be replaced");
    Settings cfg=NuclearFit::ReadSettings(settingsPath);const Shape declared=ReadShape(cfg);
    Check(declared.fraction>0,"Explicit nonzero test profile");
    const double infinity=std::numeric_limits<double>::infinity();
    Component core=Core();
    const Shape gaussian(0,0,1);Component delegated=BuildMixtureComponent(core,gaussian);
    Check(delegated.nodes.empty(),"f=0 retains exact empty-node Gaussian representation");
    for(int i=-300;i<300;++i) {
      const double lo=core.mean+.1*i,hi=lo+.1;
      Check(ComponentIntegral(delegated,gaussian,0,lo,hi)==NuclearFit::GaussianIntegral(lo,hi,core.mean,core.sigma),
            "f=0 exact Gaussian integral");
    }
    Component gl=core;gl.nodes.push_back(Node(5./18.,29,1.8));gl.nodes.push_back(Node(8./18.,30,2));gl.nodes.push_back(Node(5./18.,31,2.2));
    Check(BuildMixtureComponent(gl,gaussian).nodes.size()==3,"f=0 preserves original quadrature count");
    const Shape profiles[]={declared,Shape(.05,-1,1.5),Shape(.1,0,2),Shape(.3,1,3)};
    for(unsigned p=0;p<sizeof(profiles)/sizeof(profiles[0]);++p) {
      const Shape shape=profiles[p];Component mixed=BuildMixtureComponent(gl,shape);
      Check(RNodeCount(mixed,shape)==3&&mixed.nodes.size()==6,"Three R nodes become six weighted Gaussian nodes");
      Check(Near(KernelIntegral(mixed,-infinity,infinity),1),"Full-line mixture normalization");
      double previous=0;
      for(int i=-300;i<=300;++i) {
        double x=core.mean+.2*i,cdf=KernelIntegral(mixed,-infinity,x);
        Check(Finite(cdf)&&cdf>=previous-1.e-14&&cdf>=0&&cdf<=1+1.e-12,"Finite monotonic mixture CDF");previous=cdf;
        double bin=KernelIntegral(mixed,x,x+.2),conditional=0;
        for(unsigned q=0;q<3;++q)conditional+=RNodeWeight(mixed,shape,q)*ComponentIntegral(mixed,shape,q,x,x+.2);
        Check(Finite(bin)&&bin>=0&&Near(bin,conditional),"Positive finite bin and conditional R-weight closure");
      }
      double sum=0;for(int b=0;b<1000;++b)sum+=KernelIntegral(mixed,-100+.3*b,-100+.3*(b+1));
      Check(Near(sum,1,1.e-11),"Adjacent whole bins normalize without clipping");
      Moments moments=MixtureMoments(BuildMixtureComponent(core,shape));
      const double separation=shape.offsetSeedSigma*core.sigma;
      const double expectedMean=core.mean+shape.fraction*separation;
      const double expectedVariance=core.sigma*core.sigma*((1-shape.fraction)+shape.fraction*shape.widthRatio*shape.widthRatio)+
                                      shape.fraction*(1-shape.fraction)*separation*separation;
      Check(Near(moments.mean,expectedMean)&&Near(moments.sigma,std::sqrt(expectedVariance)),"Mixture moments differ from core parameters by the correct amount");
      Check(mixed.mean==core.mean&&mixed.sigma==core.sigma,"Scalar fit parameters retain core meaning");
      Request request;request.components.push_back(gl);request.backgroundWindowYield=10;
      request.bins.push_back(NuclearFit::Bin(20,30,100));request.bins.push_back(NuclearFit::Bin(30,40,100));
      std::vector<Shape> shapes(1,shape);Request expanded=BuildMixtureRequest(request,shapes);
      NuclearFit::Objective objective(expanded,cfg);NuclearFit::Vec point=objective.Point();
      for(unsigned k=0;k<objective.kind.size();++k) {
        if(objective.kind[k]==0)point[k]+=1.25;
        if(objective.kind[k]==1)point[k]+=std::log(1.4);
      }
      Request moved=objective.Model(point.empty()?0:&point[0]);
      for(unsigned q=0;q<3;++q) {
        Check(Near(moved.components[0].nodes[2*q].mean,gl.nodes[q].mean+1.25),"Common core mean shift follows old engine semantics");
        Check(Near(moved.components[0].nodes[2*q+1].mean-moved.components[0].nodes[2*q].mean,
                   shape.offsetSeedSigma*gl.nodes[q].sigma),"Offset stays in SEED-sigma dE/dx units while fitted width changes");
        Check(Near(moved.components[0].nodes[2*q+1].sigma/moved.components[0].nodes[2*q].sigma,shape.widthRatio),"Tail/core width ratio is preserved during width profiling");
        Check(Near(moved.components[0].nodes[2*q].sigma,gl.nodes[q].sigma*1.4),"Core width uses common fitted scale");
      }
      Check(expanded.components.size()==1,"One isotope remains one yield column, no independent tail yield");
    }
    bool rejected=false;try{ValidateShape(Shape(.5,0,2));}catch(const std::exception&){rejected=true;}Check(rejected,"Equal-weight label ambiguity rejected");
    rejected=false;try{ValidateShape(Shape(.1,0,1));}catch(const std::exception&){rejected=true;}Check(rejected,"Nonzero tail must be broader than core");
    rejected=false;try{ValidateShape(Shape(-.1,0,2));}catch(const std::exception&){rejected=true;}Check(rejected,"Negative tail weight rejected");
    rejected=false;try{ValidateShape(Shape(.1,infinity,2));}catch(const std::exception&){rejected=true;}Check(rejected,"Nonfinite tail offset rejected");
    rejected=false;try{Request request;request.components.push_back(core);std::vector<Shape> empty;BuildMixtureRequest(request,empty);}catch(const std::exception&){rejected=true;}Check(rejected,"Missing per-isotope profile rejected");
    Component wrong=BuildMixtureComponent(core,declared);rejected=false;
    try{ComponentIntegral(wrong,Shape(declared.fraction/2,declared.offsetSeedSigma,declared.widthRatio),0,10,40);}catch(const std::exception&){rejected=true;}Check(rejected,"Mismatched profile cannot silently corrupt gate integration");
    std::ofstream output(outputPath);Require(bool(output),"Cannot CREATE core+tail numeric result");
    output<<"{\"passed\":true,\"checks\":"<<checks<<",\"fixedProfilesOnly\":true,\"fitInvoked\":false,\"productionChanged\":false,\"scope\":\"Mixture mathematics and existing objective transform, not real-data fitting or PID calibration\"}\n";
    output.close();Require(output.good(),"Numeric result close failed");
    std::cout<<"NUCLEAR_CORE_TAIL_SELFTEST_PASS checks="<<checks<<std::endl;
  }catch(const std::exception& error){std::cerr<<"SelfTestCoreTail: "<<error.what()<<std::endl;gSystem->Exit(2);}
#endif
}
