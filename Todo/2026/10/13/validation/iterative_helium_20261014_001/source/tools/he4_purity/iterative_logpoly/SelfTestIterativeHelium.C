void SelfTestIterativeHelium(const char* settings,const char* output);
#ifndef __CINT__
#include "IterativeHeliumFit.C"
namespace HeliumIterativeTest {
void Require(bool ok,const char* name,int& tests){++tests;if(!ok)throw std::runtime_error(std::string("iterative selftest: ")+name);}
Data Toy(){
  Data d;d.seed.index=0;d.seed.lo=.4;d.seed.hi=.41;d.seed.anchorMean=100;d.seed.anchorSigma=6;d.seed.bg="linear";d.seed.tPresent=true;
  d.h=new TH1D("iterative_test_raw","",240,0,120);d.h->SetDirectory(0);d.raw.lo=.4;d.raw.hi=.41;d.raw.r=.405;d.raw.fitLo=0;d.raw.fitHi=120;d.raw.count=100000;
  d.shape.mu[0]=20;d.shape.mu[1]=40;d.shape.mu[2]=70;d.shape.sigma[0]=2;d.shape.sigma[1]=3;d.shape.sigma[2]=4;
  for(int b=1;b<=240;++b){GlobalHelium::Bin bin;bin.lo=d.h->GetBinLowEdge(b);bin.hi=bin.lo+.5;bin.n=0;d.raw.bins.push_back(bin);}
  d.a.push_back(.2);d.a.push_back(.1);d.a.push_back(.05);d.a.push_back(.3);d.a.push_back(.35);RNodes(d,3);return d;
}
}
#endif
void SelfTestIterativeHelium(const char* settings,const char* output){
#ifndef __CINT__
  using namespace HeliumIterative;using namespace HeliumIterativeTest;
  try{
    if(!gSystem->AccessPathName(output))throw std::runtime_error("new selftest output required");
    if(gSystem->Load("libMinuit2")<0)throw std::runtime_error("Minuit2 unavailable");
    Settings cfg=ReadSettings(settings);cfg["phase"]="free";cfg["penaltyStrength"]="0";
    Target none;Study study(cfg,none);Vec original=GlobalHelium::Legacy();int tests=0;
    Vec normal=ToNormal(original,study.mid,study.half),back=FromNormal(normal,study.mid,study.half);
    for(int j=0;j<14;++j)Require(std::fabs(back[j]-original[j])<1.e-10,"normalized polynomial roundtrip",tests);
    Vec optimized=ToOptim(original,study.mid,study.half,Value(cfg,"minimumSigma"));back=FromOptim(optimized,study.mid,study.half,Value(cfg,"minimumSigma"));
    for(int j=0;j<14;++j)Require(std::fabs(back[j]-original[j])<1.e-10,"positive width map roundtrip",tests);
    Vec alternate=normal;alternate[4]=Value(cfg,"minimumSigma")+.01;alternate[5]=0;alternate[6]=1;
    Vec family=FromNormal(alternate,study.mid,study.half);back=FromOptim(ToOptim(family,study.mid,study.half,Value(cfg,"minimumSigma")),study.mid,study.half,Value(cfg,"minimumSigma"));
    for(int j=0;j<14;++j)Require(std::fabs(back[j]-family[j])<1.e-10,"positive convex family including negative Bernstein control",tests);
    alternate[4]=1;alternate[5]=.1;alternate[6]=-.2;family=FromNormal(alternate,study.mid,study.half);back=FromOptim(ToOptim(family,study.mid,study.half,Value(cfg,"minimumSigma")),study.mid,study.half,Value(cfg,"minimumSigma"));
    for(int j=0;j<14;++j)Require(std::fabs(back[j]-family[j])<1.e-10,"positive concave quadratic family",tests);
    alternate[4]=Value(cfg,"minimumSigma");alternate[5]=alternate[6]=0;family=FromNormal(alternate,study.mid,study.half);bool reject=false;try{ToOptim(family,study.mid,study.half,Value(cfg,"minimumSigma"));}catch(const std::runtime_error&){reject=true;}
    Require(reject,"strict positive seed margin fails loudly at equality",tests);
    Data d=Toy();Vec acceptance;Matrix t=study.Templates(d,d.shape,Vec(),&acceptance);
    for(unsigned col=0;col<t[0].size();++col){double sum=0;for(unsigned b=0;b<t.size();++b)sum+=t[b][col];Require(std::fabs(sum-1)<1.e-12,"window conditional templates normalized",tests);}
    for(unsigned b=0;b<t.size();++b){double n=0;for(unsigned col=0;col<d.a.size();++col)n+=d.raw.count*d.a[col]*t[b][col];d.raw.bins[b].n=n;}
    Profile p=WindowProfile(d.raw,t,Vec(d.a.size(),.2),cfg);
    Require(p.converged,"window yields converged",tests);Require(p.d<1.e-7,"window yields recover exact synthetic mixture",tests);
    double sum=0;for(unsigned j=0;j<p.a.size();++j)sum+=p.a[j];Require(std::fabs(sum-1)<1.e-7,"unconstrained nonnegative window yields recover total",tests);
    GlobalHelium::Slice identityData;identityData.count=10000;GlobalHelium::Bin b;b.lo=0;b.hi=1;b.n=10000;identityData.bins.push_back(b);b.lo=1;b.hi=2;b.n=0;identityData.bins.push_back(b);
    Matrix identity(2,Vec(2,0));identity[0][0]=identity[1][1]=1;p=WindowProfile(identityData,identity,Vec(2,.5),cfg);
    Require(p.converged&&p.a[1]<1.e-12,"zero yield active boundary",tests);
    identityData.bins[0].n=7000;identityData.bins[1].n=3000;p=WindowProfile(identityData,identity,Vec(2,0),cfg);
    Require(p.converged&&std::fabs(p.a[0]-.7)<1.e-7,"zero start nonnegative yield optimizer",tests);
    double anchorMean=100,anchorSigma=6,height=100,uMean=-3,uSigma=1.2,binWidth=.5;
    double full=height*uSigma*std::sqrt(2*TMath::Pi())/(binWidth/anchorSigma),mu=anchorMean+anchorSigma*uMean,sigma=anchorSigma*uSigma;
    double accepted=Gaussian(60,120,mu,sigma),window=full*accepted;
    for(int j=0;j<120;++j){double lo=60+j*.5,hi=lo+.5;double old=height*uSigma*std::sqrt(2*TMath::Pi())/(binWidth/anchorSigma)*Gaussian((lo-anchorMean)/anchorSigma,(hi-anchorMean)/anchorSigma,uMean,uSigma);double current=window*Gaussian(lo,hi,mu,sigma)/accepted;Require(std::fabs(old-current)<1.e-10*std::max(1.,old),"old TF1 height/bin -> conditional-window count exact reproduction",tests);}
    Target target;target.present=true;target.coefficients.assign(14,0);target.coefficients[0]=40;target.coefficients[4]=3;target.coefficients[7]=70;target.coefficients[11]=4;
    Settings softCfg=cfg;softCfg["phase"]="regularized";softCfg["penaltyStrength"]="2";Study soft(softCfg,target);
    d.weight.info[0]=30;d.weight.meanScale[0]=d.weight.widthScale[0]=5;d.weight.info[1]=0;d.shape.mu[1]=42;d.shape.sigma[1]=4;
    Require(std::fabs(soft.Penalty(d,d.shape)-12)<1.e-12,"fixed information raw-sigma QR-consistent penalty",tests);
    d.shape.mu[1]=40;d.shape.sigma[1]=3;d.a.assign(5,.2);d.shape.mu[1]+=1;d.shape.sigma[1]*=1.1;
    FitLocal(study,d,Vec());Require(d.finite&&d.ampConverged&&d.accepted&&d.dataD<1,"free detector-support fit starts near per-slice solution",tests);
    d.shape.mu[1]=80;d.shape.mu[2]=40;Vec before=d.a;CanonicalizeFree(d);
    Require(d.labelSwapped&&d.shape.mu[1]<d.shape.mu[2]&&d.a[1]==before[2],"free label permutation preserves paired yield",tests);
    Data common=Toy();Settings commonCfg=cfg;commonCfg["phase"]="common";commonCfg["penaltyStrength"]="0";Study commonStudy(commonCfg,target);
    Matrix folded=commonStudy.Templates(common,common.shape,target.coefficients);
    common.nodes.clear();common.nodeWeights.clear();RNodes(common,5);Matrix folded5=commonStudy.Templates(common,common.shape,target.coefficients);
    for(unsigned b=0;b<folded.size();++b)for(unsigned col=0;col<folded[b].size();++col)Require(std::fabs(folded[b][col]-folded5[b][col])<1.e-12,"constant shared kernel GL3/GL5 equivalence",tests);
    delete common.h;delete d.h;
    std::ofstream file(output);if(!file)throw std::runtime_error("cannot create numerical report");file<<"{\"passed\":true,\"tests\":"<<tests<<",\"scope\":\"mathematical, seed-units and numerical closure; NOT isotope identification or production PID validation\"}\n";file.close();
    std::cout<<"ITERATIVE_HELIUM_SELFTEST_PASS tests="<<tests<<std::endl;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
