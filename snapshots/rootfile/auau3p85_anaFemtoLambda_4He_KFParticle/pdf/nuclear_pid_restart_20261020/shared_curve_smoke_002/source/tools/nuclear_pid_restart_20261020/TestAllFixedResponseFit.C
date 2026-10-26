void TestAllFixedResponseFit(const char*);
#ifndef __CINT__
#include "AllFixedResponseFit.h"
#include <iostream>
namespace RestartAllTest {
unsigned checks=0;
void Check(bool ok,const std::string& msg){++checks;if(!ok)throw std::runtime_error("ALL fixed-response test: "+msg);}
bool In(double x,const RestartAllFixed::Intervals& g){for(unsigned i=0;i<g.size();++i)if((x>g[i].lo||(x==g[i].lo&&g[i].includeLo))&&(x<g[i].hi||(x==g[i].hi&&g[i].includeHi)))return true;return false;}
void Gates(){
  NuclearFit::Vec m(4),s(4,1.);m[0]=10.;m[1]=12.;m[2]=50.;m[3]=60.;
  RestartAllFixed::Intervals g=RestartAllFixed::NodeGate(m,s,1,2.,true);
  Check(!In(11.,g)&&In(std::nextafter(11.,12.),g),"exact tie belongs to earlier species");Check(!In(14.,g),"strict positive endpoint");
  RestartAllFixed::Intervals ns=RestartAllFixed::NodeGate(m,s,1,2.,false);
  Check(!In(10.,ns)&&In(std::nextafter(10.,12.),ns)&&!In(14.,ns),"strict endpoints independent of winner");
  for(unsigned i=0;i<=2000;++i){float x=float(10.+4.*i/2000.);double q=RestartFloat16::QuantizeDefault12(x);bool got=false;
    for(unsigned k=0;k<g.size();++k){RestartAllFixed::Interval z=RestartAllFixed::EffectiveInterval(g[k],"root_default12");if(double(x)>=z.lo&&double(x)<z.hi)got=true;}
    Check(got==RestartAllFixed::Selected(q,m,s,1,2.,true),"Float16 strict winner oracle");}
  m[0]=m[1]=10.;s[1]=2.;g=RestartAllFixed::NodeGate(m,s,1,2.,true);
  Check(In(9.,g)&&In(11.,g)&&!In(10.,g)&&g.size()==2,"unequal-width punctured gate");m[0]=11.;
  g=RestartAllFixed::NodeGate(m,s,1,2.,true);
  for(unsigned i=0;i<=2000;++i){double x=6.+8.*i/2000.;Check(In(x,g)==RestartAllFixed::Selected(x,m,s,1,2.,true),"two-root gate oracle");}
  m[0]=m[1]=10.;s[0]=s[1]=1.;Check(RestartAllFixed::NodeGate(m,s,1,2.,true).empty(),"identical later curve never wins");
  Check(In(10.,RestartAllFixed::NodeGate(m,s,0,2.,true)),"identical earliest curve wins");
  RestartAllFixed::Interval a=RestartAllFixed::EffectiveInterval(RestartAllFixed::Interval(11.,12.,false,false),"root_default12"),b=RestartAllFixed::EffectiveInterval(RestartAllFixed::Interval(11.,12.,true,true),"root_default12");
  Check(a.lo>b.lo&&a.hi<b.hi,"quantized endpoint mass preserved");
}
std::vector<RestartAllFixed::Shape> Shapes(){std::vector<RestartAllFixed::Shape> a(4);const char* names[]={"d","t","he3","he4"};double w[]={5./18.,8./18.,5./18.},shift[]={-.3,0.,.3};
  for(unsigned j=0;j<4;++j){RestartAllFixed::Shape& s=a[j];s.supportProvided=true;s.component.name=names[j];s.component.mean=10.+5.*j;s.component.sigma=2.;s.component.active=true;s.component.freeMean=false;s.component.freeSigma=false;
    for(unsigned q=0;q<3;++q){s.component.nodes.push_back(NuclearFit::Node(w[q],s.component.mean+shift[q],2.+.03*(int(q)-1)));s.ratios.push_back(.8+.2*j+.02*(int(q)-1));}}return a;}
NuclearFit::Request Toy(const std::vector<RestartAllFixed::Shape>& a,bool mapped){NuclearFit::Request r;r.backgroundModel="linear";r.freeBackgroundSlope=false;r.backgroundWindowYield=90000.;r.nativeAxisMin=0.;r.nativeAxisMax=40.;
  for(unsigned j=0;j<4;++j){NuclearFit::Component c=a[j].component;c.startWindowYield=30000.+j*4000.;r.components.push_back(c);}
  for(unsigned b=0;b<300;++b){double lo=4.+.1*b,hi=4.+.1*(b+1);r.bins.push_back(NuclearFit::Bin(mapped?RestartFloat16::ContinuumEdge(lo):lo,mapped?RestartFloat16::ContinuumEdge(hi):hi,0.));}return r;}
void Synthetic(NuclearFit::Settings cfg){std::vector<RestartAllFixed::Shape> a=Shapes();
  for(unsigned mapped=0;mapped<2;++mapped){cfg["storageMapping"]=mapped?"root_default12":"continuous";NuclearFit::Request r=Toy(a,mapped);NuclearFit::Matrix t=RestartAllFixed::Templates(r,a,cfg);
    NuclearFit::Vec truth(6);truth[0]=40000.;truth[1]=50000.;truth[2]=60000.;truth[3]=70000.;truth[4]=45000.;truth[5]=55000.;
    for(unsigned b=0;b<t.size();++b)for(unsigned j=0;j<truth.size();++j)r.bins[b].count+=truth[j]*t[b][j];
    RestartAllFixed::Result f=RestartAllFixed::Fit(r,a,3,cfg);Check(f.fit.finite&&f.fit.amplitudeConverged,"Asimov convergence");Check(f.fit.dataDeviance<1.e-4,"Asimov deviance");
    for(unsigned j=0;j<truth.size();++j)Check(std::fabs(f.fit.columnWindowYields[j]-truth[j])/truth[j]<1.e-4,"Known amplitude recovery");
    Check(f.fit.yieldCovarianceValid&&f.winner.covarianceValid&&f.winner.error>0.,"positive conditional covariance/error");Check(f.winner.insideFit&&f.winner.completeSupport,"gate inside fit");
    double total=0.;for(unsigned j=0;j<truth.size();++j)total+=truth[j]*f.winner.fractions[j];Check(std::fabs(f.winner.purity-truth[3]*f.winner.fractions[3]/total)<1.e-5,"known purity recovery");
    Check(f.fit.conditionalDropDeltaD[3]>0.,"target identifiable");Check(f.winner.nominal.size()==3&&f.winner.nominal[0][0].lo!=f.winner.nominal[2][0].lo,"per-node gates retained");
    std::vector<RestartAllFixed::Shape> missing=a;missing[0].supportProvided=false;NuclearFit::Request mr=r;mr.components[0].active=false;
    RestartAllFixed::Result mf=RestartAllFixed::Fit(mr,missing,3,cfg);Check(!mf.purityAvailable&&!mf.winner.completeSupport,"missing winner curve unassessed");Check(mf.nsigma.purity>=0.,"nsigma-only diagnostic separate");
  }
  for(unsigned j=0;j<4;++j)for(unsigned q=0;q<3;++q)a[j].ratios[q]=1.;NuclearFit::Request r=Toy(a,false);
  Check(RestartAllFixed::Templates(r,a,cfg)==NuclearFit::BuildTemplates(r,cfg),"Gaussian limit exact");
  a[1].component.nodes=a[0].component.nodes;a[1].component.mean=a[0].component.mean;a[1].component.sigma=a[0].component.sigma;r=Toy(a,false);NuclearFit::Matrix t=RestartAllFixed::Templates(r,a,cfg);
  for(unsigned b=0;b<t.size();++b)for(unsigned j=0;j<t[b].size();++j)r.bins[b].count+=50000.*t[b][j];
  RestartAllFixed::Result f=RestartAllFixed::Fit(r,a,3,cfg);Check(!f.fit.yieldCovarianceValid&&!f.purityAvailable,"duplicate kernels unassessed");
}
void Run(const char* settings){NuclearFit::Settings cfg=NuclearFit::ReadSettings(settings);RestartAllFixed::ValidateSettings(cfg);Gates();Synthetic(cfg);std::cout<<"ALL_FIXED_RESPONSE_TEST_COMPLETE checks="<<checks<<std::endl;}
}
#endif
void TestAllFixedResponseFit(const char* settings){
#ifndef __CINT__
  try{RestartAllTest::Run(settings);}catch(const std::exception& e){std::cerr<<"ALL_FIXED_RESPONSE_TEST_FATAL: "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
