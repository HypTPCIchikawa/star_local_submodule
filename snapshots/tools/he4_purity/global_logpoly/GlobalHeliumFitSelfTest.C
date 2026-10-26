void GlobalHeliumFitSelfTest(const char* settings,const char* output);
#ifndef __CINT__
#include "GlobalHeliumFit.C"
#include "TH2D.h"
#include <limits>
namespace GlobalHeliumTest {
void Require(bool ok,const char* label,int& tests){++tests;if(!ok)throw std::runtime_error(std::string("selftest failed: ")+label);}
bool CountRejected(double n,double v,const GlobalHelium::Settings& cfg){try{GlobalHelium::InputAudit a;GlobalHelium::AuditCount(n,v,cfg,a);return false;}catch(const std::runtime_error&){return true;}}
bool AxesRejected(const TAxis& a,const TAxis& b,const GlobalHelium::Settings& cfg){try{GlobalHelium::AuditAxes(a,b,cfg);return false;}catch(const std::runtime_error&){return true;}}
}
#endif
void GlobalHeliumFitSelfTest(const char* settings,const char* output){
#ifndef __CINT__
  using namespace GlobalHelium;using namespace GlobalHeliumTest;
  try{
    if(!gSystem->AccessPathName(output))throw std::runtime_error("selftest output already exists");
    Settings cfg=ReadSettings(settings);Engine e(cfg);Vec g=Legacy();int tests=0;
    Require(e.Feasible(g),"legacy common coefficients feasible",tests);
    for(int sp=2;sp<4;++sp)for(int j=0;j<=10;++j){double r=Value(cfg,"rigidityMin")+(Value(cfg,"rigidityMax")-Value(cfg,"rigidityMin"))*j/10;Require(std::fabs(Mean(g,sp,r)-NuclearIdDeDxVsMom::GetMean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(sp),r))<1.e-10,"mean uses original shared polynomial",tests);Require(std::fabs(Sigma(g,sp,r)-NuclearIdDeDxVsMom::GetSigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(sp),r))<1.e-10,"sigma uses original shared polynomial",tests);}
    double q[]={1.,-2.,1.},c[]={0.,1.,0.,-1.};
    Require(std::fabs(Extreme(q,2,0,2,true))<1.e-12,"interior quadratic extremum",tests);
    Require(std::fabs(Extreme(c,3,-1,1,false)-2/(3*std::sqrt(3.)))<1.e-12,"interior cubic extremum",tests);
    Require(std::fabs(Gaussian(-12,12,0,1)-1)<1.e-12,"Gaussian full-line normalization",tests);
    Require(Gaussian(10,11,0,1)>0&&Gaussian(-11,-10,0,1)>0,"Gaussian nonzero remote tails",tests);
    Require(std::fabs(Gaussian(10,11,0,1)-Gaussian(-11,-10,0,1))<1.e-35,"Gaussian tail symmetry",tests);
    Slice s;s.lo=.40;s.hi=.42;s.r=.41;s.fitLo=0;s.fitHi=120;
    for(int j=0;j<240;++j){Bin b;b.lo=.5*j;b.hi=b.lo+.5;b.n=0;s.bins.push_back(b);}
    for(int order=1;order<=5;order+=2){Settings branch=cfg;std::ostringstream k;k<<order;branch["rQuadratureOrder"]=k.str();Engine quad(branch);quad.Nodes(s);double w=0;for(unsigned j=0;j<s.weights.size();++j)w+=s.weights[j];Require(std::fabs(w-1)<1.e-14,"R quadrature normalized average",tests);Matrix t=quad.Templates(s,g,s.shape);double bg1=0,bg2=0,he3=0;for(unsigned b=0;b<t.size();++b){bg1+=t[b][4];if(!quad.exponential)bg2+=t[b][5];he3+=t[b][2];}Require(std::fabs(bg1-1)<1.e-12,"background window normalization",tests);if(!quad.exponential)Require(std::fabs(bg2-1)<1.e-12,"second linear basis normalization",tests);Require(he3>0&&he3<=1.00000000001,"Gaussian finite-window acceptance",tests);}
    e.Nodes(s);
    Vec bad=g;bad[4]=-20;Require(!e.Feasible(bad),"reject negative width function",tests);
    bad=g;bad[7]=bad[0]-10;bad[8]=bad[1];bad[9]=bad[2];bad[10]=bad[3];Require(!e.Feasible(bad),"reject reversed helium means",tests);
    bad=g;bad[0]=std::numeric_limits<double>::quiet_NaN();Require(!e.Feasible(bad),"reject nonfinite coefficients",tests);
    Vec shape=s.shape;shape[3]=0;Require(!e.ShapeFeasible(s,g,shape),"reject zero neighbor width",tests);
    shape=s.shape;shape[2]=std::numeric_limits<double>::quiet_NaN();Require(!e.ShapeFeasible(s,g,shape),"reject nonfinite local shape",tests);
    InputAudit audit;AuditCount(17,17,cfg,audit);Require(audit.cells==1,"accept unit counts",tests);
    Require(CountRejected(1.25,1.25,cfg),"reject fractional counts",tests);
    Require(CountRejected(4,8,cfg),"reject weighted count variance",tests);
    Require(CountRejected(-1,1,cfg),"reject negative count",tests);
    TH2D all("test_all","",2,0,2,2,0,2),tag("test_tag","",2,0,2,2,0,2),shift("test_shift","",2,.01,2.01,2,0,2);
    all.SetBinContent(1,1,8);tag.SetBinContent(1,1,3);
    InputAudit cellAudit;AuditSourceCells(all,&tag,1,2,cfg,cellAudit);
    Require(cellAudit.cells==24,"audit native cells includes Y flows and complement",tests);
    AuditAxes(*all.GetXaxis(),*tag.GetXaxis(),cfg);
    Require(AxesRejected(*all.GetXaxis(),*shift.GetXaxis(),cfg),"reject axis shift",tests);
    tag.SetBinContent(1,1,9);bool reject=false;try{AuditSourceCells(all,&tag,1,2,cfg,cellAudit);}catch(const std::runtime_error&){reject=true;}
    Require(reject,"reject non-subset TAG",tests);
    // Known convex Poisson mixtures: zero start, exact zero-bound optimum,
    // nearly collinear templates, and high-statistics close-to-saturation.
    Slice toy;toy.count=10000;Bin first;first.lo=0;first.hi=1;first.n=7000;Bin second;second.lo=1;second.hi=2;second.n=3000;toy.bins.push_back(first);toy.bins.push_back(second);
    Matrix identity(2,Vec(2,0));identity[0][0]=identity[1][1]=1;
    Profile p=e.Amplitudes(toy,identity,Vec(2,0));
    Require(p.converged,"AMP zero-start KKT convergence",tests);Require(std::fabs(p.amplitude[0]-.7)<1.e-7&&std::fabs(p.amplitude[1]-.3)<1.e-7,"AMP known interior optimum",tests);
    toy.bins[0].n=10000;toy.bins[1].n=0;p=e.Amplitudes(toy,identity,Vec(2,.5));
    Require(p.converged&&p.amplitude[1]<1.e-12&&std::fabs(p.amplitude[0]-1)<1.e-7,"AMP exact active zero boundary",tests);
    Matrix close(2,Vec(2));close[0][0]=.51;close[1][0]=.49;close[0][1]=.49;close[1][1]=.51;
    toy.bins[0].n=5100;toy.bins[1].n=4900;p=e.Amplitudes(toy,close,Vec(2,.5));
    Require(p.converged&&p.kktResidual<Value(cfg,"amplitudeGradientTolerance"),"AMP near-collinear KKT",tests);
    toy.count=1.e12;toy.bins[0].n=7.e11;toy.bins[1].n=3.e11;p=e.Amplitudes(toy,identity,Vec(2,.5));
    Require(p.converged&&p.deviance<1.e-7,"AMP high-count stable deviance",tests);
    Settings limit=cfg;limit["amplitudeMaxIterations"]="1";Engine limited(limit);p=limited.Amplitudes(toy,identity,Vec(2,.01));
    Require(!p.converged,"AMP iteration limit not false convergence",tests);
    std::ofstream result(output);if(!result)throw std::runtime_error("cannot create selftest result");result<<"{\"passed\":true,\"tests\":"<<tests<<",\"scope\":\"numeric and input guards; NOT empirical isotope identification or PID validation\"}\n";result.close();
    std::cout<<"GLOBAL_HELIUM_SELFTEST_PASS tests="<<tests<<std::endl;
  }catch(const std::exception& ex){std::cerr<<ex.what()<<std::endl;gSystem->Exit(2);}
#endif
}
