void TestNoSmooth(const char*);
#ifndef __CINT__
#include "FitLocalNuclearCutoff.C"
#include <iostream>
template<class E> void CheckNoSmooth(const char* path){
 NuclearFit::Settings cfg=NuclearFit::ReadSettings(path);typename E::Request req;req.backgroundModel="none";req.backgroundWindowYield=0;req.freeBackgroundSlope=true;req.nativeAxisMin=0;req.nativeAxisMax=120;
 for(int k=0;k<4;++k){typename E::Component c;c.name=NuclearLocalCutoff::Species[k];c.mean=k==2?40:52;c.sigma=k==2?3:4;c.active=k>=2;c.freeMean=false;c.freeSigma=false;c.startWindowYield=k==2?8000:k==3?24000:0;c.nodes.push_back(typename E::Node(1,c.mean,c.sigma));req.components.push_back(c);}
 double acc2=E::Integral(20,80,40,3),acc3=E::Integral(20,80,52,4);
 for(int b=0;b<600;++b){double lo=20+.1*b,hi=20+.1*(b+1);double n=10000*E::Integral(lo,hi,40,3)/acc2+20000*E::Integral(lo,hi,52,4)/acc3;req.bins.push_back(typename E::Bin(lo,hi,n));}
 typename E::Result fit=E::Fit(req,cfg);
 NuclearFit::Require(fit.finite&&fit.amplitudeConverged&&fit.columnNames.size()==4&&fit.columnWindowYields.size()==4,"No-background fit must contain exactly four peak columns");
 NuclearFit::Require(fit.shapeNames.empty()&&fit.conditionalNdf==598,"No hidden background shape/yield degree of freedom");
 for(unsigned b=0;b<fit.smoothCounts.size();++b){NuclearFit::Require(fit.smoothCounts[b]==0,"Smooth prediction must be identically zero");double sum=0;for(int k=0;k<4;++k)sum+=fit.componentCounts[b][k];NuclearFit::Require(std::fabs(sum-fit.totalCounts[b])<1e-9*std::max(1.,sum),"Peak sum differs from total");}
 NuclearFit::Require(std::fabs(fit.columnWindowYields[2]/10000-1)<1e-4&&std::fabs(fit.columnWindowYields[3]/20000-1)<1e-4,"Failed known-yield recovery");
 std::cout<<"NO_SMOOTH_TOY_PASS D="<<fit.dataDeviance<<" ndf="<<fit.conditionalNdf<<std::endl;
}
#endif
void TestNoSmooth(const char* settings){
#ifndef __CINT__
 try{CheckNoSmooth<NuclearLocalCutoff::GaussianEngine>(settings);CheckNoSmooth<NuclearLocalCutoff::LognormalEngine>(settings);std::cout<<"NO_SMOOTH_SELFTEST_COMPLETE"<<std::endl;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
