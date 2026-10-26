void SelfTestAllNuclearPurity(const char*);
#ifndef __CINT__
#include "FitAllNuclearPurity.C"
namespace NuclearAllTest {
using namespace NuclearAll;
unsigned checks=0;
void Check(bool condition,const char* name){++checks;Require(condition,std::string("ALL purity selftest: ")+name);}
bool Contains(const Regions& regions,double x){for(unsigned j=0;j<regions.size();++j)if(x>regions[j].lo&&x<regions[j].hi)return true;return false;}
Choice C(double mean,double sigma,bool active=true){Choice c;c.mean=mean;c.sigma=sigma;c.active=active;return c;}
void TestSelectors(){
  std::vector<Choice> choices;choices.push_back(C(15,2));choices.push_back(C(26,3));choices.push_back(C(35,4));choices.push_back(C(53,5));
  std::vector<Regions> regions(4);for(int s=0;s<4;++s)regions[s]=SelectionRegions(choices,s,2,0,120);
  for(int j=0;j<3000;++j){double x=(j+.317)*120/3000;int selected=Select(choices,x,2),matches=0;
    for(int s=0;s<4;++s){bool inside=Contains(regions[s],x);matches+=inside;Check(inside==(selected==s),"partition matches original best-pull selector");}Check(matches<=1,"selector intervals are disjoint");}
  choices[1]=C(15,2);Check(Select(choices,15,2)==0,"identical-pull tie chooses first species");Check(SelectionRegions(choices,1,2,0,120).empty(),"identical later species has no positive-width selection");
  choices[1]=C(15,3);Check(Select(choices,15,2)==0,"equal-center tie is stable");Check(Contains(SelectionRegions(choices,1,2,0,120),16),"larger sigma wins away from equal center");
  choices[0].active=false;Check(Select(choices,15,2)==1,"inactive current hypothesis skipped");
  std::vector<Choice> equal(4);equal[0]=C(20,2);equal[1]=C(26,2);equal[2]=C(0,0,false);equal[3]=C(0,0,false);
  Check(Contains(SelectionRegions(equal,0,2,0,120),22.9),"equal-width lower-mean region");Check(!Contains(SelectionRegions(equal,0,2,0,120),23.1),"equal-width crossing partitions exactly");
  Check(Select(equal,16,2)==-1,"strict nSigma boundary rejected");
}
void TestIntegrals(){
  Request request;request.bins.push_back(Bin(0,120,1890));request.backgroundModel="linear";request.backgroundSlope=0;
  Result result;const double means[]={15,26,35,53},sigmas[]={2,3,4,5},yields[]={300,250,100,40,700,500};Vec weights;weights.push_back(5./18);weights.push_back(4./9);weights.push_back(5./18);
  for(int s=0;s<4;++s){Component c;c.name=Species[s];c.mean=means[s];c.sigma=sigmas[s];c.freeMean=c.freeSigma=false;
    for(int q=0;q<3;++q)c.nodes.push_back(Node(weights[q],means[s]+(q-1)*.2,sigmas[s]));result.components.push_back(c);result.windowAcceptance.push_back(KernelIntegral(c,0,120));}
  result.columnWindowYields.assign(yields,yields+6);result.backgroundSlope=0;result.yieldCovarianceValid=true;result.yieldCovarianceColumnActive.assign(6,true);result.fixedShapeYieldCovariance.assign(6,Vec(6,0));for(int j=0;j<6;++j)result.fixedShapeYieldCovariance[j][j]=yields[j];
  std::vector<Regions> selected(3);for(int q=0;q<3;++q){std::vector<Choice> c(4);for(int s=0;s<4;++s)c[s]=C(means[s]+(q-1)*.2,sigmas[s]);selected[q]=SelectionRegions(c,2,2,0,120);}
  Estimate exact=IntegrateSelection(result,request,selected,weights,2,0,1.e-9);Check(exact.available&&exact.covered&&exact.errorAvailable,"finite covered purity/error");Check(exact.purity>0&&exact.purity<1,"background-inclusive purity is inside range");
  double numericSignal=0,numericTotal=0;const int steps=240000;const double dx=120./steps;
  for(int b=0;b<steps;++b){double x=(b+.5)*dx;for(int q=0;q<3;++q)if(Contains(selected[q],x)){double total=0,signal=0;
      for(int s=0;s<4;++s){const Node& n=result.components[s].nodes[q];double z=(x-n.mean)/n.sigma,density=std::exp(-.5*z*z)/(std::sqrt(2*TMath::Pi())*n.sigma);double count=yields[s]*density/result.windowAcceptance[s];total+=count;if(s==2)signal=count;}
      total+=yields[4]*2*(120-x)/(120*120)+yields[5]*2*x/(120*120);numericSignal+=weights[q]*signal*dx;numericTotal+=weights[q]*total*dx;}}
  Check(std::fabs(exact.signal-numericSignal)<.01,"GL Gaussian signal equals dense numerical integral");Check(std::fabs(exact.total-numericTotal)<.02,"positive BG and Gaussian denominator equals dense integral");
  Check(std::fabs(exact.purity-numericSignal/numericTotal)<1.e-5,"exact best-species purity equals numerical integral");
  Vec gradient(6,0);for(int j=0;j<6;++j){Result up=result,down=result;up.columnWindowYields[j]+=.001;down.columnWindowYields[j]-=.001;gradient[j]=(IntegrateSelection(up,request,selected,weights,2,0,1.e-9).purity-IntegrateSelection(down,request,selected,weights,2,0,1.e-9).purity)/.002;}
  double variance=0;for(int j=0;j<6;++j)variance+=gradient[j]*gradient[j]*yields[j];Check(std::fabs(exact.error-std::sqrt(variance))<1.e-10,"conditional covariance error matches finite-difference gradient");
  for(int s=0;s<2;++s)Check(std::fabs(BackgroundIntegral(0,120,0,120,"linear",0,s)-1)<1.e-15,"positive pol1 basis normalizes exactly");
  const double slopes[]={-20,-.001,0,.001,20};for(unsigned j=0;j<5;++j)Check(std::fabs(BackgroundIntegral(0,120,0,120,"exponential",slopes[j],0)-1)<1.e-12,"exp basis normalizes at negative/zero/positive slope");
  request.bins[0].lo=32;request.bins[0].hi=40;Check(!IntegrateSelection(result,request,selected,weights,2,0,1.e-9).covered,"fit-window leakage marked uncovered");
  Check(GaussianIntegral(0,120,500,.05)==0,"distant Gaussian tail underflow is not fake signal");
  result.yieldCovarianceColumnActive[2]=false;Check(!IntegrateSelection(result,request,selected,weights,2,0,1.e-9).errorAvailable,"zero-boundary target does not get fake Gaussian error");
}
}
#endif
void SelfTestAllNuclearPurity(const char* outputJson){
#ifndef __CINT__
  using namespace NuclearAllTest;
  try{Require(outputJson&&*outputJson&&gSystem->AccessPathName(outputJson),"Selftest result will not be overwritten");TestSelectors();TestIntegrals();std::ofstream out(outputJson);Require(bool(out),"Cannot write new selftest result");out<<"{\"allChecksPassed\":true,\"checks\":"<<checks<<",\"physicsValidationComplete\":false}\n";out.close();std::cout<<"NUCLEAR_ALL_PURITY_SELFTEST_PASS checks="<<checks<<std::endl;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
