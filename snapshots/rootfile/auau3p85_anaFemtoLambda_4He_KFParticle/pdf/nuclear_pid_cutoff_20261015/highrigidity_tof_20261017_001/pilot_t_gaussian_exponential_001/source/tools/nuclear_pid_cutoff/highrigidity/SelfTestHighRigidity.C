void SelfTestHighRigidity();
#ifndef __CINT__
#include "FitHighRigidity.C"
namespace NuclearHighTest {
using namespace NuclearHighRigidity;
int checks=0;
void Check(bool condition,const char* message) {++checks;Require(condition,message);}
void Run() {
  Check(AllowedSample("TAG_D")&&AllowedSample("TAG_T")&&AllowedSample("TAG_HE3"),"Exact saved TAG keys rejected");
  Check(!AllowedSample("ALL")&&!AllowedSample("TAG_HE4")&&!AllowedSample(""),"Unsupported/fake high-R sample accepted");
  std::vector<Choice> choices(4);
  for(int s=0;s<4;++s) {choices[s].active=true;choices[s].mean=5+4*s;choices[s].sigma=1+.3*s;}
  for(int target=0;target<4;++target) {
    Regions regions=SelectionRegions(choices,target,2,0,25);
    for(int b=0;b<25000;++b) {
      double x=(b+.5)/1000.;bool selected=false;
      for(unsigned i=0;i<regions.size();++i) selected=selected||(x>regions[i].lo&&x<regions[i].hi);
      Check(selected==(Select(choices,x,2)==target),"High-R analytic gate disagrees with original TPC winner");
    }
  }
  Slice slice;slice.index=0;slice.rBin=201;slice.lo=2;slice.hi=2.01;slice.fitLo=0;slice.fitHi=25;
  const double weights[]={5./18.,8./18.,5./18.};
  for(int q=0;q<3;++q) for(int s=0;s<4;++s) {
    Row row;row.index=0;row.rBin=201;row.node=q;row.species=s;row.lo=2;row.hi=2.01;row.weight=weights[q];row.gate=choices[s];slice.rows.push_back(row);
  }
  TH1D raw("high_r_raw_gate_fixture","",250,0,25);raw.SetDirectory(0);raw.Sumw2();
  for(int b=1;b<=250;++b) {raw.SetBinContent(b,100);raw.SetBinError(b,10);}
  double count=0;
  for(int s=0;s<4;++s) {
    RawGate measured=Measure(slice,raw,s,2,1.e-9);
    Check(measured.count>0&&measured.lower<=measured.count+1.e-8&&measured.upper+1.e-8>=measured.count,"Raw count/full/touch bounds invalid");
    count+=measured.count;
  }
  Check(count>0&&count<=raw.Integral(),"Raw counts depend on no fit and cannot exceed total");
  choices[1]=choices[0];Check(Select(choices,choices[0].mean,2)==0,"Original tie ordering changed");
  choices[1].active=choices[2].active=choices[3].active=false;
  Check(Select(choices,choices[0].mean+2*choices[0].sigma,2)==-1,"Strict nSigma boundary changed");
}
}
#endif
void SelfTestHighRigidity() {
#ifndef __CINT__
  try {NuclearHighTest::Run();std::cout<<"NUCLEAR_HIGH_R_SELFTEST_PASS checks="<<NuclearHighTest::checks<<std::endl;}
  catch(const std::exception& e) {std::cerr<<"High-R selftest: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
