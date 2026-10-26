// Input-only exact R grouping. No fit, PID, normalization, or dE/dx rebin.
void ExtractCoarseInputs(const char*,const char*,const char*);
#ifndef __CINT__
#include "TFile.h"
#include "TH2D.h"
#include "TArrayD.h"
#include "TList.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace NuclearCoarseInputs {
const char* keys[]={"ALL","TAG_D","TAG_T","TAG_HE3","REMAINING"};
typedef std::map<std::string,double> Settings;
unsigned checks=0;
void Require(bool condition,const std::string& reason){++checks;if(!condition)throw std::runtime_error(reason);}
bool Near(double a,double b,double tolerance){return std::fabs(a-b)<=tolerance*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
double Need(const Settings& cfg,const std::string& key){Settings::const_iterator it=cfg.find(key);Require(it!=cfg.end()&&std::isfinite(it->second),"Missing/nonfinite setting: "+key);return it->second;}
int Integer(const Settings& cfg,const std::string& key){double n=Need(cfg,key);Require(n==std::floor(n),"Noninteger setting: "+key);return int(n);}
Settings Read(const char* path){std::ifstream in(path);Require(bool(in),"Missing explicit grouping settings");std::string key;double value;Settings cfg;while(in>>key>>value)Require(cfg.insert(std::make_pair(key,value)).second,"Duplicate grouping setting");Require(in.eof(),"Malformed grouping settings");return cfg;}
std::vector<double> Edges(TAxis& axis){std::vector<double> out;for(int i=1;i<=axis.GetNbins()+1;++i)out.push_back(axis.GetBinLowEdge(i));return out;}
void EqualAxis(TAxis& a,TAxis& b){Require(a.GetNbins()==b.GetNbins(),"Axis size changed");for(int i=1;i<=a.GetNbins()+1;++i)Require(a.GetBinLowEdge(i)==b.GetBinLowEdge(i),"Axis edge changed");}
double Variance(TH2D& hist,int x,int y){Require(hist.GetSumw2N()==(hist.GetNbinsX()+2)*(hist.GetNbinsY()+2),"Missing explicit native Sumw2");return hist.GetSumw2()->At(hist.GetBin(x,y));}
TH2D* Group(TH2D& source,int factor,const char* name){
  int nx=source.GetNbinsX(),ny=source.GetNbinsY();Require(factor>0&&nx%factor==0,"R grouping must exactly divide the native axis");
  std::vector<double> xedges,yedges=Edges(*source.GetYaxis());for(int x=1;x<=nx+1;x+=factor)xedges.push_back(source.GetXaxis()->GetBinLowEdge(x));
  TH2D* output=new TH2D(name,source.GetTitle(),nx/factor,&xedges[0],ny,&yedges[0]);output->SetDirectory(0);output->Sumw2();
  for(int x=0;x<=nx+1;++x){int target=x==0?0:(x==nx+1?nx/factor+1:1+(x-1)/factor);
    for(int y=0;y<=ny+1;++y){double n=source.GetBinContent(x,y),v=Variance(source,x,y);Require(std::isfinite(n)&&n>=0&&std::isfinite(v)&&v>=0,"Invalid stored source count/variance");int bin=output->GetBin(target,y);output->SetBinContent(bin,output->GetBinContent(bin)+n);output->GetSumw2()->GetArray()[bin]+=v;}}
  output->SetEntries(source.GetEntries());output->GetXaxis()->SetTitle(source.GetXaxis()->GetTitle());output->GetYaxis()->SetTitle(source.GetYaxis()->GetTitle());return output;
}
void VerifyGrouped(TH2D& source,TH2D& target,int factor,double tolerance){
  Require(source.GetNbinsX()==factor*target.GetNbinsX(),"Wrong grouped R dimensions");EqualAxis(*source.GetYaxis(),*target.GetYaxis());Require(source.GetEntries()==target.GetEntries(),"Entries changed by grouping");
  for(int x=1;x<=target.GetNbinsX()+1;++x)Require(target.GetXaxis()->GetBinLowEdge(x)==source.GetXaxis()->GetBinLowEdge(1+(x-1)*factor),"Grouped R boundary not exactly inherited");
  // Independent target-centric summation, including all corner/edge flows.
  for(int x=0;x<=target.GetNbinsX()+1;++x){int first=x==0?0:(x==target.GetNbinsX()+1?source.GetNbinsX()+1:1+(x-1)*factor);int last=(x==0||x==target.GetNbinsX()+1)?first:x*factor;
    for(int y=0;y<=target.GetNbinsY()+1;++y){double n=0,v=0;for(int q=first;q<=last;++q){n+=source.GetBinContent(q,y);v+=Variance(source,q,y);}Require(target.GetBinContent(x,y)==n,"Grouped count differs from direct native sum");Require(Variance(target,x,y)==v,"Grouped stored variance differs from direct native sum");Require(Near(std::pow(target.GetBinError(x,y),2),v,tolerance),"Grouped error squared differs beyond roundoff");}}
}
void SelfTest(double tolerance){
  unsigned firstCheck=checks;TH2D toy("coarse_input_toy","toy",10,0.,1.,3,0.,3.);toy.SetDirectory(0);toy.Sumw2();
  for(int x=0;x<=11;++x)for(int y=0;y<=4;++y){int b=toy.GetBin(x,y);toy.SetBinContent(b,(x==4&&y==2)?0:100*x+7*y+.25);toy.GetSumw2()->GetArray()[b]=.0625*(b+1);}
  toy.SetEntries(1234.5);for(int factor=2;factor<=5;factor+=3){TH2D* grouped=Group(toy,factor,"grouped_toy");VerifyGrouped(toy,*grouped,factor,tolerance);Require(grouped->GetBinContent(0,0)==toy.GetBinContent(0,0),"Underflow corner moved");Require(grouped->GetBinContent(grouped->GetNbinsX()+1,4)==toy.GetBinContent(11,4),"Overflow corner moved");delete grouped;}
  bool rejected=false;try{TH2D* bad=Group(toy,3,"invalid_toy");delete bad;}catch(const std::exception&){rejected=true;}Require(rejected,"Incomplete R grouping was accepted");
  std::cout<<"COARSE_INPUT_SELFTEST_PASS checks="<<checks-firstCheck<<std::endl;
}
void ValidateSource(TH2D& source,const Settings& cfg){
  Require(source.GetNbinsX()==Integer(cfg,"inputXBins")&&source.GetNbinsY()==Integer(cfg,"inputYBins"),"Unexpected native cache dimensions");
  Require(source.GetXaxis()->GetXmin()==Need(cfg,"inputXMin")&&source.GetXaxis()->GetXmax()==Need(cfg,"inputXMax")&&source.GetYaxis()->GetXmin()==Need(cfg,"inputYMin")&&source.GetYaxis()->GetXmax()==Need(cfg,"inputYMax"),"Native cache endpoints differ");
  for(int x=1;x<=source.GetNbinsX();++x)Require(Near(source.GetXaxis()->GetBinWidth(x),Need(cfg,"inputXWidth"),Need(cfg,"axisTolerance")),"Native R width differs");
  for(int y=1;y<=source.GetNbinsY();++y)Require(Near(source.GetYaxis()->GetBinWidth(y),Need(cfg,"inputYWidth"),Need(cfg,"axisTolerance")),"Native dEdx width differs");
}
void VerifyWritten(TH2D* expected[5],const std::string& path){
  TFile readback(path.c_str(),"READ");Require(!readback.IsZombie()&&!readback.TestBit(TFile::kRecovered),"Closed coarse cache cannot be read");Require(readback.GetListOfKeys()->GetSize()==5,"Coarse cache must have exactly five keys");
  for(int s=0;s<5;++s){TH2D* actual=dynamic_cast<TH2D*>(readback.Get(keys[s]));Require(actual!=0,"Missing/wrong-class coarse histogram");EqualAxis(*actual->GetXaxis(),*expected[s]->GetXaxis());EqualAxis(*actual->GetYaxis(),*expected[s]->GetYaxis());Require(actual->GetEntries()==expected[s]->GetEntries(),"Closed cache Entries changed");
    for(int x=0;x<=actual->GetNbinsX()+1;++x)for(int y=0;y<=actual->GetNbinsY()+1;++y)Require(actual->GetBinContent(x,y)==expected[s]->GetBinContent(x,y)&&Variance(*actual,x,y)==Variance(*expected[s],x,y)&&actual->GetBinError(x,y)==expected[s]->GetBinError(x,y),"Closed count/error/variance/flow differs");}
  readback.Close();
}
void AxisJson(std::ofstream& out,TH2D& h){out<<"{\"rBins\":"<<h.GetNbinsX()<<",\"rRange\":["<<h.GetXaxis()->GetXmin()<<','<<h.GetXaxis()->GetXmax()<<"],\"rBinWidth\":"<<(h.GetXaxis()->GetXmax()-h.GetXaxis()->GetXmin())/h.GetNbinsX()<<",\"dedxBins\":"<<h.GetNbinsY()<<",\"dedxRange\":["<<h.GetYaxis()->GetXmin()<<','<<h.GetYaxis()->GetXmax()<<"],\"dedxBinWidth\":"<<(h.GetYaxis()->GetXmax()-h.GetYaxis()->GetXmin())/h.GetNbinsY()<<'}';}
void Audit(TH2D* source[5],TH2D* grouped[5],int factor,const char* output,double tolerance){
  std::ofstream out((std::string(output)+"/native_audit.json").c_str());Require(bool(out),"Cannot create grouping audit");out<<std::setprecision(17)<<"{\"schemaVersion\":1,\"processingStage\":\"exact_r_coarsening\",\"noFit\":true,\"productionEnabled\":false,\"coarseningFactor\":"<<factor<<",\"inputNativeAxis\":";AxisJson(out,*source[0]);out<<",\"coarseAxis\":";AxisJson(out,*grouped[0]);
  out<<",\"cellsIncludingFlowsPerHistogram\":"<<(grouped[0]->GetNbinsX()+2)*(grouped[0]->GetNbinsY()+2)<<",\"exactFiveBinSumPassed\":true,\"countErrorFlowPreservedBySummation\":true,\"storedSumw2ArrayExact\":true,\"entriesPreserved\":true,\"writtenCacheReadbackExact\":true,\"axisEdgesExactlyInherited\":true,\"nativeRInputUnchanged\":true,\"errorPolicy\":\"Sum original stored Sumw2 directly; reported error is sqrt(sum), without reconstructing variances from rounded errors\",\"tagSubsetPolicy\":\"ALL/TAG are correlated samples; REMAINING is summed independently, not a new independent likelihood\",\"histograms\":[";
  for(int s=0;s<5;++s){double n=0,v=0,regular=0,sourceN=0,sourceV=0;for(int x=0;x<=source[s]->GetNbinsX()+1;++x)for(int y=0;y<=source[s]->GetNbinsY()+1;++y){sourceN+=source[s]->GetBinContent(x,y);sourceV+=Variance(*source[s],x,y);}for(int x=0;x<=grouped[s]->GetNbinsX()+1;++x)for(int y=0;y<=grouped[s]->GetNbinsY()+1;++y){double value=grouped[s]->GetBinContent(x,y);n+=value;v+=Variance(*grouped[s],x,y);if(x>0&&x<=grouped[s]->GetNbinsX()&&y>0&&y<=grouped[s]->GetNbinsY())regular+=value;}Require(n==sourceN&&Near(v,sourceV,tolerance),"Global count/variance sum changed");if(s)out<<',';out<<"{\"cacheKey\":\""<<keys[s]<<"\",\"storedCount\":"<<n<<",\"varianceSum\":"<<v<<",\"sourceVarianceSum\":"<<sourceV<<",\"varianceSumRelativeRoundoff\":"<<std::fabs(v-sourceV)/std::max(1.,sourceV)<<",\"normalCount\":"<<regular<<",\"flowCount\":"<<n-regular<<",\"entries\":"<<grouped[s]->GetEntries()<<'}';}
  out<<"],\"checks\":"<<checks<<"}\n";out.close();Require(out.good(),"Cannot finish grouping audit");
}
}
#endif
void ExtractCoarseInputs(const char* input,const char* output,const char* settings){
#ifndef __CINT__
  using namespace NuclearCoarseInputs;
  try{Require(input&&output&&settings,"Explicit input/output/settings required");Settings cfg=Read(settings);Require(Integer(cfg,"schemaVersion")==1&&Integer(cfg,"coarseningFactor")==5,"Only declared exact five-native-bin grouping is authorized");Require(Need(cfg,"axisTolerance")>0&&Need(cfg,"varianceTolerance")>0,"Invalid grouping tolerance");SelfTest(Need(cfg,"varianceTolerance"));
    std::string path=std::string(output)+"/cache.root";Require(gSystem->AccessPathName(path.c_str())&&gSystem->AccessPathName((std::string(output)+"/native_audit.json").c_str()),"Existing coarse artifacts will not be overwritten");TFile native(input,"READ");Require(!native.IsZombie()&&!native.TestBit(TFile::kRecovered),"Native input cache is unreadable/recovered");Require(native.GetListOfKeys()->GetSize()==5,"Native cache must have five flat histograms");TH2D* source[5];TH2D* grouped[5];int factor=Integer(cfg,"coarseningFactor");
    for(int s=0;s<5;++s){source[s]=dynamic_cast<TH2D*>(native.Get(keys[s]));Require(source[s]!=0,"Native key missing or wrong class");ValidateSource(*source[s],cfg);if(s){EqualAxis(*source[0]->GetXaxis(),*source[s]->GetXaxis());EqualAxis(*source[0]->GetYaxis(),*source[s]->GetYaxis());}grouped[s]=Group(*source[s],factor,keys[s]);VerifyGrouped(*source[s],*grouped[s],factor,Need(cfg,"varianceTolerance"));}
    for(int x=0;x<=grouped[0]->GetNbinsX()+1;++x)for(int y=0;y<=grouped[0]->GetNbinsY()+1;++y)Require(grouped[0]->GetBinContent(x,y)==grouped[1]->GetBinContent(x,y)+grouped[2]->GetBinContent(x,y)+grouped[3]->GetBinContent(x,y)+grouped[4]->GetBinContent(x,y),"Grouped disjoint subset sum differs from ALL");
    TFile written(path.c_str(),"NEW");Require(!written.IsZombie(),"Cannot create new coarse cache");for(int s=0;s<5;++s)Require(grouped[s]->Write()>0,"Cannot write coarse histogram");written.Close();VerifyWritten(grouped,path);Audit(source,grouped,factor,output,Need(cfg,"varianceTolerance"));native.Close();for(int s=0;s<5;++s)delete grouped[s];std::cout<<"COARSE_INPUTS_COMPLETE groupedRWidth="<<Need(cfg,"inputXWidth")*factor<<"; no fit; checks="<<checks<<std::endl;
  }catch(const std::exception& error){std::cerr<<"COARSE_INPUTS_ERROR "<<error.what()<<std::endl;gSystem->Exit(2);}
#endif
}
