// Verify saved native-to-fit rebin sums and variances without reopening input TH2.
void CheckRebinAudit(const char* fitRoot,const char* settings,const char* report);
#ifndef __CINT__
#include "../../StMaker/common/NuclearIdDeDxVsMom.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TString.h"
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace {
void RebinRequire(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
double RebinSetting(const std::map<std::string,std::string>& settings,const char* key){
  std::map<std::string,std::string>::const_iterator i=settings.find(key);
  RebinRequire(i!=settings.end(),"missing explicit rebin audit setting");
  double value=0.;std::istringstream input(i->second);input>>value;
  RebinRequire(bool(input)&&std::isfinite(value),"invalid numeric rebin audit setting");return value;
}
void RebinCell(const TH1D* native,const TH1D* coarse,int bin,int first,int last){
  double sum=0.,variance=0.;
  for(int n=first;n<=last;++n){
    const double content=native->GetBinContent(n),error=native->GetBinError(n);
    RebinRequire(std::isfinite(content)&&std::isfinite(error)&&content>=0.&&error>=0.,"finite nonnegative native bin");
    sum+=content;variance+=error*error;
  }
  RebinRequire(coarse->GetBinContent(bin)==sum,"rebinned bin is exact native count sum");
  const double error=coarse->GetBinError(bin);
  RebinRequire(std::fabs(error*error-variance)<=1.e-10*std::max(1.,variance),"rebinned squared error is native variance sum");
}
}
#endif
void CheckRebinAudit(const char* fitRoot,const char* settings,const char* report){
#ifndef __CINT__
  try{
    std::ifstream configuration(settings);RebinRequire(configuration.good(),"settings TSV readable");
    std::map<std::string,std::string> config;std::string line;
    while(std::getline(configuration,line)){const size_t tab=line.find('\t');RebinRequire(tab!=std::string::npos,"settings TSV record");config[line.substr(0,tab)]=line.substr(tab+1);}
    const double minR=RebinSetting(config,"rigidityMin"),maxR=RebinSetting(config,"rigidityMax"),step=RebinSetting(config,"rigidityStep");
    const double factorValue=RebinSetting(config,"dedxRebinFactor"),gate=RebinSetting(config,"gateNSigma");
    const double fitMin=RebinSetting(config,"fitNSigmaMin"),fitMax=RebinSetting(config,"fitNSigmaMax");
    const int factor=int(factorValue+.5),slices=int((maxR-minR)/step+.5);
    RebinRequire(step>0.&&slices>0&&factor>1&&std::fabs(factor-factorValue)<1.e-9,"valid explicit rebin audit range/factor");
    RebinRequire(gSystem->AccessPathName(report),"refusing earlier rebin audit output");
    TFile input(fitRoot,"READ");RebinRequire(!input.IsZombie()&&!input.TestBit(TFile::kRecovered),"valid saved fit ROOT");
    std::ofstream output(report);RebinRequire(output.good(),"new rebin audit report writable");
    output<<std::setprecision(17)<<"index\trigidity\tnativeCountAll\tfitCountAll\tnativeGateCount\tnativeFitRangeCount\n";
    for(int index=0;index<slices;++index){
      TH1D* native=dynamic_cast<TH1D*>(input.Get(Form("raw_native_%03d",index)));
      TH1D* coarse=dynamic_cast<TH1D*>(input.Get(Form("raw_%03d",index)));
      RebinRequire(native&&coarse,"native and coarse raw projections exist for every rigidity slice");
      RebinRequire(std::string(native->GetName())==Form("raw_native_%03d",index) &&
                   std::string(coarse->GetName())==Form("raw_%03d",index),
                   "ROOT key and internal histogram names match without native/coarse cache collisions");
      RebinRequire(native!=coarse,"native/coarse projections are distinct cached ROOT objects");
      RebinRequire(native->GetNbinsX()==factor*coarse->GetNbinsX(),"native/coarse histogram bin-count contract");
      RebinRequire(native->GetEntries()==coarse->GetEntries(),"rebinning preserves histogram entries");
      RebinCell(native,coarse,0,0,0);
      RebinCell(native,coarse,coarse->GetNbinsX()+1,native->GetNbinsX()+1,native->GetNbinsX()+1);
      for(int b=1;b<=coarse->GetNbinsX();++b){
        const int first=1+(b-1)*factor,last=b*factor;
        RebinRequire(coarse->GetBinLowEdge(b)==native->GetBinLowEdge(first),"rebinned low edge matches original edge");
        RebinCell(native,coarse,b,first,last);
      }
      const double lo=minR+index*step,hi=std::min(maxR,lo+step),r=(lo+hi)/2.;
      const double mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);
      RebinRequire(sigma>0.,"positive calibration sigma");
      double gateCount=0.,fitCount=0.;
      for(int b=1;b<=native->GetNbinsX();++b){const double u=(native->GetBinCenter(b)-mean)/sigma,count=native->GetBinContent(b);if(u>=-gate&&u<gate)gateCount+=count;if(u>=fitMin&&u<fitMax)fitCount+=count;}
      const double allNative=native->Integral(0,native->GetNbinsX()+1),allFit=coarse->Integral(0,coarse->GetNbinsX()+1);
      RebinRequire(allNative==allFit,"counts including under/overflow conserved exactly");
      output<<index<<'\t'<<r<<'\t'<<allNative<<'\t'<<allFit<<'\t'<<gateCount<<'\t'<<fitCount<<'\n';
    }
    output.close();RebinRequire(output.good(),"rebin report completed");
    std::cout<<"PASS rebin audit: "<<slices<<" projections; every bin count/variance and under/overflow conserved; native gate counts exported"<<std::endl;
  }catch(const std::exception& error){std::cerr<<"Rebin audit: "<<error.what()<<std::endl;gSystem->Exit(1);}
#endif
}
