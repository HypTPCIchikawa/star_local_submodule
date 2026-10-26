#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
void ValidateAllProduct(const char* cachePath,const char* productPath,const char* sample,
                        int firstR,int rows,double axisTolerance,double sumTolerance){
  try{
    TFile cache(cachePath,"READ"),product(productPath,"READ");
    if(cache.IsZombie()||product.IsZombie()||product.TestBit(TFile::kRecovered))
      throw std::runtime_error("unreadable/recovered ALL product");
    TH2* input=dynamic_cast<TH2*>(cache.Get(sample));if(!input)throw std::runtime_error("missing exact sample");
    const char* parts[]={"smooth","d","t","he3","he4"};
    for(int index=0;index<rows;++index){
      TH1* raw=dynamic_cast<TH1*>(product.Get(Form("raw_native_%03d",index)));
      TH1* total=dynamic_cast<TH1*>(product.Get(Form("prediction_total_%03d",index)));
      TH1* mask=dynamic_cast<TH1*>(product.Get(Form("prediction_fit_mask_%03d",index)));
      if(!raw||!total||!mask||raw->GetNbinsX()!=input->GetNbinsY())throw std::runtime_error("incomplete native/model product");
      TH1* component[5];for(int j=0;j<5;++j){component[j]=dynamic_cast<TH1*>(product.Get(Form("prediction_%s_%03d",parts[j],index)));if(!component[j])throw std::runtime_error("missing model component");}
      for(int b=0;b<=raw->GetNbinsX()+1;++b){
        double original=input->GetBinContent(firstR+index,b),saved=raw->GetBinContent(b);
        double error=input->GetBinError(firstR+index,b),savedError=raw->GetBinError(b);
        if(original!=saved||std::fabs(error-savedError)>sumTolerance*std::max(1.,std::fabs(error)))
          throw std::runtime_error("native raw count/error/flow changed");
        if(b>=1&&b<=raw->GetNbinsX()&&std::fabs(raw->GetBinLowEdge(b)-input->GetYaxis()->GetBinLowEdge(b))>axisTolerance)
          throw std::runtime_error("native dEdx axis changed");
        double sum=0;for(int j=0;j<5;++j)sum+=component[j]->GetBinContent(b);double expected=total->GetBinContent(b);
        if(!TMath::Finite(sum)||!TMath::Finite(expected)||std::fabs(sum-expected)>sumTolerance*std::max(1.,std::fabs(expected)))
          throw std::runtime_error("model component sum changed");
        if(mask->GetBinContent(b)!=0&&mask->GetBinContent(b)!=1)throw std::runtime_error("invalid fit mask");
      }
    }
    product.Close();cache.Close();std::cout<<"NUCLEAR_ALL_READBACK_PASS rows="<<rows<<std::endl;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
}
