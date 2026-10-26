void ValidateR061(const char*,int,int);
#ifndef __CINT__
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>
#endif
void ValidateR061(const char* directory,int first,int last){
#ifndef __CINT__
try{
const char* responses[2]={"gaussian","lognormal"};const char* parts[4]={"d","t","he3","he4"};
for(int r=0;r<2;++r){std::string base=std::string(directory)+"/";TFile none((base+"full_"+responses[r]+"_none/cutoff_fits.root").c_str(),"READ");if(none.IsZombie())throw std::runtime_error("No no-background ROOT");
for(int i=first;i<=last;++i){TH1* total=dynamic_cast<TH1*>(none.Get(Form("prediction_total_%03d",i)));TH1* smooth=dynamic_cast<TH1*>(none.Get(Form("prediction_smooth_%03d",i)));if(!total||!smooth)throw std::runtime_error("Missing prediction");for(int b=0;b<=total->GetNbinsX()+1;++b){if(smooth->GetBinContent(b)!=0)throw std::runtime_error("Nonzero smooth component");double sum=0;for(int s=0;s<4;++s){TH1* h=dynamic_cast<TH1*>(none.Get(Form("prediction_%s_%03d",parts[s],i)));if(!h)throw std::runtime_error("Missing nucleus");sum+=h->GetBinContent(b);}if(std::abs(sum-total->GetBinContent(b))>1e-9*std::max(1.,std::abs(sum)))throw std::runtime_error("Total differs from nuclei");}}
const char* variants[4]={"full_linear","full_exponential","warm_linear","warm_exponential"};
for(int v=0;v<4;++v){std::string key=variants[v];std::string stage=key.substr(0,key.find('_')),bg=key.substr(key.find('_')+1);TFile full((base+stage+"_"+responses[r]+"_"+bg+"/cutoff_fits.root").c_str(),"READ");if(full.IsZombie())throw std::runtime_error("No control ROOT");for(int i=first;i<=last;++i){TH1* a=dynamic_cast<TH1*>(none.Get(Form("raw_native_%03d",i)));TH1* b=dynamic_cast<TH1*>(full.Get(Form("raw_native_%03d",i)));if(!a||!b||a->GetNbinsX()!=b->GetNbinsX())throw std::runtime_error("Missing raw bins");for(int j=0;j<=a->GetNbinsX()+1;++j)if(a->GetBinContent(j)!=b->GetBinContent(j)||a->GetBinError(j)!=b->GetBinError(j))throw std::runtime_error("Raw spectrum changed");}}
}
std::cout<<"R061_ROOT_VALIDATION_PASS zero smooth; identical raw spectra in all controls"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
