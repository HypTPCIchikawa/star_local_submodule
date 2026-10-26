void InspectPeaks(const char*);
#ifndef __CINT__
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <stdexcept>
#include "TSystem.h"
#endif
void InspectPeaks(const char* input){
#ifndef __CINT__
TFile f(input,"READ");int ids[6]={20,21,23,24,25,29};for(int k=0;k<6;++k){TH1* h=dynamic_cast<TH1*>(f.Get(Form("raw_native_%03d",ids[k])));if(!h)throw std::runtime_error("Missing raw");h->Rebin(5);for(int s=0;s<2;++s){h->GetXaxis()->SetRangeUser(s?59:43,s?74:56);int b=h->GetMaximumBin();std::cout<<"RAW_LOCAL_MAX index="<<ids[k]<<" nucleus="<<(s?"he4":"he3")<<" x="<<h->GetBinCenter(b)<<" count="<<h->GetBinContent(b)<<std::endl;}}
#endif
}
