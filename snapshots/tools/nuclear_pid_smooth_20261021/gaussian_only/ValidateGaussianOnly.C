void ValidateGaussianOnly(const char*);
#ifndef __CINT__
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <stdexcept>
#endif
void ValidateGaussianOnly(const char* directory){
#ifndef __CINT__
try{
 std::ifstream list((std::string(directory)+"/validation_roots.tsv").c_str());std::string line;int count=0;std::getline(list,line);
 while(std::getline(list,line)){
  std::istringstream in(line);std::string path,oldpath;int first,last;std::getline(in,path,'\t');std::getline(in,oldpath,'\t');in>>first>>last;
  TFile file(path.c_str(),"READ"),prior(oldpath.c_str(),"READ");if(file.IsZombie()||prior.IsZombie())throw std::runtime_error("Missing root file");
  const char* parts[]={"d","t","he3","he4"};
  for(int index=first;index<=last;++index){
   TH1* total=dynamic_cast<TH1*>(file.Get(Form("prediction_total_%03d",index)));TH1* smooth=dynamic_cast<TH1*>(file.Get(Form("prediction_smooth_%03d",index)));TH1* raw=dynamic_cast<TH1*>(file.Get(Form("raw_native_%03d",index)));TH1* old=dynamic_cast<TH1*>(prior.Get(Form("raw_native_%03d",index)));
   if(!total||!smooth||!raw||!old)throw std::runtime_error("Missing spectra");
   TH1* components[4];for(int s=0;s<4;++s){components[s]=dynamic_cast<TH1*>(file.Get(Form("prediction_%s_%03d",parts[s],index)));if(!components[s])throw std::runtime_error("Missing nucleus");}
   for(int b=0;b<=total->GetNbinsX()+1;++b){
    if(smooth->GetBinContent(b)!=0)throw std::runtime_error("Nonzero smooth BG");
    double sum=0;for(int s=0;s<4;++s)sum+=components[s]->GetBinContent(b);
    if(std::abs(sum-total->GetBinContent(b))>1e-9*std::max(1.,std::abs(sum)))throw std::runtime_error("Component sum mismatch");
    if(raw->GetBinContent(b)!=old->GetBinContent(b)||raw->GetBinError(b)!=old->GetBinError(b))throw std::runtime_error("Raw histogram changed");
   }++count;
  }
 }
 std::cout<<"GAUSSIAN_ONLY_ROOT_PASS rows="<<count<<" zero smooth; exact raw; total=nuclei"<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
