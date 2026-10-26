void ValidateHe4(const char*);
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
void ValidateHe4(const char*directory){
#ifndef __CINT__
try{std::ifstream input((std::string(directory)+"/validation_roots.tsv").c_str());std::string line;std::getline(input,line);int count=0;
while(std::getline(input,line)){std::istringstream s(line);std::string region,path,oldpath;int first,last;std::getline(s,region,'\t');std::getline(s,path,'\t');std::getline(s,oldpath,'\t');s>>first>>last;TFile file(path.c_str(),"READ"),prior(oldpath.c_str(),"READ");if(file.IsZombie()||prior.IsZombie())throw std::runtime_error("Missing ROOT");
for(int i=first;i<=last;++i){const char* names[]={"smooth","d","t","he3","he4"};TH1* parts[5];for(int k=0;k<5;++k){parts[k]=dynamic_cast<TH1*>(file.Get(Form("prediction_%s_%03d",names[k],i)));if(!parts[k])throw std::runtime_error("Missing component");}TH1* total=dynamic_cast<TH1*>(file.Get(Form("prediction_total_%03d",i)));TH1* raw=dynamic_cast<TH1*>(file.Get(Form("raw_native_%03d",i)));TH1* old=dynamic_cast<TH1*>(prior.Get(Form("raw_native_%03d",i)));if(!total||!raw||!old)throw std::runtime_error("Missing spectrum");for(int b=0;b<=total->GetNbinsX()+1;++b){double sum=0;for(int k=0;k<5;++k)sum+=parts[k]->GetBinContent(b);if(std::abs(sum-total->GetBinContent(b))>1e-9*std::max(1.,std::abs(sum)))throw std::runtime_error("Sum mismatch");if((region=="low_none"||region=="low_fixed")&&parts[0]->GetBinContent(b)!=0)throw std::runtime_error("Nonzero smooth background");if(region=="high_linear")for(int k=1;k<=3;++k)if(parts[k]->GetBinContent(b)!=0)throw std::runtime_error("Nonzero forbidden nucleus");if(raw->GetBinContent(b)!=old->GetBinContent(b)||raw->GetBinError(b)!=old->GetBinError(b))throw std::runtime_error("Raw changed");}++count;}}
std::cout<<"HE4_ROOT_VALIDATION_PASS rows="<<count<<" raw unchanged; no smooth in middle; no d/t/He3 at high R"<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
