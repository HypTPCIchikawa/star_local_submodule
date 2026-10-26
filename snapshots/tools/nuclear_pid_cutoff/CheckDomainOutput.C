#include "TFile.h"
#include "TNamed.h"
#include "TH2.h"
#include "TError.h"
#include "TSystem.h"
#include <iostream>
#include <string>

void CheckDomainOutput(const char* input,const char* expectedConfigName) {
  if(!input||!*input||!expectedConfigName||!*expectedConfigName){gSystem->Exit(1);return;}
  TFile file(input,"READ");
  if(file.IsZombie()||file.TestBit(TFile::kRecovered)){gSystem->Exit(1);return;}
  TNamed* domain=dynamic_cast<TNamed*>(file.Get("FemtoNuclearPidDomain"));
  TH2* all=dynamic_cast<TH2*>(file.Get("true/hDedxP"));
  if(!domain||!all||!(all->GetEntries()>0)||
     std::string(domain->GetTitle()).find(expectedConfigName)==std::string::npos) {
    Error("CheckDomainOutput","Required domain provenance / raw QA absent");gSystem->Exit(1);return;
  }
  std::cout<<"NUCLEAR_DOMAIN_OUTPUT_READBACK_PASS rawQAEntries="<<all->GetEntries()<<"\n"<<domain->GetTitle()<<std::endl;
}
