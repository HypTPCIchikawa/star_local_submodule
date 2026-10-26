#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include <TError.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
bool nfError=false;
void nfHandler(Int_t level,Bool_t abort,const char* location,const char* msg){if(level>=kError)nfError=true;DefaultErrorHandler(level,abort,location,msg);}
std::string nfPath(std::string s){std::string::size_type p;while((p=s.find(";1"))!=std::string::npos)s.erase(p,2);return s;}
void native_float_check(const char* inputList,const char* nameList,const char* target,const char* report){
 TH1::AddDirectory(kFALSE);SetErrorHandler(nfHandler);
 std::ifstream il(inputList),nl(nameList);std::vector<std::string> files,names;std::string line;
 while(std::getline(il,line))if(!line.empty())files.push_back(line);
 while(std::getline(nl,line))if(!line.empty())names.push_back(line);
 if(files.empty()||names.empty()){gSystem->Exit(71);return;}
 std::vector<TH1*> refs(names.size(),(TH1*)0);
 for(unsigned i=0;i<files.size();++i){
  TFile f(files[i].c_str(),"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered)){gSystem->Exit(72);return;}
  for(unsigned j=0;j<names.size();++j){
   TH1* h=dynamic_cast<TH1*>(f.Get(nfPath(names[j]).c_str()));if(!h){gSystem->Exit(73);return;}
   if(!refs[j]){refs[j]=(TH1*)h->Clone("native_reference");refs[j]->SetDirectory(0);}
   else if(!refs[j]->Add(h)){gSystem->Exit(74);return;}
   delete h;
  }
  f.Close();if(nfError){gSystem->Exit(75);return;}
 }
 TFile actual(target,"READ");if(actual.IsZombie()||actual.TestBit(TFile::kRecovered)){gSystem->Exit(76);return;}
 std::ofstream out(report);out<<std::setprecision(17);int mismatches=0;
 for(unsigned j=0;j<names.size();++j){
  TH1* h=dynamic_cast<TH1*>(actual.Get(nfPath(names[j]).c_str()));TH1* ref=refs[j];
  if(!h||std::string(h->ClassName())!=ref->ClassName()||h->GetEntries()!=ref->GetEntries()){gSystem->Exit(77);return;}
  int cells=ref->GetNbinsX()+2;if(ref->GetDimension()>1)cells*=ref->GetNbinsY()+2;if(ref->GetDimension()>2)cells*=ref->GetNbinsZ()+2;
  int bad=0;double maxDiff=0;
  for(int b=0;b<cells;++b){double a=ref->GetBinContent(b),v=h->GetBinContent(b);if(!TMath::Finite(a)||!TMath::Finite(v)||a!=v){++bad;if(TMath::Abs(a-v)>maxDiff)maxDiff=TMath::Abs(a-v);}}
  out<<"HIST\t"<<names[j]<<"\t"<<cells<<"\t"<<bad<<"\t"<<maxDiff<<"\t"<<h->GetEntries()<<"\n";mismatches+=bad;
  delete h;delete ref;
 }
 out<<"DONE\t"<<files.size()<<"\t"<<names.size()<<"\t"<<mismatches<<"\n";out.close();
 if(nfError||mismatches||!out.good()){gSystem->Exit(78);return;}gSystem->Exit(0);
}
