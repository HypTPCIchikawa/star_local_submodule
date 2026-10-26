#include <TFile.h>
#include <TFileMerger.h>
#include <TKey.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TAxis.h>
#include <TObjString.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TMath.h>
#include <TError.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
bool splitError=false;
void splitHandler(Int_t level,Bool_t abort,const char* where,const char* msg){if(level>=kError)splitError=true;DefaultErrorHandler(level,abort,where,msg);}
bool sameAxis(TAxis* a,TAxis* b){
 if(a->GetNbins()!=b->GetNbins()||std::string(a->GetName())!=b->GetName()||std::string(a->GetTitle())!=b->GetTitle()||a->GetFirst()!=b->GetFirst()||a->GetLast()!=b->GetLast())return false;
 for(int i=1;i<=a->GetNbins()+1;++i)if(a->GetBinLowEdge(i)!=b->GetBinLowEdge(i))return false;
 for(int i=1;i<=a->GetNbins();++i)if(std::string(a->GetBinLabel(i))!=b->GetBinLabel(i))return false;
 return true;
}
void split_sparse(const char* mode,const char* input,const char* output,const char* sparseName,const char* reportOrManifest){
 TH1::AddDirectory(kFALSE);SetErrorHandler(splitHandler);
 if(!gSystem->AccessPathName(output)){gSystem->Exit(21);return;}
 if(std::string(mode)=="merge"){
  std::ifstream list(input);std::vector<std::string> files;std::string p;while(std::getline(list,p))if(!p.empty())files.push_back(p);
  if(files.empty()){gSystem->Exit(22);return;}
  {TFileMerger merger(kFALSE,kTRUE);merger.SetPrintLevel(1);if(!merger.OutputFile(output,"CREATE")){gSystem->Exit(23);return;}
   merger.AddObjectNames(sparseName);
   for(unsigned i=0;i<files.size();++i)if(!merger.AddFile(files[i].c_str())){gSystem->Exit(24);return;}
   if(!merger.PartialMerge(TFileMerger::kAll|TFileMerger::kSkipListed)){gSystem->Exit(25);return;}
  }
  if(splitError){gSystem->Exit(26);return;}
  TFile f(output,"UPDATE");if(f.IsZombie()||f.GetKey(sparseName)){gSystem->Exit(27);return;}
  std::ifstream in(reportOrManifest);std::ostringstream contents;contents<<in.rdbuf();
  TObjString metadata(contents.str().c_str());if(!in.good()||metadata.Write("NuclearPidPairJointSplitManifest")<=0){gSystem->Exit(28);return;}
  f.Close();if(splitError){gSystem->Exit(29);return;}gSystem->Exit(0);return;
 }
 if(std::string(mode)!="extract"){gSystem->Exit(30);return;}
 TFile source(input,"READ");if(source.IsZombie()||source.TestBit(TFile::kRecovered)){gSystem->Exit(31);return;}
 THnSparse* a=dynamic_cast<THnSparse*>(source.Get(sparseName));if(!a){gSystem->Exit(32);return;}
 {TFile f(output,"CREATE");if(f.IsZombie()){gSystem->Exit(33);return;}
  if(a->Write(sparseName)<=0){gSystem->Exit(34);return;}
  TNamed origin("SparseSplitSourceFile",input);origin.Write();
  TNamed note("SparseSplitNote","Exact copy of this source partition; combine lower-dimensional projections of all parts; keep explicit Sumw2; do not persist the full 8D sum");note.Write();f.Close();}
 if(splitError){gSystem->Exit(35);return;}
 TFile target(output,"READ");if(target.IsZombie()||target.TestBit(TFile::kRecovered)){gSystem->Exit(36);return;}
 THnSparse* b=dynamic_cast<THnSparse*>(target.Get(sparseName));
 if(!b||std::string(a->ClassName())!=b->ClassName()||a->GetNdimensions()!=b->GetNdimensions()||a->GetNbins()!=b->GetNbins()||a->GetEntries()!=b->GetEntries()||a->GetCalculateErrors()!=b->GetCalculateErrors()||a->GetSumw()!=b->GetSumw()||a->GetSumw2()!=b->GetSumw2()){gSystem->Exit(37);return;}
 int dims=a->GetNdimensions();std::vector<Int_t> ac(dims),bc(dims);
 for(int d=0;d<dims;++d)if(!sameAxis(a->GetAxis(d),b->GetAxis(d))||a->GetSumwx(d)!=b->GetSumwx(d)||a->GetSumwx2(d)!=b->GetSumwx2(d)){gSystem->Exit(38);return;}
 long double sum=0,error=0;
 for(Long64_t i=0;i<a->GetNbins();++i){
  double av=a->GetBinContent(i,&ac[0]),bv=b->GetBinContent(i,&bc[0]),ae=a->GetBinError2(i),be=b->GetBinError2(i);
  if(!TMath::Finite(av)||!TMath::Finite(bv)||!TMath::Finite(ae)||!TMath::Finite(be)||av!=bv||ae!=be||ac!=bc){gSystem->Exit(39);return;}
  sum+=av;error+=ae;
 }
 if(splitError){gSystem->Exit(40);return;}
 std::ofstream out(reportOrManifest);out<<std::setprecision(20)<<"CLASS\t"<<a->ClassName()<<"\nDIMENSIONS\t"<<dims<<"\nENTRIES\t"<<a->GetEntries()<<"\nFILLED_BINS\t"<<a->GetNbins()<<"\nSUM_CONTENT\t"<<sum<<"\nSUM_ERROR2\t"<<error<<"\nEXPLICIT_SUMW2\t"<<a->GetCalculateErrors()<<"\nDONE\t1\n";out.close();
 if(!out.good()){gSystem->Exit(41);return;}gSystem->Exit(0);
}
