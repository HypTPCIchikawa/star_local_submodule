#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TNamed.h"
#include "TSystem.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
Long64_t qaHaddChecks=0;
int qaHaddHistograms=0;
void QaHaddNeed(bool ok,const std::string& message) {
  ++qaHaddChecks;
  if(!ok) throw std::runtime_error(message);
}
bool QaHaddNear(double actual,double expected) {
  return std::fabs(actual-expected)<=1.e-11*(1.+std::fabs(expected));
}
void QaHaddAxis(const TAxis& a,const TAxis& b) {
  QaHaddNeed(a.GetNbins()==b.GetNbins(),"axis shape");
  for(int i=1;i<=a.GetNbins()+1;++i)
    QaHaddNeed(a.GetBinLowEdge(i)==b.GetBinLowEdge(i),"exact axis edges");
}
void QaHaddDirectory(TDirectory* source,TFile& before,TFile& after,const std::string& prefix) {
  TIter next(source->GetListOfKeys());
  TKey* key=0;
  while((key=(TKey*)next())) {
    const std::string name=key->GetName(),path=prefix+name;
    if(std::string(key->GetClassName())=="TDirectoryFile") {
      QaHaddDirectory(source->GetDirectory(name.c_str()),before,after,path+"/");
    } else if(name.find("hKstarNSigma")==0 || name.find("hKstarPOverQ")==0) {
      TH2* a=dynamic_cast<TH2*>(before.Get(path.c_str()));
      TH2* b=dynamic_cast<TH2*>(after.Get(path.c_str()));
      QaHaddNeed(a&&b,"missing dense QA: "+path);
      QaHaddNeed(std::string(a->ClassName())=="TH2D" && std::string(b->ClassName())=="TH2D","QA type");
      QaHaddAxis(*a->GetXaxis(),*b->GetXaxis()); QaHaddAxis(*a->GetYaxis(),*b->GetYaxis());
      QaHaddNeed(b->GetEntries()==2.*a->GetEntries(),"dense entries doubled");
      for(int x=0;x<=a->GetNbinsX()+1;++x) for(int y=0;y<=a->GetNbinsY()+1;++y) {
        const int cell=a->GetBin(x,y);
        QaHaddNeed(b->GetBinContent(cell)==2.*a->GetBinContent(cell),"dense content doubled incl flow");
        QaHaddNeed(QaHaddNear(std::pow(b->GetBinError(cell),2),2.*std::pow(a->GetBinError(cell),2)),
                   "dense variance doubled incl flow");
      }
      ++qaHaddHistograms;
      delete a; delete b;
    }
  }
}
}

void check_femto_nuclear_pair_pid_qa_hadd(const char* sourcePath,const char* mergedPath) {
  try {
    TFile source(sourcePath,"READ"),merged(mergedPath,"READ");
    QaHaddNeed(!source.IsZombie()&&!merged.IsZombie(),"closed ROOT readback");
    QaHaddNeed(!source.TestBit(TFile::kRecovered)&&!merged.TestBit(TFile::kRecovered),"not recovered ROOT");
    QaHaddDirectory(&source,source,merged,"");
    QaHaddNeed(qaHaddHistograms==1080,"all 1080 dense QA maps");
    THnSparseD* a=dynamic_cast<THnSparseD*>(source.Get("hNuclearPidPairJoint"));
    THnSparseD* b=dynamic_cast<THnSparseD*>(merged.Get("hNuclearPidPairJoint"));
    QaHaddNeed(a&&b,"joint sparse exists after hadd");
    QaHaddNeed(a->GetNdimensions()==8&&b->GetNdimensions()==8,"joint dimension");
    for(int d=0;d<8;++d)QaHaddAxis(*a->GetAxis(d),*b->GetAxis(d));
    QaHaddNeed(b->GetEntries()==2.*a->GetEntries(),"joint entries doubled");
    QaHaddNeed(b->GetNbins()==a->GetNbins(),"same occupied coordinates for duplicated input");
    for(Long64_t bin=0;bin<a->GetNbins();++bin) {
      Int_t coordinate[8];
      const double content=a->GetBinContent(bin,coordinate);
      const Long64_t target=b->GetBin(coordinate,kFALSE);
      QaHaddNeed(target>=0,"joint coordinate retained including flows");
      QaHaddNeed(b->GetBinContent(target)==2.*content,"joint coordinate content doubled");
      QaHaddNeed(b->GetBinError2(target)==2.*a->GetBinError2(bin),"joint coordinate variance doubled");
    }
    TNamed* metadata=dynamic_cast<TNamed*>(merged.Get("NuclearPidPairJointMetadata"));
    QaHaddNeed(metadata!=0,"joint definition metadata retained (not an additive counter)");
    std::cout<<"NUCLEAR_PID_PAIR_QA_HADD_PASS qaTH2="<<qaHaddHistograms
             <<" sparseEntries="<<b->GetEntries()<<" sparseFilledBins="<<b->GetNbins()
             <<" checks="<<qaHaddChecks<<std::endl;
    source.Close(); merged.Close();
    gSystem->Exit(0);
  } catch(const std::exception& error) {
    std::cerr<<"NUCLEAR_PID_PAIR_QA_HADD_FAIL "<<error.what()<<std::endl;
    gSystem->Exit(1);
  }
}
