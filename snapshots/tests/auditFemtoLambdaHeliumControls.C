// Read-only audit of real-event He4 + physical/mapped He3 output against old PID.
int auditFemtoLambdaHeliumControls(const char*,const char*);
#ifndef __CINT__
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace HeliumOutputAudit {
unsigned checks=0;
void Require(bool ok,const std::string& message) {
  ++checks;if(!ok)throw std::runtime_error(message);
}
std::string Index(int i) { std::ostringstream s;s<<i;return s.str(); }
TH1* Hist(TFile& f,const std::string& name) {
  TH1* h=dynamic_cast<TH1*>(f.Get(name.c_str()));Require(h!=0,"Missing "+name);return h;
}
TH2* Hist2(TFile& f,const std::string& name) {
  TH2* h=dynamic_cast<TH2*>(Hist(f,name));Require(h!=0,"Not TH2: "+name);return h;
}
double Sum(const TH1& h) {
  const int cells=(h.GetNbinsX()+2)*(h.GetDimension()>1?h.GetNbinsY()+2:1)*
    (h.GetDimension()>2?h.GetNbinsZ()+2:1);
  double sum=0;for(int i=0;i<cells;++i)sum+=h.GetBinContent(i);return sum;
}
void Equal(const TH1& a,const TH1& b,const std::string& label) {
  Require(a.GetDimension()==b.GetDimension()&&a.GetNbinsX()==b.GetNbinsX()&&
    a.GetNbinsY()==b.GetNbinsY()&&a.GetNbinsZ()==b.GetNbinsZ(),label+" shape");
  Require(a.GetEntries()==b.GetEntries(),label+" entries");
  const int cells=(a.GetNbinsX()+2)*(a.GetDimension()>1?a.GetNbinsY()+2:1)*
    (a.GetDimension()>2?a.GetNbinsZ()+2:1);
  for(int i=0;i<cells;++i) {
    Require(a.GetBinContent(i)==b.GetBinContent(i),label+" content");
    Require(a.GetBinError(i)==b.GetBinError(i),label+" error");
  }
}
void Ledger(TFile& current,TFile& legacy,const char* name) {
  TTree* a=dynamic_cast<TTree*>(current.Get(name));TTree* b=dynamic_cast<TTree*>(legacy.Get(name));
  Require(a&&b,std::string("Missing ledger ")+name);
  Require(a->GetEntries()==b->GetEntries(),std::string(name)+" row count");
  for(int l=0;l<a->GetListOfLeaves()->GetEntries();++l) {
    TLeaf* x=static_cast<TLeaf*>(a->GetListOfLeaves()->At(l));TLeaf* y=b->GetLeaf(x->GetName());
    Require(y&&std::string(x->GetTypeName())==y->GetTypeName(),std::string(name)+" leaf "+x->GetName());
    for(Long64_t i=0;i<a->GetEntries();++i) {
      Require(a->GetEntry(i)>0&&b->GetEntry(i)>0,"Ledger read");
      Require(x->GetValue()==y->GetValue(),std::string(name)+" exact row value "+x->GetName());
    }
  }
}
void Projection(TFile& f,const std::string& canonical,const std::string& legacy,
                const std::string& mode) {
  const std::string prefix=mode=="ME"?"Mixed_":"",dir=mode=="ME"?"mix/":"true/";
  const char* regions[]={"","_signal","_leftSB","_rightSB","_signalLow","_signalHigh"};
  const char* suffixes[]={"","","_SBNeg","_SBPos","_signalLow","_signalHigh"};
  for(int r=0;r<6;++r) {
    const std::string key="hKstar"+mode+"_lambda_"+canonical+regions[r];
    TH1* h=Hist(f,key);TH2* cent=Hist2(f,"hKstar"+mode+"VsCent_lambda_"+canonical+regions[r]);
    for(int x=0;x<=h->GetNbinsX()+1;++x) {
      double content=0,variance=0;
      for(int c=0;c<9;++c) {
        const int y=cent->GetYaxis()->FindBin(c);
        content+=cent->GetBinContent(x,y);
        variance+=std::pow(cent->GetBinError(x,y),2);
        const std::string ending=prefix+legacy+suffixes[r]+"_CentBin"+Index(c);
        if(r==0) {
          TH2* mass=Hist2(f,dir+"hKstarMass_"+prefix+legacy+"_CentBin"+Index(c));
          double massSum=0;for(int m=0;m<=mass->GetNbinsY()+1;++m)massSum+=mass->GetBinContent(x,m);
          Require(massSum==cent->GetBinContent(x,y),key+" mass projection");
        } else {
          TH2* mt=Hist2(f,dir+"hKstarMt_"+ending);
          double mtSum=0;for(int m=0;m<=mt->GetNbinsY()+1;++m)mtSum+=mt->GetBinContent(x,m);
          Require(mtSum==cent->GetBinContent(x,y),key+" mT projection");
          if(r<4)Require(Hist(f,dir+"hKstar_"+ending)->GetBinContent(x)==mtSum,key+" legacy cent alias");
        }
      }
      Require(content==h->GetBinContent(x),key+" cent sum");
      Require(std::fabs(variance-std::pow(h->GetBinError(x),2))<1.e-10*std::max(1.,variance),key+" variance");
      if(r==1)Require(Hist(f,"hKstar"+mode+"_lambda_"+canonical+"_signalLow")->GetBinContent(x)+
        Hist(f,"hKstar"+mode+"_lambda_"+canonical+"_signalHigh")->GetBinContent(x)==content,key+" low+high");
    }
    if(r>0&&r<4)Equal(*h,*Hist(f,dir+"hKstar_"+prefix+legacy+suffixes[r]),key+" legacy alias");
    if(canonical=="he3") {
      TH1* mapped=Hist(f,"hKstar"+mode+"_lambda_he3_as_he4"+regions[r]);
      Require(Sum(*h)==Sum(*mapped)&&h->GetEntries()==mapped->GetEntries(),key+" mapped same pair count");
      TH2* mappedCent=Hist2(f,"hKstar"+mode+"VsCent_lambda_he3_as_he4"+regions[r]);
      for(int c=0;c<9;++c) {
        double n=0,m=0;const int y=cent->GetYaxis()->FindBin(c);
        for(int x=0;x<=h->GetNbinsX()+1;++x){n+=cent->GetBinContent(x,y);m+=mappedCent->GetBinContent(x,y);}
        Require(n==m,key+" mapped same centrality pair count");
      }
      if(r>0&&r<4)Equal(*Hist(f,dir+"hQlab_"+prefix+legacy+suffixes[r]),
        *Hist(f,dir+"hQlab_"+prefix+"3HeAs4He"+suffixes[r]),key+" mapped unchanged qlab");
    }
  }
}
}

int auditFemtoLambdaHeliumControls(const char* currentFile,const char* legacyFile) {
  using namespace HeliumOutputAudit;
  try {
    checks=0;TFile current(currentFile,"READ"),legacy(legacyFile,"READ");
    Require(!current.IsZombie()&&!legacy.IsZombie(),"Open real-event ROOT inputs");
    TNamed* status=dynamic_cast<TNamed*>(current.Get("FemtoLambdaRunStatus"));
    Require(status&&std::string(status->GetTitle())=="completed","Completed processing");
    TParameter<Long64_t>* input=dynamic_cast<TParameter<Long64_t>*>(current.Get("inputEvents"));
    Require(input&&input->GetVal()>0,"Positive input count");
    TNamed* cal=dynamic_cast<TNamed*>(current.Get("FemtoNuclearDedxCalibration"));
    Require(cal&&std::string(cal->GetTitle()).find("exact consumed calibration YAML")!=std::string::npos,"Exact calibration metadata");
    Require(current.Get("FemtoLambdaHe3ControlDefinition")&&current.Get("FemtoLambdaPairMtDefinition_he3")&&
      current.Get("FemtoLambdaPairMtDefinition_he3_as_he4"),"Mass mapping metadata");
    TH2* sampler=Hist2(current,"hMixSamplerQA");
    Require(sampler->GetNbinsY()==8,"Sampler QA covers all eight declared channels");
    for(int x=0;x<=sampler->GetNbinsX()+1;++x)
      Require(sampler->GetBinContent(x,0)==0&&sampler->GetBinContent(x,9)==0,"No channel QA in flow bins");
    Ledger(current,legacy,"eventLedger");Ledger(current,legacy,"selectedLambda");
    const char* unchanged[]={"hLambda_InvMass","hLambda_Pt","hLambda_InvMass_vs_Cent9",
      "hLambda_PtVsYLab_signal","hLambdaProton_PtVsYLab_signal","hLambdaPion_PtVsYLab_signal","true/hDedxP"};
    for(unsigned i=0;i<sizeof(unchanged)/sizeof(unchanged[0]);++i)
      Equal(*Hist(current,unchanged[i]),*Hist(legacy,unchanged[i]),std::string("Unchanged ")+unchanged[i]);
    const char* names[]={"he4","he3","he3_as_he4"};const char* oldNames[]={"4He","3He","3HeAs4He"};
    for(int i=0;i<3;++i)for(int m=0;m<2;++m)Projection(current,names[i],oldNames[i],m?"ME":"SE");
    TH2* real=Hist2(current,"hNucleus_PtVsYLab_he3");TH2* mapped=Hist2(current,"hNucleus_PtVsYLab_he3_as_he4");
    Require(Sum(*real)==Sum(*mapped)&&real->GetEntries()==mapped->GetEntries(),"He3 mapped acceptance same tracks");
    for(int y=0;y<=real->GetNbinsY()+1;++y){double a=0,b=0;
      for(int x=0;x<=real->GetNbinsX()+1;++x){a+=real->GetBinContent(x,y);b+=mapped->GetBinContent(x,y);}
      Require(a==b,"He3 mapped unchanged physical Pt");
    }
    const int he3ChannelBin=sampler->GetYaxis()->FindBin(4);
    double he3SamplerSum=0;for(int x=0;x<=sampler->GetNbinsX()+1;++x)he3SamplerSum+=sampler->GetBinContent(x,he3ChannelBin);
    Require(he3SamplerSum>0,"Real He3 channel sampler QA is recorded in range");
    std::cout<<"{\"status\":\"passed\",\"checks\":"<<checks<<",\"inputEvents\":"<<input->GetVal()
      <<",\"he4Tracks\":"<<Hist(current,"hNucleus_PtVsYLab_he4")->GetEntries()
      <<",\"he3Tracks\":"<<real->GetEntries()<<",\"he3SignalSE\":"<<Hist(current,"hKstarSE_lambda_he3_signal")->GetEntries()
      <<",\"he3SignalME\":"<<Hist(current,"hKstarME_lambda_he3_signal")->GetEntries()
      <<",\"he3SamplerQaSum\":"<<he3SamplerSum<<"}"<<std::endl;
    return 0;
  }catch(const std::exception& error) {
    std::cerr<<"Helium output audit failed: "<<error.what()<<std::endl;gSystem->Exit(1);return 1;
  }
}
#endif
