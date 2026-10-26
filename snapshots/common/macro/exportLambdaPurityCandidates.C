// Read-only, validated export for the separate Lambda purity QA study.
// No PicoDst event loop, recuts, mass constraints, or production edits.
#include "../../tests/check_kfparticle_imp5_output.C"
#include <TObjArray.h>
#include <fstream>
#include <iomanip>

Int_t exportLambdaPurityCandidates(const char* input, const char* outputDir,
                                  Long64_t expectedInput, Long64_t expectedReco) {
  if (!input || !outputDir || !*outputDir || expectedInput<=0 || expectedReco<=0) return 1;
  if (!gSystem->AccessPathName(outputDir)) {
    std::cerr<<"Refusing existing output directory"<<std::endl; return 2;
  }
  if (KfImp5OutputQaInspect(input,kTRUE)) return 3;
  TFile file(input,"READ");
  if(file.IsZombie() || file.TestBit(TFile::kRecovered)) return 4;
  TH1* stages=dynamic_cast<TH1*>(file.Get("hKfStages"));
  TTree* tree=dynamic_cast<TTree*>(file.Get("KfLambdaCandidates"));
  TH1* selected=dynamic_cast<TH1*>(file.Get("hKfLambdaMassSelected"));
  if(!stages || !tree || !selected || stages->GetBinContent(1)!=expectedInput ||
     stages->GetBinContent(2)!=expectedReco) return 5;
  if(gSystem->mkdir(outputDir,kTRUE)!=0) return 6;
  std::ofstream csv((TString(outputDir)+"/candidates.csv").Data());
  csv<<std::setprecision(17);
  TObjArray* leaves=tree->GetListOfLeaves();
  for(Int_t j=0;j<leaves->GetEntriesFast();++j) {
    if(j) csv<<",";
    csv<<leaves->At(j)->GetName();
  }
  csv<<"\n";
  for(Long64_t i=0;i<tree->GetEntries();++i) {
    if(tree->GetEntry(i)<=0) return 7;
    for(Int_t j=0;j<leaves->GetEntriesFast();++j) {
      if(j) csv<<",";
      csv<<static_cast<TLeaf*>(leaves->At(j))->GetValue();
    }
    csv<<"\n";
  }
  csv.close(); if(!csv) return 8;
  std::ofstream hist((TString(outputDir)+"/selected_histogram.csv").Data());
  hist<<std::setprecision(17)<<"bin,low,high,count,error\n";
  for(Int_t i=0;i<=selected->GetNbinsX()+1;++i)
    hist<<i<<","<<selected->GetBinLowEdge(i)<<","
        <<selected->GetBinLowEdge(i)+selected->GetBinWidth(i)<<","
        <<selected->GetBinContent(i)<<","<<selected->GetBinError(i)<<"\n";
  hist.close(); if(!hist) return 9;
  std::ofstream meta((TString(outputDir)+"/metadata.txt").Data());
  meta<<"inputRoot: "<<input<<"\ninputEvents: "<<expectedInput
      <<"\nreconstructedEvents: "<<expectedReco
      <<"\nrawEntries: "<<tree->GetEntries()<<"\n";
  const char* keys[]={"KFParticleEffectiveConfiguration","KFEventSelectionConfiguration"};
  for(Int_t i=0;i<2;++i) {
    TObjString* obj=dynamic_cast<TObjString*>(file.Get(keys[i]));
    if(!obj) return 10;
    meta<<"\n"<<keys[i]<<"\n"<<obj->GetString().Data()<<"\n";
  }
  meta.close(); if(!meta) return 11;
  std::cout<<"PASS purity export: "<<tree->GetEntries()<<" raw candidates, "<<expectedInput
           <<" input and "<<expectedReco<<" reconstructed events"<<std::endl;
  return 0;
}
