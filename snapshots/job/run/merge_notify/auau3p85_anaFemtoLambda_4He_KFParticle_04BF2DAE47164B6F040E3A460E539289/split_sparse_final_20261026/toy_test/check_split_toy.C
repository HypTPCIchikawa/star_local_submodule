#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
void check_split_toy(const char* path){TFile f(path,"READ");if(f.IsZombie()||f.GetKey("hNuclearPidPairJoint")||!f.GetKey("NuclearPidPairJointSplitManifest")){gSystem->Exit(1);return;}TH1* h=(TH1*)f.Get("counts");TH1* nested=(TH1*)f.Get("true/map");if(!h||!nested||h->GetEntries()!=3||h->GetBinContent(1)!=6||nested->GetEntries()!=3||nested->Integral()!=6){gSystem->Exit(2);return;}gSystem->Exit(0);}
