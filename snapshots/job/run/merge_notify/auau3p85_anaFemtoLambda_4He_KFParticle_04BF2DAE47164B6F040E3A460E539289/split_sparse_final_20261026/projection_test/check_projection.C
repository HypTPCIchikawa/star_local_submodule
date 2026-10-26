#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <cmath>
void check_projection(const char* path){TFile f(path,"READ");TH1D* h=dynamic_cast<TH1D*>(f.Get("hProjectedNuclearPidPairJoint"));if(!h||h->GetBinContent(1)!=15||h->GetBinContent(4)!=3||std::fabs(h->GetBinError(1)*h->GetBinError(1)-39)>1e-10||std::fabs(h->GetBinError(4)*h->GetBinError(4)-3)>1e-10){gSystem->Exit(1);return;}gSystem->Exit(0);}
