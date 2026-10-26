#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TSystem.h>
#include <TString.h>
void make_split_toy(const char* dir){
 TH1::AddDirectory(kFALSE);
 for(int i=0;i<3;++i){TFile f(TString::Format("%s/input_%d.root",dir,i),"CREATE");if(f.IsZombie()){gSystem->Exit(1);return;}
 TH1F h("counts","counts",10,0,10);h.Sumw2();h.Fill(0.5,i+1);h.Write();
 Int_t bins[8];Double_t lo[8],hi[8],x[8];for(int j=0;j<8;++j){bins[j]=3;lo[j]=0;hi[j]=3;x[j]=j%3+0.1;}
 THnSparseD sp("hNuclearPidPairJoint","toy sparse",8,bins,lo,hi);sp.Sumw2();sp.GetAxis(4)->SetBinLabel(1,"test");
 sp.Fill(x,0.5*(i+1));x[0]=-1;sp.Fill(x,2.);x[0]=4;sp.Fill(x,-0.25);sp.Write();
 f.mkdir("true")->cd();TH2D map("map","map",2,0,2,2,0,2);map.Fill(0.5,0.5,i+1);map.Write();f.Close();}
 gSystem->Exit(0);
}
