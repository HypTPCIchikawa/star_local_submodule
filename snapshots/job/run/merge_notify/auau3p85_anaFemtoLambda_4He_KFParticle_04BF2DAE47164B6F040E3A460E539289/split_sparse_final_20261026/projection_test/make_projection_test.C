#include <TFile.h>
#include <THnSparse.h>
#include <TSystem.h>
void make_projection_test(const char* path){TFile f(path,"CREATE");if(f.IsZombie()){gSystem->Exit(1);return;}Int_t n[8];Double_t lo[8],hi[8];for(int j=0;j<8;++j){n[j]=3;lo[j]=0;hi[j]=3;}THnSparseD h("hNuclearPidPairJoint","test",8,n,lo,hi);h.Sumw2();Double_t x[8]={0.2,1.1,0.1,0,0,0,1,0};h.Fill(x,2);x[1]=4;h.Fill(x,3);x[0]=4;h.Fill(x,1);h.Write();f.Close();gSystem->Exit(0);}
