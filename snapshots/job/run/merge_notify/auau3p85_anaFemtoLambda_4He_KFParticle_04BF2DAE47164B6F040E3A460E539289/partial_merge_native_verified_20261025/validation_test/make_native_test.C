#include <TFile.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TString.h>
void make_native_test(const char* dir){
 TH1::AddDirectory(kFALSE);TH1F ref("h","h",1,0,1);bool first=true;
 for(int i=0;i<5;++i){TFile f(TString::Format("%s/input_%d.root",dir,i),"CREATE");if(f.IsZombie()){gSystem->Exit(1);return;}
 TH1F h("h","h",1,0,1);double n=i?1:16777216;h.SetBinContent(1,n);h.SetEntries(n);h.Write();f.Close();
 if(first){ref=h;first=false;}else ref.Add(&h);}
 {TFile f(TString::Format("%s/good.root",dir),"CREATE");ref.Write();f.Close();}
 double entries=ref.GetEntries();ref.SetBinContent(1,ref.GetBinContent(1)+2);ref.SetEntries(entries);
 {TFile f(TString::Format("%s/corrupt.root",dir),"CREATE");ref.Write();f.Close();}
 gSystem->Exit(0);
}
