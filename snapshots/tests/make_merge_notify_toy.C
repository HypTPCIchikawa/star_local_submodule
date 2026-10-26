#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>

// Creates only new toy files under the explicit temporary directory.
void make_merge_notify_toy(const char* directory) {
  TH1::AddDirectory(kFALSE);
  for (int i=0; i<3; ++i) {
    TString path=TString::Format("%s/toy_%d.root",directory,i);
    if (!gSystem->AccessPathName(path)) { gSystem->Exit(20); return; }
    TFile f(path,"CREATE");
    if (f.IsZombie()) { gSystem->Exit(21); return; }
    TH1D counts("counts","counts",10,0,10);
    for (int j=0; j<=i; ++j) counts.Fill(j+0.5);
    counts.Write();
    TNamed tag("definition","shared definition"); tag.Write();
    TTree tree("events","events"); int value=i;
    tree.Branch("value",&value,"value/I"); tree.Fill(); tree.Write(); tree.SetDirectory(0);
    f.mkdir("true")->cd();
    TH2D map("map","map",10,0,10,10,0,10);
    map.Fill(i+0.5,i+0.5); map.Write();
    f.Close();
  }
  gSystem->Exit(0);
}
