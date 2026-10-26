// QA/export only. Compile the existing exporter: ROOT5 CINT cannot interpret
// its TLeaf static_cast expression. Each process has a separate build directory.
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

void exportLambdaKfUpstreamStudy(const char* input, const char* outputDir,
                                Long64_t expectedInput) {
  Long64_t reconstructed = -1;
  {
    TFile file(input, "READ");
    TH1* stages = dynamic_cast<TH1*>(file.Get("hKfStages"));
    if (file.IsZombie() || file.TestBit(TFile::kRecovered) || !stages ||
        stages->GetBinContent(1) != expectedInput) { gSystem->Exit(21); return; }
    reconstructed = static_cast<Long64_t>(stages->GetBinContent(2));
  }
  TString build = TString::Format("%s/kf-upstream-export-%d", gSystem->TempDirectory(), gSystem->GetPid());
  if (gSystem->mkdir(build, kTRUE) != 0) { gSystem->Exit(22); return; }
  gSystem->SetBuildDir(build, kTRUE);
  if (!gSystem->CompileMacro("common/macro/exportLambdaPurityCandidates.C", "k")) {
    gSystem->Exit(23); return;
  }
  TString call = TString::Format("exportLambdaPurityCandidates(\"%s\",\"%s\",%lld,%lld)",
                                 input, outputDir, expectedInput, reconstructed);
  Int_t error = 0;
  Long_t result = gROOT->ProcessLine(call, &error);
  gSystem->Exit(error ? 24 : result);
}
