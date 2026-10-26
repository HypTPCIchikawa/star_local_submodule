// Run compiled, read-only study QA in a process-private ROOT5 build directory.
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
void runLambdaKfPurityQa(const char* source, const char* expression) {
  TString dir=TString::Format("%s/kf-purity-qa-%d",gSystem->TempDirectory(),gSystem->GetPid());
  if(gSystem->mkdir(dir,kTRUE)!=0) { gSystem->Exit(21); return; }
  gSystem->SetBuildDir(dir,kTRUE);
  if(!gSystem->CompileMacro(source,"k")) { gSystem->Exit(22); return; }
  Int_t error=0;
  Long_t result=gROOT->ProcessLine(expression,&error);
  gSystem->Exit(error?23:result);
}
