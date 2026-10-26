#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
TString NuclearDomainTestQuote(const char* text) {
  TString value(text?text:"");value.ReplaceAll("\\","\\\\");value.ReplaceAll("\"","\\\"");return "\""+value+"\"";
}
void run_femto_nuclear_pid_domain_test(const char* sourceconf,const char* outputDirectory) {
  if(!sourceconf || !*sourceconf || !outputDirectory || !*outputDirectory) {gSystem->Exit(1);return;}
  const TString cwd=gSystem->WorkingDirectory();
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0) {gSystem->Exit(1);return;}
  Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);
  if(error){gSystem->Exit(1);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};
  for(int i=0;i<4;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(1);return;}
  const char* local[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so"};
  for(int i=0;i<3;++i)if(gSystem->Load(cwd+"/lib/"+local[i])<0){gSystem->Exit(1);return;}
  gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");
  gInterpreter->AddIncludePath(cwd+"/StMaker/common");gInterpreter->AddIncludePath(cwd+"/src/third_party/yaml-cpp/include");
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  const TString build=TString::Format("%s/tmp/femto-nuclear-domain-test-%d",cwd.Data(),gSystem->GetPid());
  gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);
  if(!gSystem->CompileMacro("tests/test_femto_nuclear_pid_domain.C","kf")){gSystem->Exit(1);return;}
  const TString call=TString::Format("test_femto_nuclear_pid_domain(%s,%s)",
    NuclearDomainTestQuote(sourceconf).Data(),NuclearDomainTestQuote(outputDirectory).Data());
  const Long_t result=gROOT->ProcessLine(call,&error);if(error || result)gSystem->Exit(1);
}
