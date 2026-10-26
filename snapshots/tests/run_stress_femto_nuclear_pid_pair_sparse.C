#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
TString SparseStressQuote(const char* value){TString s(value?value:"");s.ReplaceAll("\\","\\\\");s.ReplaceAll("\"","\\\"");return "\""+s+"\"";}
void run_stress_femto_nuclear_pid_pair_sparse(const char* mainconf,const char* output,Long64_t count,const char* build,const char* source,const char* library){
  if(!mainconf||!*mainconf||!output||!*output||!build||!*build||!source||!*source){gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath(TString(gSystem->WorkingDirectory())+"/src/third_party/yaml-cpp/include");
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(2);return;}
  Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);if(error){gSystem->Exit(2);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};for(unsigned i=0;i<4;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(2);return;}
  const TString cwd=gSystem->WorkingDirectory();const char* libs[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so"};
  for(unsigned i=0;i<3;++i)if(gSystem->Load(cwd+"/lib/"+libs[i])<0){gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");gInterpreter->AddIncludePath(cwd+"/StMaker/common");gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  if(count<0){gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);const int ok=gSystem->CompileMacro(source,"kf");gSystem->Exit(ok?0:2);return;}
  if(!library||!*library||gSystem->Load(library)<0){gSystem->Exit(2);return;}
  const TString call=TString::Format("stress_femto_nuclear_pid_pair_sparse(%s,%s,%lld)",SparseStressQuote(mainconf).Data(),SparseStressQuote(output).Data(),count);
  const Long_t result=gROOT->ProcessLine(call,&error);gSystem->Exit(error?2:int(result));
}
