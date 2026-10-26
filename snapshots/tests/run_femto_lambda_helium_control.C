#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
#include <iostream>

TString HeliumControlQuote(const char* value){
  TString out(value?value:"");out.ReplaceAll("\\","\\\\");out.ReplaceAll("\"","\\\"");return "\""+out+"\"";
}
void run_femto_lambda_helium_control(const char* mainconf,const char* rollback){
  if(!mainconf||!*mainconf||!rollback||!*rollback){gSystem->Exit(2);return;}
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(2);return;}
  Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);
  if(error){gSystem->Exit(2);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};
  for(unsigned i=0;i<4;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(2);return;}
  const TString cwd=gSystem->WorkingDirectory();
  const char* libs[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so"};
  for(unsigned i=0;i<3;++i)if(gSystem->Load(cwd+"/lib/"+libs[i])<0){gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");
  gInterpreter->AddIncludePath(cwd+"/StMaker/common");gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  const TString build=TString::Format("%s/tmp/helium-control-test-%d",cwd.Data(),gSystem->GetPid());
  gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);
  if(!gSystem->CompileMacro("tests/femto_lambda_helium_control.C","kf")){gSystem->Exit(2);return;}
  const TString call=TString::Format("femto_lambda_helium_control(%s,%s)",HeliumControlQuote(mainconf).Data(),HeliumControlQuote(rollback).Data());
  const Long_t result=gROOT->ProcessLine(call,&error);gSystem->Exit(error?2:int(result));
}
