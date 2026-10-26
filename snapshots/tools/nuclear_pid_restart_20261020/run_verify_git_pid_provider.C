// Invoke through singularity_root_batch.sh with an explicit real mainconf.
void run_verify_git_pid_provider(const char* oracle,const char* report,
                                 double relative,double absolute) {
  if(!oracle || !oracle[0] || !report || !report[0]) {
    Error("run_verify_git_pid_provider","Explicit paths are required");
    gSystem->Exit(1);return;
  }
  Int_t error=0;
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gROOT->ProcessLine("loadSharedLibraries();",&error);
  if(error){gSystem->Exit(1);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};
  for(unsigned i=0;i<sizeof(star)/sizeof(star[0]);++i)
    if(gSystem->Load(star[i])<0){gSystem->Exit(1);return;}
  TString cwd=gSystem->WorkingDirectory();
  const char* local[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so"};
  for(unsigned i=0;i<sizeof(local)/sizeof(local[0]);++i)
    if(gSystem->Load((cwd+"/lib/"+local[i]).Data())<0){gSystem->Exit(1);return;}
  TString include="-I"+cwd+" -I"+cwd+"/include -I"+cwd+"/StMaker/common -I"+
                  cwd+"/src/third_party/yaml-cpp/include -I$STAR/StRoot";
  gSystem->AddIncludePath(include.Data());
  TString links="-L"+cwd+"/lib -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,"+cwd+"/lib";
  gSystem->AddLinkedLibs(links.Data());
  TString build=TString(report)+"_build";
  if(!gSystem->AccessPathName(build.Data()) || gSystem->mkdir(build.Data(),kTRUE)!=0) {
    Error("run_verify_git_pid_provider","Build directory must be new: %s",build.Data());
    gSystem->Exit(1);return;
  }
  gSystem->SetBuildDir(build.Data(),kFALSE);
  if(!gSystem->CompileMacro("tools/nuclear_pid_restart_20261020/verify_git_pid_provider.C","kf")) {
    gSystem->Exit(1);return;
  }
  TString call;
  call.Form("verify_git_pid_provider(\"%s\",\"%s\",%.17g,%.17g)",oracle,report,relative,absolute);
  Long_t result=gROOT->ProcessLine(call.Data(),&error);
  gSystem->Exit(error || result!=0 ? 1 : 0);
}
