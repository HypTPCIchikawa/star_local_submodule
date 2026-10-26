void RunSelfTestIterativeHelium(){
  const char* build="/tmp/iterative_helium_numeric_20261014_001";
  if(!gSystem->AccessPathName(build)){Error("iterative numeric","new build directory required");gSystem->Exit(2);return;}
  if(gSystem->mkdir(build,kTRUE)){gSystem->Exit(2);return;}
  gSystem->SetBuildDir(build,kTRUE);
  int compileError=0;
  gROOT->ProcessLine(".L /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/Todo/20261013/validation/iterative_helium_20261014_001/source/tools/he4_purity/iterative_logpoly/SelfTestIterativeHelium.C++",&compileError);
  if(compileError){gSystem->Exit(2);return;}
  int testError=0;
  gROOT->ProcessLine("SelfTestIterativeHelium(\"/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/Todo/20261013/validation/iterative_helium_20261014_001/settings.tsv\",\"/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/Todo/20261013/validation/iterative_helium_20261014_001/selftest.json\")",&testError);
  if(testError)gSystem->Exit(2);
}
