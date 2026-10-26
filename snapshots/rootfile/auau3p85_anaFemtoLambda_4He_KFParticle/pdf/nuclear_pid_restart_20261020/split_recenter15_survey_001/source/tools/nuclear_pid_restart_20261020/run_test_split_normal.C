// Explicit inputs and a fresh ACLiC build directory; no production PID changes.
void run_test_split_normal(const char* source, const char* buildDirectory,
                           const char* settingsPath, const char* reportPath) {
  if (!source || !source[0] || !buildDirectory || !buildDirectory[0] ||
      !settingsPath || !settingsPath[0] || !reportPath || !reportPath[0]) {
    Error("run_test_split_normal", "Explicit source/build/settings/report required");
    gSystem->Exit(1); return;
  }
  if (!gSystem->AccessPathName(buildDirectory) || !gSystem->AccessPathName(reportPath)) {
    Error("run_test_split_normal", "Build directory and report must be fresh");
    gSystem->Exit(1); return;
  }
  if (gSystem->mkdir(buildDirectory, kTRUE) != 0) {
    Error("run_test_split_normal", "Cannot create build directory");
    gSystem->Exit(1); return;
  }
  gSystem->SetBuildDir(buildDirectory, kTRUE);
  if (gSystem->CompileMacro(source, "k") == 0) {
    Error("run_test_split_normal", "CompileMacro failed");
    gSystem->Exit(1); return;
  }
  TString call;
  call.Form("TestSplitNormalFit(\"%s\",\"%s\");", settingsPath, reportPath);
  Int_t error = 0;
  const Long_t result = gROOT->ProcessLine(call.Data(), &error);
  gSystem->Exit(error == 0 && result == 0 ? 0 : 1);
}
