// Runtime wrapper; explicit paths only. Does not fit physical data.
void run_test_restart_slices(const char* macro, const char* build, const char* output) {
  if (!macro || !build || !output || !gSystem->AccessPathName(build) ||
      !gSystem->AccessPathName(output) || gSystem->mkdir(build, true) != 0) {
    Error("run_test_restart_slices", "Require fresh explicit build and test paths");
    gSystem->Exit(1); return;
  }
  gSystem->SetBuildDir(build, true);
  if (!gSystem->CompileMacro(macro, "k")) { gSystem->Exit(1); return; }
  TString invocation;
  invocation.Form("TestRestartSlices(\"%s\");", output);
  Int_t error = 0;
  gROOT->ProcessLine(invocation.Data(), &error);
  if (error) gSystem->Exit(1);
}
