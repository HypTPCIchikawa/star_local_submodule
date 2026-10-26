// Settings and all output paths are explicit. Synthetic data only.
void run_test_peak_response(const char* macro, const char* build,
                            const char* settings, const char* report) {
  if (!macro || !build || !settings || !report || !gSystem->AccessPathName(build) ||
      !gSystem->AccessPathName(report) || gSystem->mkdir(build, true) != 0) {
    Error("run_test_peak_response", "Require fresh explicit build and report paths");
    gSystem->Exit(1); return;
  }
  gSystem->SetBuildDir(build, true);
  if (!gSystem->CompileMacro(macro, "k")) { gSystem->Exit(1); return; }
  TString invocation;
  invocation.Form("TestPeakResponseFit(\"%s\",\"%s\");", settings, report);
  Int_t error = 0;
  Long_t status = gROOT->ProcessLine(invocation.Data(), &error);
  if (error || status != 0) gSystem->Exit(1);
}
