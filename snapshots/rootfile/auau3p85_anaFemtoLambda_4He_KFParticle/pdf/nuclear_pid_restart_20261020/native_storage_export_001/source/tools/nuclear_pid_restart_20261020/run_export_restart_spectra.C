// Explicit-path, fresh-build wrapper. Does not fit or modify the input ROOT.
void run_export_restart_spectra(const char* macro, const char* build,
                               const char* cache, const char* settings,
                               const char* slices, const char* output) {
  if (!macro || !build || !cache || !settings || !slices || !output ||
      !gSystem->AccessPathName(build) || !gSystem->AccessPathName(output) ||
      gSystem->mkdir(build, true) != 0) {
    Error("run_export_restart_spectra", "Require explicit fresh build and output paths");
    gSystem->Exit(1); return;
  }
  gSystem->SetBuildDir(build, true);
  if (!gSystem->CompileMacro(macro, "k")) { gSystem->Exit(1); return; }
  TString invocation;
  invocation.Form("ExportRestartSpectra(\"%s\",\"%s\",\"%s\",\"%s\");",
                  cache, settings, slices, output);
  Int_t error = 0;
  gROOT->ProcessLine(invocation.Data(), &error);
  if (error) gSystem->Exit(1);
}
