// The caller supplies a fresh build directory; test data stay in memory.
void run_test_float16_mapping(const char* source, const char* buildDirectory) {
  if (!source || !source[0] || !buildDirectory || !buildDirectory[0]) {
    Error("run_test_float16_mapping", "Explicit test source and fresh build directory required");
    gSystem->Exit(1); return;
  }
  if (!gSystem->AccessPathName(buildDirectory)) {
    Error("run_test_float16_mapping", "Build directory must be fresh");
    gSystem->Exit(1); return;
  }
  if (gSystem->mkdir(buildDirectory, kTRUE) != 0) {
    Error("run_test_float16_mapping", "Cannot create test build directory");
    gSystem->Exit(1); return;
  }
  gSystem->SetBuildDir(buildDirectory, kTRUE);
  if (gSystem->CompileMacro(source, "k") == 0) {
    Error("run_test_float16_mapping", "CompileMacro failed");
    gSystem->Exit(1); return;
  }
  Int_t error = 0;
  const Long_t result = gROOT->ProcessLine("TestFloat16BinMapping();", &error);
  gSystem->Exit(error == 0 && result == 0 ? 0 : 1);
}
