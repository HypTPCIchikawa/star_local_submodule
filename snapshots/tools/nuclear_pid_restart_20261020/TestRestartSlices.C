// Synthetic data/format tests only; no physical calibration values are used.
void TestRestartSlices(const char*);
#ifndef __CINT__
#include "FitRestartSlices.C"

namespace RestartPidTest {

static void Expect(bool condition, const std::string& message, int& passed) {
  RestartPid::Check(condition, "Self-test failure: " + message);
  ++passed;
}

static void Write(const std::string& path, const std::string& value) {
  std::ofstream out(path.c_str());
  RestartPid::Check(out.good(), "Cannot write synthetic test fixture");
  out << value;
  out.close();
  RestartPid::Check(out.good(), "Synthetic fixture write failed");
}

static void Run(const char* outputDir) {
  RestartPid::Check(gSystem->AccessPathName(outputDir), "Use a fresh synthetic test directory");
  RestartPid::Check(gSystem->mkdir(outputDir, true) == 0, "Cannot create test directory");
  int passed = 0;
  TH1D histogram("synthetic_native", "synthetic native counts", 10, 0, 1);
  histogram.SetDirectory(0);
  for (int bin = 1; bin <= 10; ++bin) histogram.SetBinContent(bin, bin == 5 ? 0 : bin);
  RestartPid::Slice slice;
  slice.id = "synthetic"; slice.target = "he4"; slice.mode = "TAG_D"; slice.rBin = 1;
  slice.fitLo = 0.2; slice.fitHi = 0.9; slice.gateLo = 0.3; slice.gateHi = 0.8;
  std::vector<NuclearFit::Bin> bins = RestartPid::NativeBins(histogram, slice);
  Expect(bins.size() == 7, "Native-bin count", passed);
  Expect(bins[2].count == 0, "Keep native zero bins", passed);
  Expect(std::fabs(bins[0].lo - 0.2) < 1e-12 && std::fabs(bins.back().hi - 0.9) < 1e-12,
         "Contiguous native boundaries", passed);
  slice.fitLo = 0.25;
  bool rejected = false;
  try { RestartPid::NativeBins(histogram, slice); } catch (const std::exception&) { rejected = true; }
  Expect(rejected, "Reject partial/rebinned fit edges", passed);

  const std::string slicesPath = std::string(outputDir) + "/synthetic_slices.tsv";
  const std::string header = "id\tmode\ttarget\trBin\tfitLo\tfitHi\tgateLo\tgateHi\n";
  Write(slicesPath, header + "one\tTAG_D\the4\t1\t0.2\t0.9\t0.3\t0.8\n");
  std::vector<RestartPid::Slice> parsed = RestartPid::ReadSlices(slicesPath.c_str());
  Expect(parsed.size() == 1 && parsed[0].target == "he4" && parsed[0].mode == "TAG_D",
         "He4 shares deuteron TOF tag", passed);
  Write(slicesPath, header + "one\tTAG_HE3\the4\t1\t0.2\t0.9\t0.3\t0.8\n");
  rejected = false;
  try { RestartPid::ReadSlices(slicesPath.c_str()); } catch (const std::exception&) { rejected = true; }
  Expect(rejected, "Reject different He4 TOF tag", passed);
  Write(slicesPath, header + "one\tTAG_D\td\t1\t0.2\t0.9\t0.3\t0.8\n"
                           "one\tTAG_D\td\t2\t0.2\t0.9\t0.3\t0.8\n");
  rejected = false;
  try { RestartPid::ReadSlices(slicesPath.c_str()); } catch (const std::exception&) { rejected = true; }
  Expect(rejected, "Reject duplicate slice ids", passed);
  Write(slicesPath, header + "one\tTAG_D\td\t1.5\t0.2\t0.9\t0.3\t0.8\n");
  rejected = false;
  try { RestartPid::ReadSlices(slicesPath.c_str()); } catch (const std::exception&) { rejected = true; }
  Expect(rejected, "Reject non-native R bin", passed);

  const std::string nodesPath = std::string(outputDir) + "/synthetic_nodes.tsv";
  Write(nodesPath, "id\tspecies\tweight\tmean\tsigma\n"
                   "one\td\t2\t1\t0.1\none\td\t2\t3\t0.3\n"
                   "one\tt\t1\t4\t0.4\none\the3\t1\t5\t0.5\none\the4\t1\t6\t0.6\n");
  std::map<std::string, RestartPid::Shapes> shapes = RestartPid::ReadNodes(nodesPath.c_str());
  Expect(shapes.size() == 1 && shapes["one"].size() == 4, "Require four known PID curves", passed);
  Expect(std::fabs(shapes["one"]["d"].mean - 2) < 1e-12 &&
         std::fabs(shapes["one"]["d"].sigma - 0.2) < 1e-12,
         "Weighted means and intrinsic widths", passed);
  Expect(std::fabs(shapes["one"]["d"].nodes[0].weight - 0.5) < 1e-12,
         "Normalize native-R node weights", passed);
  Expect(RestartPid::Quote("x\n\"\\") == "\"x\\n\\\"\\\\\"", "Escape JSON strings", passed);

  std::ofstream result((std::string(outputDir) + "/selftest.json").c_str());
  result << "{\"status\":\"passed\",\"checks\":" << passed
         << ",\"physicalDataFit\":false,\"productionEnabled\":false}\n";
  result.close();
  RestartPid::Check(result.good(), "Cannot write test result");
  std::cout << "RESTART_SYNTHETIC_TESTS_PASS checks=" << passed << std::endl;
}

}  // namespace RestartPidTest
#endif

void TestRestartSlices(const char* outputDir) {
#ifndef __CINT__
  try { RestartPidTest::Run(outputDir); }
  catch (const std::exception& error) {
    std::cerr << "RESTART_SYNTHETIC_TESTS_FAILED: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
