// Synthetic response tests only. The explicit fresh JSON is a test report,
// not a PID calibration or data result. No input ROOT or production file changes.
int TestSplitNormalFit(const char*, const char*);
#ifndef __CINT__
#include "SplitNormalFit.h"
#include "Float16BinMapping.h"
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
unsigned SplitChecks = 0;
void SplitCheck(bool condition, const std::string& label) {
  ++SplitChecks;
  if (!condition) throw std::runtime_error(label);
}
void SplitNear(double value, double expected, double absoluteTolerance, const std::string& label) {
  SplitCheck(NuclearFit::Finite(value) && NuclearFit::Finite(expected) &&
              std::fabs(value - expected) <= absoluteTolerance, label);
}

void SplitIntegralTests() {
  const double infinity = std::numeric_limits<double>::infinity();
  const double ratios[] = {.25, .5, .65, 1., 1.6, 2., 4.};
  for (unsigned r = 0; r < sizeof(ratios) / sizeof(ratios[0]); ++r) {
    const double ratio = ratios[r];
    SplitNear(RestartSplitNormal::Integral(-infinity, infinity, 40., 2., ratio),
               1., 2.e-14, "Full split-normal normalization");
    SplitNear(RestartSplitNormal::Cdf(-infinity, 40., 2., ratio), 0., 1.e-14, "CDF lower infinity");
    SplitNear(RestartSplitNormal::Cdf(infinity, 40., 2., ratio), 1., 2.e-14, "CDF upper infinity");
    SplitNear(RestartSplitNormal::Cdf(40., 40., 2., ratio), 1. / (1. + ratio),
               2.e-14, "The mode CDF is the left-width fraction");
    double left, right;
    RestartSplitNormal::Widths(2., ratio, left, right);
    SplitNear(std::sqrt(left * right), 2., 1.e-14, "Geometric-width convention");
    SplitNear(right / left, ratio, 1.e-14, "Side-width ratio convention");
    const double peak = RestartSplitNormal::Density(40., 40., 2., ratio);
    SplitCheck(peak > RestartSplitNormal::Density(39.99, 40., 2., ratio) &&
                peak > RestartSplitNormal::Density(40.01, 40., 2., ratio), "mu is the unique response mode");
    SplitNear(RestartSplitNormal::Density(40. - 1.e-8, 40., 2., ratio),
               RestartSplitNormal::Density(40. + 1.e-8, 40., 2., ratio),
               1.e-12, "Density is continuous at the mode");
    SplitNear(RestartSplitNormal::DensityDerivative(40. - 1.e-8, 40., 2., ratio),
               RestartSplitNormal::DensityDerivative(40. + 1.e-8, 40., 2., ratio),
               1.e-8, "First derivative is continuous and zero at the mode");
    for (int i = 0; i < 240; ++i) {
      const double lo = 10. + .25 * i, mid = lo + .11, hi = lo + .25;
      const double whole = RestartSplitNormal::Integral(lo, hi, 40., 2., ratio);
      SplitCheck(whole >= 0. && whole <= 1., "Nonnegative split-normal bin integral");
      SplitNear(RestartSplitNormal::Integral(lo, mid, 40., 2., ratio) +
                 RestartSplitNormal::Integral(mid, hi, 40., 2., ratio), whole,
                 3.e-15, "Exact integral additivity");
      SplitNear(whole, RestartSplitNormal::Integral(80. - hi, 80. - lo, 40., 2., 1. / ratio),
                 3.e-15, "Left/right widths mirror by reciprocal ratio");
      if (ratio == 1.) SplitNear(whole, NuclearFit::GaussianIntegral(lo, hi, 40., 2.),
                                  0., "Ratio one is bit-exact Gaussian integral");
    }
  }
  bool rejected = false;
  try { RestartSplitNormal::Integral(35., 45., 40., 2., 0.); }
  catch (const std::exception&) { rejected = true; }
  SplitCheck(rejected, "Nonpositive width ratio is rejected");
}

NuclearFit::Request SplitToy(double ratio, bool mapped) {
  NuclearFit::Request request;
  request.nativeAxisMin = 0.; request.nativeAxisMax = 120.;
  request.backgroundModel = "linear"; request.backgroundWindowYield = 4000.;
  const char* names[] = {"d", "t", "he3", "he4"};
  for (unsigned j = 0; j < 4; ++j) {
    NuclearFit::Component component;
    component.name = names[j]; component.mean = 10. * (j + 1); component.sigma = 2.;
    component.active = j == 3; component.freeMean = j == 3; component.freeSigma = j == 3;
    component.startWindowYield = j == 3 ? 20000. : 0.;
    request.components.push_back(component);
  }
  NuclearFit::Component& target = request.components[3];
  target.nodes.push_back(NuclearFit::Node(5. / 18., 39.5, 1.9));
  target.nodes.push_back(NuclearFit::Node(8. / 18., 40., 2.));
  target.nodes.push_back(NuclearFit::Node(5. / 18., 40.5, 2.1));
  std::vector<double> edges;
  for (unsigned i = 0; i <= 220; ++i) edges.push_back(18. + .2 * i);
  if (mapped) edges = RestartFloat16::IntegrationEdges(edges);
  for (unsigned i = 0; i < 220; ++i) request.bins.push_back(NuclearFit::Bin(edges[i], edges[i + 1], 0.));
  const double lo = request.bins.front().lo, hi = request.bins.back().hi;
  const double acceptance = RestartSplitNormal::KernelIntegral(target, lo, hi, ratio);
  for (unsigned b = 0; b < request.bins.size(); ++b) {
    const NuclearFit::Bin& bin = request.bins[b];
    const double probability = RestartSplitNormal::KernelIntegral(target, bin.lo, bin.hi, ratio) / acceptance;
    request.bins[b].count = 20000. * probability +
      3000. * NuclearFit::BackgroundIntegral(bin.lo, bin.hi, lo, hi, "linear", 0., 0) +
      1000. * NuclearFit::BackgroundIntegral(bin.lo, bin.hi, lo, hi, "linear", 0., 1);
  }
  return request;
}

void SplitShiftSeed(NuclearFit::Request& request) {
  NuclearFit::Component& target = request.components[3];
  const double oldMean = target.mean, oldSigma = target.sigma;
  target.mean = 39.7; target.sigma = 2.2; target.startWindowYield = 18000.;
  request.backgroundWindowYield = 4500.;
  for (unsigned q = 0; q < target.nodes.size(); ++q) {
    target.nodes[q].mean += target.mean - oldMean;
    target.nodes[q].sigma *= target.sigma / oldSigma;
  }
}

void SplitFitTests(const NuclearFit::Settings& cfg) {
  NuclearFit::Request gaussian = SplitToy(1., false);
  SplitCheck(RestartSplitNormal::BuildTemplates(gaussian, cfg, 1.) ==
              NuclearFit::BuildTemplates(gaussian, cfg), "Ratio-one template matrix is exactly Gaussian");
  const double ratios[] = {.65, 1., 1.6};
  for (unsigned r = 0; r < 3; ++r) for (unsigned mapped = 0; mapped < 2; ++mapped) {
    const double truthRatio = ratios[r];
    SplitCheck(NuclearFit::Number(cfg, "splitRatioMin") < truthRatio &&
                truthRatio < NuclearFit::Number(cfg, "splitRatioMax"),
                "Synthetic truth lies strictly inside explicit ratio support");
    NuclearFit::Request request = SplitToy(truthRatio, mapped != 0);
    SplitShiftSeed(request);
    const RestartSplitNormal::Result response = RestartSplitNormal::Fit(request, cfg, 1.);
    const NuclearFit::Result& fit = response.fit;
    std::cout << "SPLIT_NORMAL_ASIMOV ratio=" << truthRatio << " mapped=" << mapped
              << " status=" << fit.status << " deviance=" << fit.dataDeviance
              << " fittedRatio=" << response.ratio << std::endl;
    SplitCheck(response.optimizerDiagnosticsAvailable && response.minimizeReturned &&
                response.minimizeConverged && response.migradStatus == 0,
                "Successful MIGRAD phase is recorded");
    SplitCheck(response.hesseAttempted && response.hesseSucceeded && response.hesseStatus == 0 &&
                fit.status == 0, "Successful Hesse phase cannot conceal MIGRAD failure");
    SplitCheck(fit.finite && fit.amplitudeConverged && fit.pointAccepted && !fit.atBoundary,
                "Interior Asimov profile converges without a boundary solution");
    SplitNear(fit.components[3].mean, 40., .005, "Asimov common-mode recovery");
    SplitNear(fit.components[3].sigma, 2., .005, "Asimov geometric-width recovery");
    SplitNear(response.ratio, truthRatio, .005, "Asimov left/right ratio recovery");
    SplitNear(response.sigmaLeft, 2. / std::sqrt(truthRatio), .01, "Asimov left-width recovery");
    SplitNear(response.sigmaRight, 2. * std::sqrt(truthRatio), .01, "Asimov right-width recovery");
    SplitCheck(fit.dataDeviance < .01, "Correct split-normal Asimov model has near-zero deviance");
    SplitCheck(fit.shapeCovarianceValid && fit.yieldCovarianceValid && response.ratioError > 0. &&
                response.sigmaLeftError > 0. && response.sigmaRightError > 0.,
                "Full profile-shape covariance and conditional amplitude covariance are valid");
    SplitCheck(response.freeParameters == 3, "Two widths and a mode add exactly three shape coordinates");
    unsigned activeYields = 0;
    for (unsigned j = 0; j < fit.yieldCovarianceColumnActive.size(); ++j)
      activeYields += fit.yieldCovarianceColumnActive[j];
    SplitCheck(fit.conditionalNdf == int(request.bins.size()) - response.freeParameters - int(activeYields),
                "NDF subtracts actual free shape and interior yield dimensions");
    for (unsigned j = 0; j < 3; ++j) {
      SplitCheck(fit.components[j].mean == request.components[j].mean &&
                  fit.components[j].sigma == request.components[j].sigma &&
                  !fit.components[j].active && !fit.components[j].freeMean && !fit.components[j].freeSigma,
                  "Foreign nucleus curves remain fixed and inactive");
    }
    for (unsigned b = 0; b < fit.totalCounts.size(); ++b) {
      double sum = fit.smoothCounts[b];
      for (unsigned j = 0; j < 4; ++j) sum += fit.componentCounts[b][j];
      SplitNear(sum, fit.totalCounts[b], 1.e-10, "Saved prediction equals actual component sum");
    }
  }
  NuclearFit::Request invalid = gaussian;
  invalid.components[2].freeMean = true;
  bool rejected = false;
  try { RestartSplitNormal::Fit(invalid, cfg, 1.); } catch (const std::exception&) { rejected = true; }
  SplitCheck(rejected, "A foreign nucleus cannot become an arbitrary free component");
  invalid = gaussian;
  invalid.components[2].active = true; invalid.components[2].freeMean = true; invalid.components[2].freeSigma = true;
  rejected = false;
  try { RestartSplitNormal::Fit(invalid, cfg, 1.); } catch (const std::exception&) { rejected = true; }
  SplitCheck(rejected, "A second active target is rejected");
  NuclearFit::Settings limited = cfg;
  limited["shapeMaxCalls"] = "1";
  NuclearFit::Request difficult = SplitToy(.65, false);
  SplitShiftSeed(difficult);
  const RestartSplitNormal::Result unfinished = RestartSplitNormal::Fit(difficult, limited, 1.);
  SplitCheck(!unfinished.minimizeConverged && unfinished.fit.status != 0 &&
              !unfinished.hesseAttempted && !unfinished.fit.shapeCovarianceValid,
              "Call-limited MIGRAD cannot acquire a successful Hesse status or covariance");
}

void SplitPhaseTests() {
  RestartSplitNormal::Result result;
  result.migradStatus = 4; result.minimizeConverged = false;
  result.fit.pointAccepted = true; result.hesseAttempted = true; result.hesseSucceeded = true; result.hesseStatus = 0;
  SplitCheck(RestartSplitNormal::PhaseStatus(result) == 4, "Hesse cannot erase MIGRAD failure");
  result.migradStatus = 0;
  SplitCheck(RestartSplitNormal::PhaseStatus(result) == -10, "Nonconverged status-zero MIGRAD remains a failure");
  result.minimizeConverged = true; result.hesseSucceeded = false;
  SplitCheck(RestartSplitNormal::PhaseStatus(result) == -11, "Unsuccessful Hesse remains a failure");
  result.hesseSucceeded = true;
  SplitCheck(RestartSplitNormal::PhaseStatus(result) == 0, "Both phases must succeed");
  result.fit.pointAccepted = false;
  SplitCheck(RestartSplitNormal::PhaseStatus(result) == -12, "Rejected candidate is never a success");
}
}  // namespace
#endif

int TestSplitNormalFit(const char* settingsPath, const char* reportPath) {
#ifndef __CINT__
  try {
    SplitChecks = 0;
    SplitCheck(settingsPath && settingsPath[0] && reportPath && reportPath[0], "Explicit settings/report required");
    SplitCheck(gSystem->AccessPathName(reportPath), "Synthetic report path must be fresh");
    const NuclearFit::Settings cfg = NuclearFit::ReadSettings(settingsPath);
    RestartSplitNormal::ValidateSettings(cfg);
    SplitIntegralTests(); SplitPhaseTests(); SplitFitTests(cfg);
    std::ofstream out(reportPath);
    SplitCheck(out.good(), "Create fresh synthetic report");
    out << "{\n  \"status\": \"passed\",\n  \"checks\": " << (SplitChecks + 1)
        << ",\n  \"syntheticOnly\": true,\n  \"model\": \"split_normal\","
        << "\n  \"sigmaMeaning\": \"geometric mean of left/right widths, not RMS\","
        << "\n  \"productionEnabled\": false,\n  \"coefficientsUpdated\": false\n}\n";
    out.close();
    SplitCheck(out.good(), "Write synthetic report");
    std::cout << "SPLIT_NORMAL_TEST_PASS checks=" << SplitChecks << std::endl;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "SPLIT_NORMAL_TEST_FAIL " << error.what() << std::endl;
    return 1;
  }
#else
  return 1;
#endif
}
