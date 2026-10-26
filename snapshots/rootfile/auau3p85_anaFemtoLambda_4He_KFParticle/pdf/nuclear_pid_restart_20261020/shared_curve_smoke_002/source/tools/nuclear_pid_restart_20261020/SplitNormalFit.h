#ifndef STAR_RESTART_SPLIT_NORMAL_FIT_H
#define STAR_RESTART_SPLIT_NORMAL_FIT_H
// Study-only, single-target split-normal response. No PID curve adoption.
// mu is the common mode. Component::sigma is sqrt(sigmaLeft*sigmaRight),
// NOT the distribution RMS. ratio=sigmaRight/sigmaLeft. The density and its
// first derivative are continuous at mu; its second derivative need not be.
// Native or explicitly storage-mapped adjacent bin edges are accepted unchanged.
#include "PeakResponseFit.h"

namespace RestartSplitNormal {
using NuclearFit::Component;
using NuclearFit::Matrix;
using NuclearFit::Request;
using NuclearFit::Settings;
using NuclearFit::Vec;

inline void Widths(double sigma, double ratio, double& left, double& right) {
  NuclearFit::Require(NuclearFit::Finite(sigma) && sigma > 0. &&
                      NuclearFit::Finite(ratio) && ratio > 0.,
                      "Split-normal geometric width and ratio must be positive/finite");
  const double root = std::sqrt(ratio);
  left = sigma / root;
  right = sigma * root;
  NuclearFit::Require(NuclearFit::Finite(left) && left > 0. &&
                      NuclearFit::Finite(right) && right > 0. &&
                      NuclearFit::Finite(left + right), "Invalid split-normal side widths");
}

inline double Density(double x, double mean, double sigma, double ratio) {
  NuclearFit::Require(NuclearFit::Finite(mean) && x == x, "Invalid split-normal coordinate/mode");
  double left, right;
  Widths(sigma, ratio, left, right);
  const double side = x < mean ? left : right;
  const double z = (x - mean) / side;
  return std::sqrt(2. / TMath::Pi()) / (left + right) * std::exp(-.5 * z * z);
}

inline double DensityDerivative(double x, double mean, double sigma, double ratio) {
  NuclearFit::Require(NuclearFit::Finite(x), "Split-normal derivative needs a finite coordinate");
  double left, right;
  Widths(sigma, ratio, left, right);
  const double side = x < mean ? left : right;
  return -(x - mean) / (side * side) * Density(x, mean, sigma, ratio);
}

inline double Integral(double lo, double hi, double mean, double sigma, double ratio) {
  NuclearFit::Require(NuclearFit::Finite(mean) && lo == lo && hi == hi,
                      "Invalid split-normal integration endpoints/mode");
  double left, right;
  Widths(sigma, ratio, left, right);
  if (!(hi > lo)) return 0.;
  if (ratio == 1.) return NuclearFit::GaussianIntegral(lo, hi, mean, sigma);
  double value = 0.;
  if (lo < mean)
    value += 2. * left / (left + right) *
      NuclearFit::GaussianIntegral(lo, std::min(hi, mean), mean, left);
  if (hi > mean)
    value += 2. * right / (left + right) *
      NuclearFit::GaussianIntegral(std::max(lo, mean), hi, mean, right);
  return value;
}

inline double Cdf(double x, double mean, double sigma, double ratio) {
  return Integral(-std::numeric_limits<double>::infinity(), x, mean, sigma, ratio);
}

inline double KernelIntegral(const Component& component, double lo, double hi, double ratio) {
  if (!component.active) return 0.;
  if (component.nodes.empty())
    return Integral(lo, hi, component.mean, component.sigma, ratio);
  double value = 0.;
  for (unsigned i = 0; i < component.nodes.size(); ++i)
    value += component.nodes[i].weight *
      Integral(lo, hi, component.nodes[i].mean, component.nodes[i].sigma, ratio);
  return value;
}

inline Matrix BuildTemplates(const Request& request, const Settings& cfg,
                             double ratio, Vec* acceptance = 0) {
  if (ratio == 1.) return NuclearFit::BuildTemplates(request, cfg, acceptance);
  if (request.bins.empty()) return Matrix();
  const double lo = request.bins.front().lo, hi = request.bins.back().hi;
  const unsigned n = request.components.size(), bg = request.backgroundModel == "linear" ? 2 : 1;
  Matrix result(request.bins.size(), Vec(n + bg, 0.));
  Vec accepted(n, 0.);
  for (unsigned j = 0; j < n; ++j) if (request.components[j].active) {
    accepted[j] = KernelIntegral(request.components[j], lo, hi, ratio);
    if (!(accepted[j] > NuclearFit::Number(cfg, "minimumGaussianAcceptance"))) return Matrix();
    for (unsigned b = 0; b < result.size(); ++b)
      result[b][j] = KernelIntegral(request.components[j], request.bins[b].lo,
                                    request.bins[b].hi, ratio) / accepted[j];
  }
  for (unsigned b = 0; b < result.size(); ++b) for (unsigned j = 0; j < bg; ++j)
    result[b][n + j] = NuclearFit::BackgroundIntegral(request.bins[b].lo,
      request.bins[b].hi, lo, hi, request.backgroundModel, request.backgroundSlope, j);
  for (unsigned j = 0; j < n + bg; ++j) {
    double sum = 0.;
    for (unsigned b = 0; b < result.size(); ++b) {
      if (!NuclearFit::Finite(result[b][j]) || result[b][j] < 0.) return Matrix();
      sum += result[b][j];
    }
    if (sum > 0.) for (unsigned b = 0; b < result.size(); ++b) result[b][j] /= sum;
  }
  if (acceptance) *acceptance = accepted;
  return result;
}

struct Result {
  NuclearFit::Result fit;
  double ratio, ratioError, sigmaLeft, sigmaRight, sigmaLeftError, sigmaRightError;
  int freeParameters;
  bool optimizerDiagnosticsAvailable, minimizeReturned, minimizeConverged;
  int migradStatus, hesseStatus;
  double migradEdm;
  bool hesseAttempted, hesseSucceeded;
  Result() : ratio(1.), ratioError(-1.), sigmaLeft(-1.), sigmaRight(-1.),
    sigmaLeftError(-1.), sigmaRightError(-1.), freeParameters(0),
    optimizerDiagnosticsAvailable(false), minimizeReturned(false), minimizeConverged(false),
    migradStatus(-1), hesseStatus(-1), migradEdm(NuclearFit::invalid),
    hesseAttempted(false), hesseSucceeded(false) {}
};

inline void ValidateSettings(const Settings& cfg) {
  NuclearFit::ValidateSettings(cfg);
  const double lo = NuclearFit::Number(cfg, "splitRatioMin");
  const double hi = NuclearFit::Number(cfg, "splitRatioMax");
  NuclearFit::Require(lo > 0. && hi > lo &&
                      NuclearFit::Number(cfg, "splitRatioLogStep") > 0.,
                      "Invalid explicit split-normal ratio support/log step");
}

inline int PhaseStatus(const Result& result) {
  if (result.migradStatus != 0) return result.migradStatus;
  if (!result.minimizeConverged) return -10;
  if (!result.fit.pointAccepted) return -12;
  if (!result.hesseAttempted || !result.hesseSucceeded)
    return result.hesseStatus > 0 ? result.hesseStatus : -11;
  return result.hesseStatus;
}

class Objective {
 public:
  NuclearFit::Objective base;
  double ratioStart;
  Objective(const Request& request, const Settings& cfg, double start) :
    base(request, cfg), ratioStart(start) {}
  Vec Point() const {
    Vec point = base.Point();
    point.push_back(std::log(ratioStart));
    return point;
  }
  double Ratio(const double* point) const { return std::exp(point[base.kind.size()]); }
  Request Model(const double* point) const {
    Request request = base.Model(point);
    for (unsigned j = 0; j < request.components.size(); ++j)
      if (!base.original.components[j].active) request.components[j] = base.original.components[j];
    return request;
  }
  bool Feasible(const Request& request, double ratio) const {
    if (!base.Feasible(request) || !NuclearFit::Finite(ratio) ||
        ratio < NuclearFit::Number(base.cfg, "splitRatioMin") ||
        ratio > NuclearFit::Number(base.cfg, "splitRatioMax")) return false;
    for (unsigned j = 0; j < request.components.size(); ++j) if (request.components[j].active) {
      double left, right;
      Widths(request.components[j].sigma, ratio, left, right);
      for (unsigned k = 0; k < request.components[j].nodes.size(); ++k)
        Widths(request.components[j].nodes[k].sigma, ratio, left, right);
    }
    return true;
  }
  double operator()(const double* point) const {
    const Request request = Model(point);
    const double ratio = Ratio(point);
    if (!Feasible(request, ratio)) return NuclearFit::invalid;
    const Matrix templates = BuildTemplates(request, base.cfg, ratio);
    const NuclearFit::ProfileResult profile = NuclearFit::Profile(request.bins,
      templates, base.initial, base.cfg);
    return profile.converged ? profile.dataDeviance : NuclearFit::invalid;
  }
};

inline void ParameterLimits(const Objective& objective, unsigned index,
                            double& low, double& high, double& step, std::string& name) {
  const Settings& cfg = objective.base.cfg;
  if (index == objective.base.kind.size()) {
    low = std::log(NuclearFit::Number(cfg, "splitRatioMin"));
    high = std::log(NuclearFit::Number(cfg, "splitRatioMax"));
    step = NuclearFit::Number(cfg, "splitRatioLogStep");
    name = "split_logRatio";
    return;
  }
  const int kind = objective.base.kind[index];
  name = objective.base.names[index];
  low = kind == 0 ? objective.base.original.nativeAxisMin - NuclearFit::Number(cfg, "meanAxisPadding") :
    (kind == 1 ? std::log(NuclearFit::Number(cfg, "minimumSigma")) : -NuclearFit::Number(cfg, "exponentialSlopeLimit"));
  high = kind == 0 ? objective.base.original.nativeAxisMax + NuclearFit::Number(cfg, "meanAxisPadding") :
    (kind == 1 ? std::log(NuclearFit::Number(cfg, "maximumSigma")) : NuclearFit::Number(cfg, "exponentialSlopeLimit"));
  step = NuclearFit::Number(cfg, kind == 0 ? "meanStep" :
    (kind == 1 ? "sigmaLogStep" : "backgroundSlopeStep"));
}

inline Result Fit(const Request& request, const Settings& cfg, double ratioStart) {
  ValidateSettings(cfg);
  const unsigned target = RestartResponse::ValidateRequest(request);
  Result result;
  result.ratio = ratioStart;
  NuclearFit::Require(NuclearFit::Finite(ratioStart) &&
    ratioStart >= NuclearFit::Number(cfg, "splitRatioMin") &&
    ratioStart <= NuclearFit::Number(cfg, "splitRatioMax"), "Split-normal ratio start outside explicit support");
  if (!(NuclearFit::TotalCount(request.bins) > 0.)) { result.fit.status = -2; return result; }
  Objective objective(request, cfg, ratioStart);
  NuclearFit::Require(objective.Feasible(request, ratioStart), "Split-normal seed outside mathematical support");
  Vec point = objective.Point(), chosen = point;
  const double start = objective(&point[0]);
  NuclearFit::Require(start < NuclearFit::invalid, "Split-normal seed template/amplitude profile failed");
  NuclearFit::Require(gSystem->Load("libMinuit2") >= 0, "Cannot load split-normal ROOT Minuit2");
  std::unique_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  NuclearFit::Require(minimizer.get() != 0, "Cannot create split-normal Minuit2");
  minimizer->SetMaxFunctionCalls(NuclearFit::Integer(cfg, "shapeMaxCalls"));
  minimizer->SetTolerance(NuclearFit::Number(cfg, "shapeTolerance"));
  minimizer->SetPrintLevel(-1);
  minimizer->SetErrorDef(1.);
  ROOT::Math::Functor functor(objective, point.size());
  minimizer->SetFunction(functor);
  for (unsigned j = 0; j < point.size(); ++j) {
    double low, high, step; std::string name;
    ParameterLimits(objective, j, low, high, step, name);
    NuclearFit::Require(minimizer->SetLimitedVariable(j, name, point[j], step, low, high),
                        "Cannot register split-normal fit parameter");
  }
  result.optimizerDiagnosticsAvailable = true;
  result.minimizeReturned = minimizer->Minimize();
  result.migradStatus = minimizer->Status();
  result.migradEdm = minimizer->Edm();
  result.minimizeConverged = RestartResponse::MigradConverged(result.minimizeReturned,
    result.migradStatus, result.migradEdm);
  Vec candidate(minimizer->X(), minimizer->X() + point.size());
  const double next = objective(&candidate[0]);
  NuclearFit::Result& fit = result.fit;
  fit.pointAccepted = next < NuclearFit::invalid && next <= start;
  chosen = fit.pointAccepted ? candidate : point;
  if (result.minimizeConverged && fit.pointAccepted) {
    result.hesseAttempted = true;
    result.hesseSucceeded = minimizer->Hesse();
    result.hesseStatus = minimizer->Status();
  }
  fit.status = PhaseStatus(result);
  fit.edm = result.migradEdm;
  fit.covarianceStatus = minimizer->CovMatrixStatus();
  bool covarianceAtChosen = result.minimizeConverged && fit.pointAccepted;
  for (unsigned j = 0; j < chosen.size(); ++j)
    covarianceAtChosen = covarianceAtChosen && chosen[j] == minimizer->X()[j];
  if (covarianceAtChosen && result.hesseAttempted && result.hesseSucceeded &&
      result.hesseStatus == 0 && fit.covarianceStatus == 3) {
    fit.shapeCovariance.assign(chosen.size(), Vec(chosen.size(), 0.));
    for (unsigned j = 0; j < chosen.size(); ++j) for (unsigned k = 0; k < chosen.size(); ++k)
      fit.shapeCovariance[j][k] = minimizer->CovMatrix(j, k);
    fit.shapeCovarianceValid = RestartResponse::PositiveDefinite(fit.shapeCovariance);
  }
  if (!fit.shapeCovarianceValid) fit.shapeCovariance.clear();
  const Request fitted = objective.Model(&chosen[0]);
  result.ratio = objective.Ratio(&chosen[0]);
  Widths(fitted.components[target].sigma, result.ratio, result.sigmaLeft, result.sigmaRight);
  result.freeParameters = chosen.size();
  fit.components = fitted.components;
  fit.shapePoint = chosen;
  fit.shapeNames = objective.base.names;
  fit.shapeNames.push_back("split_logRatio");
  fit.backgroundSlope = fitted.backgroundSlope;
  const Matrix templates = BuildTemplates(fitted, cfg, result.ratio, &fit.windowAcceptance);
  const NuclearFit::ProfileResult profile = NuclearFit::Profile(fitted.bins, templates, objective.base.initial, cfg);
  fit.dataDeviance = profile.dataDeviance;
  fit.amplitudeKkt = profile.kkt;
  fit.amplitudeConverged = profile.converged;
  fit.amplitudeIterations = profile.iterations;
  fit.totalTemplates = templates;
  const double count = NuclearFit::TotalCount(fitted.bins);
  const unsigned n = fitted.components.size();
  fit.meanError.assign(n, -1.);
  fit.sigmaError.assign(n, -1.);
  for (unsigned j = 0; j < n; ++j) fit.columnNames.push_back(fitted.components[j].name);
  if (fitted.backgroundModel == "linear") {
    fit.columnNames.push_back("BG_left"); fit.columnNames.push_back("BG_right");
  } else fit.columnNames.push_back("BG_exp");
  for (unsigned j = 0; j < profile.fractions.size(); ++j)
    fit.columnWindowYields.push_back(count * profile.fractions[j]);
  if (fit.shapeCovarianceValid) {
    int sigmaIndex = -1;
    for (unsigned j = 0; j < objective.base.kind.size(); ++j) {
      const int species = objective.base.species[j];
      if (species >= 0) {
        const double error = std::sqrt(fit.shapeCovariance[j][j]);
        if (objective.base.kind[j] == 0) fit.meanError[species] = error;
        else if (objective.base.kind[j] == 1) {
          fit.sigmaError[species] = fit.components[species].sigma * error;
          if (unsigned(species) == target) sigmaIndex = j;
        }
      }
    }
    const unsigned ratioIndex = chosen.size() - 1;
    result.ratioError = result.ratio * std::sqrt(fit.shapeCovariance[ratioIndex][ratioIndex]);
    bool errorsValid = sigmaIndex >= 0;
    if (errorsValid) {
      const double vSigma = fit.shapeCovariance[sigmaIndex][sigmaIndex];
      const double vRatio = fit.shapeCovariance[ratioIndex][ratioIndex];
      const double covariance = fit.shapeCovariance[sigmaIndex][ratioIndex];
      const double leftVariance = vSigma + .25 * vRatio - covariance;
      const double rightVariance = vSigma + .25 * vRatio + covariance;
      errorsValid = leftVariance > 0. && rightVariance > 0.;
      if (errorsValid) {
        result.sigmaLeftError = result.sigmaLeft * std::sqrt(leftVariance);
        result.sigmaRightError = result.sigmaRight * std::sqrt(rightVariance);
      }
    }
    errorsValid = errorsValid && NuclearFit::Finite(result.ratioError) && result.ratioError > 0. &&
      NuclearFit::Finite(result.sigmaLeftError) && result.sigmaLeftError > 0. &&
      NuclearFit::Finite(result.sigmaRightError) && result.sigmaRightError > 0. &&
      NuclearFit::Finite(fit.meanError[target]) && fit.meanError[target] > 0. &&
      NuclearFit::Finite(fit.sigmaError[target]) && fit.sigmaError[target] > 0.;
    if (!errorsValid) {
      fit.shapeCovarianceValid = false; fit.shapeCovariance.clear();
      fit.meanError.assign(n, -1.); fit.sigmaError.assign(n, -1.);
      result.ratioError = result.sigmaLeftError = result.sigmaRightError = -1.;
    }
  }
  fit.componentCounts.assign(fitted.bins.size(), Vec(n, 0.));
  fit.totalCounts.assign(fitted.bins.size(), 0.);
  fit.smoothCounts.assign(fitted.bins.size(), 0.);
  fit.effectiveCount.assign(n, 0.);
  fit.finite = profile.converged && !templates.empty();
  for (unsigned b = 0; b < templates.size(); ++b) {
    for (unsigned j = 0; j < profile.fractions.size(); ++j) {
      const double value = fit.columnWindowYields[j] * templates[b][j];
      fit.totalCounts[b] += value;
      if (j < n) fit.componentCounts[b][j] = value;
      else fit.smoothCounts[b] += value;
    }
    const double mu = fit.totalCounts[b];
    if (!NuclearFit::Finite(mu) || mu < 0. || (mu == 0. && fitted.bins[b].count > 0.)) fit.finite = false;
    if (mu > 0.) {
      fit.pearson += std::pow(fitted.bins[b].count - mu, 2) / mu;
      for (unsigned j = 0; j < n; ++j)
        fit.effectiveCount[j] += fit.componentCounts[b][j] * fit.componentCounts[b][j] / mu;
    }
  }
  fit.finite = fit.finite && NuclearFit::Finite(fit.dataDeviance) &&
    fit.dataDeviance < NuclearFit::invalid && NuclearFit::Finite(fit.pearson);
  std::vector<unsigned> active;
  for (unsigned j = 0; j < profile.fractions.size(); ++j) if (profile.fractions[j] > 1.e-12) active.push_back(j);
  fit.yieldCovarianceColumnActive.assign(profile.fractions.size(), false);
  for (unsigned j = 0; j < active.size(); ++j) fit.yieldCovarianceColumnActive[active[j]] = true;
  Matrix hessian(active.size(), Vec(active.size(), 0.));
  for (unsigned j = 0; j < active.size(); ++j) for (unsigned k = 0; k < active.size(); ++k)
    hessian[j][k] = profile.hessian[active[j]][active[k]];
  fit.fixedShapeYieldCovariance.assign(profile.fractions.size(), Vec(profile.fractions.size(), 0.));
  fit.yieldCovarianceValid = profile.converged && !active.empty();
  for (unsigned k = 0; k < active.size() && fit.yieldCovarianceValid; ++k) {
    Vec unit(active.size(), 0.), inverse; unit[k] = 1.;
    if (!NuclearFit::Solve(hessian, unit, inverse)) { fit.yieldCovarianceValid = false; break; }
    for (unsigned j = 0; j < active.size(); ++j)
      fit.fixedShapeYieldCovariance[active[j]][active[k]] = count * inverse[j];
  }
  Matrix activeCovariance(active.size(), Vec(active.size(), 0.));
  if (fit.yieldCovarianceValid) {
    for (unsigned j = 0; j < active.size(); ++j) for (unsigned k = 0; k < active.size(); ++k)
      activeCovariance[j][k] = fit.fixedShapeYieldCovariance[active[j]][active[k]];
    fit.yieldCovarianceValid = RestartResponse::PositiveDefinite(activeCovariance);
  }
  if (!fit.yieldCovarianceValid) fit.fixedShapeYieldCovariance.clear();
  for (unsigned j = 0; j < chosen.size(); ++j) {
    double low, high, step; std::string name;
    ParameterLimits(objective, j, low, high, step, name);
    fit.atBoundary = fit.atBoundary || std::min(chosen[j] - low, high - chosen[j]) <=
      NuclearFit::Number(cfg, "boundaryTolerance") * std::max(1., high - low);
  }
  fit.conditionalDropDeltaD.assign(n, -1.);
  for (unsigned j = 0; j < n; ++j) if (fitted.components[j].active) {
    Request dropped = fitted; dropped.components[j].active = false;
    const Matrix droppedTemplates = BuildTemplates(dropped, cfg, result.ratio);
    Vec initial = profile.fractions; initial[j] = 0.;
    const NuclearFit::ProfileResult reduced = NuclearFit::Profile(dropped.bins, droppedTemplates, initial, cfg);
    if (reduced.converged) fit.conditionalDropDeltaD[j] = std::max(0., reduced.dataDeviance - fit.dataDeviance);
  }
  fit.conditionalNdf = int(fitted.bins.size()) - int(chosen.size()) - int(active.size());
  return result;
}
}  // namespace RestartSplitNormal
#endif
