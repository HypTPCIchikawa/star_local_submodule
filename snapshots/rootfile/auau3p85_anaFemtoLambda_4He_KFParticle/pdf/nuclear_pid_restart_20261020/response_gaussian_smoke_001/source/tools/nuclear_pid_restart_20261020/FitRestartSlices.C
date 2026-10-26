// Fresh Gaussian peak-location diagnostic. No historical calibration is read.
// The caller supplies original-Git PID values evaluated at native-R nodes.
void FitRestartSlices(const char*, const char*, const char*, const char*, const char*);
#ifndef __CINT__
#include "../nuclear_pid_tof/NuclearFitCommon.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TMath.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace RestartPid {

struct Slice {
  std::string id, mode, target;
  int rBin;
  double fitLo, fitHi, gateLo, gateHi;
};

struct ShapeInput {
  std::vector<NuclearFit::Node> nodes;
  double mean, sigma;
  ShapeInput() : mean(0), sigma(0) {}
};

typedef std::map<std::string, ShapeInput> Shapes;

static std::vector<std::string> Split(const std::string& text, char separator) {
  std::vector<std::string> fields;
  std::string field;
  std::istringstream stream(text);
  while (std::getline(stream, field, separator)) fields.push_back(field);
  return fields;
}

static void Check(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error(message);
}

static double Number(const std::string& value, const std::string& field) {
  char* end = 0;
  double result = std::strtod(value.c_str(), &end);
  Check(end != value.c_str() && end && *end == 0 && TMath::Finite(result),
        "Invalid finite number for " + field + ": " + value);
  return result;
}

static bool Species(const std::string& value) {
  return value == "d" || value == "t" || value == "he3" || value == "he4";
}

static bool Mode(const std::string& value) {
  return value == "ALL" || value == "TAG_D" || value == "TAG_T" || value == "TAG_HE3";
}

static void SafeId(const std::string& value) {
  Check(!value.empty(), "Empty slice id");
  for (size_t i = 0; i < value.size(); ++i) {
    char c = value[i];
    Check((c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') ||
          (c >= 'A' && c <= 'Z') || c == '_' || c == '-', "Unsafe slice id: " + value);
  }
}

static std::vector<Slice> ReadSlices(const char* path) {
  std::ifstream input(path);
  Check(input.good(), std::string("Cannot read slices: ") + path);
  std::string line;
  Check(bool(std::getline(input, line)), "Empty slices table");
  Check(line == "id\tmode\ttarget\trBin\tfitLo\tfitHi\tgateLo\tgateHi",
        "Unexpected slices header");
  std::vector<Slice> result;
  std::map<std::string, bool> ids;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::vector<std::string> row = Split(line, '\t');
    Check(row.size() == 8, "Each slice must have eight columns");
    Slice slice;
    slice.id = row[0]; slice.mode = row[1]; slice.target = row[2];
    SafeId(slice.id);
    Check(!ids[slice.id], "Duplicate slice id: " + slice.id);
    ids[slice.id] = true;
    Check(Mode(slice.mode) && Species(slice.target), "Unknown input mode or species");
    Check(slice.mode == "ALL" ||
          ((slice.target == "d" || slice.target == "he4") && slice.mode == "TAG_D") ||
          (slice.target == "t" && slice.mode == "TAG_T") ||
          (slice.target == "he3" && slice.mode == "TAG_HE3"),
          "TOF mode must match the species; He4 must use the same TAG_D sample as d");
    double rBin = Number(row[3], "rBin");
    Check(rBin >= 1 && rBin <= INT_MAX && rBin == std::floor(rBin), "rBin must be a positive integer");
    slice.rBin = int(rBin);
    slice.fitLo = Number(row[4], "fitLo"); slice.fitHi = Number(row[5], "fitHi");
    slice.gateLo = Number(row[6], "gateLo"); slice.gateHi = Number(row[7], "gateHi");
    Check(slice.fitLo < slice.fitHi && slice.gateLo < slice.gateHi, "Invalid fit or gate range");
    result.push_back(slice);
  }
  Check(!result.empty(), "No requested slices");
  return result;
}

static std::map<std::string, Shapes> ReadNodes(const char* path) {
  std::ifstream input(path);
  Check(input.good(), std::string("Cannot read shape nodes: ") + path);
  std::string line;
  Check(bool(std::getline(input, line)), "Empty node table");
  Check(line == "id\tspecies\tweight\tmean\tsigma", "Unexpected shape-node header");
  std::map<std::string, Shapes> result;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::vector<std::string> row = Split(line, '\t');
    Check(row.size() == 5, "Each shape node must have five columns");
    SafeId(row[0]);
    Check(Species(row[1]), "Unknown node species");
    double weight = Number(row[2], "weight");
    double mean = Number(row[3], "mean");
    double sigma = Number(row[4], "sigma");
    Check(weight > 0 && sigma > 0, "Node weights and widths must be positive");
    result[row[0]][row[1]].nodes.push_back(NuclearFit::Node(weight, mean, sigma));
  }
  const char* species[] = {"d", "t", "he3", "he4"};
  for (std::map<std::string, Shapes>::iterator row = result.begin(); row != result.end(); ++row) {
    Check(row->second.size() == 4, "Every slice must provide all four original PID curves");
    for (int s = 0; s < 4; ++s) {
      ShapeInput& shape = row->second[species[s]];
      Check(!shape.nodes.empty(), "Missing species nodes");
      double sumWeight = 0;
      for (size_t n = 0; n < shape.nodes.size(); ++n) {
        sumWeight += shape.nodes[n].weight;
        shape.mean += shape.nodes[n].weight * shape.nodes[n].mean;
        shape.sigma += shape.nodes[n].weight * shape.nodes[n].sigma;
      }
      Check(sumWeight > 0, "Invalid total node weight");
      shape.mean /= sumWeight; shape.sigma /= sumWeight;
      for (size_t n = 0; n < shape.nodes.size(); ++n) shape.nodes[n].weight /= sumWeight;
    }
  }
  return result;
}

static TH2* InputHistogram(TFile& input, const NuclearFit::Settings& settings,
                          const std::string& mode) {
  const std::string key = "hist." + mode;
  const std::string path = NuclearFit::Text(settings, key.c_str());
  TH2* histogram = dynamic_cast<TH2*>(input.Get(path.c_str()));
  Check(histogram != 0, "Missing native TH2: " + path);
  Check(histogram->GetXaxis()->GetNbins() == NuclearFit::Integer(settings, "expectedRBins") &&
        histogram->GetYaxis()->GetNbins() == NuclearFit::Integer(settings, "expectedDedxBins"),
        "Unexpected native histogram bin counts");
  Check(std::fabs(histogram->GetXaxis()->GetXmin() - NuclearFit::Number(settings, "expectedRMin")) < 1e-10 &&
        std::fabs(histogram->GetXaxis()->GetXmax() - NuclearFit::Number(settings, "expectedRMax")) < 1e-10 &&
        std::fabs(histogram->GetYaxis()->GetXmin() - NuclearFit::Number(settings, "expectedDedxMin")) < 1e-10 &&
        std::fabs(histogram->GetYaxis()->GetXmax() - NuclearFit::Number(settings, "expectedDedxMax")) < 1e-10,
        "Unexpected native histogram axis limits");
  return histogram;
}

static std::vector<NuclearFit::Bin> NativeBins(const TH1& spectrum, const Slice& slice) {
  const TAxis* axis = spectrum.GetXaxis();
  Check(slice.fitLo >= axis->GetXmin() && slice.fitHi <= axis->GetXmax(), "Fit range outside native axis");
  std::vector<NuclearFit::Bin> result;
  bool foundLo = false, foundHi = false;
  for (int bin = 1; bin <= axis->GetNbins(); ++bin) {
    double lo = axis->GetBinLowEdge(bin), hi = axis->GetBinUpEdge(bin);
    if (std::fabs(lo - slice.fitLo) < 1e-8) foundLo = true;
    if (std::fabs(hi - slice.fitHi) < 1e-8) foundHi = true;
    if (lo >= slice.fitLo - 1e-8 && hi <= slice.fitHi + 1e-8) {
      double count = spectrum.GetBinContent(bin);
      Check(TMath::Finite(count) && count >= 0, "Invalid native count");
      result.push_back(NuclearFit::Bin(lo, hi, count));
    }
  }
  Check(foundLo && foundHi && !result.empty(), "Fit edges must coincide with native dEdx bin edges");
  return result;
}

static std::string Quote(const std::string& value) {
  std::ostringstream output;
  output << '"';
  for (size_t i = 0; i < value.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(value[i]);
    if (c == '"' || c == '\\') output << '\\' << char(c);
    else if (c == '\n') output << "\\n";
    else if (c == '\r') output << "\\r";
    else if (c == '\t') output << "\\t";
    else if (c < 32) output << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else output << char(c);
  }
  output << '"';
  return output.str();
}

static void JsonNumber(std::ostream& output, double value) {
  if (TMath::Finite(value)) output << std::setprecision(17) << value;
  else output << "null";
}

static void JsonVector(std::ostream& output, const NuclearFit::Vec& values) {
  output << '[';
  for (size_t i = 0; i < values.size(); ++i) {
    if (i) output << ',';
    JsonNumber(output, values[i]);
  }
  output << ']';
}

static void JsonMatrix(std::ostream& output, const NuclearFit::Matrix& matrix) {
  output << '[';
  for (size_t i = 0; i < matrix.size(); ++i) {
    if (i) output << ',';
    JsonVector(output, matrix[i]);
  }
  output << ']';
}

static std::string HistogramName(const std::string& prefix, size_t index) {
  std::ostringstream name;
  name << prefix << '_' << std::setw(4) << std::setfill('0') << index;
  return name.str();
}

static bool Close(double one, double two, const NuclearFit::Settings& settings) {
  return TMath::Finite(one) && TMath::Finite(two) &&
    std::fabs(one - two) <= NuclearFit::Number(settings, "readbackToleranceAbsolute") +
    NuclearFit::Number(settings, "readbackToleranceRelative") * std::max(std::fabs(one), std::fabs(two));
}

static NuclearFit::Request Request(const TH1& native, const Slice& slice, const Shapes& shapes,
                                  const NuclearFit::Settings& settings) {
  NuclearFit::Request request;
  request.bins = NativeBins(native, slice);
  request.nativeAxisMin = native.GetXaxis()->GetXmin();
  request.nativeAxisMax = native.GetXaxis()->GetXmax();
  request.backgroundModel = NuclearFit::Text(settings, "backgroundModel");
  Check(request.backgroundModel == "linear", "Fresh peak diagnostic requires linear background");
  request.backgroundSlope = NuclearFit::Number(settings, "backgroundSlope");
  int freeSlope = NuclearFit::Integer(settings, "freeBackgroundSlope");
  Check(freeSlope == 0 || freeSlope == 1, "freeBackgroundSlope must be 0 or 1");
  request.freeBackgroundSlope = freeSlope != 0;
  double fraction = NuclearFit::Number(settings, "startTargetFraction");
  Check(fraction > 0 && fraction < 1, "startTargetFraction must be between zero and one");
  double total = 0;
  for (size_t b = 0; b < request.bins.size(); ++b) total += request.bins[b].count;
  request.backgroundWindowYield = total * (1 - fraction);
  const char* species[] = {"d", "t", "he3", "he4"};
  for (int s = 0; s < 4; ++s) {
    Shapes::const_iterator shape = shapes.find(species[s]);
    Check(shape != shapes.end(), "Missing original PID species");
    NuclearFit::Component component;
    component.name = species[s];
    component.mean = shape->second.mean; component.sigma = shape->second.sigma;
    component.nodes = shape->second.nodes;
    component.active = slice.target == species[s];
    component.freeMean = component.active;
    component.freeSigma = component.active;
    component.startWindowYield = component.active ? total * fraction : 0;
    request.components.push_back(component);
  }
  return request;
}

static void SavePrediction(TFile& output, const TH1D& native, const Slice& slice,
                           const NuclearFit::Result& result, size_t index,
                           const NuclearFit::Settings& settings) {
  TH1D total(native), smooth(native), mask(native);
  total.SetName(HistogramName("prediction_total", index).c_str());
  smooth.SetName(HistogramName("prediction_smooth", index).c_str());
  mask.SetName(HistogramName("fit_mask", index).c_str());
  total.Reset(); smooth.Reset(); mask.Reset();
  total.SetDirectory(0); smooth.SetDirectory(0); mask.SetDirectory(0);
  std::vector<TH1D*> components;
  for (size_t s = 0; s < result.components.size(); ++s) {
    TH1D* component = new TH1D(native);
    component->SetName(HistogramName("prediction_" + result.components[s].name, index).c_str());
    component->Reset(); component->SetDirectory(0);
    components.push_back(component);
  }
  std::vector<NuclearFit::Bin> bins = NativeBins(native, slice);
  Check(result.totalCounts.size() == bins.size() && result.smoothCounts.size() == bins.size() &&
        result.componentCounts.size() == bins.size(), "Unexpected prediction dimension");
  for (size_t b = 0; b < bins.size(); ++b) {
    int nativeBin = native.GetXaxis()->FindFixBin(0.5 * (bins[b].lo + bins[b].hi));
    Check(result.componentCounts[b].size() == components.size(), "Unexpected component prediction dimension");
    double sum = result.smoothCounts[b];
    Check(TMath::Finite(sum) && sum >= 0, "Invalid smooth-background prediction");
    for (size_t s = 0; s < components.size(); ++s) {
      double count = result.componentCounts[b][s];
      Check(TMath::Finite(count) && count >= 0, "Invalid isotope prediction");
      components[s]->SetBinContent(nativeBin, count);
      sum += count;
    }
    Check(Close(sum, result.totalCounts[b], settings), "Prediction components do not add to total");
    total.SetBinContent(nativeBin, result.totalCounts[b]);
    smooth.SetBinContent(nativeBin, result.smoothCounts[b]);
    mask.SetBinContent(nativeBin, 1);
  }
  output.cd(); total.Write(); smooth.Write(); mask.Write();
  for (size_t s = 0; s < components.size(); ++s) {
    components[s]->Write(); delete components[s];
  }
}

static void WriteResult(std::ostream& out, const Slice& slice, const Shapes& initial,
                        const NuclearFit::Request& request, const NuclearFit::Result& result,
                        const TH1D& native, size_t index, double rLo, double rHi) {
  double fullCounts = native.Integral(1, native.GetNbinsX());
  double fitCounts = 0;
  for (size_t b = 0; b < request.bins.size(); ++b) fitCounts += request.bins[b].count;
  bool foreignPeak = false;
  for (Shapes::const_iterator s = initial.begin(); s != initial.end(); ++s) {
    if (s->first != slice.target && s->second.mean >= slice.fitLo && s->second.mean <= slice.fitHi)
      foreignPeak = true;
  }
  out << "{\"id\":" << Quote(slice.id) << ",\"mode\":" << Quote(slice.mode)
      << ",\"target\":" << Quote(slice.target) << ",\"rBin\":" << slice.rBin
      << ",\"rLo\":" << std::setprecision(17) << rLo << ",\"rHi\":" << rHi
      << ",\"fitLo\":" << slice.fitLo << ",\"fitHi\":" << slice.fitHi
      << ",\"gateLo\":" << slice.gateLo << ",\"gateHi\":" << slice.gateHi
      << ",\"gateSource\":\"initial_git_curve_diagnostic_only\""
      << ",\"nativeHistogram\":" << Quote(HistogramName("raw_native", index))
      << ",\"totalHistogram\":" << Quote(HistogramName("prediction_total", index))
      << ",\"smoothHistogram\":" << Quote(HistogramName("prediction_smooth", index))
      << ",\"fitMaskHistogram\":" << Quote(HistogramName("fit_mask", index))
      << ",\"fitCompleted\":true,\"fullNativeCounts\":" << fullCounts
      << ",\"fitNativeCounts\":" << fitCounts << ",\"nFitBins\":" << request.bins.size()
      << ",\"status\":" << result.status << ",\"covarianceStatus\":" << result.covarianceStatus
      << ",\"finite\":" << (result.finite ? "true" : "false")
      << ",\"pointAccepted\":" << (result.pointAccepted ? "true" : "false")
      << ",\"amplitudeConverged\":" << (result.amplitudeConverged ? "true" : "false")
      << ",\"atBoundary\":" << (result.atBoundary ? "true" : "false")
      << ",\"yieldCovarianceValid\":" << (result.yieldCovarianceValid ? "true" : "false")
      << ",\"shapeCovarianceValid\":" << (result.shapeCovarianceValid ? "true" : "false")
      << ",\"foreignPeakInsideFit\":" << (foreignPeak ? "true" : "false")
      << ",\"curveAdoptionAllowed\":false,\"purityAvailable\":false,\"purity\":null"
      << ",\"engineConditionalNdf\":" << result.conditionalNdf
      << ",\"fittedShapeParameterCount\":" << result.shapePoint.size()
      << ",\"pearson\":";
  JsonNumber(out, result.pearson);
  out << ",\"pearsonPerNdf\":";
  if (result.conditionalNdf > 0) JsonNumber(out, result.pearson / result.conditionalNdf);
  else out << "null";
  out << ",\"dataDeviance\":"; JsonNumber(out, result.dataDeviance);
  out << ",\"edm\":"; JsonNumber(out, result.edm);
  out << ",\"amplitudeKkt\":"; JsonNumber(out, result.amplitudeKkt);
  out << ",\"backgroundSlope\":"; JsonNumber(out, result.backgroundSlope);
  out << ",\"components\":[";
  for (size_t s = 0; s < result.components.size(); ++s) {
    if (s) out << ',';
    const NuclearFit::Component& component = result.components[s];
    Shapes::const_iterator seed = initial.find(component.name);
    Check(seed != initial.end(), "Unknown fitted species");
    out << "{\"species\":" << Quote(component.name)
        << ",\"active\":" << (component.active ? "true" : "false")
        << ",\"freeMean\":" << (component.freeMean ? "true" : "false")
        << ",\"freeSigma\":" << (component.freeSigma ? "true" : "false")
        << ",\"initialMean\":" << seed->second.mean << ",\"initialSigma\":" << seed->second.sigma
        << ",\"mean\":";
    JsonNumber(out, component.mean);
    out << ",\"sigma\":"; JsonNumber(out, component.sigma);
    out << ",\"meanError\":";
    if (s < result.meanError.size()) JsonNumber(out, result.meanError[s]); else out << "null";
    out << ",\"sigmaError\":";
    if (s < result.sigmaError.size()) JsonNumber(out, result.sigmaError[s]); else out << "null";
    out << ",\"histogram\":" << Quote(HistogramName("prediction_" + component.name, index)) << '}';
  }
  out << "],\"columnNames\":[";
  for (size_t column = 0; column < result.columnNames.size(); ++column) {
    if (column) out << ',';
    out << Quote(result.columnNames[column]);
  }
  out << "],\"shapeNames\":[";
  for (size_t shape = 0; shape < result.shapeNames.size(); ++shape) {
    if (shape) out << ',';
    out << Quote(result.shapeNames[shape]);
  }
  out << "],\"yieldCovarianceColumnActive\":[";
  for (size_t column = 0; column < result.yieldCovarianceColumnActive.size(); ++column) {
    if (column) out << ',';
    out << (result.yieldCovarianceColumnActive[column] ? "true" : "false");
  }
  out << "],\"columnWindowYields\":"; JsonVector(out, result.columnWindowYields);
  out << ",\"windowAcceptance\":"; JsonVector(out, result.windowAcceptance);
  out << ",\"effectiveCount\":"; JsonVector(out, result.effectiveCount);
  out << ",\"shapePoint\":"; JsonVector(out, result.shapePoint);
  out << ",\"fixedShapeYieldCovariance\":"; JsonMatrix(out, result.fixedShapeYieldCovariance);
  out << ",\"shapeCovariance\":"; JsonMatrix(out, result.shapeCovariance);
  out << '}';
}

static void Readback(const char* rootPath, TFile& input, const std::vector<Slice>& slices,
                     const NuclearFit::Settings& settings, const std::vector<bool>& fitted) {
  TFile saved(rootPath, "READ");
  Check(!saved.IsZombie(), "Cannot reopen native spectra output");
  const char* species[] = {"d", "t", "he3", "he4"};
  for (size_t index = 0; index < slices.size(); ++index) {
    TH2* source = InputHistogram(input, settings, slices[index].mode);
    TH1D* raw = dynamic_cast<TH1D*>(saved.Get(HistogramName("raw_native", index).c_str()));
    Check(raw && raw->GetNbinsX() == source->GetNbinsY(), "Raw spectrum missing or rebinned");
    for (int b = 0; b <= raw->GetNbinsX() + 1; ++b) {
      Check(Close(raw->GetBinContent(b), source->GetBinContent(slices[index].rBin, b), settings),
            "Raw readback content mismatch");
      Check(Close(raw->GetBinError(b), source->GetBinError(slices[index].rBin, b), settings),
            "Raw readback error mismatch");
      if (b >= 1 && b <= raw->GetNbinsX())
        Check(Close(raw->GetXaxis()->GetBinLowEdge(b), source->GetYaxis()->GetBinLowEdge(b), settings) &&
              Close(raw->GetXaxis()->GetBinUpEdge(b), source->GetYaxis()->GetBinUpEdge(b), settings),
              "Raw readback native axis mismatch");
    }
    if (!fitted[index]) continue;
    TH1* total = dynamic_cast<TH1*>(saved.Get(HistogramName("prediction_total", index).c_str()));
    TH1* smooth = dynamic_cast<TH1*>(saved.Get(HistogramName("prediction_smooth", index).c_str()));
    TH1* mask = dynamic_cast<TH1*>(saved.Get(HistogramName("fit_mask", index).c_str()));
    Check(total && smooth && mask, "Missing saved prediction");
    std::vector<TH1*> components;
    for (int s = 0; s < 4; ++s) {
      TH1* component = dynamic_cast<TH1*>(saved.Get(HistogramName(std::string("prediction_") + species[s], index).c_str()));
      Check(component != 0, "Missing saved isotope prediction");
      components.push_back(component);
    }
    for (int b = 0; b <= raw->GetNbinsX() + 1; ++b) {
      double sum = smooth->GetBinContent(b);
      for (size_t s = 0; s < components.size(); ++s) sum += components[s]->GetBinContent(b);
      Check(Close(sum, total->GetBinContent(b), settings), "Saved prediction is not the component sum");
      double center = raw->GetXaxis()->GetBinCenter(b);
      bool inFit = b >= 1 && b <= raw->GetNbinsX() && center >= slices[index].fitLo && center < slices[index].fitHi;
      Check(mask->GetBinContent(b) == (inFit ? 1 : 0), "Saved native fit-mask mismatch");
    }
  }
}

static void Run(const char* cachePath, const char* outputDir, const char* settingsPath,
                const char* slicesPath, const char* nodesPath) {
  NuclearFit::Settings settings = NuclearFit::ReadSettings(settingsPath);
  NuclearFit::ValidateSettings(settings);
  Check(NuclearFit::Number(settings, "readbackToleranceAbsolute") >= 0 &&
        NuclearFit::Number(settings, "readbackToleranceRelative") >= 0, "Negative readback tolerance");
  std::vector<Slice> slices = ReadSlices(slicesPath);
  std::map<std::string, Shapes> shapes = ReadNodes(nodesPath);
  Check(shapes.size() == slices.size(), "Shape table must match exactly the requested slices");
  const std::string rootPath = std::string(outputDir) + "/native_spectra.root";
  const std::string jsonPath = std::string(outputDir) + "/fits.json";
  Check(gSystem->AccessPathName(rootPath.c_str()) && gSystem->AccessPathName(jsonPath.c_str()),
        "Output already exists; use a fresh study directory");
  Check(gSystem->mkdir(outputDir, true) == 0 || !gSystem->AccessPathName(outputDir), "Cannot create output directory");
  TFile input(cachePath, "READ");
  Check(!input.IsZombie(), "Cannot open immutable native input cache");
  TFile output(rootPath.c_str(), "CREATE");
  Check(!output.IsZombie(), "Cannot create fresh native spectra output");
  std::vector<std::string> resultRows;
  std::vector<bool> fitted;
  size_t failures = 0;
  for (size_t index = 0; index < slices.size(); ++index) {
    const Slice& slice = slices[index];
    Check(shapes.find(slice.id) != shapes.end(), "Missing slice nodes: " + slice.id);
    TH2* histogram = InputHistogram(input, settings, slice.mode);
    Check(slice.rBin <= histogram->GetNbinsX(), "Requested R bin outside native axis");
    TH1D* native = histogram->ProjectionY(HistogramName("raw_native", index).c_str(), slice.rBin, slice.rBin, "e");
    Check(native != 0, "Native projection failed");
    native->SetDirectory(0);
    output.cd(); native->Write();
    std::ostringstream row;
    bool complete = false;
    try {
      NuclearFit::Request request = Request(*native, slice, shapes[slice.id], settings);
      NuclearFit::Result result = NuclearFit::Fit(request, settings);
      SavePrediction(output, *native, slice, result, index, settings);
      WriteResult(row, slice, shapes[slice.id], request, result, *native, index,
                  histogram->GetXaxis()->GetBinLowEdge(slice.rBin), histogram->GetXaxis()->GetBinUpEdge(slice.rBin));
      complete = true;
      std::cout << "RESTART_PEAK_FIT id=" << slice.id << " status=" << result.status
                << " pearson=" << result.pearson << " ndf=" << result.conditionalNdf << std::endl;
    } catch (const std::exception& error) {
      ++failures;
      row.str(""); row.clear();
      row << "{\"id\":" << Quote(slice.id) << ",\"mode\":" << Quote(slice.mode)
          << ",\"target\":" << Quote(slice.target) << ",\"rBin\":" << slice.rBin
          << ",\"fitCompleted\":false,\"purityAvailable\":false,\"purity\":null"
          << ",\"curveAdoptionAllowed\":false,\"nativeHistogram\":"
          << Quote(HistogramName("raw_native", index)) << ",\"error\":" << Quote(error.what()) << '}';
      std::cerr << "RESTART_PEAK_FIT_FAILED id=" << slice.id << " error=" << error.what() << std::endl;
    }
    resultRows.push_back(row.str()); fitted.push_back(complete);
    delete native;
  }
  output.Write(); output.Close();
  Readback(rootPath.c_str(), input, slices, settings, fitted);
  std::ofstream json(jsonPath.c_str());
  Check(json.good(), "Cannot write fit JSON");
  json << "{\"status\":\"completed\",\"study\":\"fresh_gaussian_peak_diagnostic_20261020\","
       << "\"productionEnabled\":false,\"purityAvailable\":false,\"curveAdoptionAllowed\":false,"
       << "\"initialShapeSource\":\"git_log10_polynomial\",\"model\":\"gaussian_linear\","
       << "\"scope\":\"target-only peak calibration diagnostic; per-row TOF mode; not a purity fit\","
       << "\"rawReadback\":true,\"predictionReadback\":true,\"totalSlices\":" << slices.size()
       << ",\"fitFailures\":" << failures << ",\"cachePath\":" << Quote(cachePath) << ",\"results\":[";
  for (size_t index = 0; index < resultRows.size(); ++index) {
    if (index) json << ',';
    json << resultRows[index];
  }
  json << "]}\n";
  json.close();
  Check(json.good(), "Fit JSON write failed");
  std::cout << "RESTART_PEAK_READBACK_OK slices=" << slices.size() << " failures=" << failures << std::endl;
}

}  // namespace RestartPid
#endif

void FitRestartSlices(const char* cachePath, const char* outputDir, const char* settingsPath,
                      const char* slicesPath, const char* nodesPath) {
#ifndef __CINT__
  try {
    RestartPid::Run(cachePath, outputDir, settingsPath, slicesPath, nodesPath);
  } catch (const std::exception& error) {
    std::cerr << "RESTART_PEAK_FATAL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
