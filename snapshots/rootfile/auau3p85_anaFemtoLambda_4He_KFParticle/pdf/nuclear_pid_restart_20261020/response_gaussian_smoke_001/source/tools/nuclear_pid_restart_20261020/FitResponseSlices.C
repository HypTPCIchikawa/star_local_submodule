// Alternative peak-response diagnostics. Existing Gaussian I/O stays unchanged.
void FitResponseSlices(const char*, const char*, const char*, const char*, const char*);
#ifndef __CINT__
#include "FitRestartSlices.C"
#include "PeakResponseFit.h"

namespace RestartResponseSlices {

static void AppendResponse(std::ostream& output, const std::string& original,
                           const RestartResponse::Result& result) {
  RestartPid::Check(!original.empty() && original[original.size() - 1] == '}',
                    "Existing diagnostic JSON writer did not produce one complete object");
  output << original.substr(0, original.size() - 1)
         << ",\"model\":" << RestartPid::Quote(result.model)
         << ",\"responseAlpha\":";
  RestartPid::JsonNumber(output, result.alpha);
  output << ",\"responseAlphaError\":";
  RestartPid::JsonNumber(output, result.alphaError);
  output << ",\"responseFreeParameters\":" << result.freeParameters
         << ",\"responseSigmaMeaning\":\"Gaussian core width, not total RMS\"}";
}

static void ValidateComponents(const NuclearFit::Request& request,
                               const NuclearFit::Result& result,
                               const NuclearFit::Settings& settings) {
  RestartPid::Check(result.components.size() == request.components.size(),
                    "Response fit changed the known PID component list");
  for (size_t species = 0; species < request.components.size(); ++species) {
    const NuclearFit::Component& initial = request.components[species];
    const NuclearFit::Component& fitted = result.components[species];
    RestartPid::Check(initial.name == fitted.name && initial.active == fitted.active &&
                      initial.freeMean == fitted.freeMean && initial.freeSigma == fitted.freeSigma,
                      "Response fit changed an active/free component policy");
    if (!initial.active) {
      RestartPid::Check(RestartPid::Close(initial.mean, fitted.mean, settings) &&
                        RestartPid::Close(initial.sigma, fitted.sigma, settings),
                        "Response fit changed another nucleus' fixed PID curve");
    }
  }
}

static void Run(const char* cachePath, const char* outputDir, const char* settingsPath,
                const char* slicesPath, const char* nodesPath) {
  NuclearFit::Settings settings = NuclearFit::ReadSettings(settingsPath);
  NuclearFit::ValidateSettings(settings);
  const std::string model = NuclearFit::Text(settings, "responseModel");
  RestartPid::Check(model == "gaussian" || model == "gaussian_left_exp" || model == "gaussian_right_exp",
                    "Unsupported explicitly configured response model");
  const double alphaStart = NuclearFit::Number(settings, "responseAlphaStart");
  RestartPid::Check(NuclearFit::Number(settings, "readbackToleranceAbsolute") >= 0 &&
                    NuclearFit::Number(settings, "readbackToleranceRelative") >= 0,
                    "Negative readback tolerance");
  std::vector<RestartPid::Slice> slices = RestartPid::ReadSlices(slicesPath);
  std::map<std::string, RestartPid::Shapes> shapes = RestartPid::ReadNodes(nodesPath);
  RestartPid::Check(shapes.size() == slices.size(), "Shape table must match exactly the requested slices");
  const std::string rootPath = std::string(outputDir) + "/native_spectra.root";
  const std::string jsonPath = std::string(outputDir) + "/fits.json";
  RestartPid::Check(gSystem->AccessPathName(rootPath.c_str()) && gSystem->AccessPathName(jsonPath.c_str()),
                    "Output already exists; use a fresh response-study directory");
  RestartPid::Check(gSystem->mkdir(outputDir, true) == 0 || !gSystem->AccessPathName(outputDir),
                    "Cannot create response-study output directory");
  TFile input(cachePath, "READ");
  RestartPid::Check(!input.IsZombie(), "Cannot open immutable native input cache");
  TFile output(rootPath.c_str(), "CREATE");
  RestartPid::Check(!output.IsZombie(), "Cannot create fresh native spectra output");
  std::vector<std::string> resultRows;
  std::vector<bool> fitted;
  size_t failures = 0;
  for (size_t index = 0; index < slices.size(); ++index) {
    const RestartPid::Slice& slice = slices[index];
    RestartPid::Check(shapes.find(slice.id) != shapes.end(), "Missing slice nodes: " + slice.id);
    TH2* histogram = RestartPid::InputHistogram(input, settings, slice.mode);
    RestartPid::Check(slice.rBin <= histogram->GetNbinsX(), "Requested R bin outside native axis");
    TH1D* native = histogram->ProjectionY(RestartPid::HistogramName("raw_native", index).c_str(),
                                         slice.rBin, slice.rBin, "e");
    RestartPid::Check(native != 0, "Native projection failed");
    native->SetDirectory(0);
    output.cd(); native->Write();
    std::ostringstream row;
    bool complete = false;
    try {
      NuclearFit::Request request = RestartPid::Request(*native, slice, shapes[slice.id], settings);
      RestartResponse::Result response = RestartResponse::Fit(request, settings, model, alphaStart);
      RestartPid::Check(response.model == model, "Response fit returned a different model family");
      RestartPid::Check(response.freeParameters >= 0, "Invalid response free-parameter count");
      const NuclearFit::Result& result = response.fit;
      ValidateComponents(request, result, settings);
      RestartPid::SavePrediction(output, *native, slice, result, index, settings);
      std::ostringstream gaussianCompatible;
      RestartPid::WriteResult(gaussianCompatible, slice, shapes[slice.id], request, result, *native, index,
                              histogram->GetXaxis()->GetBinLowEdge(slice.rBin),
                              histogram->GetXaxis()->GetBinUpEdge(slice.rBin));
      AppendResponse(row, gaussianCompatible.str(), response);
      complete = true;
      std::cout << "RESTART_PEAK_FIT id=" << slice.id << " model=" << model << " status=" << result.status
                << " pearson=" << result.pearson << " ndf=" << result.conditionalNdf << std::endl;
    } catch (const std::exception& error) {
      ++failures;
      row.str(""); row.clear();
      row << "{\"id\":" << RestartPid::Quote(slice.id)
          << ",\"mode\":" << RestartPid::Quote(slice.mode)
          << ",\"target\":" << RestartPid::Quote(slice.target) << ",\"rBin\":" << slice.rBin
          << ",\"fitCompleted\":false,\"purityAvailable\":false,\"purity\":null"
          << ",\"curveAdoptionAllowed\":false,\"nativeHistogram\":"
          << RestartPid::Quote(RestartPid::HistogramName("raw_native", index))
          << ",\"model\":" << RestartPid::Quote(model)
          << ",\"responseAlpha\":null,\"responseAlphaError\":null,\"responseFreeParameters\":null"
          << ",\"responseSigmaMeaning\":\"Gaussian core width, not total RMS\""
          << ",\"error\":" << RestartPid::Quote(error.what()) << '}';
      std::cerr << "RESTART_PEAK_FIT_FAILED id=" << slice.id << " model=" << model
                << " error=" << error.what() << std::endl;
    }
    resultRows.push_back(row.str()); fitted.push_back(complete);
    delete native;
  }
  output.Write(); output.Close();
  RestartPid::Readback(rootPath.c_str(), input, slices, settings, fitted);
  std::ofstream json(jsonPath.c_str());
  RestartPid::Check(json.good(), "Cannot write response fit JSON");
  json << "{\"status\":\"completed\",\"study\":\"fresh_response_peak_diagnostic_20261020\","
       << "\"productionEnabled\":false,\"coefficientsUpdated\":false,\"purityAvailable\":false,"
       << "\"curveAdoptionAllowed\":false,\"initialShapeSource\":\"git_log10_polynomial\","
       << "\"model\":" << RestartPid::Quote(model) << ",\"backgroundModel\":\"linear\","
       << "\"responseSigmaMeaning\":\"Gaussian core width, not total RMS\","
       << "\"scope\":\"target-only response diagnostic; fixed original curves for other nuclei; not a purity fit\","
       << "\"rawReadback\":true,\"predictionReadback\":true,\"totalSlices\":" << slices.size()
       << ",\"fitFailures\":" << failures << ",\"cachePath\":" << RestartPid::Quote(cachePath)
       << ",\"results\":[";
  for (size_t index = 0; index < resultRows.size(); ++index) {
    if (index) json << ',';
    json << resultRows[index];
  }
  json << "]}\n";
  json.close();
  RestartPid::Check(json.good(), "Response fit JSON write failed");
  std::cout << "RESTART_PEAK_READBACK_OK slices=" << slices.size() << " failures=" << failures << std::endl;
}

}  // namespace RestartResponseSlices
#endif

void FitResponseSlices(const char* cachePath, const char* outputDir, const char* settingsPath,
                       const char* slicesPath, const char* nodesPath) {
#ifndef __CINT__
  try { RestartResponseSlices::Run(cachePath, outputDir, settingsPath, slicesPath, nodesPath); }
  catch (const std::exception& error) {
    std::cerr << "RESTART_RESPONSE_FATAL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
