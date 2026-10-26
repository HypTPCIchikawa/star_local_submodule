// Study-only split-normal response with explicit storage-response integration.
// Location is the common mode; sigma is sqrt(sigmaLeft * sigmaRight), not RMS.
void FitSplitNormalSlices(const char*, const char*, const char*, const char*, const char*);
#ifndef __CINT__
#include "FitStorageSlices.C"
#include "SplitNormalFit.h"

namespace RestartSplitNormalSlices {

static void ValidateComponents(const NuclearFit::Request& request,
                               const NuclearFit::Result& result,
                               const NuclearFit::Settings& settings) {
  RestartPid::Check(result.components.size() == request.components.size(),
                    "Split-normal fit changed the known component list");
  for (size_t species = 0; species < request.components.size(); ++species) {
    const NuclearFit::Component& initial = request.components[species];
    const NuclearFit::Component& fitted = result.components[species];
    RestartPid::Check(initial.name == fitted.name && initial.active == fitted.active &&
                      initial.freeMean == fitted.freeMean && initial.freeSigma == fitted.freeSigma,
                      "Split-normal fit changed the active/free component policy");
    if (initial.active) continue;
    RestartPid::Check(RestartPid::Close(initial.mean, fitted.mean, settings) &&
                      RestartPid::Close(initial.sigma, fitted.sigma, settings) &&
                      initial.nodes.size() == fitted.nodes.size(),
                      "Split-normal fit changed another nucleus' fixed PID curve");
    for (size_t node = 0; node < initial.nodes.size(); ++node) {
      RestartPid::Check(RestartPid::Close(initial.nodes[node].weight, fitted.nodes[node].weight, settings) &&
                        RestartPid::Close(initial.nodes[node].mean, fitted.nodes[node].mean, settings) &&
                        RestartPid::Close(initial.nodes[node].sigma, fitted.nodes[node].sigma, settings),
                        "Split-normal fit changed another nucleus' fixed R nodes");
    }
  }
}

static std::string AppendResponse(const std::string& original,
                                   const RestartSplitNormal::Result& result) {
  RestartPid::Check(!original.empty() && original[original.size() - 1] == '}',
                    "Existing writer did not produce one complete result object");
  std::ostringstream out;
  out << original.substr(0, original.size() - 1)
      << ",\"model\":\"split_normal\",\"responseWidthRatio\":";
  RestartPid::JsonNumber(out, result.ratio);
  out << ",\"responseWidthRatioError\":";
  RestartPid::JsonNumber(out, result.ratioError);
  out << ",\"responseSigmaLeft\":";
  RestartPid::JsonNumber(out, result.sigmaLeft);
  out << ",\"responseSigmaRight\":";
  RestartPid::JsonNumber(out, result.sigmaRight);
  out << ",\"responseSigmaLeftError\":";
  RestartPid::JsonNumber(out, result.sigmaLeftError);
  out << ",\"responseSigmaRightError\":";
  RestartPid::JsonNumber(out, result.sigmaRightError);
  out << ",\"responseOptimizerDiagnosticsAvailable\":"
      << (result.optimizerDiagnosticsAvailable ? "true" : "false")
      << ",\"responseMinimizeReturned\":" << (result.minimizeReturned ? "true" : "false")
      << ",\"responseMinimizeConverged\":" << (result.minimizeConverged ? "true" : "false")
      << ",\"responseMigradStatus\":" << result.migradStatus
      << ",\"responseMigradEdm\":";
  RestartPid::JsonNumber(out, result.migradEdm);
  out << ",\"responseHesseAttempted\":" << (result.hesseAttempted ? "true" : "false")
      << ",\"responseHesseSucceeded\":" << (result.hesseSucceeded ? "true" : "false")
      << ",\"responseHesseStatus\":" << result.hesseStatus
      << ",\"responseFreeParameters\":" << result.freeParameters
      << ",\"responseMeanMeaning\":\"common mode/location; not the distribution mean\""
      << ",\"responseSigmaMeaning\":\"geometric mean sqrt(sigmaLeft*sigmaRight); not RMS\""
      << ",\"responseRatioMeaning\":\"sigmaRight/sigmaLeft\"}";
  return out.str();
}

static std::string Failure(const RestartPid::Slice& slice, size_t index,
                            const std::exception& error) {
  std::ostringstream out;
  out << "{\"id\":" << RestartPid::Quote(slice.id)
      << ",\"mode\":" << RestartPid::Quote(slice.mode)
      << ",\"target\":" << RestartPid::Quote(slice.target)
      << ",\"rBin\":" << slice.rBin
      << ",\"fitLo\":" << std::setprecision(17) << slice.fitLo
      << ",\"fitHi\":" << slice.fitHi
      << ",\"gateLo\":" << slice.gateLo << ",\"gateHi\":" << slice.gateHi
      << ",\"fitCompleted\":false,\"purityAvailable\":false,\"purity\":null"
      << ",\"curveAdoptionAllowed\":false,\"nativeHistogram\":"
      << RestartPid::Quote(RestartPid::HistogramName("raw_native", index))
      << ",\"model\":\"split_normal\",\"responseWidthRatio\":null"
      << ",\"responseWidthRatioError\":null,\"responseSigmaLeft\":null,\"responseSigmaRight\":null"
      << ",\"responseSigmaLeftError\":null,\"responseSigmaRightError\":null"
      << ",\"responseOptimizerDiagnosticsAvailable\":false"
      << ",\"responseMinimizeReturned\":null,\"responseMinimizeConverged\":null"
      << ",\"responseMigradStatus\":null,\"responseMigradEdm\":null"
      << ",\"responseHesseAttempted\":null,\"responseHesseSucceeded\":null"
      << ",\"responseHesseStatus\":null,\"responseFreeParameters\":null"
      << ",\"responseMeanMeaning\":\"common mode/location; not the distribution mean\""
      << ",\"responseSigmaMeaning\":\"geometric mean sqrt(sigmaLeft*sigmaRight); not RMS\""
      << ",\"responseRatioMeaning\":\"sigmaRight/sigmaLeft\""
      << ",\"error\":" << RestartPid::Quote(error.what()) << '}';
  return out.str();
}

static void Run(const char* cachePath, const char* outputDir, const char* settingsPath,
                const char* slicesPath, const char* nodesPath) {
  RestartPid::Check(cachePath && *cachePath && outputDir && *outputDir &&
                    settingsPath && *settingsPath && slicesPath && *slicesPath &&
                    nodesPath && *nodesPath, "All split-normal diagnostic paths must be supplied");
  const NuclearFit::Settings settings = NuclearFit::ReadSettings(settingsPath);
  NuclearFit::ValidateSettings(settings);
  RestartPid::Check(NuclearFit::Text(settings, "responseModel") == "split_normal",
                    "responseModel must explicitly be split_normal");
  const double ratioStart = NuclearFit::Number(settings, "splitRatioStart");
  RestartPid::Check(ratioStart > 0, "splitRatioStart must be positive");
  RestartPid::Check(NuclearFit::Number(settings, "readbackToleranceAbsolute") >= 0 &&
                    NuclearFit::Number(settings, "readbackToleranceRelative") >= 0,
                    "Negative readback tolerance");
  const std::string mapping = RestartStorage::Mapping(settings);
  RestartPid::Check(mapping == "root_default12",
                    "This split-normal diagnostic requires explicit root_default12 mapping");
  const std::vector<RestartPid::Slice> slices = RestartPid::ReadSlices(slicesPath);
  const std::map<std::string, RestartPid::Shapes> shapes = RestartPid::ReadNodes(nodesPath);
  RestartPid::Check(shapes.size() == slices.size(), "Shape table does not match requested slices");
  const std::string rootPath = std::string(outputDir) + "/native_spectra.root";
  const std::string jsonPath = std::string(outputDir) + "/fits.json";
  const std::string mapPath = std::string(outputDir) + "/native_effective_bins.tsv";
  RestartPid::Check(gSystem->AccessPathName(rootPath.c_str()) &&
                    gSystem->AccessPathName(jsonPath.c_str()) &&
                    gSystem->AccessPathName(mapPath.c_str()),
                    "Output already exists; use a fresh split-normal directory");
  RestartPid::Check(gSystem->mkdir(outputDir, true) == 0 || !gSystem->AccessPathName(outputDir),
                    "Cannot create split-normal output directory");
  TFile input(cachePath, "READ");
  RestartPid::Check(!input.IsZombie() && !input.IsWritable(), "Cannot open read-only native input");
  TFile output(rootPath.c_str(), "CREATE");
  RestartPid::Check(!output.IsZombie() && output.IsWritable(), "Cannot create fresh split-normal ROOT");
  RestartStorage::FreshText mapFile(mapPath), jsonFile(jsonPath);
  mapFile.Write(RestartStorage::MapHeader() + "\n");
  std::vector<std::string> results;
  std::vector<bool> fitted;
  size_t failures = 0, mapRows = 0;
  for (size_t index = 0; index < slices.size(); ++index) {
    const RestartPid::Slice& slice = slices[index];
    const std::map<std::string, RestartPid::Shapes>::const_iterator shape = shapes.find(slice.id);
    RestartPid::Check(shape != shapes.end(), "Missing original Git nodes for slice: " + slice.id);
    TH2* source = RestartPid::InputHistogram(input, settings, slice.mode);
    RestartPid::Check(slice.rBin >= 1 && slice.rBin <= source->GetNbinsX(),
                      "Requested R bin outside native axis");
    TH1D* native = source->ProjectionY(RestartPid::HistogramName("raw_native", index).c_str(),
                                      slice.rBin, slice.rBin, "e");
    RestartPid::Check(native != 0, "Native split-normal projection failed");
    native->SetDirectory(0);
    output.cd();
    RestartPid::Check(native->Write() > 0, "Cannot write native split-normal spectrum");
    const NuclearFit::Request nominal = RestartPid::Request(*native, slice, shape->second, settings);
    const std::vector<double> effectiveEdges = RestartStorage::EffectiveEdges(nominal, mapping);
    for (size_t bin = 0; bin < nominal.bins.size(); ++bin) {
      mapFile.Write(RestartStorage::MapRow(slice, *native, nominal.bins[bin],
                      effectiveEdges[bin], effectiveEdges[bin + 1], mapping) + "\n");
      ++mapRows;
    }
    NuclearFit::Request effective = nominal;
    RestartStorage::ApplyEdges(effective, effectiveEdges);
    bool complete = false;
    std::string row;
    try {
      const RestartSplitNormal::Result response = RestartSplitNormal::Fit(effective, settings, ratioStart);
      RestartPid::Check(response.freeParameters >= 0, "Invalid split-normal free-parameter count");
      const NuclearFit::Result& result = response.fit;
      ValidateComponents(effective, result, settings);
      RestartPid::SavePrediction(output, *native, slice, result, index, settings);
      std::ostringstream compatible;
      RestartPid::WriteResult(compatible, slice, shape->second, effective, result, *native, index,
                              source->GetXaxis()->GetBinLowEdge(slice.rBin),
                              source->GetXaxis()->GetBinUpEdge(slice.rBin));
      row = AppendResponse(compatible.str(), response);
      complete = true;
      std::cout << "RESTART_SPLIT_NORMAL_FIT id=" << slice.id << " status=" << result.status
                << " ratio=" << response.ratio << " sigmaLeft=" << response.sigmaLeft
                << " sigmaRight=" << response.sigmaRight << " pearson=" << result.pearson
                << " ndf=" << result.conditionalNdf << std::endl;
    } catch (const std::exception& error) {
      ++failures;
      row = Failure(slice, index, error);
      std::cerr << "RESTART_SPLIT_NORMAL_FIT_FAILED id=" << slice.id
                << " error=" << error.what() << std::endl;
    }
    results.push_back(RestartStorage::AddMappingMetadata(row, mapping, effective));
    fitted.push_back(complete);
    delete native;
  }
  mapFile.Close();
  output.Write();
  RestartPid::Check(!output.TestBit(TFile::kWriteError), "Split-normal ROOT write failed");
  output.Close();
  RestartPid::Check(!output.TestBit(TFile::kWriteError), "Split-normal ROOT close/write failed");
  RestartPid::Readback(rootPath.c_str(), input, slices, settings, fitted);
  RestartPid::Check(RestartStorage::MappingReadback(mapPath.c_str(), input, settings, slices, shapes, mapping) == mapRows,
                    "Native/effective bin TSV row count changed");
  input.Close();
  std::ostringstream json;
  json << "{\"status\":\"completed\",\"study\":\"fresh_split_normal_storage_diagnostic_20261020\","
       << "\"productionEnabled\":false,\"coefficientsUpdated\":false,\"purityAvailable\":false,"
       << "\"curveAdoptionAllowed\":false,\"storageMapping\":" << RestartPid::Quote(mapping)
       << ",\"storageMappingIsHypothesis\":true,\"onFileEncodingAsserted\":false,"
       << "\"initialShapeSource\":\"git_log10_polynomial\",\"model\":\"split_normal\","
       << "\"backgroundModel\":\"linear\","
       << "\"responseMeanMeaning\":\"common mode/location; not the distribution mean\","
       << "\"responseSigmaMeaning\":\"geometric mean sqrt(sigmaLeft*sigmaRight); not RMS\","
       << "\"responseRatioMeaning\":\"sigmaRight/sigmaLeft\","
       << "\"scope\":\"target-only split-normal diagnostic; other original Git curves fixed/inactive; not a purity fit\","
       << "\"nominalHistogramAxesPreserved\":true,\"rawCountsModified\":false,"
       << "\"signalAndBackgroundUseEffectiveEdges\":true,\"rawReadback\":true,"
       << "\"predictionReadback\":true,\"mappingReadback\":true,\"totalSlices\":" << slices.size()
       << ",\"mappingRows\":" << mapRows << ",\"fitFailures\":" << failures
       << ",\"cachePath\":" << RestartPid::Quote(cachePath)
       << ",\"mappingTable\":\"native_effective_bins.tsv\",\"results\":[";
  for (size_t index = 0; index < results.size(); ++index) {
    if (index) json << ',';
    json << results[index];
  }
  json << "]}\n";
  jsonFile.Write(json.str());
  jsonFile.Close();
  std::ifstream verifyJson(jsonPath.c_str());
  RestartPid::Check(verifyJson.good(), "Cannot reopen split-normal diagnostic JSON");
  std::ostringstream savedJson;
  savedJson << verifyJson.rdbuf();
  RestartPid::Check(!verifyJson.bad() && savedJson.str() == json.str(), "Split-normal JSON readback mismatch");
  std::cout << "RESTART_SPLIT_NORMAL_READBACK_OK slices=" << slices.size()
            << " failures=" << failures << " mapping=" << mapping
            << " mappingRows=" << mapRows << " noPidUpdate=1" << std::endl;
}

}  // namespace RestartSplitNormalSlices
#endif

void FitSplitNormalSlices(const char* cachePath, const char* outputDir, const char* settingsPath,
                          const char* slicesPath, const char* nodesPath) {
#ifndef __CINT__
  try {
    RestartSplitNormalSlices::Run(cachePath, outputDir, settingsPath, slicesPath, nodesPath);
  } catch (const std::exception& error) {
    std::cerr << "RESTART_SPLIT_NORMAL_FATAL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
