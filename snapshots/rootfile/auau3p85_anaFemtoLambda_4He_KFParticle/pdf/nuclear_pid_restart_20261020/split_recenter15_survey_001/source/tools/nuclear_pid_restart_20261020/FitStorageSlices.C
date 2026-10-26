// Diagnostic storage-response comparison. Not a production PID calibration.
void FitStorageSlices(const char*, const char*, const char*, const char*, const char*);
#ifndef __CINT__
#include "FitRestartSlices.C"
#include "Float16BinMapping.h"
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>

namespace RestartStorage {

class FreshText {
 public:
  explicit FreshText(const std::string& path) : stream_(0) {
    const int descriptor = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
    RestartPid::Check(descriptor >= 0,
                      "Cannot create fresh-only text output: " + path +
                      " (" + std::strerror(errno) + ")");
    stream_ = ::fdopen(descriptor, "w");
    if (!stream_) {
      ::close(descriptor);
      throw std::runtime_error("Cannot open output stream; incomplete output retained: " + path);
    }
  }
  ~FreshText() { if (stream_) std::fclose(stream_); }
  void Write(const std::string& text) {
    RestartPid::Check(stream_ &&
                      std::fwrite(text.data(), 1, text.size(), stream_) == text.size(),
                      "Diagnostic text write failed");
  }
  void Close() {
    RestartPid::Check(stream_ != 0, "Diagnostic text stream already closed");
    const int flushed = std::fflush(stream_);
    FILE* closing = stream_;
    stream_ = 0;
    const int closed = std::fclose(closing);
    RestartPid::Check(flushed == 0 && closed == 0, "Diagnostic text flush/close failed");
  }
 private:
  FreshText(const FreshText&);
  FreshText& operator=(const FreshText&);
  FILE* stream_;
};

static std::string Mapping(const NuclearFit::Settings& settings) {
  const std::string mapping = NuclearFit::Text(settings, "storageMapping");
  RestartPid::Check(mapping == "continuous" || mapping == "root_default12",
                    "storageMapping must explicitly be continuous or root_default12");
  return mapping;
}

static std::vector<double> EffectiveEdges(const NuclearFit::Request& request,
                                         const std::string& mapping) {
  RestartPid::Check(!request.bins.empty(), "No whole native bins to map");
  std::vector<double> nominal;
  nominal.push_back(request.bins.front().lo);
  for (size_t bin = 0; bin < request.bins.size(); ++bin) {
    RestartPid::Check(bin == 0 || request.bins[bin].lo == request.bins[bin - 1].hi,
                      "Nominal fit bins are not adjacent");
    nominal.push_back(request.bins[bin].hi);
  }
  const std::vector<double> effective = mapping == "continuous" ? nominal :
      RestartFloat16::IntegrationEdges(nominal);
  RestartPid::Check(effective.size() == nominal.size(), "Effective-edge dimension changed");
  for (size_t edge = 0; edge < effective.size(); ++edge) {
    RestartPid::Check(TMath::Finite(effective[edge]) && effective[edge] >= 0,
                      "Invalid storage-response integration edge");
    if (edge) RestartPid::Check(effective[edge - 1] < effective[edge],
                                "Storage-response preimage has a collapsed/reversed bin");
  }
  return effective;
}

static void ApplyEdges(NuclearFit::Request& request, const std::vector<double>& effective) {
  RestartPid::Check(effective.size() == request.bins.size() + 1,
                    "Effective-edge count does not match the native fit bins");
  for (size_t bin = 0; bin < request.bins.size(); ++bin) {
    request.bins[bin].lo = effective[bin];
    request.bins[bin].hi = effective[bin + 1];
  }
}

static std::string MapHeader() {
  return "id\tmode\ttarget\trBin\tdedxBin\tnominalLo\tnominalHi\teffectiveLo\teffectiveHi\tcount\terror\tstorageMapping";
}

static std::string MapRow(const RestartPid::Slice& slice, const TH1D& native,
                           const NuclearFit::Bin& nominal, double effectiveLo,
                           double effectiveHi, const std::string& mapping) {
  const int bin = native.GetXaxis()->FindFixBin(0.5 * (nominal.lo + nominal.hi));
  RestartPid::Check(bin >= 1 && bin <= native.GetNbinsX(), "Invalid nominal dEdx bin index");
  const double count = native.GetBinContent(bin), error = native.GetBinError(bin);
  RestartPid::Check(TMath::Finite(count) && count >= 0 && TMath::Finite(error) && error >= 0 &&
                    count == nominal.count, "Invalid or modified native count/error");
  std::ostringstream row;
  row << std::setprecision(17) << slice.id << '\t' << slice.mode << '\t' << slice.target
      << '\t' << slice.rBin << '\t' << bin << '\t' << nominal.lo << '\t' << nominal.hi
      << '\t' << effectiveLo << '\t' << effectiveHi << '\t' << count << '\t' << error
      << '\t' << mapping;
  return row.str();
}

static size_t MappingReadback(const char* path, TFile& input,
                              const NuclearFit::Settings& settings,
                              const std::vector<RestartPid::Slice>& slices,
                              const std::map<std::string, RestartPid::Shapes>& shapes,
                              const std::string& mapping) {
  std::ifstream table(path);
  RestartPid::Check(table.good(), "Cannot reopen native/effective bin TSV");
  std::string line;
  RestartPid::Check(bool(std::getline(table, line)) && line == MapHeader(),
                    "Native/effective bin TSV header mismatch");
  size_t rows = 0;
  for (size_t index = 0; index < slices.size(); ++index) {
    const RestartPid::Slice& slice = slices[index];
    TH2* source = RestartPid::InputHistogram(input, settings, slice.mode);
    TH1D* native = source->ProjectionY(RestartPid::HistogramName("map_readback", index).c_str(),
                                      slice.rBin, slice.rBin, "e");
    RestartPid::Check(native != 0, "Native mapping readback projection failed");
    native->SetDirectory(0);
    const std::map<std::string, RestartPid::Shapes>::const_iterator shape = shapes.find(slice.id);
    RestartPid::Check(shape != shapes.end(), "Missing slice nodes during mapping readback");
    const NuclearFit::Request nominal = RestartPid::Request(*native, slice, shape->second, settings);
    const std::vector<double> effective = EffectiveEdges(nominal, mapping);
    for (size_t bin = 0; bin < nominal.bins.size(); ++bin) {
      RestartPid::Check(bool(std::getline(table, line)) &&
                        line == MapRow(slice, *native, nominal.bins[bin],
                                       effective[bin], effective[bin + 1], mapping),
                        "Native/effective bin TSV exact readback mismatch");
      ++rows;
    }
    delete native;
  }
  RestartPid::Check(!std::getline(table, line) && table.eof(),
                    "Unexpected rows or read error in native/effective bin TSV");
  return rows;
}

static std::string AddMappingMetadata(const std::string& original,
                                      const std::string& mapping,
                                      const NuclearFit::Request& effective) {
  RestartPid::Check(!original.empty() && original[original.size() - 1] == '}',
                    "Cannot attach storage metadata to fit result");
  std::ostringstream out;
  out << original.substr(0, original.size() - 1)
      << ",\"storageMapping\":" << RestartPid::Quote(mapping)
      << ",\"effectiveFitLo\":";
  RestartPid::JsonNumber(out, effective.bins.front().lo);
  out << ",\"effectiveFitHi\":";
  RestartPid::JsonNumber(out, effective.bins.back().hi);
  out << ",\"storageResponseDiagnosticOnly\":true}";
  return out.str();
}

static void Run(const char* cachePath, const char* outputDir, const char* settingsPath,
                const char* slicesPath, const char* nodesPath) {
  RestartPid::Check(cachePath && *cachePath && outputDir && *outputDir &&
                    settingsPath && *settingsPath && slicesPath && *slicesPath &&
                    nodesPath && *nodesPath, "All diagnostic paths must be supplied");
  const NuclearFit::Settings settings = NuclearFit::ReadSettings(settingsPath);
  NuclearFit::ValidateSettings(settings);
  RestartPid::Check(NuclearFit::Number(settings, "readbackToleranceAbsolute") >= 0 &&
                    NuclearFit::Number(settings, "readbackToleranceRelative") >= 0,
                    "Negative readback tolerance");
  const std::string mapping = Mapping(settings);
  const std::vector<RestartPid::Slice> slices = RestartPid::ReadSlices(slicesPath);
  const std::map<std::string, RestartPid::Shapes> shapes = RestartPid::ReadNodes(nodesPath);
  RestartPid::Check(shapes.size() == slices.size(), "Shape table does not match requested slices");
  const std::string rootPath = std::string(outputDir) + "/native_spectra.root";
  const std::string jsonPath = std::string(outputDir) + "/fits.json";
  const std::string mapPath = std::string(outputDir) + "/native_effective_bins.tsv";
  RestartPid::Check(gSystem->AccessPathName(rootPath.c_str()) &&
                    gSystem->AccessPathName(jsonPath.c_str()) &&
                    gSystem->AccessPathName(mapPath.c_str()),
                    "Output already exists; use a fresh storage diagnostic directory");
  RestartPid::Check(gSystem->mkdir(outputDir, true) == 0 || !gSystem->AccessPathName(outputDir),
                    "Cannot create storage diagnostic directory");
  TFile input(cachePath, "READ");
  RestartPid::Check(!input.IsZombie() && !input.IsWritable(), "Cannot open read-only native input");
  TFile output(rootPath.c_str(), "CREATE");
  RestartPid::Check(!output.IsZombie() && output.IsWritable(), "Cannot create fresh diagnostic ROOT");
  FreshText mapFile(mapPath), jsonFile(jsonPath);
  mapFile.Write(MapHeader() + "\n");
  std::vector<std::string> results;
  std::vector<bool> fitted;
  size_t failures = 0, mapRows = 0;
  for (size_t index = 0; index < slices.size(); ++index) {
    const RestartPid::Slice& slice = slices[index];
    const std::map<std::string, RestartPid::Shapes>::const_iterator shape = shapes.find(slice.id);
    RestartPid::Check(shape != shapes.end(), "Missing original-curve nodes for slice: " + slice.id);
    TH2* source = RestartPid::InputHistogram(input, settings, slice.mode);
    RestartPid::Check(slice.rBin >= 1 && slice.rBin <= source->GetNbinsX(),
                      "Requested R bin outside native axis");
    TH1D* native = source->ProjectionY(RestartPid::HistogramName("raw_native", index).c_str(),
                                      slice.rBin, slice.rBin, "e");
    RestartPid::Check(native != 0, "Native diagnostic projection failed");
    native->SetDirectory(0);
    output.cd();
    RestartPid::Check(native->Write() > 0, "Cannot write native diagnostic spectrum");
    const NuclearFit::Request nominal = RestartPid::Request(*native, slice, shape->second, settings);
    const std::vector<double> effectiveEdges = EffectiveEdges(nominal, mapping);
    for (size_t bin = 0; bin < nominal.bins.size(); ++bin) {
      mapFile.Write(MapRow(slice, *native, nominal.bins[bin],
                           effectiveEdges[bin], effectiveEdges[bin + 1], mapping) + "\n");
      ++mapRows;
    }
    NuclearFit::Request effective = nominal;
    ApplyEdges(effective, effectiveEdges);
    bool complete = false;
    std::ostringstream row;
    try {
      // Both isotope and smooth-background templates consume the same transformed
      // request.bins in NuclearFit::BuildTemplates; observed counts are untouched.
      const NuclearFit::Result result = NuclearFit::Fit(effective, settings);
      RestartPid::SavePrediction(output, *native, slice, result, index, settings);
      RestartPid::WriteResult(row, slice, shape->second, effective, result, *native, index,
                              source->GetXaxis()->GetBinLowEdge(slice.rBin),
                              source->GetXaxis()->GetBinUpEdge(slice.rBin));
      complete = true;
      std::cout << "RESTART_STORAGE_FIT id=" << slice.id << " mapping=" << mapping
                << " status=" << result.status << " pearson=" << result.pearson
                << " ndf=" << result.conditionalNdf << std::endl;
    } catch (const std::exception& error) {
      ++failures;
      row.str(""); row.clear();
      row << "{\"id\":" << RestartPid::Quote(slice.id)
          << ",\"mode\":" << RestartPid::Quote(slice.mode)
          << ",\"target\":" << RestartPid::Quote(slice.target)
          << ",\"rBin\":" << slice.rBin
          << ",\"fitCompleted\":false,\"purityAvailable\":false,\"purity\":null"
          << ",\"curveAdoptionAllowed\":false,\"nativeHistogram\":"
          << RestartPid::Quote(RestartPid::HistogramName("raw_native", index))
          << ",\"error\":" << RestartPid::Quote(error.what()) << '}';
      std::cerr << "RESTART_STORAGE_FIT_FAILED id=" << slice.id
                << " mapping=" << mapping << " error=" << error.what() << std::endl;
    }
    results.push_back(AddMappingMetadata(row.str(), mapping, effective));
    fitted.push_back(complete);
    delete native;
  }
  mapFile.Close();
  output.Write();
  RestartPid::Check(!output.TestBit(TFile::kWriteError), "Diagnostic ROOT write failed");
  output.Close();
  RestartPid::Check(!output.TestBit(TFile::kWriteError), "Diagnostic ROOT close/write failed");
  RestartPid::Readback(rootPath.c_str(), input, slices, settings, fitted);
  RestartPid::Check(MappingReadback(mapPath.c_str(), input, settings, slices, shapes, mapping) == mapRows,
                    "Native/effective bin TSV row count changed");
  input.Close();
  std::ostringstream json;
  json << "{\"status\":\"completed\",\"study\":\"fresh_storage_response_diagnostic_20261020\","
       << "\"productionEnabled\":false,\"purityAvailable\":false,\"curveAdoptionAllowed\":false,"
       << "\"storageMapping\":" << RestartPid::Quote(mapping)
       << ",\"storageMappingIsHypothesis\":true,\"onFileEncodingAsserted\":false,"
       << "\"initialShapeSource\":\"git_log10_polynomial\",\"model\":\"gaussian_linear\","
       << "\"scope\":\"target-only peak/storage diagnostic; not a purity fit or PID update\","
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
  RestartPid::Check(verifyJson.good(), "Cannot reopen diagnostic JSON");
  std::ostringstream savedJson;
  savedJson << verifyJson.rdbuf();
  RestartPid::Check(!verifyJson.bad() && savedJson.str() == json.str(), "Diagnostic JSON readback mismatch");
  std::cout << "RESTART_STORAGE_READBACK_OK slices=" << slices.size()
            << " failures=" << failures << " mapping=" << mapping
            << " mappingRows=" << mapRows << " noPidUpdate=1" << std::endl;
}

}  // namespace RestartStorage
#endif

void FitStorageSlices(const char* cachePath, const char* outputDir, const char* settingsPath,
                       const char* slicesPath, const char* nodesPath) {
#ifndef __CINT__
  try {
    RestartStorage::Run(cachePath, outputDir, settingsPath, slicesPath, nodesPath);
  } catch (const std::exception& error) {
    std::cerr << "RESTART_STORAGE_FATAL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
