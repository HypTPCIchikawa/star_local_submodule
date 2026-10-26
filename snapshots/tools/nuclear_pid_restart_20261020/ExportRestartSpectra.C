// Native-bin diagnostic export only. No fit, rebinning, or normalization.
void ExportRestartSpectra(const char*, const char*, const char*, const char*);
#ifndef __CINT__
#include "FitRestartSlices.C"
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace RestartSpectraExport {

static struct stat FileState(const char* path) {
  struct stat result;
  RestartPid::Check(::stat(path, &result) == 0,
                    std::string("Cannot stat input: ") + path);
  RestartPid::Check(S_ISREG(result.st_mode), "Input cache must be a regular local file");
  return result;
}

static void SameFile(const struct stat& before, const struct stat& after) {
  bool same = before.st_dev == after.st_dev && before.st_ino == after.st_ino &&
              before.st_size == after.st_size && before.st_mtime == after.st_mtime &&
              before.st_ctime == after.st_ctime;
#if defined(__linux__)
  same = same && before.st_mtim.tv_nsec == after.st_mtim.tv_nsec &&
                 before.st_ctim.tv_nsec == after.st_ctim.tv_nsec;
#endif
  RestartPid::Check(same, "Input cache identity/content metadata changed during export");
}

static TH1D* Project(TFile& input, const NuclearFit::Settings& settings,
                     const RestartPid::Slice& slice, size_t index) {
  TH2* source = RestartPid::InputHistogram(input, settings, slice.mode);
  RestartPid::Check(slice.rBin >= 1 && slice.rBin <= source->GetNbinsX(),
                    "Requested R bin outside native axis");
  TH1D* native = source->ProjectionY(
      RestartPid::HistogramName("export_native", index).c_str(),
      slice.rBin, slice.rBin, "e");
  RestartPid::Check(native != 0, "Native projection failed");
  native->SetDirectory(0);
  RestartPid::Check(native->GetNbinsX() == source->GetNbinsY(),
                    "Native dEdx bin count changed during projection");
  for (int bin = 0; bin <= native->GetNbinsX() + 1; ++bin) {
    const double lo = native->GetXaxis()->GetBinLowEdge(bin);
    const double hi = native->GetXaxis()->GetBinUpEdge(bin);
    const double count = native->GetBinContent(bin);
    const double error = native->GetBinError(bin);
    RestartPid::Check(TMath::Finite(lo) && TMath::Finite(hi) && lo < hi &&
                      TMath::Finite(count) && count >= 0 &&
                      TMath::Finite(error) && error >= 0,
                      "Invalid native bin edges, count, or error");
    RestartPid::Check(lo == source->GetYaxis()->GetBinLowEdge(bin) &&
                      hi == source->GetYaxis()->GetBinUpEdge(bin),
                      "Projection changed native dEdx bin edges");
    RestartPid::Check(count == source->GetBinContent(slice.rBin, bin) &&
                      error == source->GetBinError(slice.rBin, bin),
                      "Projection changed native count or error");
  }
  return native;
}

static std::string Row(const RestartPid::Slice& slice, const TH1D& native, int bin) {
  std::ostringstream row;
  row << std::setprecision(17) << slice.id << '\t' << slice.mode << '\t'
      << slice.target << '\t' << slice.rBin << '\t' << bin << '\t'
      << native.GetXaxis()->GetBinLowEdge(bin) << '\t'
      << native.GetXaxis()->GetBinUpEdge(bin) << '\t'
      << native.GetBinContent(bin) << '\t' << native.GetBinError(bin);
  return row.str();
}

static void WriteLine(FILE* output, const std::string& line) {
  RestartPid::Check(std::fwrite(line.data(), 1, line.size(), output) == line.size() &&
                    std::fputc('\n', output) != EOF, "Native TSV write failed");
}

static size_t Readback(const char* path, TFile& input,
                       const NuclearFit::Settings& settings,
                       const std::vector<RestartPid::Slice>& slices) {
  std::ifstream table(path);
  RestartPid::Check(table.good(), "Cannot reopen native TSV for readback");
  std::string line;
  RestartPid::Check(bool(std::getline(table, line)) &&
                    line == "id\tmode\ttarget\trBin\tdedxBin\tlo\thi\tcount\terror",
                    "Native TSV readback header mismatch");
  size_t rows = 0;
  for (size_t index = 0; index < slices.size(); ++index) {
    const RestartPid::Slice& slice = slices[index];
    TH1D* native = Project(input, settings, slice, index);
    for (int bin = 0; bin <= native->GetNbinsX() + 1; ++bin) {
      RestartPid::Check(bool(std::getline(table, line)), "Native TSV truncated");
      RestartPid::Check(line == Row(slice, *native, bin),
                        "Native TSV serialized row differs from the native projection");
      const std::vector<std::string> fields = RestartPid::Split(line, '\t');
      RestartPid::Check(fields.size() == 9 && fields[0] == slice.id &&
                        fields[1] == slice.mode && fields[2] == slice.target,
                        "Native TSV identity mismatch");
      RestartPid::Check(RestartPid::Number(fields[3], "rBin") == slice.rBin &&
                        RestartPid::Number(fields[4], "dedxBin") == bin &&
                        RestartPid::Number(fields[5], "lo") == native->GetXaxis()->GetBinLowEdge(bin) &&
                        RestartPid::Number(fields[6], "hi") == native->GetXaxis()->GetBinUpEdge(bin) &&
                        RestartPid::Number(fields[7], "count") == native->GetBinContent(bin) &&
                        RestartPid::Number(fields[8], "error") == native->GetBinError(bin),
                        "Native TSV floating-point round trip was not exact");
      ++rows;
    }
    delete native;
  }
  RestartPid::Check(!std::getline(table, line), "Unexpected extra native TSV rows");
  RestartPid::Check(table.eof(), "Native TSV readback stream error");
  return rows;
}

static void Run(const char* cachePath, const char* settingsPath,
                const char* slicesPath, const char* outputPath) {
  RestartPid::Check(cachePath && *cachePath && settingsPath && *settingsPath &&
                    slicesPath && *slicesPath && outputPath && *outputPath,
                    "Every exporter path must be explicitly supplied");
  const struct stat before = FileState(cachePath);
  const NuclearFit::Settings settings = NuclearFit::ReadSettings(settingsPath);
  const std::vector<RestartPid::Slice> slices = RestartPid::ReadSlices(slicesPath);
  TFile input(cachePath, "READ");
  RestartPid::Check(!input.IsZombie() && !input.IsWritable(), "Cannot open read-only input cache");

  // Check every request before creating a possibly incomplete output file.
  int dedxBins = -1;
  for (size_t index = 0; index < slices.size(); ++index) {
    TH1D* native = Project(input, settings, slices[index], index);
    RestartPid::Check(dedxBins < 0 || dedxBins == native->GetNbinsX(),
                      "Requested slices have inconsistent native dEdx bin counts");
    dedxBins = native->GetNbinsX();
    delete native;
  }
  const int descriptor = ::open(outputPath, O_WRONLY | O_CREAT | O_EXCL, 0666);
  RestartPid::Check(descriptor >= 0,
                    std::string("Cannot create fresh-only TSV: ") + outputPath +
                    " (" + std::strerror(errno) + ")");
  FILE* output = ::fdopen(descriptor, "w");
  if (!output) {
    ::close(descriptor);
    throw std::runtime_error("Cannot open fresh TSV stream; incomplete output retained");
  }
  size_t written = 0;
  try {
    WriteLine(output, "id\tmode\ttarget\trBin\tdedxBin\tlo\thi\tcount\terror");
    for (size_t index = 0; index < slices.size(); ++index) {
      TH1D* native = Project(input, settings, slices[index], index);
      for (int bin = 0; bin <= native->GetNbinsX() + 1; ++bin) {
        WriteLine(output, Row(slices[index], *native, bin));
        ++written;
      }
      delete native;
    }
    RestartPid::Check(std::fflush(output) == 0, "Native TSV flush failed");
  } catch (...) {
    std::fclose(output);
    throw;
  }
  RestartPid::Check(std::fclose(output) == 0, "Native TSV close failed");
  const size_t expected = slices.size() * size_t(dedxBins + 2);
  RestartPid::Check(written == expected &&
                    Readback(outputPath, input, settings, slices) == expected,
                    "Native TSV full row count mismatch");
  input.Close();
  SameFile(before, FileState(cachePath));
  std::cout << "RESTART_SPECTRA_EXPORT_COMPLETE slices=" << slices.size()
            << " rows=" << written << " noFit=1" << std::endl;
}

}  // namespace RestartSpectraExport
#endif

void ExportRestartSpectra(const char* cachePath, const char* settingsPath,
                          const char* slicesPath, const char* outputPath) {
#ifndef __CINT__
  try {
    RestartSpectraExport::Run(cachePath, settingsPath, slicesPath, outputPath);
  } catch (const std::exception& error) {
    std::cerr << "RESTART_SPECTRA_EXPORT_FATAL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
