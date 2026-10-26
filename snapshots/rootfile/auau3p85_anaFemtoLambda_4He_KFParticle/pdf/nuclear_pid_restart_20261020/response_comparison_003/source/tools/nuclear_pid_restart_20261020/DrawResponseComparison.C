// Compare saved model predictions on exactly the same native spectra. No fit.
void DrawResponseComparison(const char*, const char*, const char*);
#ifndef __CINT__
#include "DrawRestartFits.C"
#include <TMemFile.h>
#include <map>

namespace RestartResponseView {

struct Source {
  int caseIndex;
  std::string model, root, metrics;
};

static const char* Family(size_t index) {
  const char* names[] = {"gaussian", "gaussian_left_exp", "gaussian_right_exp"};
  RestartFitView::Require(index < 3, "Unknown response comparison column");
  return names[index];
}

static std::vector<Source> ReadSources(const char* path) {
  std::ifstream input(path);
  RestartFitView::Require(input.good(), std::string("Cannot read response comparison sources: ") + path);
  std::string line;
  RestartFitView::Require(bool(std::getline(input, line)), "Empty response comparison source table");
  RestartFitView::Require(line == "caseIndex\tmodel\troot\tmetrics", "Unexpected response comparison source header");
  std::vector<Source> sources;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> values = RestartFitView::Split(line);
    RestartFitView::Require(values.size() == 4, "Comparison sources require exactly four tab-separated columns");
    const double number = RestartFitView::Number(values[0], false);
    RestartFitView::Require(number >= 0 && number <= INT_MAX && number == std::floor(number),
                            "Comparison caseIndex must be a nonnegative integer");
    Source source;
    source.caseIndex = int(number);
    source.model = values[1]; source.root = values[2]; source.metrics = values[3];
    RestartFitView::Require(source.caseIndex == int(sources.size() / 3) &&
                            source.model == Family(sources.size() % 3),
                            "Comparison source groups must be contiguous cases and Gaussian/left/right columns");
    RestartFitView::Require(!source.root.empty() && !source.metrics.empty(), "Missing saved ROOT or metadata path");
    sources.push_back(source);
  }
  RestartFitView::Require(!sources.empty() && sources.size() % 3 == 0,
                          "Each comparison case requires exactly three response families");
  return sources;
}

static void SameCase(const RestartFitView::Row& reference, const RestartFitView::Row& candidate) {
  RestartFitView::Require(reference.id == candidate.id && reference.target == candidate.target &&
                          reference.mode == candidate.mode && reference.r == candidate.r &&
                          reference.fitLo == candidate.fitLo && reference.fitHi == candidate.fitHi,
                          "Response families do not represent the same slice identity/R/fit geometry");
  // The histogram index is local to its ROOT file, not a physical identity.
}

static TH1D* Raw(TFile& input, const RestartFitView::Row& row) {
  TH1D* raw = dynamic_cast<TH1D*>(input.Get(RestartFitView::Name("raw_native", row.index).c_str()));
  RestartFitView::Require(raw != 0, "Missing saved native spectrum for response comparison");
  return raw;
}

static void SameRaw(const TH1D& reference, const TH1D& candidate) {
  RestartFitView::Require(reference.GetNbinsX() == candidate.GetNbinsX(), "Response native bin counts differ");
  for (int bin = 0; bin <= reference.GetNbinsX() + 1; ++bin) {
    RestartFitView::Require(reference.GetBinContent(bin) == candidate.GetBinContent(bin) &&
                            reference.GetBinError(bin) == candidate.GetBinError(bin),
                            "Response comparison data/error/flow differs between model families");
    if (bin >= 1 && bin <= reference.GetNbinsX()) {
      RestartFitView::Require(reference.GetXaxis()->GetBinLowEdge(bin) == candidate.GetXaxis()->GetBinLowEdge(bin) &&
                              reference.GetXaxis()->GetBinUpEdge(bin) == candidate.GetXaxis()->GetBinUpEdge(bin),
                              "Response comparison native bin edges differ");
    }
  }
}

static void Run(const char* sourcesPath, const char* settingsPath, const char* pdfPath) {
  RestartFitView::Require(gSystem->AccessPathName(pdfPath), "Comparison PDF already exists; use a fresh filename");
  const NuclearFit::Settings settings = NuclearFit::ReadSettings(settingsPath);
  const std::vector<Source> sources = ReadSources(sourcesPath);
  const size_t cases = sources.size() / 3;
  std::map<std::string, std::vector<RestartFitView::Row> > metrics;
  std::map<std::string, TFile*> roots;
  for (size_t index = 0; index < sources.size(); ++index) {
    const Source& source = sources[index];
    if (metrics.find(source.metrics) == metrics.end()) {
      metrics[source.metrics] = RestartFitView::ReadRows(source.metrics.c_str());
      RestartFitView::Require(metrics[source.metrics].size() == cases,
                              "Each aggregated family metadata table must contain every comparison case");
    }
    if (roots.find(source.root) == roots.end()) {
      TFile* root = new TFile(source.root.c_str(), "READ");
      RestartFitView::Require(!root->IsZombie(), "Cannot open saved response ROOT: " + source.root);
      roots[source.root] = root;
    }
  }
  // Verify every raw spectrum before producing any comparison page. Models may
  // choose different alpha starts, but the observed data must be identical.
  for (size_t index = 0; index < cases; ++index) {
    const Source& first = sources[3 * index];
    const RestartFitView::Row& reference = metrics[first.metrics][index];
    TH1D* referenceRaw = Raw(*roots[first.root], reference);
    for (size_t family = 1; family < 3; ++family) {
      const Source& source = sources[3 * index + family];
      const RestartFitView::Row& candidate = metrics[source.metrics][index];
      SameCase(reference, candidate);
      SameRaw(*referenceRaw, *Raw(*roots[source.root], candidate));
    }
    std::cout << "RESPONSE_COMPARISON_RAW_EQUAL case=" << index << " id=" << reference.id
              << " families=3 counts=exact errors=exact flows=exact edges=exact" << std::endl;
  }
  gStyle->SetOptStat(0);
  // Match the wide six-panel layout; avoid a portrait page with tiny plots.
  gStyle->SetPaperSize(28., 16.);
  TCanvas canvas("restart_response_comparison", "Same-data peak-response comparison", 2100, 1150);
  canvas.Print((std::string(pdfPath) + "[").c_str());
  size_t pages = 0;
  int panels = 0;
  for (size_t index = 0; index < cases; ++index) {
    canvas.Clear(); canvas.Divide(3, 2, 0.0001, 0.0001);
    std::vector<TObject*> owned;
    for (int logarithmic = 0; logarithmic < 2; ++logarithmic) {
      for (size_t family = 0; family < 3; ++family) {
        const Source& source = sources[3 * index + family];
        const RestartFitView::Row& row = metrics[source.metrics][index];
        canvas.cd(logarithmic * 3 + int(family) + 1);
        const int serial = panels++;
        if (row.reason.find("RAW_ONLY_FALLBACK=1") != std::string::npos) {
          // A failed selection must not accidentally display an unrelated saved
          // model. Use a RAM-only view containing the exact native histogram.
          // No ROOT file is regenerated on disk and no data are transformed.
          std::ostringstream memoryName;
          memoryName << "response_raw_only_" << serial << ".root";
          TMemFile rawOnly(memoryName.str().c_str(), "RECREATE");
          TH1D* native = Raw(*roots[source.root], row);
          TH1D* rawCopy = dynamic_cast<TH1D*>(native->Clone(native->GetName()));
          RestartFitView::Require(rawCopy != 0, "Cannot make raw-only display view");
          rawCopy->SetDirectory(&rawOnly); rawOnly.cd(); rawCopy->Write();
          RestartFitView::Panel(rawOnly, row, logarithmic != 0, settings, owned, serial);
          std::cout << "RESPONSE_COMPARISON_RAW_ONLY case=" << index << " model=" << source.model
                    << " log=" << logarithmic << " serial=" << serial << std::endl;
        } else {
          RestartFitView::Panel(*roots[source.root], row, logarithmic != 0, settings, owned, serial);
        }
        std::cout << "RESPONSE_COMPARISON_PANEL case=" << index << " model=" << source.model
                  << " log=" << logarithmic << " serial=" << serial << " id=" << row.id << std::endl;
      }
    }
    canvas.Print(pdfPath); ++pages;
    canvas.Clear();
    for (size_t object = 0; object < owned.size(); ++object) delete owned[object];
  }
  canvas.Print((std::string(pdfPath) + "]").c_str());
  RestartFitView::Require(pages == cases && panels == int(6 * cases), "Comparison panel/page count mismatch");
  for (std::map<std::string, TFile*>::iterator input = roots.begin(); input != roots.end(); ++input)
    delete input->second;
  std::cout << "RESPONSE_COMPARISON_COMPLETE cases=" << cases << " pages=" << pages
            << " panels=" << panels << " refit=0 normalization=0" << std::endl;
}

}  // namespace RestartResponseView
#endif

void DrawResponseComparison(const char* sourcesPath, const char* settingsPath, const char* pdfPath) {
#ifndef __CINT__
  try { RestartResponseView::Run(sourcesPath, settingsPath, pdfPath); }
  catch (const std::exception& error) {
    std::cerr << "RESPONSE_COMPARISON_FATAL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
