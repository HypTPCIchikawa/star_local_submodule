// ROOT5-compatible rendering only: all CF values and correlated ratio errors
// are supplied by the analysis TSV. No ROOT input or physics fit is performed.
// Usage (all arguments are required):
// DrawMassHalfCF(rowsPath, outputDir, energyLabel, normMin, normMax, xMax)
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TNamed.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace MassHalfPlot {

struct Row {
  std::string method, group, label, region;
  int bin;
  double low, high, value, error, ratio, ratioError;
};

typedef std::map<std::string, std::vector<Row> > Table;

const char* const kMethods[] = {"raw", "sbsub", "cfsub", "method3", "method5"};
const char* const kMethodTitles[] = {
    "Raw SE / ME", "Count-level sideband subtraction",
    "CF-level purity correction", "Independent k* mass fits",
    "Sliding-window mass fits"};
const char* const kGroups[] = {
    "all_0_80", "cent20_80", "cent0_20", "cent0_10", "cent10_20", "cent20_60"};
const char* const kGroupTitles[] = {"0-80%", "20-80%", "0-20%", "0-10%", "10-20%", "20-60%"};
const char* const kRegions[] = {"full", "low", "high"};
const char* const kRegionTitles[] = {
    "Full: #pm3#sigma", "Low: -3#sigma to 0", "High: 0 to +3#sigma"};
const int kColors[] = {kBlack, kRed + 1, kBlue + 1};
const int kMarkers[] = {20, 24, 25};

bool Known(const std::string& value, const char* const* values, unsigned count) {
  for (unsigned i = 0; i < count; ++i) if (value == values[i]) return true;
  return false;
}

std::string Key(const std::string& method, const std::string& group,
                const std::string& region) {
  return method + "/" + group + "/" + region;
}

double Number(const std::string& token, unsigned line) {
  char* end = 0;
  errno = 0;
  const double value = std::strtod(token.c_str(), &end);
  if (end == token.c_str() || !end || *end != '\0' || errno == ERANGE) {
    std::ostringstream message;
    message << "Malformed numeric field at TSV line " << line << ": " << token;
    throw std::runtime_error(message.str());
  }
  return value;
}

Table Read(const char* path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(std::string("Cannot read TSV: ") + path);
  const char* columns[] = {"method", "group", "label", "region", "bin", "low", "high",
                           "value", "error", "ratio", "ratio_error"};
  Table table;
  std::set<std::string> seen;
  std::string line;
  bool header = false;
  unsigned number = 0;
  while (std::getline(input, line)) {
    ++number;
    if (line.empty() || line[0] == '#') continue;
    std::istringstream parser(line);
    std::vector<std::string> fields;
    std::string field;
    while (parser >> field) fields.push_back(field);
    if (fields.empty()) continue;
    if (fields.size() != 11) {
      std::ostringstream message;
      message << "Expected 11 TSV fields at line " << number << ", found " << fields.size();
      throw std::runtime_error(message.str());
    }
    if (!header) {
      for (unsigned i = 0; i < 11; ++i)
        if (fields[i] != columns[i]) throw std::runtime_error("Unexpected TSV header/schema");
      header = true;
      continue;
    }
    Row row;
    row.method = fields[0]; row.group = fields[1]; row.label = fields[2]; row.region = fields[3];
    if (!Known(row.method, kMethods, 5) || !Known(row.group, kGroups, 6) ||
        !Known(row.region, kRegions, 3)) {
      throw std::runtime_error("Unknown method/group/region at TSV line " + fields[4]);
    }
    const double bin = Number(fields[4], number);
    if (!TMath::Finite(bin) || bin < 0 || bin > 1000000 || bin != TMath::Floor(bin))
      throw std::runtime_error("Invalid integer bin in TSV");
    row.bin = static_cast<int>(bin);
    row.low = Number(fields[5], number); row.high = Number(fields[6], number);
    row.value = Number(fields[7], number); row.error = Number(fields[8], number);
    row.ratio = Number(fields[9], number); row.ratioError = Number(fields[10], number);
    if (!TMath::Finite(row.low) || !TMath::Finite(row.high) || row.high <= row.low)
      throw std::runtime_error("Invalid finite bin edges in TSV");
    if ((TMath::Finite(row.error) && row.error < 0) ||
        (TMath::Finite(row.ratioError) && row.ratioError < 0))
      throw std::runtime_error("Negative uncertainty in TSV");
    std::ostringstream unique;
    unique << Key(row.method, row.group, row.region) << "/" << row.bin;
    if (!seen.insert(unique.str()).second) throw std::runtime_error("Duplicate TSV bin: " + unique.str());
    table[Key(row.method, row.group, row.region)].push_back(row);
  }
  if (!header) throw std::runtime_error("TSV header is missing");
  return table;
}

TGraphErrors* Graph(const Table& table, const std::string& method, const std::string& group,
                    int region, bool ratio, double xMax) {
  TGraphErrors* graph = new TGraphErrors();
  const std::string name = Key(method, group, kRegions[region]) + (ratio ? "/ratio" : "/cf");
  graph->SetName(name.c_str());
  graph->SetMarkerStyle(kMarkers[region]); graph->SetMarkerColor(kColors[region]);
  graph->SetLineColor(kColors[region]); graph->SetLineWidth(1);
  graph->SetMarkerSize(region == 0 ? 0.65 : 0.95);
  Table::const_iterator found = table.find(Key(method, group, kRegions[region]));
  if (found == table.end()) return graph;
  for (unsigned i = 0; i < found->second.size(); ++i) {
    const Row& row = found->second[i];
    const double x = 0.5 * (row.low + row.high);
    const double value = ratio ? row.ratio : row.value;
    const double error = ratio ? row.ratioError : row.error;
    if (!TMath::Finite(value) || !TMath::Finite(error) || x < 0 || x > xMax) continue;
    const int point = graph->GetN();
    graph->SetPoint(point, x, value); graph->SetPointError(point, 0, error);
  }
  return graph;
}

void Bounds(TGraphErrors** graphs, int first, double& low, double& high) {
  low = high = 1.;
  for (int region = first; region < 3; ++region) {
    for (int i = 0; i < graphs[region]->GetN(); ++i) {
      double x = 0, y = 0;
      graphs[region]->GetPoint(i, x, y);
      const double error = graphs[region]->GetErrorY(i);
      low = std::min(low, y - error); high = std::max(high, y + error);
    }
  }
  double span = high - low;
  if (span < 1.e-8) span = 0.2;
  // Leave room for labels without concealing any finite data or uncertainty.
  low -= 0.10 * span; high += 0.22 * span;
  const double magnitude = TMath::Power(10., TMath::Floor(TMath::Log10((high - low) / 5.)));
  const double relative = (high - low) / (5. * magnitude);
  const double step = magnitude * (relative > 5. ? 10. : relative > 2. ? 5. : relative > 1. ? 2. : 1.);
  low = step * TMath::Floor(low / step);
  high = step * TMath::Ceil(high / step);
}

void Label(double x, double y, const char* text, double size, int align = 12, int color = kBlack) {
  TLatex* label = new TLatex(x, y, text);
  label->SetNDC(); label->SetTextFont(42); label->SetTextSize(size);
  label->SetTextAlign(align); label->SetTextColor(color); label->Draw();
}

void Frame(TPad* pad, const char* name, TGraphErrors** graphs, bool ratio,
           double normMin, double normMax, double xMax) {
  pad->cd();
  double low = 0, high = 0;
  Bounds(graphs, ratio ? 1 : 0, low, high);
  TH1F* frame = pad->DrawFrame(0, low, xMax, high);
  frame->SetName(name); frame->SetTitle("");
  frame->GetXaxis()->SetTitle("k* [GeV/c]");
  frame->GetYaxis()->SetTitle(ratio ? "Half / full" : "CF");
  frame->GetXaxis()->SetLabelFont(42); frame->GetYaxis()->SetLabelFont(42);
  frame->GetXaxis()->SetTitleFont(42); frame->GetYaxis()->SetTitleFont(42);
  frame->GetYaxis()->SetNdivisions(ratio ? 504 : 505);
  frame->GetXaxis()->SetNdivisions(505);
  if (ratio) {
    frame->GetXaxis()->SetLabelSize(0.09); frame->GetXaxis()->SetTitleSize(0.10);
    frame->GetXaxis()->SetTitleOffset(1.08);
    frame->GetYaxis()->SetLabelSize(0.085); frame->GetYaxis()->SetTitleSize(0.09);
    frame->GetYaxis()->SetTitleOffset(0.90);
  } else {
    frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTitleSize(0);
    frame->GetYaxis()->SetLabelSize(0.039); frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.45);
  }
  const double left = std::max(0., normMin), right = std::min(xMax, normMax);
  if (right > left) {
    TBox* band = new TBox(left, low, right, high);
    band->SetFillColor(kGray); band->SetFillStyle(3004);
    band->SetLineColor(kGray + 1); band->SetLineWidth(1); band->Draw();
  }
  TLine* unity = new TLine(0, 1, xMax, 1);
  unity->SetLineStyle(2); unity->SetLineColor(kGray + 2); unity->Draw();
  int count = 0;
  for (int region = ratio ? 1 : 0; region < 3; ++region) {
    if (graphs[region]->GetN()) graphs[region]->Draw("PZ SAME");
    count += graphs[region]->GetN();
  }
  if (!count) Label(0.55, 0.54, "No valid bins", ratio ? 0.085 : 0.045, 22);
  pad->RedrawAxis();
}

TCanvas* Page(const Table& table, int method, int page, const char* energy,
              double normMin, double normMax, double xMax) {
  const TString name = TString::Format("canvas_%s_page%d", kMethods[method], page + 1);
  TCanvas* canvas = new TCanvas(name, name, 1500, 860);
  canvas->SetFillColor(kWhite);
  Label(0.5, 0.967, TString::Format("%s  |  #Lambda-^{4}He  |  %s", energy,
        kMethodTitles[method]), 0.029, 22);
  Label(0.5, 0.929, TString::Format("Mass-window comparison; each CF independently normalized in %.2f-%.2f GeV/c",
        normMin, normMax), 0.0195, 22);
  for (int column = 0; column < 3; ++column) {
    const int group = page * 3 + column;
    const double left = column / 3., right = (column + 1) / 3.;
    TPad* top = new TPad(TString::Format("top_%s_%d", name.Data(), column), "", left, 0.315, right, 0.902);
    TPad* bottom = new TPad(TString::Format("bottom_%s_%d", name.Data(), column), "", left, 0.063, right, 0.315);
    top->SetLeftMargin(0.165); top->SetRightMargin(0.035); top->SetBottomMargin(0.025); top->SetTopMargin(0.082);
    bottom->SetLeftMargin(0.165); bottom->SetRightMargin(0.035); bottom->SetBottomMargin(0.29); bottom->SetTopMargin(0.035);
    top->SetFillColor(kWhite); bottom->SetFillColor(kWhite);
    canvas->cd(); top->Draw(); bottom->Draw();
    TGraphErrors* cf[3]; TGraphErrors* ratios[3];
    for (int region = 0; region < 3; ++region) {
      cf[region] = Graph(table, kMethods[method], kGroups[group], region, false, xMax);
      ratios[region] = Graph(table, kMethods[method], kGroups[group], region, true, xMax);
    }
    Frame(top, TString::Format("cf_frame_%s_%d", name.Data(), column), cf, false, normMin, normMax, xMax);
    Label(0.565, 0.971, TString::Format("Centrality %s", kGroupTitles[group]), 0.041, 22);
    TLegend* legend = new TLegend(0.48, 0.724, 0.946, 0.905);
    legend->SetTextFont(42); legend->SetTextSize(0.032);
    legend->SetFillStyle(0); legend->SetBorderSize(0);
    for (int region = 0; region < 3; ++region) {
      TString title(kRegionTitles[region]);
      if (!cf[region]->GetN()) title += " (missing)";
      legend->AddEntry(cf[region], title, "pe");
    }
    legend->Draw();
    Frame(bottom, TString::Format("ratio_frame_%s_%d", name.Data(), column), ratios, true, normMin, normMax, xMax);
    canvas->cd();
  }
  Label(0.018, 0.025, "Normalization range shaded where visible. Lower panels: low/full (red), high/full (blue).", 0.017, 12);
  Label(0.982, 0.025, "Invalid / failed bins omitted; ratios use supplied covariance errors.", 0.017, 32);
  canvas->Modified(); canvas->Update();
  return canvas;
}

void RequireAbsent(const TString& path) {
  if (!gSystem->AccessPathName(path.Data()))
    throw std::runtime_error(std::string("Refusing to overwrite existing output: ") + path.Data());
}

void Run(const char* rowsPath, const char* outputDir, const char* energyLabel,
         double normMin, double normMax, double xMax) {
  if (!rowsPath || !*rowsPath || !outputDir || !*outputDir || !energyLabel || !*energyLabel)
    throw std::runtime_error("All path and label arguments are required");
  if (!TMath::Finite(normMin) || !TMath::Finite(normMax) || !TMath::Finite(xMax) ||
      normMin < 0 || normMax <= normMin || xMax <= 0)
    throw std::runtime_error("Explicit finite normalization and display limits are required");
  const Table table = Read(rowsPath);
  const TString directory(outputDir), canvasPath = directory + "/plotcanvases.root";
  RequireAbsent(canvasPath);
  for (int method = 0; method < 5; ++method) {
    RequireAbsent(directory + "/mass_halves_" + kMethods[method] + ".pdf");
    RequireAbsent(directory + "/mass_halves_" + kMethods[method] + ".png");
  }
  if (gSystem->AccessPathName(directory.Data()) && gSystem->mkdir(directory.Data(), kTRUE) != 0)
    throw std::runtime_error("Cannot create output directory");
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(3);
  gStyle->SetCanvasColor(kWhite); gStyle->SetPadColor(kWhite);
  TH1::AddDirectory(kFALSE);
  TFile output(canvasPath, "NEW");
  if (output.IsZombie()) throw std::runtime_error("Cannot create canvas ROOT output in NEW mode");
  TNamed provenance("rows_path", rowsPath); provenance.Write();
  const TString parameters = TString::Format("energy=%s; norm=[%.17g,%.17g]; xMax=%.17g; values_and_errors_supplied_in_TSV",
                                             energyLabel, normMin, normMax, xMax);
  TNamed definition("plot_parameters", parameters.Data()); definition.Write();
  for (int method = 0; method < 5; ++method) {
    const TString pdf = directory + "/mass_halves_" + kMethods[method] + ".pdf";
    for (int page = 0; page < 2; ++page) {
      TCanvas* canvas = Page(table, method, page, energyLabel, normMin, normMax, xMax);
      if (page == 0) canvas->Print(pdf + "[");
      canvas->Print(pdf);
      if (page == 0) canvas->Print(directory + "/mass_halves_" + kMethods[method] + ".png");
      if (page == 1) canvas->Print(pdf + "]");
      output.cd(); canvas->Write();
      delete canvas;
    }
  }
  output.Close();
  std::cout << "[mass-half-plot] Wrote 5 two-page PDFs, 5 PNGs and 10 saved canvases to "
            << outputDir << std::endl;
}

} // namespace MassHalfPlot

void DrawMassHalfCF(const char* rowsPath, const char* outputDir, const char* energyLabel,
                    double normMin, double normMax, double xMax) {
  try {
    MassHalfPlot::Run(rowsPath, outputDir, energyLabel, normMin, normMax, xMax);
  } catch (const std::exception& error) {
    std::cerr << "[mass-half-plot] ERROR: " << error.what() << std::endl;
    gSystem->Exit(2);
  }
}
