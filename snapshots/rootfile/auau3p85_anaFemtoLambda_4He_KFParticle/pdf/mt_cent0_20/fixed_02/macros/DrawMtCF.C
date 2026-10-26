// Rendering only, ROOT 5 / ACLiC compatible. No physics calculation or ROOT input.
// TSV: method slice label low high value error
// low/high are k* bin edges; label is the supplied mT range without whitespace.
// The last centrality label is optional only for legacy 0-10% calls:
// DrawMtCF(tsvPath, outputDir, normMin, normMax, xMax, gridPath, edgesPath, centralityLabel)
// Grid TSV: sample kbin mtbin klow khigh mtlow mthigh count (complete SE/ME grids).
// Edge text: exactly five strictly increasing, saved-axis-aligned mT edges.
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2D.h"
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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace MtCFPlot {

const char* const kMethods[] = {"constant_purity", "method5_common_purity"};
const int kSliceColors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
const int kSliceMarkers[] = {20, 21, 24, 25};

struct Row {
  int method, slice;
  std::string label;
  double low, high, value, error;
};
struct Input {
  std::vector<Row> rows;
  std::string labels[4];
  bool methods[2];
  Input() { methods[0] = methods[1] = false; }
};
struct GridRow {
  int sample, kbin, mtbin;
  double klow, khigh, mtlow, mthigh, count;
};
struct GridInput {
  std::vector<GridRow> rows;
  std::map<int, std::pair<double, double> > kaxis, mtaxis;
  double edges[5];
};

double Number(const std::string& token, unsigned line) {
  char* end = 0;
  errno = 0;
  const double value = std::strtod(token.c_str(), &end);
  if (end == token.c_str() || !end || *end != '\0' || errno == ERANGE) {
    std::ostringstream message;
    message << "Malformed numeric field on TSV line " << line << ": " << token;
    throw std::runtime_error(message.str());
  }
  return value;
}

Input Read(const char* path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error(std::string("Cannot read TSV: ") + path);
  const char* const header[] = {"method", "slice", "label", "low", "high", "value", "error"};
  Input input;
  bool headerSeen = false;
  std::string line;
  unsigned lineNo = 0;
  std::set<std::string> seen;
  while (std::getline(file, line)) {
    ++lineNo;
    if (line.empty() || line[0] == '#') continue;
    std::istringstream parser(line);
    std::vector<std::string> fields;
    std::string field;
    while (parser >> field) fields.push_back(field);
    if (fields.empty()) continue;
    if (fields.size() != 7) throw std::runtime_error("Expected exactly 7 TSV columns");
    if (!headerSeen) {
      for (unsigned i = 0; i < 7; ++i)
        if (fields[i] != header[i]) throw std::runtime_error("Unexpected TSV header/schema");
      headerSeen = true;
      continue;
    }
    Row row;
    if (fields[0] == kMethods[0]) row.method = 0;
    else if (fields[0] == kMethods[1]) row.method = 1;
    else throw std::runtime_error("Unknown method: " + fields[0]);
    const double slice = Number(fields[1], lineNo);
    if (!TMath::Finite(slice) || slice < 0 || slice > 3 || slice != TMath::Floor(slice))
      throw std::runtime_error("Slice must be an integer from 0 to 3");
    row.slice = static_cast<int>(slice);
    row.label = fields[2];
    if (!input.labels[row.slice].empty() && input.labels[row.slice] != row.label)
      throw std::runtime_error("Conflicting mT labels for the same slice");
    input.labels[row.slice] = row.label;
    row.low = Number(fields[3], lineNo); row.high = Number(fields[4], lineNo);
    row.value = Number(fields[5], lineNo); row.error = Number(fields[6], lineNo);
    if (!TMath::Finite(row.low) || !TMath::Finite(row.high) || row.low >= row.high)
      throw std::runtime_error("Invalid finite k* bin edges");
    if (TMath::Finite(row.error) && row.error < 0)
      throw std::runtime_error("Negative uncertainty in TSV");
    std::ostringstream key;
    key.precision(17);
    key << row.method << '/' << row.slice << '/' << row.low << '/' << row.high;
    if (!seen.insert(key.str()).second) throw std::runtime_error("Duplicate TSV bin: " + key.str());
    input.methods[row.method] = true;
    input.rows.push_back(row);
  }
  if (!headerSeen || input.rows.empty()) throw std::runtime_error("TSV contains no data rows");
  return input;
}

void AxisBin(std::map<int, std::pair<double, double> >& axis, int bin, double low, double high) {
  if (!TMath::Finite(low) || !TMath::Finite(high) || low >= high)
    throw std::runtime_error("Invalid grid axis edges");
  const std::pair<double, double> edges(low, high);
  if (axis.count(bin) && axis[bin] != edges)
    throw std::runtime_error("Conflicting grid edges for the same bin index");
  axis[bin] = edges;
}

std::vector<double> AxisEdges(const std::map<int, std::pair<double, double> >& axis) {
  if (axis.empty()) throw std::runtime_error("Empty grid axis");
  std::vector<double> result;
  int previous = axis.begin()->first - 1;
  for (std::map<int, std::pair<double, double> >::const_iterator bin = axis.begin(); bin != axis.end(); ++bin) {
    if (bin->first != previous + 1) throw std::runtime_error("Grid axis has missing bin indices");
    if (!result.empty() && std::fabs(result.back() - bin->second.first) > 1.e-9)
      throw std::runtime_error("Grid axis is not contiguous");
    if (result.empty()) result.push_back(bin->second.first);
    result.push_back(bin->second.second); previous = bin->first;
  }
  return result;
}

GridInput ReadGrid(const char* path, const char* edgesPath) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error(std::string("Cannot read grid TSV: ") + path);
  const char* const header[] = {"sample", "kbin", "mtbin", "klow", "khigh", "mtlow", "mthigh", "count"};
  GridInput grid;
  bool headerSeen = false;
  unsigned lineNo = 0, totals[2] = {0, 0};
  std::set<std::string> seen;
  std::string line;
  while (std::getline(file, line)) {
    ++lineNo;
    if (line.empty() || line[0] == '#') continue;
    std::istringstream parser(line);
    std::vector<std::string> fields; std::string field;
    while (parser >> field) fields.push_back(field);
    if (fields.empty()) continue;
    if (fields.size() != 8) throw std::runtime_error("Expected exactly 8 grid TSV columns");
    if (!headerSeen) {
      for (unsigned i = 0; i < 8; ++i)
        if (fields[i] != header[i]) throw std::runtime_error("Unexpected grid TSV header/schema");
      headerSeen = true; continue;
    }
    GridRow row;
    if (fields[0] == "SE") row.sample = 0;
    else if (fields[0] == "ME") row.sample = 1;
    else throw std::runtime_error("Unknown grid sample: " + fields[0]);
    const double kb = Number(fields[1], lineNo), mb = Number(fields[2], lineNo);
    if (!TMath::Finite(kb) || !TMath::Finite(mb) || kb < 0 || mb < 0 ||
        kb > 1000000 || mb > 1000000 || kb != TMath::Floor(kb) || mb != TMath::Floor(mb))
      throw std::runtime_error("Invalid grid bin index");
    row.kbin = static_cast<int>(kb); row.mtbin = static_cast<int>(mb);
    row.klow = Number(fields[3], lineNo); row.khigh = Number(fields[4], lineNo);
    row.mtlow = Number(fields[5], lineNo); row.mthigh = Number(fields[6], lineNo);
    row.count = Number(fields[7], lineNo);
    if (!TMath::Finite(row.count) || row.count < 0) throw std::runtime_error("Nonfinite or negative grid count");
    AxisBin(grid.kaxis, row.kbin, row.klow, row.khigh);
    AxisBin(grid.mtaxis, row.mtbin, row.mtlow, row.mthigh);
    std::ostringstream key; key << row.sample << '/' << row.kbin << '/' << row.mtbin;
    if (!seen.insert(key.str()).second) throw std::runtime_error("Duplicate grid cell");
    ++totals[row.sample]; grid.rows.push_back(row);
  }
  if (!headerSeen || grid.rows.empty()) throw std::runtime_error("Empty grid TSV");
  const std::vector<double> kEdges = AxisEdges(grid.kaxis), mtEdges = AxisEdges(grid.mtaxis);
  const unsigned expected = grid.kaxis.size() * grid.mtaxis.size();
  if (totals[0] != expected || totals[1] != expected)
    throw std::runtime_error("SE and ME must both contain complete grids, including zero cells");
  std::ifstream edgesFile(edgesPath);
  if (!edgesFile) throw std::runtime_error("Cannot read mT slice edges");
  for (int i = 0; i < 5; ++i) {
    std::string token;
    if (!(edgesFile >> token)) throw std::runtime_error("Expected exactly five mT slice edges");
    grid.edges[i] = Number(token, 1);
    if (!TMath::Finite(grid.edges[i]) || (i && grid.edges[i] <= grid.edges[i-1]))
      throw std::runtime_error("mT slice edges must be finite and increasing");
    bool aligned = false;
    for (unsigned j = 0; j < mtEdges.size(); ++j)
      if (std::fabs(grid.edges[i] - mtEdges[j]) < 1.e-9) aligned = true;
    if (!aligned) throw std::runtime_error("mT slice edge is not aligned with the saved grid axis");
  }
  std::string extra;
  if (edgesFile >> extra) throw std::runtime_error("Extra fields in mT slice edge file");
  return grid;
}

void Text(double x, double y, const char* text, double size, int align = 12, int color = kBlack) {
  TLatex* label = new TLatex(x, y, text);
  label->SetNDC(); label->SetTextFont(42); label->SetTextSize(size);
  label->SetTextAlign(align); label->SetTextColor(color); label->Draw();
}

TString SliceTitle(const Input& input, int slice) {
  if (input.labels[slice].empty()) return TString::Format("m_{T} slice %d (range missing)", slice + 1);
  return TString::Format("m_{T} = %s GeV/c^{2}", input.labels[slice].c_str());
}

TGraphErrors* Graph(const Input& input, int method, int slice, double xmax, bool overlay) {
  TGraphErrors* graph = new TGraphErrors();
  graph->SetName(TString::Format("mtcf_%s_slice%d_%s", kMethods[method], slice, overlay ? "overlay" : "compare"));
  graph->SetMarkerColor(overlay ? kSliceColors[slice] : method == 0 ? kBlack : kRed + 1);
  graph->SetLineColor(graph->GetMarkerColor());
  graph->SetMarkerStyle(overlay ? kSliceMarkers[slice] : method == 0 ? 20 : 24);
  graph->SetMarkerSize(overlay ? 0.85 : method == 0 ? 0.8 : 1.05);
  graph->SetLineWidth(1);
  for (unsigned i = 0; i < input.rows.size(); ++i) {
    const Row& row = input.rows[i];
    const double x = 0.5 * (row.low + row.high);
    if (row.method != method || row.slice != slice || x < 0 || x > xmax ||
        !TMath::Finite(row.value) || !TMath::Finite(row.error)) continue;
    const int point = graph->GetN();
    graph->SetPoint(point, x, row.value); graph->SetPointError(point, 0, row.error);
  }
  return graph;
}

void Bounds(const std::vector<TGraphErrors*>& graphs, double& low, double& high) {
  low = high = 1.;
  for (unsigned graph = 0; graph < graphs.size(); ++graph) {
    for (int i = 0; i < graphs[graph]->GetN(); ++i) {
      double x = 0, y = 0;
      graphs[graph]->GetPoint(i, x, y);
      const double error = graphs[graph]->GetErrorY(i);
      if (!TMath::Finite(y - error) || !TMath::Finite(y + error))
        throw std::runtime_error("Finite value/error overflowed the display bounds");
      low = std::min(low, y - error); high = std::max(high, y + error);
    }
  }
  double span = high - low;
  if (span < 1.e-8) span = 0.2;
  low -= 0.10 * span; high += 0.26 * span;
  const double magnitude = TMath::Power(10., TMath::Floor(TMath::Log10((high - low) / 5.)));
  const double relative = (high - low) / (5. * magnitude);
  const double step = magnitude * (relative > 5. ? 10. : relative > 2. ? 5. : relative > 1. ? 2. : 1.);
  low = step * TMath::Floor(low / step);
  high = step * TMath::Ceil(high / step);
}

void Plot(TPad* pad, const char* name, const std::vector<TGraphErrors*>& graphs,
          double normMin, double normMax, double xmax) {
  pad->cd();
  double low = 0, high = 0;
  Bounds(graphs, low, high);
  TH1F* frame = pad->DrawFrame(0, low, xmax, high);
  frame->SetName(name); frame->SetTitle("");
  frame->GetXaxis()->SetTitle("k* [GeV/c]"); frame->GetYaxis()->SetTitle("CF");
  frame->GetXaxis()->SetLabelFont(42); frame->GetYaxis()->SetLabelFont(42);
  frame->GetXaxis()->SetTitleFont(42); frame->GetYaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetLabelSize(0.035); frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.041); frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetTitleOffset(1.10); frame->GetYaxis()->SetTitleOffset(1.10);
  frame->GetXaxis()->SetNdivisions(505); frame->GetYaxis()->SetNdivisions(505);
  const double left = std::max(0., normMin), right = std::min(xmax, normMax);
  if (right > left) {
    TBox* shade = new TBox(left, low, right, high);
    shade->SetFillStyle(3004); shade->SetFillColor(kGray); shade->SetLineColor(kGray + 1); shade->Draw();
  }
  TLine* unity = new TLine(0, 1, xmax, 1);
  unity->SetLineStyle(2); unity->SetLineColor(kGray + 2); unity->Draw();
  int points = 0;
  for (unsigned i = 0; i < graphs.size(); ++i) {
    points += graphs[i]->GetN();
    if (graphs[i]->GetN()) graphs[i]->Draw("PZ SAME");
  }
  if (!points) Text(0.55, 0.52, "No valid bins", 0.047, 22);
  pad->RedrawAxis();
}

TPad* Pad(TCanvas* canvas, const char* name, double x1, double y1, double x2, double y2) {
  canvas->cd();
  TPad* pad = new TPad(name, "", x1, y1, x2, y2);
  pad->SetFillColor(kWhite); pad->SetLeftMargin(0.115); pad->SetRightMargin(0.035);
  pad->SetTopMargin(0.14); pad->SetBottomMargin(0.14); pad->Draw();
  return pad;
}

void Header(TCanvas* canvas, bool hasMethod5, double normMin, double normMax, const char* centralityLabel) {
  canvas->cd();
  Text(0.5, 0.973, TString::Format("#Lambda-^{4}He: m_{T}-selected correlation functions | centrality %s",
                                centralityLabel), 0.027, 22);
  Text(0.5, 0.939, TString::Format("Signal: #pm3#sigma; normalization %.2f-%.2f GeV/c; no hidden clipping of finite errors",
                                                normMin, normMax), 0.018, 22);
  if (hasMethod5)
    Text(0.5, 0.910, "Shared inclusive P(k*) (Method 5 approximation): purity is NOT fitted separately in each m_{T} slice",
         0.018, 22, kRed + 1);
  else
    Text(0.5, 0.910, "Constant-purity result only; Method 5 approximation has not been calculated", 0.018, 22);
}

TCanvas* Comparison(const Input& input, double normMin, double normMax, double xmax, const char* centralityLabel) {
  TCanvas* canvas = new TCanvas("mt_cf_comparison", "mT CF method comparison", 1400, 1050);
  canvas->SetFillColor(kWhite);
  Header(canvas, input.methods[1], normMin, normMax, centralityLabel);
  for (int slice = 0; slice < 4; ++slice) {
    const int column = slice % 2, row = slice / 2;
    const double x1 = 0.5 * column, x2 = x1 + 0.5;
    const double y2 = row == 0 ? 0.88 : 0.46, y1 = row == 0 ? 0.46 : 0.04;
    TPad* pad = Pad(canvas, TString::Format("mt_compare_slice%d", slice), x1, y1, x2, y2);
    std::vector<TGraphErrors*> graphs;
    int methodIndex[2] = {-1, -1};
    for (int method = 0; method < 2; ++method) {
      if (!input.methods[method]) continue;
      methodIndex[method] = graphs.size();
      graphs.push_back(Graph(input, method, slice, xmax, false));
    }
    Plot(pad, TString::Format("mt_compare_frame%d", slice), graphs, normMin, normMax, xmax);
    Text(0.54, 0.955, SliceTitle(input, slice), 0.041, 22);
    TLegend* legend = new TLegend(0.48, 0.69, 0.946, 0.84);
    legend->SetTextFont(42); legend->SetTextSize(0.031); legend->SetFillStyle(0); legend->SetBorderSize(0);
    for (int method = 0; method < 2; ++method) {
      if (methodIndex[method] < 0) continue;
      TGraphErrors* graph = graphs[methodIndex[method]];
      TString label = method == 0 ? "Constant purity" : "Shared inclusive P(k*)";
      if (!graph->GetN()) label += " (no valid bins)";
      legend->AddEntry(graph, label.Data(), "pe");
    }
    legend->Draw();
  }
  canvas->cd();
  Text(0.018, 0.019, "Invalid / missing bins are omitted, not set to zero. Normalization is shaded where visible.", 0.0165, 12);
  canvas->Modified(); canvas->Update();
  return canvas;
}

TCanvas* Overlays(const Input& input, double normMin, double normMax, double xmax, const char* centralityLabel) {
  TCanvas* canvas = new TCanvas("mt_cf_overlays", "mT slice overlays by method", 1500, 860);
  canvas->SetFillColor(kWhite);
  Header(canvas, input.methods[1], normMin, normMax, centralityLabel);
  for (int method = 0; method < 2; ++method) {
    TPad* pad = Pad(canvas, TString::Format("mt_overlay_method%d", method), method * 0.5, 0.06,
                   (method + 1) * 0.5, 0.88);
    std::vector<TGraphErrors*> graphs;
    for (int slice = 0; slice < 4; ++slice) graphs.push_back(Graph(input, method, slice, xmax, true));
    Plot(pad, TString::Format("mt_overlay_frame%d", method), graphs, normMin, normMax, xmax);
    Text(0.54, 0.966, method == 0 ? "Constant purity" : "Shared inclusive P(k*)", 0.038, 22);
    if (method == 1) Text(0.54, 0.925, "(Method 5 approximation)", 0.030, 22, kRed + 1);
    if (!input.methods[method]) {
      Text(0.55, 0.44, "Not calculated", 0.046, 22);
      continue;
    }
    TLegend* legend = new TLegend(0.49, 0.655, 0.94, 0.842);
    legend->SetTextFont(42); legend->SetTextSize(0.028); legend->SetFillStyle(0); legend->SetBorderSize(0);
    for (int slice = 0; slice < 4; ++slice) {
      TString label = SliceTitle(input, slice);
      if (!graphs[slice]->GetN()) label += " (missing)";
      legend->AddEntry(graphs[slice], label.Data(), "pe");
    }
    legend->Draw();
  }
  canvas->cd();
  Text(0.018, 0.023, "The same supplied CF values are reused; this renderer performs no fit or background correction.", 0.0165, 12);
  canvas->Modified(); canvas->Update();
  return canvas;
}

TCanvas* Distribution(const GridInput& grid, double normMin, double normMax, const char* centralityLabel) {
  const std::vector<double> kEdges = AxisEdges(grid.kaxis), mtEdges = AxisEdges(grid.mtaxis);
  if (kEdges.front() > 0 || kEdges.back() < normMax)
    throw std::runtime_error("Grid k* range does not cover zero through normalization maximum");
  TH2D* histograms[2];
  for (int sample = 0; sample < 2; ++sample) {
    histograms[sample] = new TH2D(sample ? "mt_distribution_ME" : "mt_distribution_SE", "",
         kEdges.size()-1, &kEdges[0], mtEdges.size()-1, &mtEdges[0]);
    histograms[sample]->SetDirectory(0);
  }
  double minimum[2] = {0, 0}, maximum[2] = {0, 0};
  for (unsigned i = 0; i < grid.rows.size(); ++i) {
    const GridRow& row = grid.rows[i];
    const int bx = row.kbin-grid.kaxis.begin()->first+1;
    const int by = row.mtbin-grid.mtaxis.begin()->first+1;
    histograms[row.sample]->SetBinContent(bx, by, row.count);
    const double k = 0.5*(row.klow+row.khigh), mt = 0.5*(row.mtlow+row.mthigh);
    if (k >= 0 && k < normMax && mt >= grid.edges[0] && mt < grid.edges[4] && row.count > 0) {
      maximum[row.sample] = std::max(maximum[row.sample], row.count);
      if (minimum[row.sample] == 0 || row.count < minimum[row.sample]) minimum[row.sample] = row.count;
    }
  }
  TCanvas* canvas = new TCanvas("mt_distribution", "mT split map", 1500, 850);
  canvas->SetFillColor(kWhite);
  Text(0.5, 0.973, TString::Format("#Lambda-^{4}He: m_{T} split map | centrality %s", centralityLabel), 0.029, 22);
  Text(0.5, 0.934, TString::Format("Signal: #pm3#sigma; original bins; horizontal lines: slice boundaries; normalization %.2f-%.2f GeV/c",
                                                normMin, normMax), 0.018, 22);
  for (int sample = 0; sample < 2; ++sample) {
    TPad* pad = Pad(canvas, sample ? "mt_distribution_pad_ME" : "mt_distribution_pad_SE",
                   sample*0.5, 0.065, (sample+1)*0.5, 0.90);
    pad->SetRightMargin(0.17); pad->SetLeftMargin(0.135); pad->SetTopMargin(0.09);
    pad->cd();
    TH2D* hist = histograms[sample];
    hist->GetXaxis()->SetTitle("k* [GeV/c]"); hist->GetYaxis()->SetTitle("m_{T} [GeV/c^{2}]");
    hist->GetZaxis()->SetTitle("Pair counts");
    hist->GetXaxis()->SetRangeUser(0, normMax);
    hist->GetYaxis()->SetRangeUser(grid.edges[0], grid.edges[4]);
    hist->GetXaxis()->SetTitleSize(0.038); hist->GetYaxis()->SetTitleSize(0.041);
    hist->GetZaxis()->SetTitleSize(0.034); hist->GetZaxis()->SetTitleOffset(1.25);
    hist->GetXaxis()->SetLabelSize(0.032); hist->GetYaxis()->SetLabelSize(0.032); hist->GetZaxis()->SetLabelSize(0.030);
    hist->GetYaxis()->SetTitleOffset(1.35); hist->GetXaxis()->SetTitleOffset(1.20);
    if (maximum[sample] > 0) {
      hist->SetMinimum(0.5*minimum[sample]); hist->SetMaximum(1.05*maximum[sample]); pad->SetLogz();
    }
    hist->Draw("COLZ");
    for (int edge = 0; edge < 5; ++edge) {
      TLine* line = new TLine(0, grid.edges[edge], normMax, grid.edges[edge]);
      line->SetLineColor(kBlack); line->SetLineStyle(edge == 0 || edge == 4 ? 1 : 2);
      line->SetLineWidth(2); line->Draw();
    }
    TLine* norm = new TLine(normMin, grid.edges[0], normMin, grid.edges[4]);
    norm->SetLineColor(kRed + 1); norm->SetLineStyle(3); norm->SetLineWidth(2); norm->Draw();
    Text(0.49, 0.956, sample ? "Mixed event (ME), log z" : "Same event (SE), log z", 0.040, 22);
    if (maximum[sample] <= 0) Text(0.49, 0.50, "No positive counts; log z unavailable", 0.032, 22);
    pad->RedrawAxis();
  }
  canvas->cd();
  Text(0.018, 0.026, "Black horizontal lines: supplied m_{T} edges. Red vertical line: normalization start. No cut or reweighting is applied by this renderer.", 0.016, 12);
  canvas->Modified(); canvas->Update();
  return canvas;
}

void Absent(const TString& path) {
  if (!gSystem->AccessPathName(path))
    throw std::runtime_error(std::string("Refusing to overwrite existing output: ") + path.Data());
}

void Run(const char* path, const char* outputDir, double normMin, double normMax, double xmax,
         const char* gridPath, const char* edgesPath, const char* centralityLabel) {
  if (!path || !*path || !outputDir || !*outputDir || !gridPath || !*gridPath || !edgesPath || !*edgesPath)
    throw std::runtime_error("Explicit TSV, grid, slice-edge paths and output directory are required");
  if (!centralityLabel || !*centralityLabel)
    throw std::runtime_error("An explicit nonempty centrality label is required");
  if (!TMath::Finite(normMin) || !TMath::Finite(normMax) || !TMath::Finite(xmax) ||
      normMin < 0 || normMax <= normMin || xmax <= 0)
    throw std::runtime_error("Invalid normalization/display arguments");
  const Input input = Read(path);
  const GridInput grid = ReadGrid(gridPath, edgesPath);
  const TString out(outputDir), rootPath = out + "/mt_cf_canvases.root";
  Absent(rootPath);
  const char* const stems[] = {"mt_cf_comparison", "mt_cf_overlays", "mt_distribution"};
  for (int i = 0; i < 3; ++i) {
    Absent(out + "/" + stems[i] + ".pdf"); Absent(out + "/" + stems[i] + ".png");
  }
  if (gSystem->AccessPathName(out) && gSystem->mkdir(out, kTRUE) != 0)
    throw std::runtime_error("Cannot create output directory");
  gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(3);
  gStyle->SetCanvasColor(kWhite); gStyle->SetPadColor(kWhite); TH1::AddDirectory(kFALSE);
  TFile output(rootPath, "NEW");
  if (output.IsZombie()) throw std::runtime_error("Cannot create canvas ROOT output in NEW mode");
  TNamed source("tsv_path", path); source.Write();
  TNamed gridSource("grid_path", gridPath); gridSource.Write();
  TNamed edgeSource("edges_path", edgesPath); edgeSource.Write();
  TNamed centrality("centrality_label", centralityLabel); centrality.Write();
  const TString definition = TString::Format("normalization=[%.17g,%.17g]; xMax=%.17g; constant_purity=%s; method5_common_purity=%s; renderer_only",
                   normMin, normMax, xmax, input.methods[0] ? "present" : "not_calculated",
                   input.methods[1] ? "present" : "not_calculated");
  TNamed config("plot_parameters", definition.Data()); config.Write();
  TNamed approximation("method5_definition",
      "Shared inclusive P(k*) (Method 5 approximation); not an mT-dependent purity determination");
  approximation.Write();
  TCanvas* canvases[3] = {Comparison(input, normMin, normMax, xmax, centralityLabel),
                         Overlays(input, normMin, normMax, xmax, centralityLabel),
                         Distribution(grid, normMin, normMax, centralityLabel)};
  for (int i = 0; i < 3; ++i) {
    canvases[i]->Print(out + "/" + stems[i] + ".pdf");
    canvases[i]->Print(out + "/" + stems[i] + ".png");
    output.cd(); canvases[i]->Write(); delete canvases[i];
  }
  output.Close();
  std::cout << "[mt-cf-plot] Created three one-page PDFs, three PNGs and saved canvases in " << outputDir << std::endl;
}
} // namespace MtCFPlot

void DrawMtCF(const char* tsvPath, const char* outputDir, double normMin, double normMax, double xMax,
              const char* gridPath, const char* edgesPath, const char* centralityLabel = "0-10%") {
  try { MtCFPlot::Run(tsvPath, outputDir, normMin, normMax, xMax, gridPath, edgesPath, centralityLabel); }
  catch (const std::exception& error) {
    std::cerr << "[mt-cf-plot] ERROR: " << error.what() << std::endl;
    gSystem->Exit(2);
  }
}
