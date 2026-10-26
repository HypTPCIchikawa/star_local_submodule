// Display saved native spectra and saved predictions only. Never refits.
void DrawRestartFits(const char*, const char*, const char*, const char*);
#ifndef __CINT__
#include "../nuclear_pid_tof/NuclearFitCommon.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace RestartFitView {

struct Row {
  int index;
  std::string id, target, mode, numerical, reason;
  double r, fitLo, fitHi, mean, sigma, chi2Ndf;
};

static void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

static std::vector<std::string> Split(const std::string& value) {
  std::vector<std::string> fields;
  std::string field;
  std::istringstream input(value);
  while (std::getline(input, field, '\t')) fields.push_back(field);
  if (!value.empty() && value[value.size() - 1] == '\t') fields.push_back("");
  return fields;
}

static double Number(const std::string& value, bool optional) {
  if (optional && (value.empty() || value == "null" || value == "nan" || value == "NA"))
    return std::numeric_limits<double>::quiet_NaN();
  char* end = 0;
  double number = std::strtod(value.c_str(), &end);
  Require(end != value.c_str() && end && *end == 0 && TMath::Finite(number), "Invalid numerical plot metadata: " + value);
  return number;
}

static std::vector<Row> ReadRows(const char* path) {
  std::ifstream input(path);
  Require(input.good(), std::string("Cannot read plot metadata: ") + path);
  std::string line;
  Require(bool(std::getline(input, line)), "Empty plot metadata");
  Require(line == "index\tid\ttarget\tmode\tr\tfitLo\tfitHi\tmean\tsigma\tchi2Ndf\tnumerical\treason",
          "Unexpected plot metadata header");
  std::vector<Row> rows;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::vector<std::string> values = Split(line);
    Require(values.size() == 12, "Plot metadata must have 12 tab-separated columns");
    Row row;
    double index = Number(values[0], false);
    Require(index >= 0 && index <= INT_MAX && index == std::floor(index), "Invalid saved spectrum index");
    row.index = int(index); row.id = values[1]; row.target = values[2]; row.mode = values[3];
    Require(row.target == "d" || row.target == "t" || row.target == "he3" || row.target == "he4", "Unknown target label");
    Require(row.mode == "ALL" || row.mode == "TAG_D" || row.mode == "TAG_T" || row.mode == "TAG_HE3", "Unknown sample label");
    row.r = Number(values[4], false);
    row.fitLo = Number(values[5], false); row.fitHi = Number(values[6], false);
    row.mean = Number(values[7], true); row.sigma = Number(values[8], true);
    row.chi2Ndf = Number(values[9], true); row.numerical = values[10]; row.reason = values[11];
    Require(row.fitLo < row.fitHi, "Invalid local fit range");
    rows.push_back(row);
  }
  Require(!rows.empty(), "No saved spectra requested for plotting");
  return rows;
}

static std::string Name(const std::string& prefix, int index) {
  std::ostringstream name;
  name << prefix << '_' << std::setw(4) << std::setfill('0') << index;
  return name.str();
}

static std::string Format(double value, int precision) {
  if (!TMath::Finite(value)) return "not available";
  std::ostringstream out;
  out << std::fixed << std::setprecision(precision) << value;
  return out.str();
}

static const char* Label(const std::string& target) {
  if (target == "d") return "d";
  if (target == "t") return "t";
  if (target == "he3") return "^{3}He";
  return "^{4}He";
}

static int Color(const std::string& target) {
  if (target == "d") return kBlue + 1;
  if (target == "t") return kOrange + 7;
  if (target == "he3") return kGreen + 2;
  return kMagenta + 1;
}

static void Latex(double x, double y, const std::string& value, double size) {
  TLatex text;
  text.SetNDC(); text.SetTextFont(42); text.SetTextSize(size);
  TLatex* drawn = text.DrawLatex(x, y, value.c_str());
  drawn->SetBit(TObject::kCanDelete);
}

static void Plain(double x, double y, const std::string& value, double size) {
  TText text;
  text.SetNDC(); text.SetTextFont(42); text.SetTextSize(size);
  TText* drawn = text.DrawTextNDC(x, y, value.c_str());
  drawn->SetBit(TObject::kCanDelete);
}

static TH1D* Clone(TFile& input, const std::string& key, const std::string& unique,
                   int rebin, std::vector<TObject*>& owned, bool required) {
  TH1D* source = dynamic_cast<TH1D*>(input.Get(key.c_str()));
  if (!source) {
    Require(!required, "Missing saved native spectrum: " + key);
    return 0;
  }
  Require(source->GetNbinsX() % rebin == 0, "displayRebin must divide native bin count");
  TH1D* histogram = dynamic_cast<TH1D*>(source->Clone(unique.c_str()));
  Require(histogram != 0, "Cannot clone native display spectrum");
  histogram->SetDirectory(0); histogram->SetBit(TObject::kCanDelete, false);
  if (rebin != 1) histogram->Rebin(rebin);
  owned.push_back(histogram);
  return histogram;
}

static TGraph* MaskedGraphs(const TH1D& histogram, const TH1D& mask, int rebin,
                            int color, int style, std::vector<TObject*>& owned) {
  TGraph* first = 0;
  TGraph* current = 0;
  for (int b = 1; b <= histogram.GetNbinsX(); ++b) {
    // Partial display bins at a fit boundary are not extrapolated or filled.
    bool active = std::fabs(mask.GetBinContent(b) - rebin) < 1e-10;
    double value = histogram.GetBinContent(b);
    if (!active || !TMath::Finite(value) || value < 0) { current = 0; continue; }
    if (!current) {
      current = new TGraph();
      current->SetBit(TObject::kCanDelete, false);
      current->SetLineColor(color); current->SetLineStyle(style); current->SetLineWidth(2);
      owned.push_back(current);
      if (!first) first = current;
    }
    current->SetPoint(current->GetN(), histogram.GetBinCenter(b), value);
  }
  // Draw only completed contiguous segments; never bridge an unobserved gap.
  for (size_t item = 0; item < owned.size(); ++item) {
    TGraph* graph = dynamic_cast<TGraph*>(owned[item]);
    if (graph && graph->GetLineColor() == color && graph->GetLineStyle() == style && graph->GetN())
      graph->Draw("L SAME");
  }
  return first;
}

static void Panel(TFile& input, const Row& row, bool logarithmic,
                   const NuclearFit::Settings& settings, std::vector<TObject*>& owned,
                   int serial) {
  const int rebin = NuclearFit::Integer(settings, "displayRebin");
  const double padding = NuclearFit::Number(settings, "viewPaddingSigma");
  const double headroom = NuclearFit::Number(settings, logarithmic ? "logHeadroom" : "linearHeadroom");
  const double logMinimum = NuclearFit::Number(settings, "logMinimum");
  Require(rebin >= 1 && padding >= 0 && headroom > 1 && logMinimum > 0, "Invalid explicit display settings");
  std::ostringstream unique;
  unique << "restart_display_" << serial;
  TH1D* raw = Clone(input, Name("raw_native", row.index), unique.str() + "_raw", rebin, owned, true);
  TH1D* total = Clone(input, Name("prediction_total", row.index), unique.str() + "_total", rebin, owned, false);
  TH1D* smooth = Clone(input, Name("prediction_smooth", row.index), unique.str() + "_smooth", rebin, owned, false);
  TH1D* target = Clone(input, Name("prediction_" + row.target, row.index), unique.str() + "_target", rebin, owned, false);
  TH1D* mask = Clone(input, Name("fit_mask", row.index), unique.str() + "_mask", rebin, owned, false);
  bool predictions = total && smooth && target && mask;
  Require((!total && !smooth && !target && !mask) || predictions, "Incomplete saved fit prediction set");
  double width = TMath::Finite(row.sigma) && row.sigma > 0 ? row.sigma : 0;
  double xLo = std::max(raw->GetXaxis()->GetXmin(), row.fitLo - padding * width);
  double xHi = std::min(raw->GetXaxis()->GetXmax(), row.fitHi + padding * width);
  Require(xLo < xHi, "Invalid native local display range");
  double maximum = 0;
  for (int b = 1; b <= raw->GetNbinsX(); ++b) {
    if (raw->GetBinCenter(b) < xLo || raw->GetBinCenter(b) > xHi) continue;
    maximum = std::max(maximum, raw->GetBinContent(b) + raw->GetBinError(b));
    if (predictions) maximum = std::max(maximum, total->GetBinContent(b));
  }
  if (!(maximum > 0)) maximum = 1;
  gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.035);
  gPad->SetTopMargin(0.225); gPad->SetBottomMargin(0.13);
  gPad->SetLogy(logarithmic);
  raw->SetTitle(""); raw->SetStats(false);
  raw->SetMarkerStyle(20); raw->SetMarkerSize(0.28); raw->SetLineWidth(1);
  raw->SetLineColor(kBlack); raw->SetMarkerColor(kBlack);
  raw->GetXaxis()->SetTitle("dE/dx [keV/cm]");
  raw->GetXaxis()->SetRangeUser(xLo, xHi);
  std::ostringstream yTitle;
  yTitle << "Counts / " << std::setprecision(4) << raw->GetXaxis()->GetBinWidth(1) << " keV/cm";
  raw->GetYaxis()->SetTitle(yTitle.str().c_str());
  raw->GetXaxis()->SetTitleSize(0.052); raw->GetYaxis()->SetTitleSize(0.048);
  raw->GetXaxis()->SetLabelSize(0.044); raw->GetYaxis()->SetLabelSize(0.043);
  raw->GetYaxis()->SetTitleOffset(1.20);
  raw->SetMinimum(logarithmic ? logMinimum : 0);
  raw->SetMaximum(logarithmic ? std::max(maximum * headroom, logMinimum * 10) : maximum * headroom);
  raw->Draw("E1");
  TGraph* totalLine = 0; TGraph* targetLine = 0; TGraph* smoothLine = 0;
  if (predictions) {
    // Retain only objects created for this panel while drawing its segments.
    std::vector<TObject*> graphs;
    smoothLine = MaskedGraphs(*smooth, *mask, rebin, kGray + 2, 3, graphs);
    targetLine = MaskedGraphs(*target, *mask, rebin, Color(row.target), 2, graphs);
    totalLine = MaskedGraphs(*total, *mask, rebin, kRed + 1, 1, graphs);
    owned.insert(owned.end(), graphs.begin(), graphs.end());
    raw->Draw("E1 SAME");
  }
  TLegend* legend = new TLegend(0.68, 0.55, 0.955, 0.765);
  legend->SetBit(TObject::kCanDelete); legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.032);
  legend->AddEntry(raw, "Native data", "lep");
  if (totalLine) legend->AddEntry(totalLine, "Saved total", "l");
  if (targetLine) legend->AddEntry(targetLine, (std::string(Label(row.target)) + " signal response").c_str(), "l");
  if (smoothLine) legend->AddEntry(smoothLine, "Linear background", "l");
  legend->Draw();
  Latex(0.13, 0.958, std::string(Label(row.target)) + "; R = " + Format(row.r, 3) + " GeV/c; " + row.mode +
        (logarithmic ? "; log" : "; linear"), 0.038);
  Plain(0.13, 0.914, "Numerical: " + row.numerical + "; chi2/ndf = " + Format(row.chi2Ndf, 3), 0.031);
  Plain(0.13, 0.873, "Saved fit [" + Format(row.fitLo, 2) + ", " + Format(row.fitHi, 2) +
        "]; mean = " + Format(row.mean, 3) + "; sigma = " + Format(row.sigma, 3), 0.029);
  std::string reason = row.reason.empty() ? "No additional numerical note" : row.reason;
  // Split long diagnostics, without hiding or truncating their content.
  for (size_t start = 0, line = 0; start < reason.size(); start += 108, ++line)
    Plain(0.13, 0.831 - line * 0.031, reason.substr(start, 108), 0.026);
  Plain(0.13, 0.025, "Diagnostic only: purity not assessed; no refit or normalization", 0.027);
  std::cout << "RESTART_VIEW_PANEL serial=" << serial << " id=" << row.id
            << " target=" << row.target << " mode=" << row.mode << " log=" << logarithmic
            << " xLo=" << xLo << " xHi=" << xHi << " displayRebin=" << rebin
            << " prediction=" << predictions << std::endl;
}

static void Run(const char* rootPath, const char* metricsPath, const char* settingsPath, const char* pdfPath) {
  Require(gSystem->AccessPathName(pdfPath), "Plot output already exists; use a fresh filename");
  NuclearFit::Settings settings = NuclearFit::ReadSettings(settingsPath);
  std::vector<Row> rows = ReadRows(metricsPath);
  TFile input(rootPath, "READ");
  Require(!input.IsZombie(), "Cannot open saved native fit spectra");
  gStyle->SetOptStat(0);
  TCanvas canvas("restart_peak_fits", "Fresh PID peak-response diagnostics", 1400, 1000);
  canvas.Print((std::string(pdfPath) + "[").c_str());
  int pages = 0, panels = 0;
  for (size_t start = 0; start < rows.size(); start += 2) {
    canvas.Clear(); canvas.Divide(2, 2, 0.0001, 0.0001);
    std::vector<TObject*> owned;
    for (size_t offset = 0; offset < 2 && start + offset < rows.size(); ++offset) {
      for (int logarithmic = 0; logarithmic < 2; ++logarithmic) {
        canvas.cd(int(offset * 2) + logarithmic + 1);
        Panel(input, rows[start + offset], logarithmic != 0, settings, owned, panels++);
      }
    }
    canvas.Print(pdfPath); ++pages;
    canvas.Clear();
    for (size_t object = 0; object < owned.size(); ++object) delete owned[object];
  }
  canvas.Print((std::string(pdfPath) + "]").c_str());
  Require(pages == int((rows.size() + 1) / 2) && panels == int(rows.size() * 2), "Unexpected PDF panel/page count");
  std::cout << "RESTART_VIEW_COMPLETE rows=" << rows.size() << " pages=" << pages
            << " panels=" << panels << " expectedPages=" << (rows.size() + 1) / 2
            << " refit=0 normalization=0 purityAssessed=0" << std::endl;
}

}  // namespace RestartFitView
#endif

void DrawRestartFits(const char* rootPath, const char* metricsPath, const char* settingsPath, const char* pdfPath) {
#ifndef __CINT__
  try { RestartFitView::Run(rootPath, metricsPath, settingsPath, pdfPath); }
  catch (const std::exception& error) {
    std::cerr << "RESTART_VIEW_FATAL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
#endif
}
