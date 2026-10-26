#ifndef KF_SUMMARY_SUPPORT_H
#define KF_SUMMARY_SUPPORT_H
// Downstream-only safety/audit support. Input ROOT files remain read-only.
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TSystem.h"
#include "TPaveText.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
namespace KFReport {
static std::ofstream audit;
static double meanMass = 0, sigmaMass = 0;
static double massFitSettings[9] = {0};
static std::map<std::string, std::string> llWarnings;
inline void Init(const char* path, double mean, double sigma, const double* fitSettings) {
  audit.open(path, std::ios::out | std::ios::trunc);
  if (!audit.good()) { std::cerr << "Cannot create audit " << path << std::endl; gSystem->Exit(31); }
  audit << "type\tobject\tstatus\tcovariance\tboundary\treason\tvalue\terror\tnormalization\n" << std::setprecision(17);
  meanMass = mean; sigmaMass = sigma;
  for (int i=0; i<9; ++i) massFitSettings[i] = fitSettings[i];
}
inline void Record(const char* type, const char* object, int status, int covariance,
                   bool boundary, const char* reason, double value, double error, double normalization) {
  audit << type << '\t' << object << '\t' << status << '\t' << covariance << '\t'
        << boundary << '\t' << reason << '\t' << value << '\t' << error << '\t' << normalization << std::endl;
}
inline bool Boundary(TF1* f) {
  for (int p = 0; p < f->GetNpar(); ++p) {
    double lo = 0, hi = 0; f->GetParLimits(p, lo, hi);
    if (hi > lo) {
      const double eps = (hi-lo)*1.e-4;
      if (f->GetParameter(p) <= lo+eps || f->GetParameter(p) >= hi-eps) return true;
    }
  }
  return false;
}
// Legacy gaus+pol2 model and interval. Total/background independence is an
// approximation; sliding-window cross-bin covariance is not propagated.
inline void Purity(TH1D* h, double minM, double maxM, double& p, double& ep,
                   double& signal, double& esignal, double& total, double& etotal) {
  p = ep = signal = esignal = total = etotal = 0;
  if (!h || h->GetEntries() < massFitSettings[7]) {
    Record("mass_fit", h ? h->GetName() : "missing", -1, -1, false,
           "insufficient_entries", 0, 0, h ? h->GetEntries() : 0); return;
  }
  const int lo = h->FindBin(minM), hi = h->FindBin(maxM);
  total = h->IntegralAndError(lo, hi, etotal);
  TF1 f("kf_mass_total", "gaus(0)+pol2(3)", massFitSettings[0], massFitSettings[1]);
  f.SetParameters(h->GetMaximum(), meanMass, massFitSettings[6], 0, 0, 0);
  f.SetParLimits(1, massFitSettings[2], massFitSettings[3]);
  f.SetParLimits(2, massFitSettings[4], massFitSettings[5]);
  TFitResultPtr r = h->Fit(&f, "SQN", "", massFitSettings[0], massFitSettings[1]);
  const int status = int(r), covStatus = r.Get() ? r->CovMatrixStatus() : -1;
  const bool boundary = Boundary(&f);
  if (!r.Get() || !r->IsValid() || status != 0 || covStatus < massFitSettings[8] || total <= 0) {
    Record("mass_fit", h->GetName(), status, covStatus, boundary,
           "invalid_fit_or_covariance", 0, 0, total); return;
  }
  TF1 bg("kf_mass_bg", "pol2", massFitSettings[0], massFitSettings[1]);
  bg.SetParameters(f.GetParameter(3), f.GetParameter(4), f.GetParameter(5));
  const TMatrixDSym cov = r->GetCovarianceMatrix(); double c[9];
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) c[3*i+j] = cov(i+3,j+3);
  const double low = h->GetBinLowEdge(lo), high = h->GetBinLowEdge(hi)+h->GetBinWidth(hi);
  const double width = h->GetBinWidth(lo);
  const double background = bg.Integral(low, high)/width;
  const double ebackground = bg.IntegralError(low, high, bg.GetParameters(), c)/width;
  if (!TMath::Finite(background) || !TMath::Finite(ebackground) || background < 0 || background >= total) {
    Record("mass_fit", h->GetName(), status, covStatus, boundary,
           "nonphysical_background", background, ebackground, total); return;
  }
  signal = total-background; esignal = TMath::Sqrt(etotal*etotal+ebackground*ebackground);
  p = signal/total;
  ep = TMath::Sqrt(TMath::Power(ebackground/total,2)+TMath::Power(background*etotal/(total*total),2));
  Record("mass_fit", h->GetName(), status, covStatus, boundary,
         boundary ? "valid_boundary_warning" : "valid", p, ep, total);
}
// Zero-error failed-fit bins are absent, not plotted as physical zero points.
inline void DrawMeasured(TH1* h, const char* option) {
  TString opt(option); opt.ToUpper();
  if (!opt.Contains("SAME")) {
    h->Draw("AXIS");
    TString title(h->GetTitle()); const Ssiz_t separator = title.Index(";");
    if (separator >= 0) title.Resize(separator);
    TPaveText* label = new TPaveText(.12,.91,.91,.99,"NDC");
    label->SetFillColor(kWhite); label->SetBorderSize(0); label->SetTextSize(.04);
    label->AddText(title); label->Draw();
  }
  TGraphErrors* graph = new TGraphErrors(); graph->SetName(Form("%s_valid_points", h->GetName()));
  graph->SetMarkerStyle(h->GetMarkerStyle()); graph->SetMarkerSize(h->GetMarkerSize());
  graph->SetMarkerColor(h->GetMarkerColor()); graph->SetLineColor(h->GetLineColor());
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double y = h->GetBinContent(b), e = h->GetBinError(b);
    if (!TMath::Finite(y) || !TMath::Finite(e) || e <= 0) {
      Record("masked_bin", h->GetName(), b, -1, false, "no_valid_error", y, e, h->GetBinCenter(b)); continue;
    }
    const int n = graph->GetN(); graph->SetPoint(n, h->GetBinCenter(b), y); graph->SetPointError(n, 0, e);
  }
  graph->Draw("P SAME");
}
// Legacy repeated fits are all recorded; use the last row for each function.
inline void FitLL(TH1* h, TF1* f, const char* options) {
  int points = 0;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double x = h->GetBinCenter(b);
    if (x >= f->GetXmin() && x <= f->GetXmax() && h->GetBinError(b) > 0 &&
        TMath::Finite(h->GetBinContent(b)) && TMath::Finite(h->GetBinError(b))) ++points;
  }
  if (points <= f->GetNpar()) {
    llWarnings[f->GetName()] = "unavailable: insufficient valid points";
    Record("ll_fit", f->GetName(), -1, -1, false, "insufficient_valid_points", 0, 0, points); return;
  }
  TString opts(options); opts += "SQN"; TFitResultPtr r = h->Fit(f, opts);
  const int status = int(r), covariance = r.Get() ? r->CovMatrixStatus() : -1;
  const bool boundary = Boundary(f);
  llWarnings[f->GetName()] = status == 0 && r.Get() && r->IsValid() && covariance >= 2 && !boundary
    ? "" : "unreliable: fit/covariance/boundary check";
  if (llWarnings[f->GetName()].empty()) f->Draw("SAME");
  Record("ll_fit", f->GetName(), status, covariance, boundary,
         status == 0 && r.Get() && r->IsValid() && covariance >= 2 && !boundary ? "valid_exploratory" : "unreliable_fit",
         f->GetParameter(0), f->GetParError(0), f->GetNDF());
  for (int p = 0; p < f->GetNpar(); ++p)
    Record("ll_parameter", Form("%s:%s", f->GetName(), f->GetParName(p)), status, covariance, boundary,
           "conditional_fit_parameter", f->GetParameter(p), f->GetParError(p), f->GetChisquare());
}
// Replaces the legacy parameter text; auto-fit text size avoids sub-pad clipping.
inline void DrawLLWarnings(TF1* fixed, TF1* free) {
  TF1* functions[2] = {fixed, free};
  for (int i=0; i<2; ++i) {
    TF1* f = functions[i]; const char* label = i ? "Free" : "Fixed";
    const std::string warning = llWarnings[f->GetName()];
    TPaveText* box = new TPaveText(.43, i ? .36 : .57, .88, i ? .56 : .74, "NDC");
    box->SetFillColor(kWhite); box->SetBorderSize(0); box->SetTextFont(42);
    box->SetTextColor(i ? kRed+1 : kBlue); box->SetTextSize(.018); box->SetTextAlign(12);
    if (!warning.empty()) {
      box->AddText(Form("%s: %s", label, warning.find("unavailable") == 0 ? "no fit" : "unreliable"));
      box->AddText("See audit");
    } else {
      box->AddText(label);
      box->AddText(Form("R = %.2f +/- %.2f fm", f->GetParameter(0), f->GetParError(0)));
      box->AddText(Form("chi2/ndf = %.1f/%d", f->GetChisquare(), f->GetNDF()));
      if (i) {
        box->AddText(Form("a0 = %.2f +/- %.2f fm", f->GetParameter(1), f->GetParError(1)));
        box->AddText(Form("reff = %.2f +/- %.2f fm", f->GetParameter(2), f->GetParError(2)));
      }
    }
    box->Draw();
  }
}
} // namespace KFReport
#endif
