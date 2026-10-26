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
#include <vector>
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
// Keep the original purity-estimation grid independent of the output CF grid.
// Apply the purity belonging to each original k* bin BEFORE summing those bins
// into the requested CF bin. This avoids assigning an entire 20-MeV bin to one
// 50-MeV cell when it straddles that cell boundary. The same purity error is
// shared inside an output bin, not counted as independent for every fine bin.
// Count--fit and neighboring-fit covariance remain unrecorded approximations.
inline void Method5IndependentGrid(TH2* source, TH1D* signal, TH1D* purity,
                                  int cfRebin, int purityRebin, int halfWidth,
                                  double minM, double maxM) {
  const int nOriginal=source ? source->GetNbinsX() : 0;
  if (!source || !signal || !purity || cfRebin<1 || purityRebin<1 || halfWidth<0 ||
      nOriginal%cfRebin || nOriginal%purityRebin ||
      signal->GetNbinsX()!=nOriginal/cfRebin || purity->GetNbinsX()!=signal->GetNbinsX()) {
    std::cerr<<"Invalid independent Method5 grid"<<std::endl;gSystem->Exit(32);return;
  }
  const int nCoarse=nOriginal/purityRebin;
  std::vector<double> p(nCoarse,0),vP(nCoarse,0);
  for (int c=0;c<nCoarse;++c) {
    const int first=TMath::Max(0,c-halfWidth)*purityRebin+1;
    const int last=(TMath::Min(nCoarse-1,c+halfWidth)+1)*purityRebin;
    TH1D* projection=source->ProjectionY(Form("%s_coarse_purity_%d",signal->GetName(),c),first,last,"e");
    projection->SetDirectory(0);
    double ep=0,ns=0,es=0,nt=0,et=0;
    Purity(projection,minM,maxM,p[c],ep,ns,es,nt,et);vP[c]=ep*ep;
    delete projection;
  }
  const int massLo=source->GetYaxis()->FindBin(minM),massHi=source->GetYaxis()->FindBin(maxM);
  for (int b=1;b<=signal->GetNbinsX();++b) {
    std::map<int,double> n,vN;
    for (int x=(b-1)*cfRebin+1;x<=b*cfRebin;++x) {
      const int c=(x-1)/purityRebin;
      for (int y=massLo;y<=massHi;++y) {
        n[c]+=source->GetBinContent(x,y);
        const double e=source->GetBinError(source->GetBin(x,y));vN[c]+=e*e;
      }
    }
    double total=0,value=0,variance=0;bool valid=true;
    for (std::map<int,double>::const_iterator it=n.begin();it!=n.end();++it) {
      const int c=it->first;const double ni=it->second;
      if ((ni!=0 || vN[c]>0) && (!(p[c]>0 && p[c]<=1) || !TMath::Finite(vP[c]))) valid=false;
      total+=ni;value+=ni*p[c];variance+=p[c]*p[c]*vN[c]+ni*ni*vP[c];
    }
    if (!valid || total<=0 || !TMath::Finite(value) || !TMath::Finite(variance)) {
      Record("method5_output_bin",signal->GetName(),b,-1,false,"unavailable_purity_component",0,0,total);
      continue;
    }
    const double effective=value/total;double vEffective=0;
    for (std::map<int,double>::const_iterator it=n.begin();it!=n.end();++it) {
      const int c=it->first;const double ni=it->second;
      vEffective+=(p[c]-effective)*(p[c]-effective)*vN[c]+ni*ni*vP[c];
    }
    signal->SetBinContent(b,value);signal->SetBinError(b,TMath::Sqrt(variance));
    purity->SetBinContent(b,effective);purity->SetBinError(b,TMath::Sqrt(vEffective)/total);
    Record("method5_output_bin",signal->GetName(),b,-1,false,
           "original_bin_purity_then_cf_sum",effective,TMath::Sqrt(vEffective)/total,total);
  }
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
  // The plotting window must not truncate the independently configured fit
  // interval (e.g. show 0--0.3 but retain the requested LL fit through 0.4).
  TAxis* axis=h->GetXaxis();const int first=axis->GetFirst(),last=axis->GetLast();
  const bool hadRange=axis->TestBit(TAxis::kAxisRange);
  axis->SetRange(0,0);
  TString opts(options); opts += "SQN"; TFitResultPtr r = h->Fit(f, opts);
  if (hadRange) axis->SetRange(first,last); else axis->SetRange(0,0);
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
// Replaces legacy parameter text; short lines and a fixed size avoid ROOT5 PDF clipping.
inline void DrawLLWarnings(TF1* fixed, TF1* free) {
  TF1* functions[2] = {fixed, free};
  for (int i=0; i<2; ++i) {
    TF1* f = functions[i]; const char* label = i ? "Free" : "Fixed";
    const std::string warning = llWarnings[f->GetName()];
    TPaveText* box = new TPaveText(.43, i ? .36 : .57, .88, i ? .56 : .74, "NDC");
    box->SetFillColor(kWhite); box->SetBorderSize(0); box->SetTextFont(42);
    box->SetTextColor(i ? kRed+1 : kBlue); box->SetTextSize(.028); box->SetTextAlign(12);
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
