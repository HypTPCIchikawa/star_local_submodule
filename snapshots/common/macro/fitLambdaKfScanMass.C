// Read-only mass-fit cross-check for scanLambdaKfImp5Cuts.C outputs (ROOT 5).
// All mass windows, fit ranges, resolution bounds and minimum statistics are
// explicit QA arguments; this macro never changes a production selection.
//
// fitLambdaKfScanMass(inputRoot, outputPrefix, histogramPathsCsv, fitRangesCsv,
//                    signalLow, signalHigh, sigmaMin, sigmaMax, minFitEntries)
// histogramPathsCsv: comma-separated ROOT histogram paths (first = reference).
// fitRangesCsv: comma-separated lower:upper pairs in GeV/c^2.
// Outputs: .csv, .root, and _rangeNN.png/.pdf. Existing outputs are not replaced.
//
// Fit model: Gaussian + quadratic polynomial. The polynomial is evaluated about
// the fixed signal-window midpoint for numerical conditioning; this spans the
// same functions as gaus(0)+pol2(3). Fits use unweighted Poisson likelihood with
// bin integration (ROOT options RSQLI0N). S and B are analytic model integrals
// over [signalLow,signalHigh), divided by the common histogram bin width.
// Errors use the full fit covariance (including S-B correlations), conditional
// on this model/range and ignoring correlations between selections. Neither
// S/B nor counts in the window constitute a reconstruction-efficiency estimate.

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TString.h>
#include <TMath.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace LambdaKfMassFit {

struct Range { Double_t low, high; };

struct Result {
  TF1* total;
  TF1* background;
  TF1* signal;
  Int_t status, covarianceQuality;
  Bool_t accepted;
  Double_t nWindow, nWindowError, nFit, s, b, sError, bError;
  Double_t sb, sbError, purity, purityError, covarianceSB;
  Double_t chi2, ndf, mean, sigma;
  std::string flags;
  Result() : total(0), background(0), signal(0), status(-999),
    covarianceQuality(-1), accepted(kFALSE), nWindow(0), nWindowError(0), nFit(0) {
    const Double_t nan = std::numeric_limits<Double_t>::quiet_NaN();
    s = b = sError = bError = sb = sbError = purity = purityError = nan;
    covarianceSB = chi2 = ndf = mean = sigma = nan;
  }
};

static Int_t Fail(const char* message) {
  std::cerr << "[fitLambdaKfScanMass] ERROR: " << message << std::endl;
  return 1;
}

static std::vector<std::string> Split(const char* value) {
  std::vector<std::string> result;
  std::stringstream input(value ? value : "");
  std::string token;
  while (std::getline(input, token, ',')) {
    const std::string::size_type first = token.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) { result.clear(); return result; }
    const std::string::size_type last = token.find_last_not_of(" \t\r\n");
    result.push_back(token.substr(first, last - first + 1));
  }
  return result;
}

static std::string Quote(const std::string& value) {
  std::string result("\"");
  for (std::string::size_type i = 0; i < value.size(); ++i) {
    if (value[i] == '"') result += '"';
    result += value[i];
  }
  return result + "\"";
}

static void Flag(Result& result, const char* flag) {
  if (!result.flags.empty()) result.flags += ";";
  result.flags += flag;
}

// Return the first bin of an exactly aligned edge. Nbins+1 is the upper edge.
static Int_t Edge(const TH1* histogram, Double_t value) {
  for (Int_t bin = 1; bin <= histogram->GetNbinsX() + 1; ++bin)
    if (TMath::Abs(histogram->GetXaxis()->GetBinLowEdge(bin) - value) < 1.e-9)
      return bin;
  return -1;
}

static Double_t Model(Double_t* x, Double_t* p) {
  const Double_t z = (x[0] - p[1]) / p[2];
  const Double_t dx = x[0] - p[6];
  return p[0] * TMath::Exp(-0.5 * z * z) + p[3] + p[4] * dx + p[5] * dx * dx;
}

static Double_t Quadratic(Double_t* x, Double_t* p) {
  const Double_t dx = x[0] - p[3];
  return p[0] + p[1] * dx + p[2] * dx * dx;
}

static Double_t Product(const Double_t* left, const TMatrixDSym& covariance,
                        const Double_t* right) {
  Double_t result = 0.;
  for (Int_t i = 0; i < 7; ++i)
    for (Int_t j = 0; j < 7; ++j)
      result += left[i] * covariance(i, j) * right[j];
  return result;
}

static void FitOne(TH1* histogram, const Range& range, Double_t low,
                   Double_t high, Double_t sigmaMin, Double_t sigmaMax,
                   Double_t minEntries, Int_t index, TFile* output, Result& r) {
  const Double_t width = histogram->GetXaxis()->GetBinWidth(1);
  const Int_t first = Edge(histogram, low), last = Edge(histogram, high) - 1;
  r.nWindow = histogram->IntegralAndError(first, last, r.nWindowError);
  const Int_t fitFirst = Edge(histogram, range.low);
  const Int_t fitLast = Edge(histogram, range.high) - 1;
  r.nFit = histogram->Integral(fitFirst, fitLast);
  if (r.nFit < minEntries || r.nFit <= 0. || fitLast - fitFirst + 1 <= 6) {
    Flag(r, "insufficient_fit_statistics");
    return;
  }

  const Double_t center = 0.5 * (low + high);
  Double_t sideCount = 0.;
  Int_t sideBins = 0, peakBin = first;
  for (Int_t bin = fitFirst; bin <= fitLast; ++bin) {
    if (bin < first || bin > last) { sideCount += histogram->GetBinContent(bin); ++sideBins; }
    if (bin >= first && bin <= last &&
        histogram->GetBinContent(bin) > histogram->GetBinContent(peakBin)) peakBin = bin;
  }
  const Double_t backgroundSeed = sideBins > 0 ? sideCount / sideBins : 0.;
  const Double_t amplitudeSeed = TMath::Max(1., histogram->GetBinContent(peakBin) - backgroundSeed);
  r.total = new TF1(Form("KfMassTotal_%d", index), Model, range.low, range.high, 7);
  r.total->SetParNames("GaussianHeight", "Mean", "Sigma", "B0", "B1", "B2", "FixedCenter");
  r.total->SetParameters(amplitudeSeed, histogram->GetBinCenter(peakBin),
                        TMath::Sqrt(sigmaMin * sigmaMax), backgroundSeed, 0., 0., center);
  // The amplitude upper bound is deliberately remote; it only protects the
  // minimizer. The physical mean and width bounds are explicit caller choices.
  r.total->SetParLimits(0, 0., 1000. * TMath::Max(1., r.nFit));
  r.total->SetParLimits(1, low, high);
  r.total->SetParLimits(2, sigmaMin, sigmaMax);
  r.total->FixParameter(6, center);
  TFitResultPtr fit = histogram->Fit(r.total, "RSQLI0N");
  r.status = Int_t(fit);
  if (!fit.Get()) { Flag(r, "missing_fit_result"); return; }
  r.covarianceQuality = fit->CovMatrixStatus();
  r.chi2 = fit->Chi2();
  r.ndf = fit->Ndf();
  r.mean = r.total->GetParameter(1);
  r.sigma = r.total->GetParameter(2);
  output->cd();
  fit->Write(Form("KfMassFitResult_%d", index));
  r.total->Write();
  if (r.status != 0) Flag(r, "fit_status_nonzero");
  if (r.covarianceQuality != 3) Flag(r, "covariance_not_accurate");

  Double_t p[7];
  r.total->GetParameters(p);
  for (Int_t i = 0; i < 7; ++i)
    if (!TMath::Finite(p[i])) { Flag(r, "nonfinite_parameters"); return; }
  if (p[0] <= 0. || p[2] <= 0.) { Flag(r, "nonpositive_signal_parameters"); return; }
  if (p[1] - low <= 1.e-5 * (high - low) || high - p[1] <= 1.e-5 * (high - low))
    Flag(r, "mean_at_bound");
  if (p[2] - sigmaMin <= 1.e-5 * (sigmaMax - sigmaMin) ||
      sigmaMax - p[2] <= 1.e-5 * (sigmaMax - sigmaMin)) Flag(r, "sigma_at_bound");

  r.background = new TF1(Form("KfMassBackground_%d", index), Quadratic, range.low, range.high, 4);
  r.background->SetParameters(p[3], p[4], p[5], center);
  r.signal = new TF1(Form("KfMassSignal_%d", index), "gaus", range.low, range.high);
  r.signal->SetParameters(p[0], p[1], p[2]);
  Double_t minimumBackground = TMath::Min(r.background->Eval(range.low), r.background->Eval(range.high));
  if (p[5] > 0.) {
    const Double_t minimumX = center - p[4] / (2. * p[5]);
    if (minimumX > range.low && minimumX < range.high)
      minimumBackground = TMath::Min(minimumBackground, r.background->Eval(minimumX));
  }
  if (!TMath::Finite(minimumBackground) || minimumBackground < 0.)
    Flag(r, "negative_background_in_fit_range");

  const Double_t a = (low - p[1]) / p[2], b = (high - p[1]) / p[2];
  const Double_t ea = TMath::Exp(-0.5 * a * a), eb = TMath::Exp(-0.5 * b * b);
  const Double_t base = p[2] * TMath::Sqrt(TMath::Pi() / 2.) *
    (TMath::Erf(b / TMath::Sqrt(2.)) - TMath::Erf(a / TMath::Sqrt(2.))) / width;
  const Double_t dl = low - center, dh = high - center;
  Double_t gradientS[7] = {0., 0., 0., 0., 0., 0., 0.};
  Double_t gradientB[7] = {0., 0., 0., 0., 0., 0., 0.};
  gradientS[0] = base;
  gradientS[1] = p[0] * (ea - eb) / width;
  r.s = p[0] * base;
  gradientS[2] = r.s / p[2] + p[0] * (a * ea - b * eb) / width;
  gradientB[3] = (high - low) / width;
  gradientB[4] = (dh * dh - dl * dl) / (2. * width);
  gradientB[5] = (dh * dh * dh - dl * dl * dl) / (3. * width);
  r.b = gradientB[3] * p[3] + gradientB[4] * p[4] + gradientB[5] * p[5];
  if (!TMath::Finite(r.s) || !TMath::Finite(r.b) || r.s <= 0. || r.b <= 0.)
    Flag(r, "nonpositive_or_nonfinite_window_components");
  r.background->Write();
  r.signal->Write();
  if (!r.flags.empty()) return;

  const TMatrixDSym covariance = fit->GetCovarianceMatrix();
  if (covariance.GetNrows() != 7) { Flag(r, "unexpected_covariance_dimension"); return; }
  Double_t gradientRatio[7], gradientPurity[7];
  for (Int_t i = 0; i < 7; ++i) {
    gradientRatio[i] = gradientS[i] / r.b - r.s * gradientB[i] / (r.b * r.b);
    gradientPurity[i] = (r.b * gradientS[i] - r.s * gradientB[i]) / ((r.s + r.b) * (r.s + r.b));
  }
  const Double_t varianceS = Product(gradientS, covariance, gradientS);
  const Double_t varianceB = Product(gradientB, covariance, gradientB);
  const Double_t varianceRatio = Product(gradientRatio, covariance, gradientRatio);
  const Double_t variancePurity = Product(gradientPurity, covariance, gradientPurity);
  if (!TMath::Finite(varianceS) || !TMath::Finite(varianceB) ||
      !TMath::Finite(varianceRatio) || !TMath::Finite(variancePurity) ||
      varianceS < 0. || varianceB < 0. || varianceRatio < 0. || variancePurity < 0.) {
    Flag(r, "invalid_propagated_covariance");
    return;
  }
  r.sError = TMath::Sqrt(varianceS);
  r.bError = TMath::Sqrt(varianceB);
  r.covarianceSB = Product(gradientS, covariance, gradientB);
  r.sb = r.s / r.b;
  r.sbError = TMath::Sqrt(varianceRatio);
  r.purity = r.s / (r.s + r.b);
  r.purityError = TMath::Sqrt(variancePurity);
  r.accepted = kTRUE;
  r.flags = "ok";
}

static void Style(TH1* histogram, Int_t color, Double_t width) {
  histogram->SetStats(kFALSE);
  histogram->SetLineColor(color);
  histogram->SetMarkerColor(color);
  histogram->SetLineWidth(2);
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitle("M_{p#pi^{-}} [GeV/c^{2}]");
  histogram->GetYaxis()->SetTitle(Form("Candidates / %.3g MeV/c^{2}", width * 1000.));
}

} // namespace LambdaKfMassFit

Int_t fitLambdaKfScanMass(const char* inputRoot, const char* outputPrefix,
                        const char* histogramPathsCsv, const char* fitRangesCsv,
                        Double_t signalLow, Double_t signalHigh,
                        Double_t sigmaMin, Double_t sigmaMax, Double_t minFitEntries) {
  using namespace LambdaKfMassFit;
  if (!inputRoot || !outputPrefix || !outputPrefix[0] || !TMath::Finite(signalLow) ||
      !TMath::Finite(signalHigh) || signalHigh <= signalLow ||
      !TMath::Finite(sigmaMin) || !TMath::Finite(sigmaMax) || sigmaMin <= 0. || sigmaMax <= sigmaMin ||
      !TMath::Finite(minFitEntries) || minFitEntries < 0.) return Fail("invalid explicit QA arguments");
  const std::vector<std::string> paths = Split(histogramPathsCsv);
  const std::vector<std::string> rangeTokens = Split(fitRangesCsv);
  if (paths.empty() || rangeTokens.empty()) return Fail("empty histogram or fit-range list");
  std::vector<Range> ranges;
  for (size_t i = 0; i < rangeTokens.size(); ++i) {
    std::stringstream parser(rangeTokens[i]);
    Range range;
    char separator;
    if (!(parser >> range.low >> separator >> range.high) || separator != ':' ||
        !TMath::Finite(range.low) || !TMath::Finite(range.high) ||
        range.low >= signalLow || range.high <= signalHigh) return Fail("invalid fit range or no sidebands");
    parser >> std::ws;
    if (!parser.eof()) return Fail("extra text in fit-range token");
    ranges.push_back(range);
  }
  const TString rootPath = TString(outputPrefix) + ".root";
  const TString csvPath = TString(outputPrefix) + ".csv";
  if (!gSystem->AccessPathName(rootPath) || !gSystem->AccessPathName(csvPath))
    return Fail("refusing to overwrite an existing fit output");
  for (size_t i = 0; i < ranges.size(); ++i)
    if (!gSystem->AccessPathName(Form("%s_range%02u.png", outputPrefix, (unsigned)i)) ||
        !gSystem->AccessPathName(Form("%s_range%02u.pdf", outputPrefix, (unsigned)i)))
      return Fail("refusing to overwrite an existing fit plot");
  TFile* input = TFile::Open(inputRoot, "READ");
  if (!input || input->IsZombie() || input->TestBit(TFile::kRecovered))
    return Fail("cannot read a complete input ROOT file");
  std::vector<TH1*> histograms;
  Double_t binWidth = 0.;
  for (size_t i = 0; i < paths.size(); ++i) {
    TH1* source = dynamic_cast<TH1*>(input->Get(paths[i].c_str()));
    if (!source || source->GetDimension() != 1) return Fail("requested 1D mass histogram missing");
    if (i == 0) binWidth = source->GetXaxis()->GetBinWidth(1);
    if (binWidth <= 0. || Edge(source, signalLow) < 0 || Edge(source, signalHigh) < 0)
      return Fail("signal-window edges must match whole source bins");
    for (Int_t bin = 1; bin <= source->GetNbinsX(); ++bin) {
      const Double_t n = source->GetBinContent(bin), error = source->GetBinError(bin);
      if (!TMath::Finite(n) || !TMath::Finite(error) || n < 0. || error < 0. ||
          TMath::Abs(n - TMath::Nint(n)) > 1.e-7 ||
          TMath::Abs(error * error - n) > 1.e-6 * TMath::Max(1., n) ||
          TMath::Abs(source->GetXaxis()->GetBinWidth(bin) - binWidth) > 1.e-9)
        return Fail("Poisson fitting requires raw unweighted counts and equal uniform bin widths");
    }
    for (size_t j = 0; j < ranges.size(); ++j)
      if (Edge(source, ranges[j].low) < 0 || Edge(source, ranges[j].high) < 0)
        return Fail("fit-range edges must match whole source bins");
    TH1* clone = dynamic_cast<TH1*>(source->Clone(Form("FitInput_%u", (unsigned)i)));
    clone->SetDirectory(0);
    histograms.push_back(clone);
  }
  const TString directory = gSystem->DirName(outputPrefix);
  if (gSystem->AccessPathName(directory) && gSystem->mkdir(directory, kTRUE) != 0)
    return Fail("cannot create output directory");
  TFile* output = TFile::Open(rootPath, "CREATE");
  if (!output || output->IsZombie()) return Fail("cannot create fit ROOT output");
  std::ofstream csv(csvPath.Data());
  if (!csv) return Fail("cannot create fit CSV output");
  csv << std::setprecision(12);
  csv << "histogram,range_index,fit_low,fit_high,window_low,window_high,bin_width,"
      << "fit_entries,Nwin,Nwin_error,accepted,flags,fit_status,covariance_quality,"
      << "mean,sigma,chi2,ndf,S_estimate,B_estimate,S_error,B_error,covariance_SB,"
      << "S_over_B,S_over_B_error,S_over_SplusB,S_over_SplusB_error\n";
  output->cd();
  TNamed sources("FitInputSource", inputRoot); sources.Write();
  TNamed pathsMetadata("FitHistogramPaths", histogramPathsCsv); pathsMetadata.Write();
  TNamed choices("FitChoices", Form("model=gaus+pol2 (centered polynomial); ROOT fit options=RSQLI0N; ranges=%s; signal=[%.17g,%.17g); sigma=[%.17g,%.17g]; minFitEntries=%.17g; binWidth=%.17g; mean restricted to signal window; accepted iff status=0,covarianceQuality=3,positive components,nonnegative quadratic over fit range,and no mean/sigma bound contact", fitRangesCsv, signalLow, signalHigh, sigmaMin, sigmaMax, minFitEntries, binWidth));
  choices.Write();
  TNamed caveat("FitInterpretation", "S/B is Gaussian-window integral divided by quadratic-window integral. Errors are first-order full fit covariance, conditional on a fixed model/range; no scan-selection, between-selection correlation, fit-range, or background-model systematic is included. NaN ratios/errors indicate rejected/unavailable fit, not zero background. S_estimate and B_estimate can exist for rejected fits and must not be used as accepted measurements. Nwin counts candidates including background, not events or true-particle efficiency. No histogram is normalized or mass constrained here.");
  caveat.Write();
  TObject* representatives = input->Get("KfRepresentativeSelections");
  if (representatives) representatives->Write("InputKfRepresentativeSelections");
  for (size_t i = 0; i < histograms.size(); ++i) histograms[i]->Write();
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  const Int_t colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
  Int_t acceptedCount = 0;
  for (size_t rangeIndex = 0; rangeIndex < ranges.size(); ++rangeIndex) {
    const Range& range = ranges[rangeIndex];
    std::vector<Result> results(histograms.size());
    for (size_t i = 0; i < histograms.size(); ++i) {
      Result& r = results[i];
      FitOne(histograms[i], range, signalLow, signalHigh, sigmaMin, sigmaMax,
             minFitEntries, Int_t(rangeIndex * histograms.size() + i), output, r);
      if (r.accepted) ++acceptedCount;
      csv << Quote(paths[i]) << ',' << rangeIndex << ',' << range.low << ',' << range.high << ','
          << signalLow << ',' << signalHigh << ',' << binWidth << ',' << r.nFit << ','
          << r.nWindow << ',' << r.nWindowError << ',' << Int_t(r.accepted) << ',' << Quote(r.flags) << ','
          << r.status << ',' << r.covarianceQuality << ',' << r.mean << ',' << r.sigma << ','
          << r.chi2 << ',' << r.ndf << ',' << r.s << ',' << r.b << ',' << r.sError << ','
          << r.bError << ',' << r.covarianceSB << ',' << r.sb << ',' << r.sbError << ','
          << r.purity << ',' << r.purityError << '\n';
    }
    const Int_t panels = Int_t(histograms.size()) + 1;
    const Int_t rows = (panels + 1) / 2;
    TCanvas* canvas = new TCanvas(Form("KfMassFitCanvas_%u", (unsigned)rangeIndex), "Absolute counts and gaus+pol2 cross-check", 1400, 420 * rows);
    canvas->Divide(2, rows);
    Double_t maximum = 0.;
    for (size_t i = 0; i < histograms.size(); ++i) {
      histograms[i]->GetXaxis()->SetRangeUser(range.low, range.high);
      maximum = TMath::Max(maximum, histograms[i]->GetMaximum());
    }
    canvas->cd(1);
    gPad->SetLeftMargin(0.13);
    TLegend* legend = new TLegend(0.43, 0.65, 0.94, 0.91);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.027);
    for (size_t i = 0; i < histograms.size(); ++i) {
      TH1* drawn = dynamic_cast<TH1*>(histograms[i]->Clone(Form("Overlay_%u_%u", (unsigned)rangeIndex, (unsigned)i)));
      drawn->SetDirectory(0);
      Style(drawn, colors[i % 6], binWidth);
      drawn->SetMinimum(0.); drawn->SetMaximum(TMath::Max(1., maximum) * 1.5);
      drawn->Draw(i == 0 ? "HIST" : "HIST SAME");
      legend->AddEntry(drawn, Form("%s: Nwin=%.0f", paths[i].c_str(), results[i].nWindow), "l");
    }
    legend->Draw();
    TLatex text; text.SetNDC(kTRUE); text.SetTextSize(0.029);
    text.DrawLatex(0.16, 0.92, "Absolute counts (no scaling)");
    text.DrawLatex(0.16, 0.86, Form("Window [%.3f, %.3f) GeV/c^{2}", signalLow, signalHigh));
    for (size_t i = 0; i < histograms.size(); ++i) {
      canvas->cd(Int_t(i) + 2); gPad->SetLeftMargin(0.13);
      TH1* drawn = dynamic_cast<TH1*>(histograms[i]->Clone(Form("Panel_%u_%u", (unsigned)rangeIndex, (unsigned)i)));
      drawn->SetDirectory(0);
      const Int_t color = colors[i % 6];
      Style(drawn, color, binWidth);
      drawn->SetMinimum(0.); drawn->SetMaximum(TMath::Max(1., drawn->GetMaximum()) * 1.5);
      drawn->Draw("E");
      Result& r = results[i];
      if (r.total) {
        r.total->SetLineColor(kBlack); r.total->SetLineStyle(r.accepted ? 1 : 3); r.total->Draw("SAME");
      }
      if (r.background) { r.background->SetLineColor(kOrange + 7); r.background->SetLineStyle(2); r.background->Draw("SAME"); }
      if (r.signal) { r.signal->SetLineColor(kGreen + 2); r.signal->SetLineStyle(3); r.signal->Draw("SAME"); }
      text.SetTextColor(kBlack); text.SetTextSize(0.028);
      text.DrawLatex(0.16, 0.92, paths[i].c_str());
      text.DrawLatex(0.16, 0.86, Form("Range %.3f-%.3f; status=%d, cov=%d", range.low, range.high, r.status, r.covarianceQuality));
      if (r.accepted) text.DrawLatex(0.16, 0.80, Form("S/B=%.3g #pm %.2g; #sigma=%.3g MeV/c^{2}", r.sb, r.sbError, 1000. * r.sigma));
      else { text.SetTextColor(kRed + 1); text.DrawLatex(0.16, 0.80, "REJECTED / UNAVAILABLE fit (see CSV flags)"); }
      text.SetTextColor(kBlack); text.SetTextSize(0.024);
      text.DrawLatex(0.16, 0.74, "Black: total; orange: pol2; green: Gaussian");
      TLine* left = new TLine(signalLow, 0., signalLow, drawn->GetMaximum() * 0.6);
      TLine* right = new TLine(signalHigh, 0., signalHigh, drawn->GetMaximum() * 0.6);
      left->SetLineStyle(2); right->SetLineStyle(2); left->Draw(); right->Draw();
      gPad->RedrawAxis();
    }
    canvas->Modified(); canvas->Update(); output->cd(); canvas->Write();
    canvas->SaveAs(Form("%s_range%02u.png", outputPrefix, (unsigned)rangeIndex));
    canvas->SaveAs(Form("%s_range%02u.pdf", outputPrefix, (unsigned)rangeIndex));
  }
  csv.close();
  if (!csv) return Fail("failed while writing fit CSV");
  TNamed completed("FitRunStatus", "completed; individual fits may be rejected, inspect accepted/flags in CSV");
  output->cd(); completed.Write(); output->Write(); output->Close(); input->Close();
  std::cout << "[fitLambdaKfScanMass] Accepted " << acceptedCount << " / "
            << ranges.size() * histograms.size() << " fits. CSV: " << csvPath << std::endl;
  return 0;
}
