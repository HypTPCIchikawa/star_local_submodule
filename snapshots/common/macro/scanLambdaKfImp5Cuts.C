// Stored-output cut scan, compiled with ROOT5 ACLiC (not interpreted CINT).
// First prepare a manifest with prepareLambdaKfImp5Scan.py, then:
// root4star -b -q 'common/macro/scanLambdaKfImp5Cuts.C+("/tmp/scan.tsv","new-output-directory")'
// No PicoDst event loop or production configuration modification occurs here.
// All thresholds are read from the manifest or source ROOT effective metadata.
#include "../../tests/check_kfparticle_imp5_output.C"
#include "compareLambdaHelixKF.C"
#include <TLine.h>
#include <TDirectory.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#include <limits>

struct LksPoint {
  Double_t pv, daughter, cosine;
  Int_t changed;
  TH1D* mass;
};

struct LksStats {
  Double_t n, left, right, alpha, background, signal, sn, purity;
  Double_t en, eb, es, esn, epurity;
};

static Int_t LksFail(const char* reason) {
  std::cerr << "[scanLambdaKfImp5Cuts] ERROR: " << reason << std::endl;
  return 1;
}

static Bool_t LksNumbers(const std::map<std::string,std::string>& config,
                        const char* key, std::vector<Double_t>& values) {
  std::map<std::string,std::string>::const_iterator found = config.find(key);
  if (found == config.end()) return kFALSE;
  std::istringstream stream(found->second);
  Double_t value = 0.;
  values.clear();
  while (stream >> value) {
    if (!TMath::Finite(value)) return kFALSE;
    values.push_back(value);
  }
  return stream.eof() && !values.empty();
}

static Bool_t LksScalar(const std::map<std::string,std::string>& config,
                       const char* key, Double_t& value) {
  std::vector<Double_t> values;
  if (!LksNumbers(config, key, values) || values.size() != 1) return kFALSE;
  value = values[0];
  return kTRUE;
}

static Int_t LksFind(const std::vector<LksPoint>& points, Double_t pv,
                     Double_t daughter, Double_t cosine) {
  for (UInt_t i = 0; i < points.size(); ++i)
    if (points[i].pv == pv && points[i].daughter == daughter && points[i].cosine == cosine)
      return i;
  return -1;
}

static void LksAdd(std::vector<LksPoint>& points, Double_t pv, Double_t daughter,
                   Double_t cosine, const Double_t* cut, Int_t bins,
                   Double_t low, Double_t high) {
  if (LksFind(points, pv, daughter, cosine) >= 0) return;
  LksPoint point;
  point.pv = pv; point.daughter = daughter; point.cosine = cosine;
  point.changed = (pv != cut[6]) + (daughter != cut[5]) + (cosine != cut[7]);
  point.mass = new TH1D(Form("kf_point_%04u", (UInt_t)points.size()), "", bins, low, high);
  point.mass->SetDirectory(0);
  point.mass->Sumw2();
  points.push_back(point);
}

static Double_t LksCount(const TH1* histogram, const std::vector<Double_t>& bounds) {
  const Double_t low = histogram->GetXaxis()->GetXmin();
  const Double_t width = histogram->GetXaxis()->GetBinWidth(1);
  const Int_t first = TMath::Nint((bounds[0] - low) / width) + 1;
  const Int_t last = TMath::Nint((bounds[1] - low) / width);
  return histogram->Integral(first, last);
}

static LksStats LksMeasure(const TH1* mass, const std::vector<Double_t>& window,
                          const std::vector<Double_t>& left,
                          const std::vector<Double_t>& right) {
  LksStats s;
  s.n = LksCount(mass, window);
  s.left = LksCount(mass, left); s.right = LksCount(mass, right);
  s.alpha = (window[1] - window[0]) / (left[1] - left[0] + right[1] - right[0]);
  s.background = s.alpha * (s.left + s.right);
  s.signal = s.n - s.background;
  s.en = TMath::Sqrt(s.n);
  const Double_t vb = s.alpha * s.alpha * (s.left + s.right);
  s.eb = TMath::Sqrt(vb); s.es = TMath::Sqrt(s.n + vb);
  s.sn = s.purity = s.esn = s.epurity = std::numeric_limits<Double_t>::quiet_NaN();
  if (s.background > 0.) {
    s.sn = s.signal / s.background;
    s.esn = TMath::Sqrt(s.n / (s.background * s.background) +
      s.n * s.n * vb / TMath::Power(s.background, 4));
  }
  if (s.n > 0.) {
    s.purity = s.signal / s.n;
    s.epurity = TMath::Sqrt(vb / (s.n * s.n) +
      s.background * s.background / TMath::Power(s.n, 3));
  }
  return s;
}

static void LksRow(std::ostream& stream, const char* name, Int_t id,
                   const LksPoint* point, const char* windowName,
                   const char* sidebandName, const TH1* mass,
                   const std::vector<Double_t>& window,
                   const std::vector<Double_t>& left,
                   const std::vector<Double_t>& right, Double_t target) {
  const LksStats s = LksMeasure(mass, window, left, right);
  stream << name << ',' << id << ',';
  if (point) stream << point->pv << ',' << point->daughter << ',' << point->cosine << ',' << point->changed;
  else stream << "nan,nan,nan,0";
  stream << ',' << windowName << ',' << sidebandName << ',' << window[0] << ',' << window[1]
         << ',' << left[0] << ',' << left[1] << ',' << right[0] << ',' << right[1]
         << ',' << mass->GetEntries() << ',' << s.n << ',' << s.en
         << ',' << s.left << ',' << s.right << ',' << s.alpha
         << ',' << s.background << ',' << s.eb << ',' << s.signal << ',' << s.es
         << ',' << s.sn << ',' << s.esn << ',' << s.purity << ',' << s.epurity
         << ',' << target << ',' << s.n - target << '\n';
}

static Int_t LksNearest(const std::vector<LksPoint>& points,
                        const std::vector<Double_t>& window,
                        Double_t target, Int_t family, const Double_t* cut) {
  Int_t best = -1;
  Double_t bestDifference = 0.;
  for (UInt_t i = 0; i < points.size(); ++i) {
    if (family == 1 && (points[i].daughter != cut[5] || points[i].cosine != cut[7])) continue;
    if (family == 2 && points[i].cosine != cut[7]) continue;
    const Double_t difference = TMath::Abs(LksCount(points[i].mass, window) - target);
    if (best < 0 || difference < bestDifference ||
        (difference == bestDifference && points[i].changed < points[best].changed)) {
      best = i; bestDifference = difference;
    }
  }
  return best;
}

static TCanvas* LksOverlay(const char* name, TH1D* helix, TH1D* baseline,
                           TH1D* match, const char* label, Double_t width,
                           const std::vector<Double_t>& window) {
  // Scale clones only: all count histograms and cut-scan statistics stay raw.
  // As in compareLambdaHelixKF, normalize by the full displayed-bin sum,
  // excluding under/overflow, with no bin-width factor and no peak scaling.
  TCanvas* canvas = new TCanvas(name, "Absolute counts and unit-area shapes", 1500, 850);
  TH1D* histograms[3] = {helix, baseline, match};
  Double_t integrals[3], maxima[2] = {0., 0.};
  Bool_t emptyShape = kFALSE;
  for (Int_t j = 0; j < 3; ++j) {
    integrals[j] = histograms[j]->Integral(1, histograms[j]->GetNbinsX());
    maxima[0] = TMath::Max(maxima[0], histograms[j]->GetMaximum());
    if (integrals[j] > 0.)
      maxima[1] = TMath::Max(maxima[1], histograms[j]->GetMaximum() / integrals[j]);
    else emptyShape = kTRUE;
  }
  for (Int_t panel = 0; panel < 2; ++panel) {
    canvas->cd();
    TPad* pad = new TPad(Form("%s_panel%d",name,panel),
                        panel == 0 ? "Raw counts" : "Unit-area shapes",
                        panel * 0.5, 0.01, (panel + 1) * 0.5, 0.76);
    pad->Draw(); pad->cd();
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.035);
    for (Int_t j = 0; j < 3; ++j) {
      TH1D* draw = dynamic_cast<TH1D*>(histograms[j]->Clone(Form("%s_panel%d_%d",name,panel,j)));
      draw->SetDirectory(0);
      draw->GetXaxis()->SetRange(1, draw->GetNbinsX());
      if (panel == 1 && integrals[j] > 0.) draw->Scale(1. / integrals[j]);
      LambdaComparisonStyle(draw, j == 0 ? kBlue + 1 : (j == 1 ? kRed + 1 : kGreen + 2),
                            j + 1, Form(panel == 0 ? "Candidates / (%.3g MeV/c^{2})" :
                            "Unit-area fraction / (%.3g MeV/c^{2})", width * 1000.));
      draw->GetXaxis()->SetLabelFont(42); draw->GetXaxis()->SetTitleFont(42);
      draw->GetYaxis()->SetLabelFont(42); draw->GetYaxis()->SetTitleFont(42);
      if (j == 0) draw->SetMaximum(maxima[panel] > 0. ? 1.12 * maxima[panel] : 1.);
      if (panel == 1 && integrals[j] <= 0.) {
        if (j == 0) draw->Draw("AXIS");
      } else draw->Draw(j == 0 ? "HIST" : "HIST SAME");
    }
    if (panel == 1 && emptyShape) {
      TLatex* empty = new TLatex(); empty->SetNDC(); empty->SetTextFont(42);
      empty->SetTextSize(0.027);
      empty->DrawLatex(0.20,0.90,"N=0: corresponding normalized shape unavailable");
    }
    gPad->RedrawAxis();
  }
  canvas->cd();
  helix->SetLineColor(kBlue+1); helix->SetLineStyle(1);
  baseline->SetLineColor(kRed+1); baseline->SetLineStyle(2);
  match->SetLineColor(kGreen+2); match->SetLineStyle(3);
  helix->SetLineWidth(2); baseline->SetLineWidth(2); match->SetLineWidth(2);
  TLegend* legend = new TLegend(0.045, 0.825, 0.97, 0.925);
  legend->SetNColumns(3); legend->SetBorderSize(0); legend->SetFillStyle(0);
  legend->SetTextFont(42); legend->SetTextSize(0.027);
  legend->AddEntry(helix, Form("Helix Imp5: Nwin=%.0f",LksCount(helix,window)), "l");
  legend->AddEntry(baseline, Form("KF baseline: Nwin=%.0f",LksCount(baseline,window)), "l");
  legend->AddEntry(match, Form("%s: Nwin=%.0f",label,LksCount(match,window)), "l");
  legend->Draw();
  TLatex* caption = new TLatex(); caption->SetNDC(); caption->SetTextFont(42);
  caption->SetTextSize(0.032);
  caption->DrawLatex(0.045,0.953,"#Lambda #rightarrow p#pi^{-}: left = absolute counts  |  right = unit area (shape only)");
  caption->SetTextSize(0.026);
  caption->DrawLatex(0.045,0.785,Form("Nwin: %.3f #leq M < %.3f GeV/c^{2}, including background. Normalize each curve over %.3f #leq M < %.3f.",
      window[0],window[1],helix->GetXaxis()->GetXmin(),helix->GetXaxis()->GetXmax()));
  return canvas;
}

static Int_t LksRun(const char* manifestPath, const char* outputDirectory) {
  if (!manifestPath || !outputDirectory || !*outputDirectory)
    return LksFail("manifest and a new output directory are required");
  if (!gSystem->AccessPathName(outputDirectory)) return LksFail("output directory already exists; refusing overwrite");
  std::ifstream manifest(manifestPath);
  if (!manifest) return LksFail("cannot open manifest");
  std::map<std::string,std::string> config;
  std::string line, manifestText;
  while (std::getline(manifest,line)) {
    manifestText += line + "\n";
    const std::string::size_type split = line.find('\t');
    if (split == std::string::npos || split == 0) return LksFail("malformed manifest row");
    const std::string key = line.substr(0,split);
    if (config.find(key) != config.end()) return LksFail("duplicate manifest key");
    config[key] = line.substr(split+1);
  }
  if (config["representativeRule"] != "nearest_window_count_then_fewest_changed_parameters")
    return LksFail("unsupported representative selection rule");
  Double_t low=0., high=0., width=0., events=0., grid=0., schema=0.;
  if (!LksScalar(config,"schemaVersion",schema) || schema != 1. ||
      !LksScalar(config,"massMin",low) || !LksScalar(config,"massMax",high) ||
      !LksScalar(config,"massBinWidth",width) || !LksScalar(config,"inputEvents",events) ||
      !LksScalar(config,"includeThreeDimensionalGrid",grid) || width<=0. || high<=low ||
      events<=0. || events!=TMath::Floor(events) || (grid!=0. && grid!=1.))
    return LksFail("invalid manifest scalar parameters");
  const Int_t bins = TMath::Nint((high-low)/width);
  if (bins<=0 || TMath::Abs(low+bins*width-high)>1.e-9) return LksFail("invalid histogram edges");
  std::vector<Double_t> windows[2], left[3], right[3], pv, daughter, cosine;
  if (!LksNumbers(config,"window",windows[0]) || !LksNumbers(config,"narrowWindow",windows[1]) ||
      !LksNumbers(config,"leftSideband",left[0]) || !LksNumbers(config,"rightSideband",right[0]) ||
      !LksNumbers(config,"nearLeftSideband",left[1]) || !LksNumbers(config,"nearRightSideband",right[1]) ||
      !LksNumbers(config,"farLeftSideband",left[2]) || !LksNumbers(config,"farRightSideband",right[2]) ||
      !LksNumbers(config,"maxDistanceToPv",pv) || !LksNumbers(config,"maxDaughterDistance",daughter) ||
      !LksNumbers(config,"minCosPointing",cosine)) return LksFail("missing scan arrays/windows");
  for (Int_t w=0; w<2; ++w) {
    if (windows[w].size()!=2 || windows[w][0]<low || windows[w][1]>high || windows[w][1]<=windows[w][0])
      return LksFail("invalid count window");
    for (Int_t edge=0; edge<2; ++edge)
      if (TMath::Abs((windows[w][edge]-low)/width-TMath::Nint((windows[w][edge]-low)/width))>1.e-7)
        return LksFail("count window not aligned with bins");
    for (Int_t s=0; s<3; ++s) {
      if (left[s].size()!=2 || right[s].size()!=2 || left[s][0]<low || right[s][1]>high ||
          left[s][1]<=left[s][0] || right[s][1]<=right[s][0] ||
          left[s][1]>windows[w][0] || right[s][0]<windows[w][1]) return LksFail("invalid/overlapping sidebands");
      for (Int_t edge=0; edge<2; ++edge)
        if (TMath::Abs((left[s][edge]-low)/width-TMath::Nint((left[s][edge]-low)/width))>1.e-7 ||
            TMath::Abs((right[s][edge]-low)/width-TMath::Nint((right[s][edge]-low)/width))>1.e-7)
          return LksFail("sideband not aligned with bins");
    }
  }
  const TString helixPath(config["helixRoot"].c_str()), kfPath(config["kfRoot"].c_str());
  if (helixPath.IsNull() || kfPath.IsNull()) return LksFail("missing input paths");
  if (KfImp5OutputQaInspect(kfPath,kTRUE)!=0) return LksFail("input failed existing structural/Imp5 QA");
  TFile* helixFile=TFile::Open(helixPath,"READ");
  TFile* kfFile=TFile::Open(kfPath,"READ");
  if (!helixFile || helixFile->IsZombie() || !kfFile || kfFile->IsZombie() ||
      helixFile->TestBit(TFile::kRecovered) || kfFile->TestBit(TFile::kRecovered))
    return LksFail("missing, incomplete or recovered input file");
  TH1* helixSource=dynamic_cast<TH1*>(helixFile->Get("hLambda_InvMass"));
  TH1* kfSource=dynamic_cast<TH1*>(kfFile->Get("hKfLambdaMassSelected"));
  TH1* oldRead=dynamic_cast<TH1*>(helixFile->Get("hRefMultVsNTOFMatch"));
  TH1* newRead=dynamic_cast<TH1*>(kfFile->Get("hRefMultVsNTOFMatch"));
  TH1* oldAccepted=dynamic_cast<TH1*>(helixFile->Get("hN"));
  TH1* newAccepted=dynamic_cast<TH1*>(kfFile->Get("hN"));
  TTree* tree=dynamic_cast<TTree*>(kfFile->Get("KfLambdaCandidates"));
  TObjString* metadata=dynamic_cast<TObjString*>(kfFile->Get("KFParticleEffectiveConfiguration"));
  if (!helixSource || !kfSource || !oldRead || !newRead || !oldAccepted || !newAccepted || !tree || !metadata)
    return LksFail("missing required input objects");
  if (oldRead->GetEntries()!=events || newRead->GetEntries()!=events ||
      oldAccepted->GetEntries()!=newAccepted->GetEntries()) return LksFail("read/accepted event counters differ");
  TH1D* helix=LambdaComparisonRebin(helixSource,"HelixLambdaCounts",low,high,width);
  TH1D* saved=LambdaComparisonRebin(kfSource,"KfSourceSelectedCounts",low,high,width);
  if (!helix || !saved) return LksFail("common histogram cannot use whole source bins");
  const char* keys[16]={"nSigmaProton","nSigmaPion","imp5MinDCAProton","imp5MinDCAPion",
    "imp5MaxPathLength","maxDaughterDistance","maxDistanceToPv","minCosPointing",
    "minMass","maxMass","maxMassError","maxChi2Ndf","maxTopoChi2Ndf",
    "minDecayLength","minDecayLengthSignificance","minVertexLineSignificance"};
  Double_t cut[16];
  for (Int_t i=0;i<16;++i)
    if (!KfOutputMetadataDouble(metadata->GetString(),keys[i],cut[i])) return LksFail("missing source final-cut metadata");
  Bool_t anti=kFALSE;
  if (!KfImp5MetadataBool(metadata->GetString(),"reconstructAntiLambda",anti) || anti)
    return LksFail("this scan requires Lambda-only Imp5 input");
  const char* names[19]={"protonPidPull","pionPidPull","protonDcaToPv","pionDcaToPv",
    "protonHelixPathLength","pionHelixPathLength","daughterDistance","distanceToPv","cosPointing",
    "mass","massError","chi2Ndf","topoChi2Ndf","decayLength","decayLengthSignificance",
    "vertexLineLengthSignificance","pdg","selected","imp5PathValid"};
  TLeaf* leaves[19];
  for (Int_t i=0;i<19;++i) {
    leaves[i]=tree->GetLeaf(names[i]);
    if (!leaves[i]) return LksFail("missing candidate leaf");
  }
  std::vector<LksPoint> points;
  LksAdd(points,cut[6],cut[5],cut[7],cut,bins,low,high);
  for (UInt_t p=0;p<pv.size();++p) LksAdd(points,pv[p],cut[5],cut[7],cut,bins,low,high);
  for (UInt_t d=0;d<daughter.size();++d) {
    if (daughter[d]<0.) return LksFail("negative daughter-distance scan threshold");
    LksAdd(points,cut[6],daughter[d],cut[7],cut,bins,low,high);
  }
  for (UInt_t c=0;c<cosine.size();++c) {
    if (cosine[c]<-1. || cosine[c]>1.) return LksFail("invalid cosine scan threshold");
    LksAdd(points,cut[6],cut[5],cosine[c],cut,bins,low,high);
  }
  for (UInt_t p=0;p<pv.size();++p)
    for (UInt_t d=0;d<daughter.size();++d) {
      LksAdd(points,pv[p],daughter[d],cut[7],cut,bins,low,high);
      if (grid!=0.) for (UInt_t c=0;c<cosine.size();++c)
        LksAdd(points,pv[p],daughter[d],cosine[c],cut,bins,low,high);
    }
  TH1* audit=dynamic_cast<TH1*>(kfSource->Clone("KfBaselineNativeBinningAudit"));
  audit->SetDirectory(0); audit->Reset();
  Long64_t storedSelected=0, fixedSurvived=0, mismatches=0;
  for (Long64_t entry=0;entry<tree->GetEntries();++entry) {
    if (tree->GetEntry(entry)<=0) return LksFail("unreadable candidate entry");
    Double_t v[19];
    for (Int_t i=0;i<19;++i) v[i]=leaves[i]->GetValue();
    if (v[16]!=3122.) return LksFail("unexpected non-Lambda candidate");
    for (Int_t i=0;i<16;++i)
      if (i!=4 && i!=5 && !TMath::Finite(v[i])) return LksFail("nonfinite candidate observable");
    const Bool_t stored=v[17]==1.;
    if (stored) ++storedSelected;
    // Original PID and DCA were already enforced before Finder. Existing QA
    // validates their Float_t serialization tolerance; do not retune those here.
    Bool_t fixed=v[18]==1. && TMath::Finite(v[4]) && TMath::Finite(v[5]) &&
      TMath::Abs(v[4])<=cut[4] && TMath::Abs(v[5])<=cut[4] && v[9]>=cut[8] && v[9]<=cut[9];
    if (cut[10]>=0. && v[10]>cut[10]) fixed=kFALSE;
    if (cut[11]>=0. && v[11]>cut[11]) fixed=kFALSE;
    if (cut[12]>=0. && v[12]>cut[12]) fixed=kFALSE;
    if (cut[13]>=0. && v[13]<cut[13]) fixed=kFALSE;
    if (cut[14]>=0. && v[14]<cut[14]) fixed=kFALSE;
    if (cut[15]>=0. && v[15]<cut[15]) fixed=kFALSE;
    const Bool_t baseline=fixed && (cut[5]<0. || v[6]<=cut[5]) &&
      (cut[6]<0. || v[7]<=cut[6]) && v[8]>=cut[7];
    if (baseline!=stored) ++mismatches;
    if (baseline) audit->Fill(v[9]);
    if (!fixed) continue;
    ++fixedSurvived;
    for (UInt_t j=0;j<points.size();++j) {
      if (v[6]>points[j].daughter || (points[j].pv>=0. && v[7]>points[j].pv) ||
          v[8]<points[j].cosine) continue;
      points[j].mass->Fill(v[9]);
    }
  }
  if (mismatches || storedSelected!=audit->GetEntries() || audit->GetEntries()!=kfSource->GetEntries())
    return LksFail("baseline selected bits/count do not exactly reproduce source");
  for (Int_t bin=0;bin<=kfSource->GetNbinsX()+1;++bin)
    if (audit->GetBinContent(bin)!=kfSource->GetBinContent(bin)) return LksFail("native baseline mass bins differ");
  for (Int_t bin=1;bin<=bins;++bin)
    if (points[0].mass->GetBinContent(bin)!=saved->GetBinContent(bin)) return LksFail("rebinned baseline mass differs");
  TH1D* baseline=dynamic_cast<TH1D*>(points[0].mass->Clone("KfBaselineCounts")); baseline->SetDirectory(0);
  const Double_t target=LksCount(helix,windows[0]);
  const Int_t bestPv=LksNearest(points,windows[0],target,1,cut);
  const Int_t bestPd=LksNearest(points,windows[0],target,2,cut);
  const Int_t best=LksNearest(points,windows[0],target,0,cut);
  if (bestPv<0 || bestPd<0 || best<0) return LksFail("no representative scan point");
  const Int_t indices[3]={bestPv,bestPd,best};
  const char* aliases[3]={"KfNearestPvCounts","KfNearestPvDaughterCounts","KfNearestCounts"};
  TH1D* representatives[3];
  std::ostringstream chosen;
  chosen << std::setprecision(17) << "Rule: absolute primary-window count difference, then fewest changed parameters; remaining ties retain manifest grid order. S/B is not used.\n";
  for (Int_t r=0;r<3;++r) {
    representatives[r]=dynamic_cast<TH1D*>(points[indices[r]].mass->Clone(aliases[r]));
    representatives[r]->SetDirectory(0);
    chosen << aliases[r] << ": point=" << indices[r] << ", maxDistanceToPv=" << points[indices[r]].pv
      << ", maxDaughterDistance=" << points[indices[r]].daughter << ", minCosPointing=" << points[indices[r]].cosine
      << ", changed=" << points[indices[r]].changed << ", Nwin=" << LksCount(representatives[r],windows[0]) << '\n';
  }
  if (gSystem->mkdir(outputDirectory,kTRUE)!=0) return LksFail("cannot exclusively create output directory");
  const TString prefix=TString(outputDirectory)+"/";
  std::ofstream csv((prefix+"scan_metrics.csv").Data());
  std::ofstream summary((prefix+"scan_summary.txt").Data());
  if (!csv || !summary) return LksFail("cannot create CSV/summary");
  csv << std::setprecision(17);
  csv << "sample,point,maxDistanceToPv,maxDaughterDistance,minCosPointing,changedParameters,window,sidebandVariation,windowLow,windowHigh,leftLow,leftHigh,rightLow,rightHigh,totalSelected,Nwin,NwinError,Nleft,Nright,alpha,B,BError,S,SError,sidebandSoverB,sidebandSoverBError,purity,purityError,HelixNwin,NwinMinusHelix\n";
  const char* windowNames[2]={"primary","narrow"};
  const char* sidebandNames[3]={"default","near","far"};
  for (Int_t w=0;w<2;++w) for (Int_t s=0;s<3;++s) {
    const Double_t reference=LksCount(helix,windows[w]);
    LksRow(csv,"helix",-1,0,windowNames[w],sidebandNames[s],helix,windows[w],left[s],right[s],reference);
    for (UInt_t j=0;j<points.size();++j)
      LksRow(csv,j==0?"kf_baseline":"kf_scan",j,&points[j],windowNames[w],sidebandNames[s],
             points[j].mass,windows[w],left[s],right[s],reference);
  }
  const TString scope="Stored raw KF candidates only: input PID/gDCA, covariance, Finder and Topo cuts cannot be recovered or loosened by this QA. Path/mass/other final cuts remain from source metadata. Nwin counts candidates INCLUDING background, not unique events, true Lambda yield or efficiency. Sideband B uses constant local-density interpolation; sideband S/B is NOT the earlier gaus+pol2 fit S/N. Negative background-subtracted S is not clipped. Errors use independent unweighted Poisson counts in disjoint windows, with no sideband-shape systematic or candidate/event correlations. Both methods share events, and scan points share candidates: errors on method differences/ratios are NOT evaluated. Same-data cut tuning is exploratory, not an optimal-cut or significance claim. Counter equality alone does not prove event identity; use the prior ordered-input run manifest.";
  summary << std::setprecision(17) << "Baseline replay: exact selected decisions and native/rebinned mass bins PASS\n"
    << "raw=" << tree->GetEntries() << ", fixedCutsSurvived=" << fixedSurvived << ", baselineSelected=" << storedSelected
    << ", points=" << points.size() << ", inputEvents=" << events << ", acceptedEach=" << oldAccepted->GetEntries()
    << ", HelixNwin=" << target << "\n" << chosen.str() << "\n" << scope.Data() << "\n";
  summary << "Manifest (input/config/code SHA256 recorded by preparation, not recomputed inside ROOT):\n" << manifestText;
  csv.close(); summary.close();
  if (!csv || !summary) return LksFail("CSV/summary write failed");
  gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetTextFont(42);
  TCanvas* overlay=LksOverlay("LambdaCountMatchedOverlay",helix,baseline,representatives[2],"KF nearest count",width,windows[0]);
  TCanvas* pvOverlay=LksOverlay("LambdaPvOnlyOverlay",helix,baseline,representatives[0],"KF PV-only nearest",width,windows[0]);
  TCanvas* pairedOverlay=LksOverlay("LambdaPairedCutOverlay",helix,baseline,representatives[1],"KF PV+daughter",width,windows[0]);
  TH1D* pvN=new TH1D("PvOnlyWindowCounts",";maxDistanceToPv [cm];Candidates in mass window",pv.size(),0,pv.size());
  TH1D* pvSn=new TH1D("PvOnlySidebandSoverB",";maxDistanceToPv [cm];Sideband-estimated S/B",pv.size(),0,pv.size());
  pvN->SetDirectory(0); pvSn->SetDirectory(0);
  const LksStats helixStats=LksMeasure(helix,windows[0],left[0],right[0]);
  for (UInt_t p=0;p<pv.size();++p) {
    const Int_t id=LksFind(points,pv[p],cut[5],cut[7]);
    const LksStats s=LksMeasure(points[id].mass,windows[0],left[0],right[0]);
    const TString label=pv[p]<0.?TString("disabled"):TString(Form("%.3g",pv[p]));
    pvN->GetXaxis()->SetBinLabel(p+1,label); pvSn->GetXaxis()->SetBinLabel(p+1,label);
    pvN->SetBinContent(p+1,s.n); pvN->SetBinError(p+1,s.en);
    if (TMath::Finite(s.sn)) { pvSn->SetBinContent(p+1,s.sn); pvSn->SetBinError(p+1,s.esn); }
  }
  TCanvas* curve=new TCanvas("LambdaPvOnlyScan","PV-distance-only scan",1500,650); curve->Divide(2,1);
  for (Int_t panel=1;panel<=2;++panel) {
    curve->cd(panel); gPad->SetLeftMargin(0.20); gPad->SetRightMargin(0.045);
    gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.28);
    TH1D* h=panel==1?pvN:pvSn; h->SetStats(0); h->SetLineColor(kRed+1); h->SetMarkerStyle(20);
    h->GetYaxis()->SetTitleOffset(1.85); h->GetYaxis()->SetTitleSize(0.042);
    h->GetYaxis()->SetLabelSize(0.037); h->GetYaxis()->SetTitleFont(42); h->GetYaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleOffset(1.25); h->GetXaxis()->SetTitleSize(0.042);
    h->GetXaxis()->SetLabelSize(0.037); h->GetXaxis()->SetTitleFont(42); h->GetXaxis()->SetLabelFont(42);
    const Double_t reference=panel==1?helixStats.n:helixStats.sn;
    if (TMath::Finite(reference)) h->SetMaximum(1.25*TMath::Max(h->GetMaximum(),reference));
    h->Draw("E1");
    if (TMath::Finite(reference)) { TLine* line=new TLine(0,reference,pv.size(),reference); line->SetLineColor(kBlue+1); line->SetLineStyle(2); line->Draw(); }
    TLatex* note=new TLatex(); note->SetNDC(); note->SetTextFont(42); note->SetTextSize(0.029);
    note->DrawLatex(0.22,0.92,"Red: KF scan; blue dashed: Helix Imp5");
    note->DrawLatex(0.22,0.86,"Other cuts fixed to KF Imp5 baseline");
    note->SetTextSize(0.027);
    note->DrawLatex(0.22,0.80,"Categorical x bins (unequal numeric spacing)");
  }
  TFile* output=TFile::Open(prefix+"scan.root","CREATE");
  if (!output || output->IsZombie()) return LksFail("cannot create scan ROOT file");
  output->cd();
  Int_t bytes=helix->Write()+baseline->Write()+saved->Write()+audit->Write();
  for (Int_t r=0;r<3;++r) bytes+=representatives[r]->Write();
  TObjString manifestObject(manifestText.c_str()); bytes+=manifestObject.Write("KfScanManifest");
  TObjString choices(chosen.str().c_str()); bytes+=choices.Write("KfRepresentativeSelections");
  TObjString scopeObject(scope); bytes+=scopeObject.Write("KfScanScope");
  TObjString plotScope("Mass overlays: left raw counts; right clones scaled by 1/Integral(1,Nbins) over the full massMin/massMax range. Under/overflow excluded; no width factor or peak normalization. Count histograms, candidate selection and S/B remain unscaled.");
  bytes+=plotScope.Write("KfOverlayDisplayConvention");
  bytes+=metadata->Write("KfSourceEffectiveConfiguration");
  bytes+=pvN->Write()+pvSn->Write()+overlay->Write()+pvOverlay->Write()+pairedOverlay->Write()+curve->Write();
  TDirectory* directory=output->mkdir("points");
  if (!directory) return LksFail("cannot create points ROOT directory");
  directory->cd();
  for (UInt_t j=0;j<points.size();++j) bytes+=points[j].mass->Write();
  output->Flush(); const Bool_t badWrite=output->TestBit(TFile::kWriteError); output->Close();
  if (bytes<=0 || badWrite || output->TestBit(TFile::kWriteError)) return LksFail("ROOT write error");
  overlay->Print(prefix+"count_matched_overlay.png"); overlay->Print(prefix+"count_matched_overlay.pdf");
  pvOverlay->Print(prefix+"pv_only_overlay.png"); pvOverlay->Print(prefix+"pv_only_overlay.pdf");
  pairedOverlay->Print(prefix+"paired_cut_overlay.png"); pairedOverlay->Print(prefix+"paired_cut_overlay.pdf");
  curve->Print(prefix+"pv_only_curve.png"); curve->Print(prefix+"pv_only_curve.pdf");
  const char* products[8]={"count_matched_overlay.png","count_matched_overlay.pdf","pv_only_overlay.png","pv_only_overlay.pdf","paired_cut_overlay.png","paired_cut_overlay.pdf","pv_only_curve.png","pv_only_curve.pdf"};
  for (Int_t f=0;f<8;++f) { FileStat_t info; if (gSystem->GetPathInfo(prefix+products[f],info)!=0 || info.fSize<=0) return LksFail("missing/empty plot output"); }
  std::cout << "[scanLambdaKfImp5Cuts] PASS baseline exact replay; " << points.size() << " unique scan points\n"
            << chosen.str() << "Output: " << outputDirectory << std::endl;
  helixFile->Close(); kfFile->Close();
  return 0;
}

void scanLambdaKfImp5Cuts(const char* manifestPath, const char* outputDirectory) {
  gSystem->Exit(LksRun(manifestPath,outputDirectory));
}
