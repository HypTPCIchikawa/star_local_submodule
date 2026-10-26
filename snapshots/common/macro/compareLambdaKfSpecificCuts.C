// ROOT5 ACLiC QA of actual reruns with KF-specific thresholds relaxed and an
// optional Lambda Finder geometry bypass; ordinary Maker cuts stay fixed.
// Uses the validated TSV from prepareLambdaKfImp5Scan.py for plotting/counting.
// runsManifest: id<TAB>ROOT path<TAB>English legend label; # comments allowed.
// No event loop, refitting, source modification, or source-file overwriting.
#include "scanLambdaKfImp5Cuts.C"
#include <set>
#include <cctype>

namespace {
const char* LkspKeys[5] = {"interfaceChiPrimaryCut", "finderChiPrimary2D",
  "finderChi2Ndf2D", "finderLdL2D", "topoChi2NdfCut"};

struct LkspRun {
  std::string id, path, label, configuration, eventConfiguration;
  Double_t cuts[5], raw, selected, accepted;
  Bool_t applyGeometry;
  Int_t changed;
  TH1D* mass;
};

static Int_t LkspFail(const std::string& message) {
  std::cerr << "[compareLambdaKfSpecificCuts] ERROR: " << message << std::endl;
  return 1;
}

static std::string LkspTrim(const std::string& value) {
  const std::string::size_type begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) return "";
  return value.substr(begin, value.find_last_not_of(" \t\r\n") - begin + 1);
}

static Bool_t LkspMap(const std::string& content, char separator,
                     std::map<std::string,std::string>& result) {
  result.clear();
  std::istringstream input(content);
  std::string line;
  while (std::getline(input,line)) {
    line = LkspTrim(line);
    if (line.empty() || line[0] == '#') continue;
    const std::string::size_type split = line.find(separator);
    if (split == std::string::npos) return kFALSE;
    const std::string key = LkspTrim(line.substr(0,split));
    const std::string value = LkspTrim(line.substr(split+1));
    if (key.empty() || value.empty() || result.find(key) != result.end()) return kFALSE;
    result[key] = value;
  }
  return !result.empty() && input.eof();
}

static Bool_t LkspRead(const char* path, std::string& content) {
  if (!path || !*path) return kFALSE;
  std::ifstream stream(path);
  if (!stream) return kFALSE;
  std::ostringstream copy; copy << stream.rdbuf();
  if (stream.bad()) return kFALSE;
  content = copy.str();
  return !content.empty();
}

static std::string LkspCsv(const std::string& value) {
  std::string quoted("\"");
  for (std::string::size_type i=0; i<value.size(); ++i) {
    if (value[i] == '"') quoted += '"';
    quoted += value[i];
  }
  return quoted + '"';
}

static Bool_t LkspSafeId(const std::string& value) {
  if (value.empty() || value == "." || value == "..") return kFALSE;
  for (std::string::size_type i=0; i<value.size(); ++i)
    if (!std::isalnum(static_cast<unsigned char>(value[i])) &&
        value[i]!='_' && value[i]!='-' && value[i]!='.') return kFALSE;
  return kTRUE;
}

static Bool_t LkspSameValue(const std::string& a, const std::string& b) {
  if (a == b) return kTRUE;
  const TString left(a.c_str()), right(b.c_str());
  return left.IsFloat() && right.IsFloat() && TMath::Finite(left.Atof()) &&
         TMath::Finite(right.Atof()) && left.Atof() == right.Atof();
}

static Bool_t LkspFixedEqual(const std::string& baseline,
                            const std::string& candidate, Bool_t events,
                            std::string& reason) {
  std::map<std::string,std::string> a,b;
  if (!LkspMap(baseline,':',a) || !LkspMap(candidate,':',b)) {
    reason = "malformed or duplicate effective metadata"; return kFALSE;
  }
  if (events) {
    // Exact legacy-to-sparse diagnostic-source spelling normalization only.
    // Mode, policy, flags, numerical center, and active maxNTr remain compared.
    std::map<std::string,std::string>* metadata[2]={&a,&b};
    for (Int_t m=0; m<2; ++m) {
      std::map<std::string,std::string>::iterator source=metadata[m]->find("KFEventSelection.vertexSource");
      const std::map<std::string,std::string>::const_iterator mode=metadata[m]->find("KFEventSelection.mode");
      if (source!=metadata[m]->end() && mode!=metadata[m]->end() &&
          source->second=="event.vertexByMode."+mode->second)
        source->second="event.qaVertexByMode."+mode->second;
    }
    const char* paths[3] = {"KFEventSelection.mainconf", "KFEventSelection.analysisInfo",
                          "KFEventSelection.eventConfig"};
    for (Int_t i=0; i<3; ++i) { a.erase(paths[i]); b.erase(paths[i]); }
    // Imp5 explicitly disables these event cuts; old metadata nevertheless
    // printed their inactive values before the 2026-09-08 cleanup.
    const char* inactive[7] = {"KFEventSelection.minVz", "KFEventSelection.maxVz",
      "KFEventSelection.maxVr", "KFEventSelection.minRefMult", "KFEventSelection.maxRefMult",
      "KFEventSelection.maxVzDiff", "KFEventSelection.maxAbsVzVpd"};
    for (Int_t i=0; i<7; ++i) { a.erase(inactive[i]); b.erase(inactive[i]); }
  } else {
    for (Int_t i=0; i<5; ++i) { a.erase(LkspKeys[i]); b.erase(LkspKeys[i]); }
    a.erase("applyLambdaGeometryCuts"); b.erase("applyLambdaGeometryCuts");
    // KfImp5OutputQaInspect has already proved lambda_imp5/TPC-only for both.
    // These keys do not select that profile and were omitted after cleanup.
    const char* inactive[16] = {"minNHitsDedx", "minPt", "maxPt", "minEta", "maxEta",
      "minDedxError", "maxDedxError", "nSigmaKaon", "tofCalibrationProfile", "tofNSigma",
      "tofPMax", "tofKaonPMin", "tofKaonPMax", "primaryProbCut", "tofMean", "tofSigma"};
    for (Int_t i=0; i<16; ++i) { a.erase(inactive[i]); b.erase(inactive[i]); }
  }
  if (a.size() != b.size()) {
    reason = "fixed metadata key sets differ"; return kFALSE;
  }
  for (std::map<std::string,std::string>::const_iterator i=a.begin(); i!=a.end(); ++i) {
    const std::map<std::string,std::string>::const_iterator j=b.find(i->first);
    if (j==b.end() || !LkspSameValue(i->second,j->second)) {
      reason = "fixed metadata differ at " + i->first; return kFALSE;
    }
  }
  return kTRUE;
}

static Bool_t LkspLoadRun(LkspRun& run, Double_t events, Double_t low,
                          Double_t high, Double_t width, std::string& reason) {
  // A completed zero-candidate trial is a valid scan point, not a successful match.
  if (KfImp5OutputQaInspect(run.path.c_str(),kFALSE) != 0) {
    reason = "run failed structural/Imp5 QA: " + run.id; return kFALSE;
  }
  TFile* source = TFile::Open(run.path.c_str(),"READ");
  if (!source || source->IsZombie() || source->TestBit(TFile::kRecovered)) {
    reason = "cannot open complete run: " + run.id; return kFALSE;
  }
  TObjString* configuration = dynamic_cast<TObjString*>(source->Get("KFParticleEffectiveConfiguration"));
  TObjString* eventConfiguration = dynamic_cast<TObjString*>(source->Get("KFEventSelectionConfiguration"));
  TH1* read = dynamic_cast<TH1*>(source->Get("hRefMultVsNTOFMatch"));
  TH1* accepted = dynamic_cast<TH1*>(source->Get("hN"));
  TH1* selected = dynamic_cast<TH1*>(source->Get("hKfLambdaMassSelected"));
  TTree* raw = dynamic_cast<TTree*>(source->Get("KfLambdaCandidates"));
  if (!configuration || !eventConfiguration || !read || !accepted || !selected || !raw ||
      read->GetEntries()!=events) {
    reason = "run objects/read event counter mismatch: " + run.id;
    source->Close(); delete source; return kFALSE;
  }
  Bool_t anti=kTRUE;
  if (!KfImp5MetadataBool(configuration->GetString(),"reconstructAntiLambda",anti) || anti) {
    reason = "Lambda-only reconstruction required: " + run.id;
    source->Close(); delete source; return kFALSE;
  }
  run.configuration = configuration->GetString().Data();
  run.eventConfiguration = eventConfiguration->GetString().Data();
  run.applyGeometry=kTRUE;
  if (!KfOutputMetadataValue(configuration->GetString(),"applyLambdaGeometryCuts").IsNull() &&
      !KfImp5MetadataBool(configuration->GetString(),"applyLambdaGeometryCuts",run.applyGeometry)) {
    reason = "invalid applyLambdaGeometryCuts flag: " + run.id;
    source->Close(); delete source; return kFALSE;
  }
  for (Int_t i=0; i<5; ++i) {
    // Historical upstream SelectParticleCandidates used chi2/NDF < 3 without a
    // configurable field. A present field must always be finite and parseable.
    if (i==4 && KfOutputMetadataValue(configuration->GetString(),LkspKeys[i]).IsNull())
      run.cuts[i] = 3.;
    else if (!KfOutputMetadataDouble(configuration->GetString(),LkspKeys[i],run.cuts[i])) {
      reason = "missing/nonfinite KF threshold " + std::string(LkspKeys[i]) + ": " + run.id;
      source->Close(); delete source; return kFALSE;
    }
    if (run.cuts[i]<0. || ((i==2 || i==4) && run.cuts[i]<=0.)) {
      reason = "invalid KF threshold: " + run.id;
      source->Close(); delete source; return kFALSE;
    }
  }
  run.raw = raw->GetEntries(); run.selected = selected->GetEntries();
  run.accepted = accepted->GetEntries();
  run.mass = LambdaComparisonRebin(selected,Form("KfSpecific_%s_Counts",run.id.c_str()),low,high,width);
  source->Close(); delete source;
  if (!run.mass) { reason = "mass bins do not match QA domain: " + run.id; return kFALSE; }
  return kTRUE;
}

static void LkspMetric(std::ostream& out, const std::string& id,
                        const LkspRun* run, const char* windowName,
                        const char* sidebandName, TH1D* mass,
                        const std::vector<Double_t>& window,
                        const std::vector<Double_t>& left,
                        const std::vector<Double_t>& right, Double_t target) {
  const LksStats s = LksMeasure(mass,window,left,right);
  out << LkspCsv(id) << ',' << LkspCsv(run ? run->label : "Helix Imp5") << ','
      << LkspCsv(run ? run->path : "") << ',';
  for (Int_t i=0; i<5; ++i) {
    if (run) out << run->cuts[i]; else out << "nan";
    out << ',';
  }
  out << (run ? (run->applyGeometry ? "true" : "false") : "N/A") << ','
      << (run ? run->changed : 0) << ',' << windowName << ',' << sidebandName
      << ',' << window[0] << ',' << window[1] << ',' << left[0] << ',' << left[1]
      << ',' << right[0] << ',' << right[1] << ',' << (run ? run->raw : -1.)
      << ',' << (run ? run->selected : mass->GetEntries())
      << ',' << mass->Integral(1,mass->GetNbinsX()) << ',' << s.n << ',' << s.en
      << ',' << s.left << ',' << s.right << ',' << s.alpha << ',' << s.background
      << ',' << s.eb << ',' << s.signal << ',' << s.es << ',' << s.sn << ',' << s.esn
      << ',' << s.purity << ',' << s.epurity << ',' << target << ',' << s.n-target << '\n';
}
} // anonymous namespace

Int_t compareLambdaKfSpecificCuts(const char* qaManifest, const char* runsManifest,
                                const char* outputDirectory) {
  if (!outputDirectory || !*outputDirectory || !gSystem->AccessPathName(outputDirectory))
    return LkspFail("a new output directory is required; refusing overwrite");
  std::string qaText, runText;
  std::map<std::string,std::string> config;
  if (!LkspRead(qaManifest,qaText) || !LkspMap(qaText,'\t',config) ||
      !LkspRead(runsManifest,runText)) return LkspFail("cannot read valid manifests");
  Double_t schema=0.,events=0.,low=0.,high=0.,width=0.;
  if (!LksScalar(config,"schemaVersion",schema) || schema!=1. ||
      !LksScalar(config,"inputEvents",events) || events<=0. || events!=TMath::Floor(events) ||
      !LksScalar(config,"massMin",low) || !LksScalar(config,"massMax",high) ||
      !LksScalar(config,"massBinWidth",width) || width<=0. || high<=low ||
      config["representativeRule"]!="nearest_window_count_then_fewest_changed_parameters")
    return LkspFail("invalid QA scalars or representative rule");
  const Int_t bins=TMath::Nint((high-low)/width);
  if (bins<=0 || TMath::Abs(low+bins*width-high)>1.e-9) return LkspFail("unaligned mass edges");
  std::vector<Double_t> windows[2],left[3],right[3];
  if (!LksNumbers(config,"window",windows[0]) || !LksNumbers(config,"narrowWindow",windows[1]) ||
      !LksNumbers(config,"leftSideband",left[0]) || !LksNumbers(config,"rightSideband",right[0]) ||
      !LksNumbers(config,"nearLeftSideband",left[1]) || !LksNumbers(config,"nearRightSideband",right[1]) ||
      !LksNumbers(config,"farLeftSideband",left[2]) || !LksNumbers(config,"farRightSideband",right[2]))
    return LkspFail("missing QA windows/sidebands");
  for (Int_t w=0; w<2; ++w) {
    if (windows[w].size()!=2 || windows[w][0]<low || windows[w][1]>high || windows[w][0]>=windows[w][1])
      return LkspFail("invalid count window");
    for (Int_t s=0; s<3; ++s) {
      if (left[s].size()!=2 || right[s].size()!=2 || left[s][0]<low || right[s][1]>high ||
          left[s][0]>=left[s][1] || right[s][0]>=right[s][1] ||
          left[s][1]>windows[w][0] || right[s][0]<windows[w][1])
        return LkspFail("invalid or overlapping sidebands");
      for (Int_t e=0; e<2; ++e)
        if (TMath::Abs((windows[w][e]-low)/width-TMath::Nint((windows[w][e]-low)/width))>1.e-7 ||
            TMath::Abs((left[s][e]-low)/width-TMath::Nint((left[s][e]-low)/width))>1.e-7 ||
            TMath::Abs((right[s][e]-low)/width-TMath::Nint((right[s][e]-low)/width))>1.e-7)
          return LkspFail("window/sideband edges must align with mass bins");
    }
  }
  LkspRun baseline; baseline.id="baseline"; baseline.path=config["kfRoot"];
  baseline.label="KF baseline"; baseline.changed=0; baseline.mass=0;
  std::string reason;
  if (!LkspLoadRun(baseline,events,low,high,width,reason)) return LkspFail(reason);
  TFile* helixFile=TFile::Open(config["helixRoot"].c_str(),"READ");
  if (!helixFile || helixFile->IsZombie() || helixFile->TestBit(TFile::kRecovered))
    return LkspFail("cannot open complete Helix reference");
  TH1* helixMass=dynamic_cast<TH1*>(helixFile->Get("hLambda_InvMass"));
  TH1* helixRead=dynamic_cast<TH1*>(helixFile->Get("hRefMultVsNTOFMatch"));
  TH1* helixAccepted=dynamic_cast<TH1*>(helixFile->Get("hN"));
  if (!helixMass || !helixRead || !helixAccepted || helixRead->GetEntries()!=events ||
      helixAccepted->GetEntries()!=baseline.accepted) return LkspFail("Helix event counters differ");
  TH1D* helix=LambdaComparisonRebin(helixMass,"HelixLambdaCounts",low,high,width);
  helixFile->Close(); delete helixFile;
  if (!helix) return LkspFail("Helix source bins do not match QA domain");
  const Double_t target=LksCount(helix,windows[0]);
  std::vector<LkspRun> runs;
  std::set<std::string> ids,paths;
  ids.insert("baseline"); ids.insert("helix"); paths.insert(baseline.path);
  std::istringstream list(runText); std::string line;
  while (std::getline(list,line)) {
    line=LkspTrim(line);
    if (line.empty() || line[0]=='#') continue;
    const std::string::size_type first=line.find('\t');
    const std::string::size_type second=first==std::string::npos ? first : line.find('\t',first+1);
    if (first==std::string::npos || second==std::string::npos || line.find('\t',second+1)!=std::string::npos)
      return LkspFail("runs manifest needs exactly three tab-separated columns");
    LkspRun run; run.id=LkspTrim(line.substr(0,first));
    run.path=LkspTrim(line.substr(first+1,second-first-1));
    run.label=LkspTrim(line.substr(second+1)); run.mass=0; run.changed=0;
    if (!LkspSafeId(run.id) || run.path.empty() || run.label.empty() ||
        !ids.insert(run.id).second || !paths.insert(run.path).second)
      return LkspFail("invalid/duplicate run identifier or path");
    for (std::string::size_type c=0; c<run.label.size(); ++c)
      if (static_cast<unsigned char>(run.label[c])<32 || static_cast<unsigned char>(run.label[c])>126)
        return LkspFail("legend labels must use printable English/ASCII text");
    if (!LkspLoadRun(run,events,low,high,width,reason)) return LkspFail(reason);
    if (run.accepted!=baseline.accepted) return LkspFail("accepted event counters differ: "+run.id);
    if (!LkspFixedEqual(baseline.configuration,run.configuration,kFALSE,reason) ||
        !LkspFixedEqual(baseline.eventConfiguration,run.eventConfiguration,kTRUE,reason))
      return LkspFail(run.id+": "+reason);
    if (!baseline.applyGeometry && run.applyGeometry)
      return LkspFail("trial re-enables Finder geometry absent from baseline: "+run.id);
    for (Int_t c=0; c<5; ++c) {
      const Bool_t lowerBound=(c==0 || c==1 || c==3);
      if ((lowerBound && run.cuts[c]>baseline.cuts[c]) || (!lowerBound && run.cuts[c]<baseline.cuts[c]))
        return LkspFail("trial tightens a KF-specific threshold: "+run.id+" / "+LkspKeys[c]);
      if (run.cuts[c]!=baseline.cuts[c]) ++run.changed;
    }
    runs.push_back(run);
  }
  if (!list.eof() || runs.empty()) return LkspFail("no valid trial runs");
  Int_t best=-1; Double_t bestDifference=0.;
  for (UInt_t i=0; i<runs.size(); ++i) {
    const Double_t difference=TMath::Abs(LksCount(runs[i].mass,windows[0])-target);
    if (best<0 || difference<bestDifference ||
        (difference==bestDifference && runs[i].changed<runs[best].changed)) {
      best=i; bestDifference=difference;
    }
  }
  const LksStats hs=LksMeasure(helix,windows[0],left[0],right[0]);
  const LksStats bs=LksMeasure(baseline.mass,windows[0],left[0],right[0]);
  const LksStats ms=LksMeasure(runs[best].mass,windows[0],left[0],right[0]);
  std::ostringstream summary;
  summary << std::setprecision(17)
      << "KF-specific scan with optional Lambda Finder geometry disabled.\n"
      << "PASS structural/Imp5 QA; active fixed-cut/event metadata equal (config paths, equivalent legacy diagnostic-source spelling, and documented profile-inactive legacy fields excepted).\n"
      << "inputEvents=" << events << ", acceptedEach=" << baseline.accepted << ", trials=" << runs.size() << '\n'
      << "Rule: nearest primary-window count, then fewest changed KF-specific parameters, then listed order. S/B is not optimized.\n"
      << "HelixNwin=" << hs.n << ", baselineNwin=" << bs.n << ", selectedId=" << runs[best].id
      << ", selectedNwin=" << ms.n << ", difference=" << ms.n-target << ", exactMatch=" << (ms.n==target ? "true" : "false") << '\n'
      << "Helix S/B=" << hs.sn << " +/- " << hs.esn << "; baseline S/B=" << bs.sn << " +/- " << bs.esn
      << "; selected S/B=" << ms.sn << " +/- " << ms.esn << '\n';
  summary << "applyLambdaGeometryCuts: baseline=" << (baseline.applyGeometry ? "true" : "false")
      << ", selected=" << (runs[best].applyGeometry ? "true" : "false")
      << "; this flag is reported separately, not counted as a changed KF-specific threshold.\n";
  for (Int_t c=0; c<5; ++c) summary << LkspKeys[c] << ": baseline=" << baseline.cuts[c] << ", selected=" << runs[best].cuts[c] << '\n';
  const std::string scope="Actual rerun selected histograms; no offline recovery of upstream-rejected candidates. Only the five listed KF-specific thresholds and applyLambdaGeometryCuts may differ. The latter optionally disables Lambda Finder geometry; ordinary Maker geometry and input covariance/PID/track quality remain fixed. Other active effective cut metadata and event-selection settings are equal; only documented profile-inactive historical metadata fields and config paths are excluded from equality checks. The equivalent event.vertexByMode.<mode>/event.qaVertexByMode.<mode> diagnostic-source spelling is normalized without changing mode, policy or numerical center. Historical missing topoChi2NdfCut means upstream chi2/NDF < 3; missing applyLambdaGeometryCuts means true. Nwin counts candidates INCLUDING background, not unique events, true Lambda yield or efficiency. Sideband background assumes locally constant density; S/B is not the earlier gaus+pol2 fit metric. Errors are unweighted independent-Poisson window errors; event/candidate correlations, method-difference covariance and sideband-shape systematics are not included. Same-data tuning is exploratory. Matching read/accepted counters does not establish ordered event identity; verify run provenance separately. Manifest hashes are recorded, not recomputed by this macro.";
  summary << scope << '\n';
  if (gSystem->mkdir(outputDirectory,kTRUE)!=0) return LkspFail("cannot exclusively create output directory");
  const TString prefix=TString(outputDirectory)+"/";
  std::ofstream csv((prefix+"kf_specific_cut_scan_metrics.csv").Data());
  std::ofstream report((prefix+"kf_specific_cut_scan_summary.txt").Data());
  if (!csv || !report) return LkspFail("cannot create metrics/summary");
  csv << std::setprecision(17) << "sample,label,rootPath";
  for (Int_t c=0; c<5; ++c) csv << ',' << LkspKeys[c];
  csv << ",applyLambdaGeometryCuts,changedKfParameters,window,sidebandVariation,windowLow,windowHigh,leftLow,leftHigh,rightLow,rightHigh,rawCandidates,totalSelected,displayedCandidates,Nwin,NwinError,Nleft,Nright,alpha,B,BError,S,SError,sidebandSoverB,sidebandSoverBError,purity,purityError,HelixNwin,NwinMinusHelix\n";
  const char* wn[2]={"primary","narrow"}; const char* sn[3]={"default","near","far"};
  for (Int_t w=0; w<2; ++w) for (Int_t s=0; s<3; ++s) {
    const Double_t reference=LksCount(helix,windows[w]);
    LkspMetric(csv,"helix",0,wn[w],sn[s],helix,windows[w],left[s],right[s],reference);
    LkspMetric(csv,baseline.id,&baseline,wn[w],sn[s],baseline.mass,windows[w],left[s],right[s],reference);
    for (UInt_t r=0; r<runs.size(); ++r)
      LkspMetric(csv,runs[r].id,&runs[r],wn[w],sn[s],runs[r].mass,windows[w],left[s],right[s],reference);
  }
  report << summary.str() << "\nQA manifest:\n" << qaText << "\nRuns manifest:\n" << runText;
  csv.close(); report.close();
  if (!csv || !report) return LkspFail("metrics/summary write error");
  gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetTextFont(42);
  TCanvas* overlay=LksOverlay("LambdaKfSpecificCutOverlay",helix,baseline.mass,runs[best].mass,
                              runs[best].label.c_str(),width,windows[0]);
  // This study uses long parameter-name labels. Keep the shared overlay
  // helper untouched and reflow only this canvas legend into three rows.
  TIter legendObjects(overlay->GetListOfPrimitives());
  TObject* legendObject=0;
  Bool_t legendFound=kFALSE;
  while ((legendObject=legendObjects())) {
    TLegend* legend=dynamic_cast<TLegend*>(legendObject);
    if (!legend) continue;
    legend->SetNColumns(1);
    legend->SetX1NDC(0.045); legend->SetX2NDC(0.97);
    legend->SetY1NDC(0.825); legend->SetY2NDC(0.925);
    legend->SetMargin(0.055);
    legend->SetTextSize(0.027);
    legendFound=kTRUE;
  }
  if (!legendFound) return LkspFail("shared overlay did not contain its expected legend");
  overlay->Modified(); overlay->Update();
  TFile* output=TFile::Open(prefix+"kf_specific_cut_scan.root","CREATE");
  if (!output || output->IsZombie()) return LkspFail("cannot create output ROOT");
  output->cd(); Int_t bytes=helix->Write()+baseline.mass->Write()+overlay->Write();
  TH1D* chosen=dynamic_cast<TH1D*>(runs[best].mass->Clone("KfSpecificNearestCounts"));
  chosen->SetDirectory(0); bytes+=chosen->Write();
  TObjString qaObject(qaText.c_str()),runObject(runText.c_str()),summaryObject(summary.str().c_str()),scopeObject(scope.c_str());
  TObjString baselineConfig(baseline.configuration.c_str()),baselineEvent(baseline.eventConfiguration.c_str());
  TObjString display("Left: raw counts. Right: clones scaled by 1/Integral(1,Nbins) over full massMin/massMax, excluding under/overflow, no bin-width or peak normalization. Count histograms and metrics remain raw.");
  bytes+=qaObject.Write("KfSpecificQaManifest")+runObject.Write("KfSpecificRunsManifest")+
    summaryObject.Write("KfSpecificSummary")+scopeObject.Write("KfSpecificScope")+
    baselineConfig.Write("KfBaselineEffectiveConfiguration")+baselineEvent.Write("KfBaselineEventSelectionConfiguration")+
    display.Write("KfOverlayDisplayConvention");
  TDirectory* directory=output->mkdir("runs");
  if (!directory) return LkspFail("cannot create ROOT runs directory");
  for (UInt_t r=0; r<runs.size(); ++r) {
    TDirectory* entry=directory->mkdir(runs[r].id.c_str());
    if (!entry) return LkspFail("cannot create ROOT run directory");
    entry->cd();
    TObjString cutMetadata(runs[r].configuration.c_str()),eventMetadata(runs[r].eventConfiguration.c_str());
    bytes+=runs[r].mass->Write("Counts")+cutMetadata.Write("KFParticleEffectiveConfiguration")+
      eventMetadata.Write("KFEventSelectionConfiguration");
  }
  output->Flush(); const Bool_t failed=output->TestBit(TFile::kWriteError); output->Close();
  if (bytes<=0 || failed || output->TestBit(TFile::kWriteError)) return LkspFail("ROOT output write error");
  overlay->Print(prefix+"kf_specific_cut_overlay.png"); overlay->Print(prefix+"kf_specific_cut_overlay.pdf");
  const char* products[5]={"kf_specific_cut_overlay.png","kf_specific_cut_overlay.pdf","kf_specific_cut_scan.root",
    "kf_specific_cut_scan_metrics.csv","kf_specific_cut_scan_summary.txt"};
  for (Int_t f=0; f<5; ++f) {
    FileStat_t info;
    if (gSystem->GetPathInfo(prefix+products[f],info)!=0 || info.fSize<=0) return LkspFail("missing/empty QA product");
  }
  std::cout << "[compareLambdaKfSpecificCuts] " << summary.str() << "Output: " << outputDirectory << std::endl;
  return 0;
}
