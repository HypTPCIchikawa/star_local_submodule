#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TTree.h>
#include <TError.h>
#include <TSystem.h>
#include <TMath.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <cstdio>

bool mergeAuditError = false;
void mergeAuditHandler(Int_t level, Bool_t abort, const char* location, const char* msg) {
  if (level >= kError) mergeAuditError = true;
  DefaultErrorHandler(level, abort, location, msg);
}
struct MergeAuditHist {
  std::string schema;
  double entries, sum, weighted;
  MergeAuditHist(): entries(0), sum(0), weighted(0) {}
};
bool mergeAuditDir(TDirectory* dir, const std::string& prefix,
                   std::map<std::string, MergeAuditHist>& totals, int& objects) {
  TIter next(dir->GetListOfKeys());
  TKey* key = 0;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (!obj) return false;
    ++objects;
    std::string name = prefix + key->GetName() + ";";
    char cycle[32]; sprintf(cycle, "%d", key->GetCycle()); name += cycle;
    if (name.find('\t') != std::string::npos || name.find('\n') != std::string::npos) return false;
    if (obj->InheritsFrom(TDirectory::Class())) {
      // Directories belong to the TFile; Close() releases them.
      if (!mergeAuditDir((TDirectory*)obj, name + "/", totals, objects)) return false;
    } else {
      if (obj->InheritsFrom(TH1::Class())) {
        TH1* h = (TH1*)obj;
        char shape[128]; sprintf(shape, ":%d:%d:%d:%d", h->GetDimension(),
            h->GetNbinsX(), h->GetNbinsY(), h->GetNbinsZ());
        std::string schema = std::string(h->ClassName()) + shape;
        MergeAuditHist& sum = totals[name];
        if (!sum.schema.empty() && sum.schema != schema) { delete obj; return false; }
        sum.schema = schema;
        sum.entries += h->GetEntries();
        const bool profile = h->InheritsFrom("TProfile") || h->InheritsFrom("TProfile2D")
                             || h->InheritsFrom("TProfile3D");
        int cells = h->GetNbinsX() + 2;
        if (h->GetDimension() > 1) cells *= h->GetNbinsY() + 2;
        if (h->GetDimension() > 2) cells *= h->GetNbinsZ() + 2;
        for (int bin = 0; bin < cells; ++bin) {
          double v = h->GetBinContent(bin);
          if (!TMath::Finite(v)) { delete obj; return false; }
          // Profile bin means are not additive; compare their entries only.
          if (!profile) { sum.sum += v; sum.weighted += v * (bin + 1.0); }
        }
        if (!TMath::Finite(sum.entries) || !TMath::Finite(sum.sum)
            || !TMath::Finite(sum.weighted)) { delete obj; return false; }
      }
      if (obj->InheritsFrom(TTree::Class())) {
        TTree* t = (TTree*)obj;
        for (Long64_t e = 0; e < t->GetEntries(); ++e)
          if (t->GetEntry(e) < 0) { delete obj; return false; }
      }
      delete obj;
    }
    if (mergeAuditError) return false;
  }
  return true;
}
void merge_notify_audit(const char* listPath, const char* reportPath) {
  TH1::AddDirectory(kFALSE);
  SetErrorHandler(mergeAuditHandler);
  std::ifstream list(listPath);
  std::ofstream report(reportPath);
  if (!list.good() || !report.good()) { gSystem->Exit(61); return; }
  std::map<std::string, MergeAuditHist> totals;
  std::string path;
  int files = 0;
  while (std::getline(list, path)) {
    if (path.empty()) { gSystem->Exit(62); return; }
    std::cout << "AUDIT " << path << std::endl;
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered) || mergeAuditError) {
      std::cerr << "AUDIT FAILED " << path << std::endl; gSystem->Exit(63); return;
    }
    int objects = 0;
    bool ok = mergeAuditDir(f, "", totals, objects);
    f->Close(); delete f;
    if (!ok || !objects || mergeAuditError) {
      std::cerr << "AUDIT FAILED " << path << std::endl; gSystem->Exit(64); return;
    }
    report << "FILE\t" << path << "\t" << objects << "\n";
    ++files;
  }
  if (!files || totals.empty()) { gSystem->Exit(65); return; }
  report << std::setprecision(17);
  for (std::map<std::string, MergeAuditHist>::const_iterator i = totals.begin(); i != totals.end(); ++i)
    report << "HIST\t" << i->first << "\t" << i->second.schema << "\t"
           << i->second.entries << "\t" << i->second.sum << "\t" << i->second.weighted << "\n";
  report << "DONE\t" << files << "\n";
  report.close();
  if (!report.good()) { gSystem->Exit(66); return; }
  gSystem->Exit(0);
}
