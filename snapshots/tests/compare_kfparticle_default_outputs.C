// ROOT5 ACLiC, read-only default-behavior regression on actual PicoDst outputs.
// Metadata strings may legitimately gain defaulted configuration keys. This
// test instead requires exact histogram values and ordered candidate values.
// Return value is zero on PASS; a runner may call gSystem->Exit(returnValue).
#include "check_kfparticle_imp5_output.C"
#include "TDirectory.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TAxis.h"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <climits>

namespace {
typedef std::map<std::string,TH1*> KfdHistogramMap;
typedef std::map<std::string,TLeaf*> KfdLeafMap;

static Bool_t KfdEqual(Double_t a, Double_t b) {
  // Compare values, not bit patterns: both NaN is the only unequal-by-== pair
  // accepted here. No tolerance, rounding, rebinning, or normalization is used.
  return a==b || (a!=a && b!=b);
}

static Bool_t KfdCollect(TDirectory* directory, const std::string& prefix,
                         KfdHistogramMap& histograms, std::string& reason) {
  if (!directory || !directory->GetListOfKeys()) {
    reason="missing ROOT directory/key list"; return kFALSE;
  }
  std::set<std::string> names;
  TIter next(directory->GetListOfKeys());
  TKey* key=0;
  while ((key=dynamic_cast<TKey*>(next()))) {
    const std::string name=key->GetName();
    // Get(name) consistently reads the newest stored cycle for each key.
    if (!names.insert(name).second) continue;
    TObject* object=directory->Get(name.c_str());
    if (!object) { reason="cannot read ROOT object: "+prefix+name; return kFALSE; }
    TH1* histogram=dynamic_cast<TH1*>(object);
    if (histogram) {
      if (!histograms.insert(std::make_pair(prefix+name,histogram)).second) {
        reason="duplicate histogram path: "+prefix+name; return kFALSE;
      }
      continue;
    }
    TDirectory* child=dynamic_cast<TDirectory*>(object);
    if (child && !KfdCollect(child,prefix+name+"/",histograms,reason)) return kFALSE;
  }
  return kTRUE;
}

static Bool_t KfdAxis(const TAxis* oldAxis, const TAxis* newAxis,
                      const std::string& location, std::string& reason) {
  if (!oldAxis || !newAxis || oldAxis->GetNbins()!=newAxis->GetNbins()) {
    reason="axis bin count differs: "+location; return kFALSE;
  }
  for (Int_t edge=1; edge<=oldAxis->GetNbins()+1; ++edge)
    if (!KfdEqual(oldAxis->GetBinLowEdge(edge),newAxis->GetBinLowEdge(edge))) {
      std::ostringstream out; out << "axis edge differs: " << location << " edge=" << edge;
      reason=out.str(); return kFALSE;
    }
  for (Int_t bin=0; bin<=oldAxis->GetNbins()+1; ++bin)
    if (TString(oldAxis->GetBinLabel(bin))!=TString(newAxis->GetBinLabel(bin))) {
      std::ostringstream out; out << "axis bin label differs: " << location << " bin=" << bin;
      reason=out.str(); return kFALSE;
    }
  return kTRUE;
}

static Bool_t KfdHistogram(TH1* oldHistogram, TH1* newHistogram,
                           const std::string& path, Long64_t& comparedCells,
                           std::string& reason) {
  const Int_t dimensions=oldHistogram->GetDimension();
  if (TString(oldHistogram->ClassName())!=TString(newHistogram->ClassName()) ||
      dimensions!=newHistogram->GetDimension() || dimensions<1 || dimensions>3) {
    reason="histogram class/dimension differs or is unsupported: "+path; return kFALSE;
  }
  const TAxis* oldAxes[3]={oldHistogram->GetXaxis(),oldHistogram->GetYaxis(),oldHistogram->GetZaxis()};
  const TAxis* newAxes[3]={newHistogram->GetXaxis(),newHistogram->GetYaxis(),newHistogram->GetZaxis()};
  Long64_t cells=1;
  const char* axisNames[3]={"X","Y","Z"};
  for (Int_t axis=0; axis<dimensions; ++axis) {
    if (!KfdAxis(oldAxes[axis],newAxes[axis],path+"/"+axisNames[axis],reason)) return kFALSE;
    cells*=Long64_t(oldAxes[axis]->GetNbins())+2;
  }
  if (cells<=0 || cells>INT_MAX) { reason="unsupported histogram cell count: "+path; return kFALSE; }
  if (!KfdEqual(oldHistogram->GetEntries(),newHistogram->GetEntries())) {
    std::ostringstream out;
    out << std::setprecision(17) << "histogram entries differ: " << path
        << " old=" << oldHistogram->GetEntries() << " new=" << newHistogram->GetEntries();
    reason=out.str(); return kFALSE;
  }
  // ROOT5 has no public GetNcells. This is its full global-bin layout,
  // including every combination of underflow and overflow coordinates.
  for (Int_t cell=0; cell<Int_t(cells); ++cell) {
    const Double_t oldValue=oldHistogram->GetBinContent(cell),newValue=newHistogram->GetBinContent(cell);
    const Double_t oldError=oldHistogram->GetBinError(cell),newError=newHistogram->GetBinError(cell);
    if (!KfdEqual(oldValue,newValue) || !KfdEqual(oldError,newError)) {
      std::ostringstream out;
      out << std::setprecision(17) << "histogram cell differs: " << path << " globalBin=" << cell
          << " oldContent=" << oldValue << " newContent=" << newValue
          << " oldError=" << oldError << " newError=" << newError;
      reason=out.str(); return kFALSE;
    }
  }
  comparedCells+=cells;
  return kTRUE;
}

static Bool_t KfdLeaves(TTree* tree, KfdLeafMap& leaves, std::string& reason) {
  if (!tree || !tree->GetListOfLeaves()) { reason="candidate tree/leaves missing"; return kFALSE; }
  TObjArray* list=tree->GetListOfLeaves();
  for (Int_t i=0; i<list->GetEntriesFast(); ++i) {
    TLeaf* leaf=dynamic_cast<TLeaf*>(list->At(i));
    if (!leaf || !leaf->GetBranch()) { reason="invalid candidate leaf/branch"; return kFALSE; }
    const std::string name=std::string(leaf->GetBranch()->GetName())+":"+leaf->GetName();
    if (!leaves.insert(std::make_pair(name,leaf)).second) {
      reason="duplicate candidate leaf identity: "+name; return kFALSE;
    }
  }
  if (leaves.empty()) { reason="candidate leaf set is empty"; return kFALSE; }
  return kTRUE;
}

static Bool_t KfdInteger(const TString& type) {
  return type=="Bool_t" || type=="Char_t" || type=="UChar_t" || type=="Short_t" ||
    type=="UShort_t" || type=="Int_t" || type=="UInt_t" || type=="Long_t" ||
    type=="ULong_t" || type=="Long64_t" || type=="ULong64_t";
}

static Int_t KfdCompare(TFile* oldFile, TFile* newFile) {
  std::string reason;
  KfdHistogramMap oldHistograms,newHistograms;
  if (!KfdCollect(oldFile,"",oldHistograms,reason) || !KfdCollect(newFile,"",newHistograms,reason)) {
    std::cerr << "FAIL " << reason << std::endl; return 3;
  }
  if (oldHistograms.empty() || oldHistograms.size()!=newHistograms.size()) {
    std::cerr << "FAIL histogram key sets differ: old=" << oldHistograms.size()
              << " new=" << newHistograms.size() << std::endl; return 3;
  }
  Long64_t cells=0;
  for (KfdHistogramMap::const_iterator i=oldHistograms.begin(); i!=oldHistograms.end(); ++i) {
    const KfdHistogramMap::const_iterator j=newHistograms.find(i->first);
    if (j==newHistograms.end()) {
      std::cerr << "FAIL missing new histogram: " << i->first << std::endl; return 3;
    }
    if (!KfdHistogram(i->second,j->second,i->first,cells,reason)) {
      std::cerr << "FAIL " << reason << std::endl; return 4;
    }
  }
  TTree* oldTree=dynamic_cast<TTree*>(oldFile->Get("KfLambdaCandidates"));
  TTree* newTree=dynamic_cast<TTree*>(newFile->Get("KfLambdaCandidates"));
  KfdLeafMap oldLeaves,newLeaves;
  if (!KfdLeaves(oldTree,oldLeaves,reason) || !KfdLeaves(newTree,newLeaves,reason)) {
    std::cerr << "FAIL " << reason << std::endl; return 5;
  }
  if (oldLeaves.size()!=newLeaves.size() || oldTree->GetEntries()!=newTree->GetEntries()) {
    std::cerr << "FAIL candidate leaf/entry counts differ: oldLeaves=" << oldLeaves.size()
        << " newLeaves=" << newLeaves.size() << " oldEntries=" << oldTree->GetEntries()
        << " newEntries=" << newTree->GetEntries() << std::endl; return 5;
  }
  std::vector<TLeaf*> before,after;
  std::vector<std::string> names;
  std::vector<Bool_t> integral;
  for (KfdLeafMap::const_iterator i=oldLeaves.begin(); i!=oldLeaves.end(); ++i) {
    const KfdLeafMap::const_iterator j=newLeaves.find(i->first);
    if (j==newLeaves.end()) {
      std::cerr << "FAIL missing new leaf: " << i->first << std::endl; return 5;
    }
    TLeaf* a=i->second; TLeaf* b=j->second;
    const TString type(a->GetTypeName());
    const TString countA=a->GetLeafCount() ? a->GetLeafCount()->GetName() : "";
    const TString countB=b->GetLeafCount() ? b->GetLeafCount()->GetName() : "";
    if (type!=b->GetTypeName() || a->GetLenStatic()!=b->GetLenStatic() ||
        a->IsUnsigned()!=b->IsUnsigned() || countA!=countB ||
        TString(a->ClassName())!=TString(b->ClassName())) {
      std::cerr << "FAIL leaf schema differs: " << i->first << std::endl; return 5;
    }
    const Bool_t integer=KfdInteger(type);
    if (!integer && type!="Float_t" && type!="Double_t" && type!="Float16_t" && type!="Double32_t") {
      std::cerr << "FAIL unsupported candidate leaf type: " << i->first << " / " << type << std::endl;
      return 5;
    }
    before.push_back(a); after.push_back(b); names.push_back(i->first); integral.push_back(integer);
  }
  Long64_t values=0;
  for (Long64_t entry=0; entry<oldTree->GetEntries(); ++entry) {
    if (oldTree->GetEntry(entry)<=0 || newTree->GetEntry(entry)<=0) {
      std::cerr << "FAIL cannot read ordered candidate entry " << entry << std::endl; return 6;
    }
    for (UInt_t leaf=0; leaf<before.size(); ++leaf) {
      const Int_t length=before[leaf]->GetLen();
      if (length!=after[leaf]->GetLen() || length<0) {
        std::cerr << "FAIL candidate leaf length differs at entry " << entry << " / " << names[leaf] << std::endl;
        return 6;
      }
      for (Int_t element=0; element<length; ++element) {
        // Preserve integer ID precision, including 64-bit leaves, rather than
        // comparing their potentially rounded Double_t GetValue conversions.
        const Bool_t equal=integral[leaf]
          ? before[leaf]->GetValueLong64(element)==after[leaf]->GetValueLong64(element)
          : KfdEqual(before[leaf]->GetValue(element),after[leaf]->GetValue(element));
        if (!equal) {
          std::cerr << std::setprecision(17) << "FAIL candidate value differs: entry=" << entry
              << " leaf=" << names[leaf] << " element=" << element;
          if (integral[leaf]) std::cerr << " old=" << before[leaf]->GetValueLong64(element)
                                      << " new=" << after[leaf]->GetValueLong64(element);
          else std::cerr << " old=" << before[leaf]->GetValue(element)
                         << " new=" << after[leaf]->GetValue(element);
          std::cerr << std::endl; return 6;
        }
        ++values;
      }
    }
  }
  std::cout << "PASS default KF output regression: " << oldHistograms.size()
      << " TH1/TH2/TH3 objects, " << cells << " global cells including all flow bins; "
      << "entries, axes, labels, contents and errors agree exactly.\n"
      << "PASS KfLambdaCandidates: " << oldTree->GetEntries() << " ordered entries, "
      << oldLeaves.size() << " leaves, " << values << " values agree exactly (both-NaN accepted).\n"
      << "Both files passed existing structural/Imp5 QA. Metadata strings are intentionally not compared; "
      << "new defaulted configuration fields are permitted. Source files were opened READ only." << std::endl;
  return 0;
}
} // anonymous namespace

Int_t compare_kfparticle_default_outputs(const char* oldRoot, const char* newRoot) {
  if (!oldRoot || !*oldRoot || !newRoot || !*newRoot) {
    std::cerr << "FAIL two source ROOT paths are required" << std::endl; return 1;
  }
  if (KfImp5OutputQaInspect(oldRoot,kFALSE)!=0 || KfImp5OutputQaInspect(newRoot,kFALSE)!=0) {
    std::cerr << "FAIL source structural/Imp5 QA" << std::endl; return 1;
  }
  TFile* oldFile=TFile::Open(oldRoot,"READ");
  TFile* newFile=TFile::Open(newRoot,"READ");
  if (!oldFile || !newFile || oldFile->IsZombie() || newFile->IsZombie() ||
      oldFile->TestBit(TFile::kRecovered) || newFile->TestBit(TFile::kRecovered)) {
    std::cerr << "FAIL cannot open complete source ROOT files" << std::endl;
    if (oldFile) { oldFile->Close(); delete oldFile; }
    if (newFile) { newFile->Close(); delete newFile; }
    return 2;
  }
  const Int_t result=KfdCompare(oldFile,newFile);
  oldFile->Close(); newFile->Close(); delete oldFile; delete newFile;
  return result;
}
