// Read-only ROOT5 PicoDst recovery audit, matching FemtoLambdaRunSupport branches.
// Compile once with .C+ before concurrent audits; ROOT5 CINT cannot install the
// error handler. Reuse the prebuilt _C.so through .C+ in subsequent processes.
// Exit: 0 = valid/nonempty; 10 = structurally valid empty; 20 = invalid.
#include "TROOT.h"
#include "TSystem.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TClass.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

Long64_t recoveryPicoRootErrors = 0;
Long64_t recoveryPicoCurrentEntry = -1;
std::vector<std::string> recoveryPicoErrors;

void recoveryPicoRecord(const char* text) {
  if (recoveryPicoErrors.size() < 100)
    recoveryPicoErrors.push_back(text ? text : "unspecified audit error");
}

void recoveryPicoErrorHandler(Int_t level, Bool_t abortFlag,
                             const char* location, const char* message) {
  if (level >= kError) {
    ++recoveryPicoRootErrors;
    TString text = TString::Format("ROOT level=%d entry=%lld %s: %s", level,
        recoveryPicoCurrentEntry, location ? location : "", message ? message : "");
    recoveryPicoRecord(text.Data());
  }
  DefaultErrorHandler(level, abortFlag, location, message);
}

std::string recoveryPicoJson(const char* raw) {
  std::string result = "\"";
  const unsigned char* p = reinterpret_cast<const unsigned char*>(raw ? raw : "");
  const char* hex = "0123456789abcdef";
  for (; *p; ++p) {
    if (*p == '\\') result += "\\\\";
    else if (*p == '"') result += "\\\"";
    else if (*p == '\n') result += "\\n";
    else if (*p == '\r') result += "\\r";
    else if (*p == '\t') result += "\\t";
    else if (*p < 32) {
      result += "\\u00";
      result += hex[*p >> 4];
      result += hex[*p & 15];
    } else result += static_cast<char>(*p);
  }
  result += "\"";
  return result;
}

void checkRecoveryPico(const char* input, const char* report) {
  recoveryPicoRootErrors = 0;
  recoveryPicoCurrentEntry = -1;
  recoveryPicoErrors.clear();
  if (!input || !*input || !report || !*report || !gSystem->AccessPathName(report)) {
    std::cerr << "ERROR: explicit input and new report path are required" << std::endl;
    gSystem->Exit(20); return;
  }
  std::ofstream output(report);
  if (!output) {
    std::cerr << "ERROR: cannot create audit report " << report << std::endl;
    gSystem->Exit(20); return;
  }
  ErrorHandlerFunc_t previousHandler = SetErrorHandler(recoveryPicoErrorHandler);
  Bool_t valid = kTRUE;
  Bool_t recovered = kFALSE;
  Long64_t entries = -1;
  Long64_t entriesRead = 0;
  Long64_t bytesRead = 0;
  TFile* file = 0;
  TTree* tree = 0;
  const char* names[] = {"Event", "Track", "TrackCovMatrix", "BTofPidTraits", "BTofHit"};
  const char* classes[] = {"StPicoEvent", "StPicoTrack", "StPicoTrackCovMatrix",
                           "StPicoBTofPidTraits", "StPicoBTofHit"};
  TBranch* branches[5] = {0, 0, 0, 0, 0};
  TClonesArray* arrays[5] = {0, 0, 0, 0, 0};
  std::vector<std::string> active;
  if (gSystem->Load("StarRoot") < 0 || gSystem->Load("St_base") < 0
      || gSystem->Load("StarClassLibrary") < 0 || gSystem->Load("StPicoEvent") < 0) {
    recoveryPicoRecord("Cannot load STAR PicoDst dictionaries");
    valid = kFALSE;
  }
  if (valid) {
    file = TFile::Open(input, "READ");
    if (!file || file->IsZombie()) {
      recoveryPicoRecord("Input is missing, unreadable, or a zombie ROOT file");
      valid = kFALSE;
    } else {
      recovered = file->TestBit(TFile::kRecovered);
      if (recovered) {
        recoveryPicoRecord("Recovered ROOT file is not accepted as a valid source");
        valid = kFALSE;
      }
      tree = dynamic_cast<TTree*>(file->Get("PicoDst"));
      if (!tree) {
        recoveryPicoRecord("Missing PicoDst TTree");
        valid = kFALSE;
      }
    }
  }
  if (valid) {
    entries = tree->GetEntries();
    if (entries < 0) {
      recoveryPicoRecord("PicoDst TTree has invalid entry count");
      valid = kFALSE;
    }
    tree->SetBranchStatus("*", 0);
    for (Int_t i = 0; i < 5; ++i) {
      branches[i] = tree->GetBranch(names[i]);
      if (!branches[i]) branches[i] = tree->GetBranch(TString(names[i]) + ".");
      if (!branches[i]) {
        if (i < 4) {
          recoveryPicoRecord((TString("Missing required branch ") + names[i]).Data());
          valid = kFALSE;
        }
        continue;
      }
      if (!TClass::GetClass(classes[i])) {
        recoveryPicoRecord((TString("Missing class dictionary ") + classes[i]).Data());
        valid = kFALSE; continue;
      }
      if (branches[i]->GetEntries() != entries) {
        recoveryPicoRecord((TString("Branch entry count differs from TTree: ") + names[i]).Data());
        valid = kFALSE; continue;
      }
      tree->SetBranchStatus(TString(names[i]) + "*", 1);
      arrays[i] = new TClonesArray(classes[i]);
      const Int_t addressResult = tree->SetBranchAddress(branches[i]->GetName(), &arrays[i]);
      if (addressResult < 0) {
        recoveryPicoRecord((TString("Cannot bind TClonesArray branch ") + names[i]).Data());
        valid = kFALSE; continue;
      }
      active.push_back(names[i]);
    }
  }
  if (recoveryPicoRootErrors) valid = kFALSE;
  for (Long64_t entry = 0; valid && entry < entries; ++entry) {
    recoveryPicoCurrentEntry = entry;
    for (Int_t b = 0; b < 5; ++b) {
      if (!branches[b]) continue;
      arrays[b]->Clear();
      const Int_t bytes = branches[b]->GetEntry(entry);
      if (bytes <= 0 || recoveryPicoRootErrors) {
        recoveryPicoRecord(TString::Format("Read failed branch=%s entry=%lld bytes=%d",
                                          names[b], entry, bytes).Data());
        valid = kFALSE; break;
      }
      bytesRead += bytes;
    }
    if (valid && arrays[0]->GetEntriesFast() != 1) {
      recoveryPicoRecord(TString::Format("Expected one StPicoEvent at entry=%lld; got=%d",
                                        entry, arrays[0]->GetEntriesFast()).Data());
      valid = kFALSE;
    }
    if (valid) ++entriesRead;
    if (entry % 10000 == 0)
      std::cout << "[recovery-audit] " << input << " read=" << entriesRead
                << "/" << entries << std::endl;
  }
  recoveryPicoCurrentEntry = -1;
  if (tree) tree->ResetBranchAddresses();
  for (Int_t i = 0; i < 5; ++i) { delete arrays[i]; arrays[i] = 0; }
  if (file) { file->Close(); delete file; file = 0; }
  if (recoveryPicoRootErrors || entriesRead != entries) valid = kFALSE;
  const char* status = valid ? (entries == 0 ? "empty" : "valid") : "invalid";
  const Int_t exitCode = valid ? (entries == 0 ? 10 : 0) : 20;
  output << "{\n  \"schemaVersion\": 1,\n  \"status\": " << recoveryPicoJson(status)
         << ",\n  \"input\": " << recoveryPicoJson(input)
         << ",\n  \"entries\": " << entries
         << ",\n  \"entriesRead\": " << entriesRead
         << ",\n  \"activeBranchBytesRead\": " << bytesRead
         << ",\n  \"recovered\": " << (recovered ? "true" : "false")
         << ",\n  \"rootErrors\": " << recoveryPicoRootErrors
         << ",\n  \"activeBranches\": [";
  for (unsigned int a = 0; a < active.size(); ++a) {
    if (a) output << ", ";
    output << recoveryPicoJson(active[a].c_str());
  }
  output << "],\n  \"errors\": [";
  for (unsigned int e = 0; e < recoveryPicoErrors.size(); ++e) {
    if (e) output << ", ";
    output << recoveryPicoJson(recoveryPicoErrors[e].c_str());
  }
  output << "],\n  \"exitCode\": " << exitCode << "\n}\n";
  output.close();
  SetErrorHandler(previousHandler);
  if (!output) {
    std::cerr << "ERROR: could not finish audit report " << report << std::endl;
    gSystem->Exit(20); return;
  }
  std::cout << "[recovery-audit] status=" << status << " entries=" << entries
            << " entriesRead=" << entriesRead << " report=" << report << std::endl;
  gSystem->Exit(exitCode);
}
