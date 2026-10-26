// Independent scalar ROOT-tree replay of chosen, additional purity cuts.
// Compile/run with ROOT5 ACLiC. The public entry point exits nonzero on failure.
// Input cuts TSV: selection<TAB>key<TAB>operator<TAB>value
// Input histogram CSV: selection,bin,low,high,count[,error], including bins 0..201.
// Every selection is intersected with stored selected==1 and pdg==3122.
// This validates saved-candidate selection, NOT upstream PicoDst reconstruction.
#include <TFile.h>
#include <TH1.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TTree.h>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

void check_lambda_purity_selection(const char* inputRoot, const char* cutsTsv,
                                  const char* expectedHistogramsCsv, const char* outputRoot);

#ifndef __CINT__
namespace LambdaPurityIndependent {
struct Cut {
  std::string key, op, leaf, secondLeaf;
  double threshold;
  bool absolute;
};
struct Choice {
  std::string name;
  std::vector<Cut> cuts;
  std::vector<double> expected;
  std::vector<Long64_t> entries;
  TH1* histogram;
  Choice() : expected(202,-1.), histogram(0) {}
};

int Fail(const std::string& message) {
  std::cerr << "FAIL independent purity selection: " << message << std::endl;
  return 1;
}
bool ReadLine(std::istream& stream, std::string& line) {
  if(!std::getline(stream,line)) return false;
  if(!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1);
  return true;
}
std::vector<std::string> Split(const std::string& s, char delimiter) {
  std::vector<std::string> fields;
  std::string field;
  std::istringstream stream(s);
  while (std::getline(stream,field,delimiter)) fields.push_back(field);
  if (!s.empty() && s[s.size()-1]==delimiter) fields.push_back("");
  return fields;
}
bool Number(const std::string& s, double& value) {
  if(s.empty()) return false;
  char* end=0;
  errno=0;
  value=std::strtod(s.c_str(),&end);
  return errno==0 && end && *end=='\0' && TMath::Finite(value);
}
bool Identifier(const std::string& s) {
  if(s.empty() || s.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")!=std::string::npos)
    return false;
  return s[0]<'0' || s[0]>'9';
}
// Mapping independently transcribed from StLambdaKFParticleMaker::PassCandidateCuts
// and the Imp5 input PID/gDCA definitions; upper/lower bounds are inclusive.
bool MapCut(Cut& cut) {
  if(cut.op=="eq") cut.op="==";
  else if(cut.op=="le") cut.op="<=";
  else if(cut.op=="ge") cut.op=">=";
  else if(cut.op=="abs_le") cut.op="abs<=";
  cut.absolute=false;
  if(cut.key=="selected") { cut.leaf="selected"; return cut.op=="==" && cut.threshold==1.; }
  if(cut.key=="maxMassError") cut.leaf="massError";
  else if(cut.key=="maxChi2Ndf") cut.leaf="chi2Ndf";
  else if(cut.key=="maxTopoChi2Ndf") cut.leaf="topoChi2Ndf";
  else if(cut.key=="maxDaughterDistance") cut.leaf="daughterDistance";
  else if(cut.key=="maxDistanceToPv") cut.leaf="distanceToPv";
  else if(cut.key=="minDecayLength") cut.leaf="decayLength";
  else if(cut.key=="minDecayLengthSignificance") cut.leaf="decayLengthSignificance";
  else if(cut.key=="minVertexLineSignificance") cut.leaf="vertexLineLengthSignificance";
  else if(cut.key=="minCosPointing") cut.leaf="cosPointing";
  else if(cut.key=="imp5MinDCAProton") cut.leaf="protonDcaToPv";
  else if(cut.key=="imp5MinDCAPion") cut.leaf="pionDcaToPv";
  else if(cut.key=="nSigmaProton") { cut.leaf="protonPidPull"; cut.absolute=true; }
  else if(cut.key=="nSigmaPion") { cut.leaf="pionPidPull"; cut.absolute=true; }
  else if(cut.key=="imp5MaxPathLength") {
    cut.leaf="protonHelixPathLength"; cut.secondLeaf="pionHelixPathLength"; cut.absolute=true;
  } else return false;
  if(cut.threshold<0.) return false; // Disabled bounds are not additional cuts.
  if(cut.key=="minCosPointing" && cut.threshold>1.) return false;
  if(cut.absolute) return cut.op=="abs<=";
  const bool minimum=cut.key.find("min")==0 || cut.key.find("imp5Min")==0;
  return cut.op==(minimum ? ">=" : "<=");
}
TLeaf* Scalar(TTree* tree, const std::string& name) {
  TLeaf* leaf=tree->GetLeaf(name.c_str());
  if(!leaf || leaf->GetLeafCount() || leaf->GetLenStatic()!=1) return 0;
  return leaf;
}
bool Pass(double value,const Cut& cut) {
  if(cut.absolute) value=TMath::Abs(value);
  if(cut.op=="==") return value==cut.threshold;
  if(cut.op==">=") return value>=cut.threshold;
  return value<=cut.threshold;
}

int Verify(const char* inputRoot, const char* cutsPath, const char* histogramPath,
           const char* outputRoot) {
  if(!inputRoot || !cutsPath || !histogramPath || !outputRoot || !*outputRoot)
    return Fail("missing file argument");
  if(!gSystem->AccessPathName(outputRoot)) return Fail("refusing existing output file");
  std::ifstream cutsFile(cutsPath), histogramFile(histogramPath);
  if(!cutsFile || !histogramFile) return Fail("cannot open cut/histogram manifest");
  std::string line;
  if(!ReadLine(cutsFile,line) || line!="selection\tkey\toperator\tvalue")
    return Fail("unexpected cuts TSV header");
  std::vector<Choice> choices;
  std::map<std::string,unsigned int> choiceIndex;
  std::set<std::string> uniqueCuts;
  std::ostringstream cutText;
  cutText<<line<<"\n";
  while(ReadLine(cutsFile,line)) {
    if(line.empty()) return Fail("empty cut row");
    cutText<<line<<"\n";
    const std::vector<std::string> f=Split(line,'\t');
    Cut cut;
    if(f.size()!=4 || !Identifier(f[0]) || !Number(f[3],cut.threshold))
      return Fail("malformed cut row");
    cut.key=f[1]; cut.op=f[2];
    if(!MapCut(cut)) return Fail("unsupported cut key/operator/value: "+f[1]+" "+f[2]);
    if(!uniqueCuts.insert(f[0]+"/"+f[1]).second) return Fail("duplicate cut key in selection");
    if(choiceIndex.find(f[0])==choiceIndex.end()) {
      choiceIndex[f[0]]=choices.size();
      Choice choice; choice.name=f[0]; choices.push_back(choice);
    }
    choices[choiceIndex[f[0]]].cuts.push_back(cut);
  }
  if(!cutsFile.eof() || choices.empty()) return Fail("incomplete/empty cut manifest");

  TFile source(inputRoot,"READ");
  if(source.IsZombie() || source.TestBit(TFile::kRecovered)) return Fail("invalid/recovered source ROOT");
  TTree* tree=dynamic_cast<TTree*>(source.Get("KfLambdaCandidates"));
  TH1* native=dynamic_cast<TH1*>(source.Get("hKfLambdaMassSelected"));
  if(!tree || !native || native->GetDimension()!=1 || native->GetNbinsX()!=200 ||
     TMath::Abs(native->GetXaxis()->GetXmin()-1.05)>1.e-12 ||
     TMath::Abs(native->GetXaxis()->GetXmax()-1.25)>1.e-12)
    return Fail("missing tree or unexpected native 200-bin mass axis");
  if(!ReadLine(histogramFile,line)) return Fail("empty histogram CSV");
  const bool withErrors=line=="selection,bin,low,high,count,error";
  if(!withErrors && line!="selection,bin,low,high,count") return Fail("unexpected histogram CSV header");
  while(ReadLine(histogramFile,line)) {
    const std::vector<std::string> f=Split(line,',');
    double binValue=0.,low=0.,high=0.,count=0.,error=0.;
    if(f.size()!=(withErrors ? 6u : 5u) || choiceIndex.find(f[0])==choiceIndex.end() ||
       !Number(f[1],binValue) || binValue<0. || binValue>201. || binValue!=TMath::Floor(binValue) ||
       !Number(f[2],low) || !Number(f[3],high) || !Number(f[4],count) ||
       count<0. || count!=TMath::Floor(count) || (withErrors && (!Number(f[5],error) || error<0.)))
      return Fail("malformed/unknown histogram row");
    const int bin=static_cast<int>(binValue);
    Choice& choice=choices[choiceIndex[f[0]]];
    if(choice.expected[bin]>=0.) return Fail("duplicate expected histogram bin");
    if(TMath::Abs(low-native->GetBinLowEdge(bin))>1.e-12 ||
       TMath::Abs(high-native->GetBinLowEdge(bin)-native->GetBinWidth(bin))>1.e-12)
      return Fail("histogram bin edge mismatch");
    if(withErrors && TMath::Abs(error*error-count)>1.e-7) return Fail("unexpected count error");
    choice.expected[bin]=count;
  }
  if(!histogramFile.eof()) return Fail("cannot read complete histogram manifest");

  std::map<std::string,TLeaf*> leaves;
  const char* required[]={"selected","pdg","mass","runId","eventId","protonId","pionId","protonIndex","pionIndex"};
  for(unsigned int i=0;i<sizeof(required)/sizeof(required[0]);++i) {
    leaves[required[i]]=Scalar(tree,required[i]);
    if(!leaves[required[i]]) return Fail(std::string("missing/non-scalar leaf: ")+required[i]);
  }
  bool hasBaseline=false;
  for(unsigned int c=0;c<choices.size();++c) {
    Choice& choice=choices[c];
    for(int b=0;b<=201;++b) if(choice.expected[b]<0.) return Fail("missing native/flow bin");
    choice.histogram=static_cast<TH1*>(native->Clone(("hMass_"+choice.name).c_str()));
    choice.histogram->SetDirectory(0); choice.histogram->Reset("ICES");
    if(choice.cuts.size()==1 && choice.cuts[0].key=="selected") hasBaseline=true;
    for(unsigned int k=0;k<choice.cuts.size();++k) {
      const Cut& cut=choice.cuts[k];
      leaves[cut.leaf]=Scalar(tree,cut.leaf);
      if(!leaves[cut.leaf]) return Fail("missing/non-scalar cut leaf: "+cut.leaf);
      if(!cut.secondLeaf.empty()) {
        leaves[cut.secondLeaf]=Scalar(tree,cut.secondLeaf);
        if(!leaves[cut.secondLeaf]) return Fail("missing/non-scalar second cut leaf");
      }
    }
  }
  if(!hasBaseline) return Fail("require one selected==1-only baseline selection");
  for(Long64_t entry=0;entry<tree->GetEntries();++entry) {
    if(tree->GetEntry(entry)<=0) return Fail("unreadable candidate entry");
    const double selected=leaves["selected"]->GetValue(), pdg=leaves["pdg"]->GetValue();
    const double mass=leaves["mass"]->GetValue();
    if((selected!=0. && selected!=1.) || (pdg!=3122. && pdg!=-3122.) || !TMath::Finite(mass))
      return Fail("invalid selection flag/pdg/mass");
    if(selected!=1. || pdg!=3122.) continue;
    // Only fields participating in the selected-candidate study must be finite;
    // raw unselected candidates may legitimately carry invalid path diagnostics.
    for(std::map<std::string,TLeaf*>::const_iterator it=leaves.begin();it!=leaves.end();++it)
      if(!TMath::Finite(it->second->GetValue())) return Fail("nonfinite selected-candidate field: "+it->first);
    for(unsigned int i=3;i<sizeof(required)/sizeof(required[0]);++i) {
      const double value=leaves[required[i]]->GetValue();
      if(value<INT_MIN || value>INT_MAX || value!=TMath::Floor(value)) return Fail("non-integral identity field");
    }
    for(unsigned int c=0;c<choices.size();++c) {
      Choice& choice=choices[c]; bool pass=true;
      for(unsigned int k=0;k<choice.cuts.size();++k) {
        const Cut& cut=choice.cuts[k];
        if(!Pass(leaves[cut.leaf]->GetValue(),cut) ||
           (!cut.secondLeaf.empty() && !Pass(leaves[cut.secondLeaf]->GetValue(),cut))) pass=false;
      }
      if(pass) { choice.histogram->Fill(mass); choice.entries.push_back(entry); }
    }
  }
  for(unsigned int c=0;c<choices.size();++c) {
    Choice& choice=choices[c];
    for(int b=0;b<=201;++b) {
      if(choice.histogram->GetBinContent(b)!=choice.expected[b]) {
        std::ostringstream error; error<<choice.name<<" bin "<<b<<": ROOT="
          <<choice.histogram->GetBinContent(b)<<", scanner="<<choice.expected[b];
        return Fail(error.str());
      }
      if(choice.cuts.size()==1 && choice.cuts[0].key=="selected" &&
         choice.histogram->GetBinContent(b)!=native->GetBinContent(b))
        return Fail("baseline tree does not reproduce original selected histogram");
    }
  }
  // All comparisons succeeded before creating the output. CREATE never overwrites.
  TFile output(outputRoot,"CREATE");
  if(output.IsZombie()) return Fail("cannot create new output ROOT");
  TObjString provenance((std::string("inputRoot: ")+inputRoot+"\ncuts: "+cutsPath+
     "\nexpectedHistograms: "+histogramPath+"\nselection: stored selected==1 && pdg==3122\n"
     "validation: exact 200 native bins plus underflow/overflow; independent scalar leaf comparisons\n"
     "scope: saved-candidate replay, not a PicoDst rerun\n").c_str());
  if(provenance.Write("IndependentSelectionProvenance")<=0) return Fail("metadata write failed");
  TObjString chosen(cutText.str().c_str());
  if(chosen.Write("ChosenCutsTSV")<=0) return Fail("cut metadata write failed");
  for(unsigned int c=0;c<choices.size();++c) {
    Choice& choice=choices[c];
    if(choice.histogram->Write()<=0) return Fail("histogram write failed");
    Long64_t sourceEntry=0; Int_t runId=0,eventId=0,protonId=0,pionId=0,protonIndex=0,pionIndex=0,pdg=0;
    Double_t mass=0.;
    TTree identities(("Candidates_"+choice.name).c_str(),"Selected candidate identities in original tree order");
    identities.Branch("sourceEntry",&sourceEntry,"sourceEntry/L");
    identities.Branch("runId",&runId,"runId/I"); identities.Branch("eventId",&eventId,"eventId/I");
    identities.Branch("protonId",&protonId,"protonId/I"); identities.Branch("pionId",&pionId,"pionId/I");
    identities.Branch("protonIndex",&protonIndex,"protonIndex/I"); identities.Branch("pionIndex",&pionIndex,"pionIndex/I");
    identities.Branch("pdg",&pdg,"pdg/I"); identities.Branch("mass",&mass,"mass/D");
    for(unsigned int i=0;i<choice.entries.size();++i) {
      sourceEntry=choice.entries[i];
      if(tree->GetEntry(sourceEntry)<=0) return Fail("candidate reread failed");
      runId=static_cast<Int_t>(leaves["runId"]->GetValue()); eventId=static_cast<Int_t>(leaves["eventId"]->GetValue());
      protonId=static_cast<Int_t>(leaves["protonId"]->GetValue()); pionId=static_cast<Int_t>(leaves["pionId"]->GetValue());
      protonIndex=static_cast<Int_t>(leaves["protonIndex"]->GetValue()); pionIndex=static_cast<Int_t>(leaves["pionIndex"]->GetValue());
      pdg=static_cast<Int_t>(leaves["pdg"]->GetValue()); mass=leaves["mass"]->GetValue();
      if(identities.Fill()<0) return Fail("candidate identity write failed");
    }
    if(identities.Write()<=0) return Fail("candidate tree write failed");
    identities.SetDirectory(0);
    std::cout<<"PASS independent selection "<<choice.name<<": candidates="<<choice.entries.size()
             <<", exact native+flow bins=202"<<std::endl;
  }
  output.Close(); source.Close();
  for(unsigned int c=0;c<choices.size();++c) delete choices[c].histogram;
  return 0;
}
} // namespace LambdaPurityIndependent

void check_lambda_purity_selection(const char* inputRoot, const char* cutsTsv,
                                  const char* expectedHistogramsCsv, const char* outputRoot) {
  gSystem->Exit(LambdaPurityIndependent::Verify(inputRoot,cutsTsv,expectedHistogramsCsv,outputRoot));
}

#endif // ROOT5 dictionary sees only the public entry point
