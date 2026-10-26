// Read-only stored-histogram extraction; no CF calculation, fits or event cuts.
// ROOT5 / ACLiC: ExtractMtCF(inputRoot, outputDir, centMin, centMax)
#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TList.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TNamed.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "THnBase.h"
#include "TSystem.h"
#include "TMath.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace MtCFExtract {
void Require(bool yes,const std::string &message) {
  if(!yes)throw std::runtime_error(message);
}
std::string Quote(const std::string &s) {
  std::ostringstream o;o<<'"';
  for(std::string::const_iterator i=s.begin();i!=s.end();++i) {
    if(*i=='"' || *i=='\\')o<<'\\'<<*i;
    else if(*i=='\n')o<<"\\n";
    else if(*i=='\r')o<<"\\r";
    else if(*i=='\t')o<<"\\t";
    else o<<*i;
  }
  o<<'"';return o.str();
}
void Num(std::ostream &o,double v) {
  Require(TMath::Finite(v),"nonfinite numeric output");o<<std::setprecision(17)<<v;
}
void Array(std::ostream &o,const std::vector<double> &v) {
  o<<'[';for(unsigned i=0;i<v.size();++i) {if(i)o<<',';Num(o,v[i]);}o<<']';
}
bool Close(double a,double b) {
  return std::fabs(a-b)<=1.e-12*std::max(1.,std::max(std::fabs(a),std::fabs(b)));
}
void SameAxis(const TAxis *a,const TAxis *b,const std::string &context) {
  Require(a->GetNbins()==b->GetNbins(),"bin count mismatch: "+context);
  for(int i=1;i<=a->GetNbins()+1;++i)
    Require(Close(a->GetBinLowEdge(i),b->GetBinLowEdge(i)),"axis edge mismatch: "+context);
}
void Uniform(const TAxis *a,const std::string &context) {
  Require(a->GetNbins()>0 && a->GetXmax()>a->GetXmin(),"invalid axis: "+context);
  for(int i=2;i<=a->GetNbins();++i)
    Require(Close(a->GetBinWidth(1),a->GetBinWidth(i)),"nonuniform axis unsupported: "+context);
}
std::vector<double> Edges(const TAxis *a) {
  std::vector<double> v;for(int i=1;i<=a->GetNbins()+1;++i)v.push_back(a->GetBinLowEdge(i));return v;
}
double Variance(const TH1 *h,int bin) {
  const double n=h->GetBinContent(bin);
  Require(TMath::Finite(n) && n>=0,std::string("negative/nonfinite stored count: ")+h->GetName());
  if(h->GetSumw2N()==0)return n;
  const double e=h->GetBinError(bin);Require(TMath::Finite(e),"nonfinite stored bin error");return e*e;
}
std::string VariancePolicy(const TH1 *h) {
  return h->GetSumw2N()==0?"unit_fill_no_Sumw2_error2_equals_content":"stored_Sumw2";
}
TH1 *Histogram(TFile *f,const std::string &name,int dimension) {
  TH1 *h=dynamic_cast<TH1*>(f->Get(name.c_str()));
  Require(h && h->GetDimension()==dimension,"missing or wrong histogram dimension: "+name);
  return h;
}
TH2D *DoubleCopy(const TH2 *h,const std::string &name) {
  Uniform(h->GetXaxis(),h->GetName());Uniform(h->GetYaxis(),h->GetName());
  TH2D *d=new TH2D(name.c_str(),h->GetTitle(),h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
    h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  d->SetDirectory(0);d->Sumw2();
  for(int x=0;x<=h->GetNbinsX()+1;++x)for(int y=0;y<=h->GetNbinsY()+1;++y) {
    const int b=h->GetBin(x,y);d->SetBinContent(x,y,h->GetBinContent(b));d->SetBinError(x,y,std::sqrt(Variance(h,b)));
  }
  d->SetEntries(h->GetEntries());return d;
}

struct KeyInfo {
  std::string path,klass;int cycle,dimension;bool directory,jointCandidate;
  Long64_t entries;std::vector<std::string> branches;
  KeyInfo():cycle(0),dimension(0),directory(false),jointCandidate(false),entries(-1) {}
};
struct NamedInfo {std::string path,name,title;int cycle;};
// Fully qualified template arguments are required by ROOT5's CINT dictionary.
void Inventory(TDirectory *dir,const std::string &prefix,
               std::vector<MtCFExtract::KeyInfo> &keys,
               std::vector<MtCFExtract::NamedInfo> &named) {
  TIter next(dir->GetListOfKeys());TKey *key=0;std::set<std::string> traversed;
  while((key=dynamic_cast<TKey*>(next()))) {
    KeyInfo k;k.path=prefix+key->GetName();k.klass=key->GetClassName();k.cycle=key->GetCycle();
    TClass *cl=TClass::GetClass(k.klass.c_str());
    k.directory=(cl && cl->InheritsFrom(TDirectory::Class())) || k.klass=="TDirectoryFile" || k.klass=="TDirectory";
    if(cl && cl->InheritsFrom("TH3"))k.dimension=3;
    else if(cl && cl->InheritsFrom(TH2::Class()))k.dimension=2;
    else if(cl && cl->InheritsFrom(TH1::Class()))k.dimension=1;
    const bool tree=(cl && cl->InheritsFrom(TTree::Class())) || k.klass=="TTree";
    const bool thn=(cl && cl->InheritsFrom(THnBase::Class())) || k.klass.find("THn")==0;
    k.jointCandidate=tree || thn || k.dimension>=3;
    // Only potential joint-data objects and exact TNamed metadata are read.
    // In particular, ordinary TH1/TH2 QA maps are never ReadObj'd for inventory.
    if(tree || thn || k.klass=="TNamed") {
      TObject *obj=key->ReadObj();Require(obj!=0,"cannot read metadata/candidate object: "+k.path);
      if(tree) {
        TTree *t=dynamic_cast<TTree*>(obj);Require(t!=0,"TTree dictionary/type mismatch: "+k.path);
        k.dimension=-1;k.entries=t->GetEntries();
        TIter branches(t->GetListOfBranches());TBranch *br=0;
        while((br=dynamic_cast<TBranch*>(branches())))k.branches.push_back(br->GetName());
      } else if(thn) {
        THnBase *h=dynamic_cast<THnBase*>(obj);Require(h!=0,"THn dictionary/type mismatch: "+k.path);
        k.dimension=h->GetNdimensions();
      } else {
        TNamed *n=dynamic_cast<TNamed*>(obj);Require(n!=0,"TNamed type mismatch: "+k.path);
        NamedInfo item;item.path=k.path;item.name=n->GetName();item.title=n->GetTitle();item.cycle=k.cycle;named.push_back(item);
      }
      delete obj;
    }
    keys.push_back(k);
    if(k.directory && traversed.insert(key->GetName()).second) {
      TDirectory *child=dir->GetDirectory(key->GetName());Require(child!=0,"cannot read directory: "+k.path);
      Inventory(child,k.path+"/",keys,named);
    }
  }
}
void WriteKey(std::ostream &o,const KeyInfo &k) {
  o<<"{\"path\":"<<Quote(k.path)<<",\"class\":"<<Quote(k.klass)<<",\"cycle\":"<<k.cycle<<",\"dimension\":";
  if(k.dimension<0)o<<"null";else o<<k.dimension;
  o<<",\"isDirectory\":"<<(k.directory?"true":"false");
  if(k.entries>=0)o<<",\"entries\":"<<k.entries;
  if(!k.branches.empty()) {o<<",\"branches\":[";for(unsigned i=0;i<k.branches.size();++i) {if(i)o<<',';o<<Quote(k.branches[i]);}o<<']';}
  o<<'}';
}

struct Summed {
  TH2D *map;
  std::vector<std::string> sources,sourceClasses,policies,referenceSources;
  std::vector<double> sourceEntries,reference,referenceError2;
  std::string referenceKind;
  Summed():map(0) {}
};
Summed SumChannel(TFile *f,bool mixed,int region,int centMin,int centMax) {
  const char *suffix[]={"","_SBNeg","_SBPos","_signalLow","_signalHigh"};
  const char *channel[]={"signal","leftSB","rightSB","signalLow","signalHigh"};
  const std::string folder=mixed?"mix/":"true/",prefix=mixed?"Mixed_":"",mode=mixed?"ME":"SE";
  Summed out;
  for(int cent=centMin;cent<=centMax;++cent) {
    const std::string path=folder+"hKstarMt_"+prefix+"4He"+suffix[region]+Form("_CentBin%d",cent);
    TH2 *h=static_cast<TH2*>(Histogram(f,path,2));
    Require(h->IsA()->InheritsFrom(TH2D::Class()) || h->IsA()->InheritsFrom(TH2F::Class()),"unsupported mT histogram storage class: "+path);
    if(out.map) {SameAxis(out.map->GetXaxis(),h->GetXaxis(),path);SameAxis(out.map->GetYaxis(),h->GetYaxis(),path);}
    TH2D *copy=DoubleCopy(h,"mt_extract_temporary");
    if(!out.map) {
      out.map=copy;out.map->SetName(Form("mt_extract_%s_%s",mode.c_str(),channel[region]));
      out.reference.assign(h->GetNbinsX()+2,0.);out.referenceError2.assign(h->GetNbinsX()+2,0.);
    } else {out.map->Add(copy);delete copy;}
    out.sources.push_back(path);out.sourceClasses.push_back(h->ClassName());out.policies.push_back(VariancePolicy(h));out.sourceEntries.push_back(h->GetEntries());
    if(region<3) {
      const std::string ref=folder+"hKstar_"+prefix+"4He"+suffix[region]+Form("_CentBin%d",cent);
      TH1 *r=Histogram(f,ref,1);SameAxis(h->GetXaxis(),r->GetXaxis(),ref);
      for(int b=0;b<=r->GetNbinsX()+1;++b) {out.reference[b]+=r->GetBinContent(b);out.referenceError2[b]+=Variance(r,b);}
      out.referenceSources.push_back(ref);out.referenceKind="legacy_1D_histogram_sum";
    }
  }
  if(region>=3) {
    const std::string ref="hKstar"+mode+"VsCent_lambda_he4_"+channel[region];
    TH2 *r=static_cast<TH2*>(Histogram(f,ref,2));SameAxis(out.map->GetXaxis(),r->GetXaxis(),ref);
    Require(r->GetNbinsY()==9 && Close(r->GetYaxis()->GetXmin(),-.5) && Close(r->GetYaxis()->GetXmax(),8.5),
      "half reference must have canonical cent9 Y axis: "+ref);
    for(int b=0;b<=r->GetNbinsX()+1;++b)for(int cent=centMin;cent<=centMax;++cent) {
      const int cell=r->GetBin(b,cent+1);out.reference[b]+=r->GetBinContent(cell);out.referenceError2[b]+=Variance(r,cell);
    }
    out.referenceSources.push_back(ref);out.referenceKind="top_level_kstar_vs_cent_projection";
  }
  return out;
}
void WriteMatrix(std::ostream &o,TH2D *h,bool flow,bool variance) {
  const int first=flow?0:1,endX=h->GetNbinsX()+(flow?1:0),endY=h->GetNbinsY()+(flow?1:0);o<<'[';
  for(int x=first;x<=endX;++x) {if(x>first)o<<',';o<<'[';
    for(int y=first;y<=endY;++y) {if(y>first)o<<',';const int b=h->GetBin(x,y);Num(o,variance?Variance(h,b):h->GetBinContent(b));}o<<']';}
  o<<']';
}
void StringArray(std::ostream &o,const std::vector<std::string> &v) {
  o<<'[';for(unsigned i=0;i<v.size();++i) {if(i)o<<',';o<<Quote(v[i]);}o<<']';
}
void WriteChannel(std::ostream &o,const Summed &s) {
  TH2D *h=s.map;o<<"{\"sourcePaths\":";StringArray(o,s.sources);o<<",\"sourceClasses\":";StringArray(o,s.sourceClasses);
  o<<",\"sourceVariancePolicies\":";StringArray(o,s.policies);o<<",\"sourceEntries\":";Array(o,s.sourceEntries);
  o<<",\"summedClass\":\"TH2D\",\"kstarEdges\":";Array(o,Edges(h->GetXaxis()));o<<",\"mtEdges\":";Array(o,Edges(h->GetYaxis()));
  o<<",\"counts\":";WriteMatrix(o,h,false,false);o<<",\"error2\":";WriteMatrix(o,h,false,true);
  o<<",\"flowCounts\":";WriteMatrix(o,h,true,false);o<<",\"flowError2\":";WriteMatrix(o,h,true,true);
  std::vector<double> pr(h->GetNbinsX()+2,0.),pv(pr.size(),0.),diff(pr.size(),0.),vdiff(pr.size(),0.),inside(pr.size(),0.);
  double maxAbs=0,maxRel=0,total=0,referenceTotal=0,flowTotal=0;
  for(int x=0;x<=h->GetNbinsX()+1;++x) {
    for(int y=0;y<=h->GetNbinsY()+1;++y) {
      const int b=h->GetBin(x,y);pr[x]+=h->GetBinContent(b);pv[x]+=Variance(h,b);
      if(y>=1 && y<=h->GetNbinsY())inside[x]+=h->GetBinContent(b);else flowTotal+=h->GetBinContent(b);
    }
    diff[x]=pr[x]-s.reference[x];vdiff[x]=pv[x]-s.referenceError2[x];maxAbs=std::max(maxAbs,std::fabs(diff[x]));
    maxRel=std::max(maxRel,std::fabs(diff[x])/std::max(1.,std::fabs(s.reference[x])));total+=pr[x];referenceTotal+=s.reference[x];
  }
  o<<",\"projectionClosure\":{\"includesKstarFlow\":true,\"includesMtFlow\":true,\"referenceKind\":"<<Quote(s.referenceKind)
   <<",\"referencePaths\":";StringArray(o,s.referenceSources);
  o<<",\"projectedIncludingMtFlow\":";Array(o,pr);o<<",\"projectedExcludingMtFlow\":";Array(o,inside);
  o<<",\"projectedError2\":";Array(o,pv);o<<",\"referenceKstar\":";Array(o,s.reference);
  o<<",\"referenceError2\":";Array(o,s.referenceError2);o<<",\"difference\":";Array(o,diff);
  o<<",\"varianceDifference\":";Array(o,vdiff);o<<",\"maxAbsDifference\":";Num(o,maxAbs);
  o<<",\"maxRelativeDifference\":";Num(o,maxRel);o<<",\"projectedTotal\":";Num(o,total);
  o<<",\"referenceTotal\":";Num(o,referenceTotal);o<<",\"mtFlowTotal\":";Num(o,flowTotal);
  o<<",\"policy\":\"Raw differences are retained; no exact-equality requirement and no replacement of either stored histogram. Float storage/merge rounding must be assessed by the caller.\"}}";
}

void Run(const char *inputRoot,const char *outputDir,int centMin,int centMax) {
  Require(inputRoot && *inputRoot && outputDir && *outputDir,"explicit input ROOT and output directory required");
  Require(0<=centMin && centMin<=centMax && centMax<=8,"invalid explicit cent9 range");
  const std::string target=std::string(outputDir)+"/mt_extract.json",temporary=target+".tmp";
  Require(gSystem->AccessPathName(target.c_str()) && gSystem->AccessPathName(temporary.c_str()),"refusing to overwrite mt_extract.json or .tmp");
  if(gSystem->AccessPathName(outputDir))Require(gSystem->mkdir(outputDir,kTRUE)==0,"cannot create output directory");
  std::cout<<"MT_EXTRACT_OPEN "<<inputRoot<<std::endl;
  TFile *f=TFile::Open(inputRoot,"READ");Require(f && !f->IsZombie(),"cannot open input ROOT read-only");
  std::vector<KeyInfo> keys;std::vector<NamedInfo> named;Inventory(f,"",keys,named);
  TNamed *definition=dynamic_cast<TNamed*>(f->Get("FemtoLambdaPairMtDefinition"));
  Require(definition!=0,"required stored FemtoLambdaPairMtDefinition is missing");
  std::cout<<"MT_DEFINITION "<<definition->GetTitle()<<std::endl;
  std::ofstream o(temporary.c_str());Require(o.good(),"cannot create JSON output");
  o<<"{\"schemaVersion\":1,\"inputRoot\":"<<Quote(inputRoot)<<",\"centMin\":"<<centMin<<",\"centMax\":"<<centMax
   <<",\"mtDefinition\":"<<Quote(definition->GetTitle())
   <<",\"arrayLayout\":\"x-major: counts[kstarBin-1][mtBin-1]; flowCounts[kstarROOTBin][mtROOTBin], index 0 underflow and N+1 overflow\","
    "\"selectionPolicy\":\"Stored histograms only; no event cuts, mass selection or mT splits are changed/chosen by this extractor.\",\"inventory\":[";
  for(unsigned i=0;i<keys.size();++i) {if(i)o<<',';WriteKey(o,keys[i]);}o<<"],\"namedMetadata\":[";
  for(unsigned i=0;i<named.size();++i) {if(i)o<<',';o<<"{\"path\":"<<Quote(named[i].path)<<",\"name\":"<<Quote(named[i].name)
    <<",\"title\":"<<Quote(named[i].title)<<",\"cycle\":"<<named[i].cycle<<'}';}
  o<<"],\"jointStorageCandidates\":[";unsigned candidates=0;
  for(unsigned i=0;i<keys.size();++i)if(keys[i].jointCandidate) {if(candidates++)o<<',';WriteKey(o,keys[i]);}
  o<<"],\"jointStorageInventoryStatus\":"<<Quote(candidates?"TH3_THn_or_TTree_candidates_require_inspection":"no_TH3_THn_or_TTree_keys_in_recursive_inventory")
   <<",\"channels\":{";
  const char *modes[]={"SE","ME"},*regions[]={"signal","leftSB","rightSB","signalLow","signalHigh"};
  TH2D *axisReference=0;
  for(int m=0;m<2;++m) {
    if(m)o<<',';o<<Quote(modes[m])<<":{";
    for(int r=0;r<5;++r) {
      if(r)o<<',';std::cout<<"MT_EXTRACT_CHANNEL "<<modes[m]<<' '<<regions[r]<<std::endl;
      Summed s=SumChannel(f,m!=0,r,centMin,centMax);
      if(axisReference) {SameAxis(axisReference->GetXaxis(),s.map->GetXaxis(),regions[r]);SameAxis(axisReference->GetYaxis(),s.map->GetYaxis(),regions[r]);}
      else axisReference=static_cast<TH2D*>(s.map->Clone("mt_axis_reference"));
      o<<Quote(regions[r])<<':';WriteChannel(o,s);delete s.map;o.flush();Require(o.good(),"JSON write failed");
    }
    o<<'}';
  }
  delete axisReference;o<<"}}\n";o.close();Require(!o.fail(),"cannot finish JSON output");f->Close();delete f;
  Require(gSystem->Rename(temporary.c_str(),target.c_str())==0,"cannot finalize JSON output");
  std::cout<<"MT_EXTRACT_COMPLETE "<<target<<" inventory="<<keys.size()<<" jointCandidates="<<candidates<<std::endl;
}
} // namespace MtCFExtract

void ExtractMtCF(const char *inputRoot,const char *outputDir,int centMin,int centMax) {
  try {MtCFExtract::Run(inputRoot,outputDir,centMin,centMax);}
  catch(const std::exception &e) {std::cerr<<"ExtractMtCF: "<<e.what()<<std::endl;gSystem->Exit(2);}
}
