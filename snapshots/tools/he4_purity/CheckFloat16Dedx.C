// Empirical ROOT5 Float16_t storage / native dEdx-bin alias diagnostic.
void CheckFloat16Dedx(const char* savedFits,const char* cachedPico,const char* output);
#ifndef __CINT__
#include "TBufferFile.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
#include <iomanip>

namespace {
void AliasRequire(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
}
#endif
void CheckFloat16Dedx(const char* savedFits,const char* cachedPico,const char* output) {
#ifndef __CINT__
  try {
    const TString dir(output),jsonPath=dir+"/float16_alias.json",rootPath=dir+"/float16_alias.root";
    AliasRequire(gSystem->AccessPathName(jsonPath)&&gSystem->AccessPathName(rootPath),"refusing existing alias diagnostics");
    gSystem->mkdir(dir,kTRUE);
    const int count=500000;std::vector<float> original(count),decoded(count);
    for(int i=0;i<count;++i)original[i]=70.f+20.f*(float(i)+.5f)/float(count);
    TBufferFile packed(TBuffer::kWrite);packed.WriteFastArrayFloat16(&original[0],count,0);
    TBufferFile unpack(TBuffer::kRead,packed.Length(),packed.Buffer(),kFALSE);
    unpack.ReadFastArrayFloat16(&decoded[0],count,0);
    std::set<float> values(decoded.begin(),decoded.end());
    double minimumStep=1.e10,maximumStep=0.,previous=0.;bool first=true;
    for(std::set<float>::const_iterator i=values.begin();i!=values.end();++i) {
      if(!first){minimumStep=std::min(minimumStep,double(*i)-previous);maximumStep=std::max(maximumStep,double(*i)-previous);}
      first=false;previous=*i;
    }
    AliasRequire(minimumStep==maximumStep && minimumStep>0.,"uniform 70-90 input decodes to one Float16 lattice");
    TH1D before("uniform_uncompressed",";dE/dx (keV/cm);entries / 0.1",1200,0.,120.);
    TH1D after("uniform_float16",";dE/dx (keV/cm);entries / 0.1",1200,0.,120.);
    TH1D capacity("float16_lattice_capacity",";dE/dx (keV/cm);representable Float16 values / 0.1",1200,0.,120.);
    before.SetDirectory(0);after.SetDirectory(0);capacity.SetDirectory(0);
    for(int i=0;i<count;++i){before.Fill(original[i]);after.Fill(decoded[i]);}
    for(double value=64.;value<128.;value+=minimumStep)capacity.Fill(value);
    std::map<int,double> uniformSum,actualRatioSum;std::map<int,int> uniformBins,actualBins;
    for(int b=before.FindBin(70.2);b<=before.FindBin(89.8);++b){int n=int(capacity.GetBinContent(b));uniformSum[n]+=after.GetBinContent(b);++uniformBins[n];}
    TFile source(savedFits,"READ");AliasRequire(!source.IsZombie(),"saved fit ROOT readable");
    TH1D* raw=dynamic_cast<TH1D*>(source.Get("raw_010"));AliasRequire(raw!=0,"saved R=.505 raw projection exists");
    double rippleBefore=0.,rippleAfter=0.;int comparisons=0;
    for(int b=raw->FindBin(70.2);b<=raw->FindBin(89.8);++b) {
      const int n=int(capacity.GetBinContent(b));double local=0.,correctedLocal=0.;
      for(int k=-2;k<=2;++k){local+=raw->GetBinContent(b+k);correctedLocal+=raw->GetBinContent(b+k)/capacity.GetBinContent(b+k);}
      local/=5.;correctedLocal/=5.;
      if(local<=0.||correctedLocal<=0.)continue;
      const double uncorrected=raw->GetBinContent(b)/local,corrected=raw->GetBinContent(b)/n/correctedLocal;
      actualRatioSum[n]+=uncorrected;++actualBins[n];
      rippleBefore+=(uncorrected-1.)*(uncorrected-1.);rippleAfter+=(corrected-1.)*(corrected-1.);++comparisons;
    }
    const char* libs[]={"StarRoot","St_base","StarClassLibrary","StPicoEvent"};
    for(int i=0;i<4;++i)AliasRequire(gSystem->Load(libs[i])>=0,"STAR dictionaries load");
    TFile pico(cachedPico,"READ");AliasRequire(!pico.IsZombie(),"cached PicoDst readable");
    TStreamerInfo* info=dynamic_cast<TStreamerInfo*>(pico.GetStreamerInfoList()->FindObject("StPicoTrack"));
    AliasRequire(info!=0,"cached PicoDst stores StPicoTrack StreamerInfo");
    TStreamerElement* member=dynamic_cast<TStreamerElement*>(info->GetElements()->FindObject("mDedx"));
    AliasRequire(member!=0,"stored track has mDedx metadata");
    TTree* tree=dynamic_cast<TTree*>(pico.Get("PicoDst"));AliasRequire(tree!=0,"cached PicoDst tree exists");
    TLeaf* leaf=tree->GetLeaf("Track.mDedx");if(!leaf)leaf=tree->GetLeaf("mDedx");
    AliasRequire(leaf!=0,"cached track dEdx leaf exists");
    tree->SetBranchStatus("*",0);tree->SetBranchStatus("Track*",1);
    long measurements=0,offGrid=0;std::set<double> cachedValues;
    for(int event=0;event<std::min<Long64_t>(1000,tree->GetEntries());++event) {
      AliasRequire(tree->GetEntry(event)>0,"cached event read succeeds");
      for(int t=0;t<leaf->GetNdata();++t) {
        const double value=leaf->GetValue(t);if(value<70.||value>=90.)continue;
        ++measurements;cachedValues.insert(value);
        if(std::fabs(value/minimumStep-std::floor(value/minimumStep+.5))>1.e-5)++offGrid;
      }
    }
    std::ofstream json(jsonPath.Data());json<<std::setprecision(17);
    json<<"{\"uniformInputCount\":"<<count<<",\"serializedBytes\":"<<packed.Length()<<",\"decodedLatticeStep\":"<<minimumStep
        <<",\"uniqueDecodedValues\":"<<values.size()<<",\"trackMemberType\":\""<<member->GetTypeName()<<"\",\"trackRangeMin\":"<<member->GetXmin()
        <<",\"trackRangeMax\":"<<member->GetXmax()<<",\"trackRangeFactor\":"<<member->GetFactor()<<",\"cachedMeasurements70to90\":"<<measurements
        <<",\"cachedUniqueValues\":"<<cachedValues.size()<<",\"cachedOffDefaultLattice\":"<<offGrid<<",\"capacityGroups\":[";
    bool comma=false;for(std::map<int,double>::const_iterator i=uniformSum.begin();i!=uniformSum.end();++i){if(comma)json<<",";comma=true;json<<"{\"valuesPerBin\":"<<i->first<<",\"meanUniformDecodedCount\":"<<i->second/uniformBins[i->first]<<",\"actualDetrendedMean\":"<<actualRatioSum[i->first]/actualBins[i->first]<<"}";}
    json<<"],\"rawFractionalRippleRMS\":"<<std::sqrt(rippleBefore/comparisons)<<",\"capacityCorrectedRippleRMS\":"<<std::sqrt(rippleAfter/comparisons)<<",\"uniformRebins\":[";
    TFile roots(rootPath,"NEW");before.Write();after.Write();capacity.Write();
    for(int factor=5;factor<=10;factor+=5) {
      TH1D* coarse=static_cast<TH1D*>(after.Clone(Form("uniform_float16_rebin%d",factor)));coarse->SetDirectory(0);coarse->Rebin(factor);
      // Exclude both sampled interval edges: rounding moves half a Float16 cell
      // through the endpoints, independently of the interior bin-capacity effect.
      double lo=1.e10,hi=0.;for(int b=coarse->FindBin(71.);b<=coarse->FindBin(88.);++b){lo=std::min(lo,coarse->GetBinContent(b));hi=std::max(hi,coarse->GetBinContent(b));}
      json<<(factor>5?",":"")<<"{\"binWidth\":"<<coarse->GetBinWidth(1)<<",\"minimumCount\":"<<lo<<",\"maximumCount\":"<<hi<<"}";
      coarse->Write();delete coarse;
    }
    json<<"]}\n";json.close();roots.Close();
    std::cout<<"Float16 step="<<minimumStep<<"; stored type="<<member->GetTypeName()<<" factor="<<member->GetFactor()
             <<"; cached values="<<measurements<<" off-grid="<<offGrid<<"; raw/corrected ripple="
             <<std::sqrt(rippleBefore/comparisons)<<"/"<<std::sqrt(rippleAfter/comparisons)<<std::endl;
    AliasRequire(measurements>0 && offGrid==0,"cached track values agree with the measured default Float16 lattice");
    std::cout<<"PASS Float16 dEdx alias diagnostic"<<std::endl;
  } catch(const std::exception& error){std::cerr<<"Float16 alias diagnostic: "<<error.what()<<std::endl;gSystem->Exit(1);}
#endif
}
