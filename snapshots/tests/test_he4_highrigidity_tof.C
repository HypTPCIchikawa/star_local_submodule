// ROOT5 ACLiC fixture. Cut values below are test data matching the YAML profile.
int test_he4_highrigidity_tof(const char* legacyconf,const char* extendedconf);
#ifndef __CINT__
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "StMaker/common/NuclearIdDeDxVsMom.h"
#include "FemtoCandidate.h"
#include "HistManager.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "yaml-cpp/yaml.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
unsigned He4Checks=0;
void He4Require(bool ok,const char* label) {
  ++He4Checks;
  if(!ok) throw std::runtime_error(label);
}
int LegacyTpcOracle(double rigidity,const double pulls[4]) {
  if(!(rigidity<2.5)) return -1;
  int best=-1;double smallest=999.;
  for(int i=0;i<4;++i) {
    const double pull=std::fabs(pulls[i]);
    if(pull<2.0 && pull<smallest) { best=i;smallest=pull; }
  }
  return best;
}
void CheckLegacyProfile(FemtoLambdaLegacy& helper) {
  const double inputs[][4]={{.1,.1,999.,999.},{999.,999.,999.,0.},
    {999.,999.,.1,.2},{999.,999.,-.2,.1},{2.,2.,2.,2.},
    {-1.999,999.,999.,999.},{999.,999.,999.,1.999}};
  const double rigidity[]={.5,1.9999,2.,2.4999,2.5,3.};
  const double mass2[]={-999.,1.95,3.5,7.8};
  for(unsigned p=0;p<sizeof(inputs)/sizeof(inputs[0]);++p)
    for(unsigned r=0;r<sizeof(rigidity)/sizeof(rigidity[0]);++r)
      for(unsigned m=0;m<sizeof(mass2)/sizeof(mass2[0]);++m)
        for(int tof=0;tof<2;++tof)
          He4Require(helper.SelectNuclearType(rigidity[r],inputs[p],tof!=0,mass2[m])==
                     LegacyTpcOracle(rigidity[r],inputs[p]),"absent block preserves legacy finite-input selection");
}
void CheckExtendedProfile(FemtoLambdaLegacy& helper) {
  double pulls[4]={999.,999.,999.,.1};
  const double nan=std::numeric_limits<double>::quiet_NaN();
  const double inf=std::numeric_limits<double>::infinity();
  He4Require(helper.SelectNuclearType(std::nextafter(2.,0.),pulls,false,-999.)==3,
             "He4 immediately below threshold needs no TOF");
  He4Require(helper.SelectNuclearType(2.,pulls,false,3.5)==-1,"exact threshold requires TOF");
  He4Require(helper.SelectNuclearType(std::nextafter(2.,inf),pulls,false,3.5)==-1,
             "He4 immediately above threshold requires TOF");
  He4Require(helper.SelectNuclearType(2.,pulls,true,3.5)==3,"exact threshold accepts valid in-window TOF");
  He4Require(helper.SelectNuclearType(2.5,pulls,true,3.5)==3,"He4 bypasses exact old cap");
  He4Require(helper.SelectNuclearType(3.,pulls,true,3.5)==3,"He4 extended above old cap");
  const double badM2[]={2.9,3.,4.,4.1,nan,inf,-inf};
  for(unsigned i=0;i<sizeof(badM2)/sizeof(badM2[0]);++i)
    He4Require(helper.SelectNuclearType(3.,pulls,true,badM2[i])==-1,"strict TOF window rejects boundary/outside/nonfinite");
  He4Require(helper.SelectNuclearType(3.,pulls,true,std::nextafter(3.,4.))==3,
             "TOF immediately inside lower boundary passes");
  He4Require(helper.SelectNuclearType(3.,pulls,true,std::nextafter(4.,3.))==3,
             "TOF immediately inside upper boundary passes");
  const double badR[]={nan,inf,-inf,-1.};
  for(unsigned i=0;i<sizeof(badR)/sizeof(badR[0]);++i)
    He4Require(helper.SelectNuclearType(badR[i],pulls,true,3.5)==-1,"nonfinite/negative rigidity rejected");
  pulls[3]=nan;
  He4Require(helper.SelectNuclearType(3.,pulls,true,3.5)==-1,"nonfinite He4 pull never passes");
  for(int i=0;i<4;++i)pulls[i]=nan;
  He4Require(helper.SelectNuclearType(1.,pulls,false,-999.)==-1,"all nonfinite hypotheses rejected");
  for(int species=0;species<3;++species) {
    for(int i=0;i<4;++i)pulls[i]=999.;
    pulls[species]=.1;pulls[3]=.2;
    He4Require(helper.SelectNuclearType(2.4,pulls,true,3.5)==species,
               "TOF window cannot promote He4 over a better TPC hypothesis");
    He4Require(helper.SelectNuclearType(2.5,pulls,true,3.5)==-1,"other species keep exact legacy rigidity cap");
    He4Require(helper.SelectNuclearType(3.,pulls,true,3.5)==-1,"other species cannot promote He4 above cap");
  }
  for(int i=0;i<4;++i)pulls[i]=999.;
  pulls[2]=pulls[3]=.1;
  He4Require(helper.SelectNuclearType(2.4,pulls,false,-999.)==2,"He3-He4 TPC tie retains lower-index species");
  pulls[2]=999.;pulls[3]=2.;
  He4Require(helper.SelectNuclearType(3.,pulls,true,3.5)==-1,"TPC pull boundary remains strict");
}

struct He4PicoFixture {
  TClonesArray events,tracks,tofs;
  TClonesArray* arrays[StPicoArrays::NAllPicoArrays];
  StPicoDst dst;
  He4PicoFixture():events("StPicoEvent"),tracks("StPicoTrack"),tofs("StPicoBTofPidTraits") {
    for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)arrays[i]=0;
    arrays[StPicoArrays::Event]=&events;arrays[StPicoArrays::Track]=&tracks;
    arrays[StPicoArrays::BTofPidTraits]=&tofs;
    new(events[0])StPicoEvent();StPicoDst::set(arrays);
  }
  ~He4PicoFixture() { StPicoDst::unset(); }
  StPicoTrack* Add(double rigidity,bool tof,double mass2,int species=3) {
    const int index=tracks.GetEntriesFast();
    StPicoTrack* t=new(tracks[index])StPicoTrack();
    t->setId(index+100);t->setNHitsFit(25);t->setNHitsDedx(20);
    t->setPrimaryMomentum(rigidity,0.,0.);t->setGlobalMomentum(rigidity,0.,0.);
    t->setDedx(1.e-6*NuclearIdDeDxVsMom::GetMean(
      static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(species),t->pMom().Mag()));
    if(tof) {
      const int ti=tofs.GetEntriesFast();
      StPicoBTofPidTraits* b=new(tofs[ti])StPicoBTofPidTraits();
      b->setTrackIndex(index);b->setBTofMatchFlag(1);
      b->setBeta(static_cast<float>(1./std::sqrt(1.+mass2/(rigidity*rigidity))));
      t->setBTofPidTraitsIndex(ti);
    }
    return t;
  }
};
void CheckCollectedTracks(FemtoLambdaLegacy& helper,bool extended) {
  He4PicoFixture fixture;
  fixture.Add(1.99,false,-999.);                 // 0: TPC-only low rigidity.
  fixture.Add(2.,false,-999.);                   // 1: exact TOF threshold.
  fixture.Add(2.,true,3.5);                      // 2: valid compressed beta.
  fixture.Add(2.5,true,3.5);                     // 3: exact old rigidity cap.
  fixture.Add(3.,true,3.5);                      // 4: extended population.
  fixture.Add(3.,true,1.97);                     // 5: He3-like TOF m2.
  StPicoTrack* zero=fixture.Add(3.,true,3.5);     // 6: unusable beta.
  fixture.dst.btofPidTraits(zero->bTofPidTraitsIndex())->setBeta(0.f);
  fixture.Add(3.,false,-999.)->setBTofPidTraitsIndex(30000); // 7: missing TOF traits.
  fixture.Add(5.1,true,3.5);                     // 8: outside calibrated He4 range.
  fixture.Add(3.,true,3.5)->setNHitsDedx(14);     // 9: track-quality failure.
  fixture.Add(3.,true,3.5)->setDedx(std::numeric_limits<float>::quiet_NaN()); // 10
  StPicoTrack* nonfinite=fixture.Add(3.,true,3.5); // 11: invalid primary momentum.
  nonfinite->setPrimaryMomentum(std::numeric_limits<double>::quiet_NaN(),0.,0.);
  fixture.Add(2.4,false,-999.,2);                // 12: retained He3, no TOF promotion.
  fixture.Add(2.5,true,3.5,2);                   // 13: He3 still capped.
  FemtoCandidateStore store;
  helper.CollectNuclei(&fixture.dst,7,8,store);
  const int extendedIndices[]={0,2,3,4},legacyIndices[]={0,1,2};
  const int* expected=extended?extendedIndices:legacyIndices;
  const unsigned count=extended?4:3;
  He4Require(store["he4"].size()==count,"real CollectNuclei has exact expected He4 population");
  He4Require(store["he3"].size()==1 && store["he3"][0].trk.trackIndex==12,
             "real CollectNuclei preserves He3 population and cap");
  for(unsigned i=0;i<count;++i) {
    const FemtoCandidate& c=store["he4"][i];
    He4Require(c.trk.trackIndex==expected[i],"accepted He4 track identity");
    const StPicoTrack* track=fixture.dst.track(c.trk.trackIndex);
    const double rigidity=track->pMom().Mag();
    He4Require(std::fabs(c.P4().P()-2.*rigidity)<1.e-6,"physical He4 momentum doubles rigidity exactly once");
    He4Require(c.charge==2,"physical He4 charge is two");
    if(track->isTofTrack()) {
      const double beta=fixture.dst.btofPidTraits(track->bTofPidTraitsIndex())->btofBeta();
      const double expectedM2=rigidity*rigidity*(1./(beta*beta)-1.);
      He4Require(c.trk.tofMatch && std::fabs(c.trk.mass2-expectedM2)<1.e-5,
                 "TOF stores rigidity m2, not the physical Z2 mass squared");
      He4Require(c.trk.mass2>3. && c.trk.mass2<4.,"compressed beta fixture remains within He4 TOF window");
    }
  }
  TH2* acceptance=dynamic_cast<TH2*>(helper.RootHistograms()->Get("hNucleus_PtVsYLab_he4"));
  He4Require(acceptance && acceptance->GetEntries()==count,"one acceptance entry per selected He4 track");
  if(extended) {
    TH1* stages=helper.NuclearHistograms()->Get("hHe4PidStage");
    const int expectedStages[]={8,1,7,4,3,4};
    He4Require(stages!=0,"He4 cut-stage histogram exists");
    for(int stage=0;stage<6;++stage)
      He4Require(stages->GetBinContent(stages->FindBin(stage))==expectedStages[stage],
                 "He4 cut-stage counts follow TPC identity and TOF gates");
    TH1* preTof=helper.NuclearHistograms()->Get("hM2P_4He_TPC");
    TH1* postTof=helper.NuclearHistograms()->Get("hM2P_4He_selected");
    TH1* selected=helper.NuclearHistograms()->Get("hDedxP_4He_selected");
    TH1* control=helper.NuclearHistograms()->Get("hDedxP_4He_TOFControl");
    TH1* pulls=helper.NuclearHistograms()->Get("hNSigmaHe4VsP_TOFControl");
    He4Require(preTof && preTof->GetEntries()==4,"pre-TOF-window QA retains the rejected He3-like mass2 track");
    He4Require(postTof && postTof->GetEntries()==3,"selected TOF QA excludes low-rigidity no-TOF track");
    He4Require(selected && selected->GetEntries()==count,"selected dEdx QA matches final He4 population");
    He4Require(control && pulls && control->GetEntries()==5 && pulls->GetEntries()==5,
               "TOF control sample includes a competing TPC identity without truncating He4 pulls");
  }
}

std::string ModifiedConfig(const char* source,const YAML::Node& replacement) {
  std::string mainPath(source);
  if(mainPath[0]!='/')mainPath=std::string(gSystem->WorkingDirectory())+"/"+mainPath;
  const size_t pos=mainPath.rfind("/config/");
  He4Require(pos!=std::string::npos,"fixture mainconf is beneath config/");
  const std::string configBase=mainPath.substr(0,pos)+"/config/";
  YAML::Node main=YAML::LoadFile(mainPath);
  // Absolute references keep all unmodified physics/hist files tied to the source profile.
  for(YAML::const_iterator i=main.begin();i!=main.end();++i) {
    const std::string key=i->first.as<std::string>(),value=i->second.as<std::string>();
    if(!value.empty() && value[0]!='/')main[key]=configBase+value;
  }
  YAML::Node nuclear=YAML::LoadFile(main["nuclearid"].as<std::string>());
  nuclear["he4HighRigidityTof"]=replacement;
  char temporary[]="/tmp/star-he4-tof-test-XXXXXX";
  const char* directory=mkdtemp(temporary);
  He4Require(directory!=0,"unique runtime config directory");
  const std::string config=std::string(directory)+"/config";
  He4Require(gSystem->mkdir(config.c_str(),kTRUE)==0,"runtime config subdirectory");
  const std::string nuclearPath=config+"/nuclear.yaml",path=config+"/main.yaml";
  std::ofstream nuclearFile(nuclearPath.c_str());nuclearFile<<nuclear; nuclearFile.close();
  He4Require(nuclearFile.good(),"write nuclear runtime fixture");
  main["nuclearid"]=nuclearPath;
  std::ofstream mainFile(path.c_str());mainFile<<main; mainFile.close();
  He4Require(mainFile.good(),"write mainconf runtime fixture");
  return path;
}
YAML::Node BaseBlock() {
  YAML::Node block;block["enabled"]=true;block["thresholdPOverQ"]=2.;
  block["mass2Min"]=3.;block["mass2Max"]=4.;return block;
}
void CheckConfiguration(const char* extendedconf) {
  for(int variant=0;variant<11;++variant) {
    YAML::Node block=BaseBlock();
    if(variant==0)block.remove("thresholdPOverQ");
    if(variant==1)block.remove("mass2Min");
    if(variant==2)block.remove("mass2Max");
    if(variant==3)block["thresholdPOverQ"]=-1.;
    if(variant==4)block["thresholdPOverQ"]=std::numeric_limits<double>::quiet_NaN();
    if(variant==5)block["mass2Max"]=3.;
    if(variant==6)block["mass2Min"]=std::numeric_limits<double>::infinity();
    if(variant==7)block="not-a-map";
    if(variant==8)block["thresholdPOverq"]=2.;
    if(variant==9)block.remove("enabled");
    if(variant==10)block["mass2Min"]=-1.;
    const std::string path=ModifiedConfig(extendedconf,block);
    FemtoLambdaLegacy helper;std::ostringstream diagnostics;
    He4Require(!helper.Init(path.c_str(),diagnostics),"malformed/missing/nonfinite enabled TOF config fails Init");
    He4Require(diagnostics.str().find("Init failed")!=std::string::npos,"invalid config diagnostic reports failure");
  }
  {
    YAML::Node block=BaseBlock();block["enabled"]=false;
    const std::string path=ModifiedConfig(extendedconf,block);
    FemtoLambdaLegacy helper;
    He4Require(helper.Init(path.c_str(),std::cerr),"explicit disabled TOF block Init");
    CheckLegacyProfile(helper);
  }
  {
    YAML::Node block=BaseBlock();block["thresholdPOverQ"]=3.;
    block["mass2Min"]=2.8;block["mass2Max"]=4.2;
    const std::string path=ModifiedConfig(extendedconf,block);
    FemtoLambdaLegacy helper;
    He4Require(helper.Init(path.c_str(),std::cerr),"custom TOF thresholds Init");
    const double pulls[]={999.,999.,999.,.1};
    He4Require(helper.SelectNuclearType(2.8,pulls,false,-999.)==3,"threshold is configurable, not hardcoded at two");
    He4Require(helper.SelectNuclearType(3.,pulls,false,3.5)==-1,"custom exact threshold requires TOF");
    He4Require(helper.SelectNuclearType(3.2,pulls,true,2.9)==3,"custom lower TOF window takes effect");
    He4Require(helper.SelectNuclearType(3.2,pulls,true,4.1)==3,"custom upper TOF window takes effect");
  }
}
}

int test_he4_highrigidity_tof(const char* legacyconf,const char* extendedconf) {
  try {
    {
      FemtoLambdaLegacy legacy;
      He4Require(legacy.Init(legacyconf,std::cerr),"legacy production Init");
      CheckLegacyProfile(legacy);CheckCollectedTracks(legacy,false);
    }
    {
      FemtoLambdaLegacy extended;
      He4Require(extended.Init(extendedconf,std::cerr),"extended production Init");
      CheckExtendedProfile(extended);CheckCollectedTracks(extended,true);
    }
    CheckConfiguration(extendedconf);
    std::cout<<"PASS He4 high-rigidity TOF: "<<He4Checks<<" checks; legacy/new selection, "
             <<"PicoDst TOF, Z=2, invalid configs and configurable thresholds"<<std::endl;
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<"FAIL He4 high-rigidity TOF after "<<He4Checks<<" checks: "<<e.what()<<std::endl;
    return 1;
  }
}
#endif
