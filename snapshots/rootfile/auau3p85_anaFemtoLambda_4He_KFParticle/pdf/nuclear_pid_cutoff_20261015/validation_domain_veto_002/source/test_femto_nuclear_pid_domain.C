// ROOT5 ACLiC tests. All numerical domains/windows below are artificial test data,
// not recommended physics cuts. The caller supplies the source mainconf explicitly.
int test_femto_nuclear_pid_domain(const char*,const char*);
#ifndef __CINT__
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "FemtoCandidate.h"
#include "HistManager.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "yaml-cpp/yaml.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TNamed.h"
#include "TH1.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
unsigned DomainChecks=0;
const char* const DomainSpecies[]={"d","t","he3","he4"};
const double DomainMinimum[]={.55,.65,.75,.85};
void DomainRequire(bool pass,const char* label) {
  ++DomainChecks;if(!pass)throw std::runtime_error(label);
}
YAML::Node MinimumBlock() {
  YAML::Node out;
  for(int i=0;i<4;++i)out[DomainSpecies[i]]=DomainMinimum[i];
  return out;
}
YAML::Node HighBlock(bool enabled=true,bool extend=true) {
  YAML::Node out;out["enabled"]=enabled;out["thresholdPOverQ"]=1.8;
  out["extendHe4BeyondMaxPOverQ"]=extend;
  for(int i=0;i<4;++i) {
    out["species"][DomainSpecies[i]]["mass2Min"]=double(2*i+1);
    out["species"][DomainSpecies[i]]["mass2Max"]=double(2*i+2);
  }
  return out;
}
YAML::Node EmptyVetoBlock() {
  YAML::Node out;
  for(int i=0;i<4;++i)out[DomainSpecies[i]]=YAML::Node(YAML::NodeType::Sequence);
  return out;
}
void AddVeto(YAML::Node& block,int species,double minimum,double maximum) {
  YAML::Node range;range["min"]=minimum;range["max"]=maximum;
  block[DomainSpecies[species]].push_back(range);
}
YAML::Node VetoBlock() {
  YAML::Node out=EmptyVetoBlock();
  AddVeto(out,0,.95,1.05);AddVeto(out,0,1.75,2.);
  AddVeto(out,1,1.15,1.25);AddVeto(out,2,1.35,1.45);AddVeto(out,3,.95,1.05);
  return out;
}
std::string Fixture(const char* source,const std::string& directory,
                    const YAML::Node& minimum,const YAML::Node& high,
                    bool oldTof=false,bool globalTof=false,
                    const YAML::Node& veto=YAML::Node()) {
  std::string sourcePath(source?source:"");
  DomainRequire(!sourcePath.empty(),"explicit source mainconf required");
  if(sourcePath[0]!='/')sourcePath=std::string(gSystem->WorkingDirectory())+"/"+sourcePath;
  const size_t marker=sourcePath.rfind("/config/");
  DomainRequire(marker!=std::string::npos,"source mainconf lies beneath config/");
  const std::string base=sourcePath.substr(0,marker)+"/config/";
  YAML::Node main=YAML::LoadFile(sourcePath);
  for(YAML::const_iterator i=main.begin();i!=main.end();++i) {
    const std::string key=i->first.as<std::string>(),value=i->second.as<std::string>();
    if(!value.empty() && value[0]!='/')main[key]=base+value;
  }
  YAML::Node nuclear=YAML::LoadFile(main["nuclearid"].as<std::string>());
  nuclear.remove("speciesMinPOverQ");nuclear.remove("highRigidityTof");
  nuclear.remove("speciesPOverQVeto");
  nuclear.remove("he4HighRigidityTof");nuclear["m2_selection"]=globalTof;
  if(nuclear["dedxCalibration"] && nuclear["dedxCalibration"]["calibrationFile"]) {
    const std::string path=nuclear["dedxCalibration"]["calibrationFile"].as<std::string>();
    if(!path.empty() && path[0]!='/')nuclear["dedxCalibration"]["calibrationFile"]=base+path;
  }
  if(minimum.IsDefined() && !minimum.IsNull())nuclear["speciesMinPOverQ"]=minimum;
  if(high.IsDefined() && !high.IsNull())nuclear["highRigidityTof"]=high;
  if(veto.IsDefined() && !veto.IsNull())nuclear["speciesPOverQVeto"]=veto;
  if(oldTof) {
    nuclear["he4HighRigidityTof"]["enabled"]=true;
    nuclear["he4HighRigidityTof"]["thresholdPOverQ"]=1.8;
    nuclear["he4HighRigidityTof"]["mass2Min"]=7.;
    nuclear["he4HighRigidityTof"]["mass2Max"]=8.;
  }
  const std::string config=directory+"/config";
  DomainRequire(gSystem->AccessPathName(directory.c_str()),"fixture destination must be new");
  DomainRequire(gSystem->mkdir(config.c_str(),kTRUE)==0,"create fixture config directory");
  const std::string nuclearPath=config+"/nuclear.yaml",mainPath=config+"/main.yaml";
  std::ofstream cuts(nuclearPath.c_str());cuts<<nuclear;cuts.close();
  DomainRequire(cuts.good(),"write fixture nuclear config");
  main["nuclearid"]=nuclearPath;
  std::ofstream entry(mainPath.c_str());entry<<main;entry.close();
  DomainRequire(entry.good(),"write fixture mainconf");
  return mainPath;
}
void CheckMinimum(FemtoLambdaLegacy& helper) {
  const double infinity=std::numeric_limits<double>::infinity();
  for(int species=0;species<4;++species) {
    double pulls[]={999.,999.,999.,999.};pulls[species]=.1;
    DomainRequire(helper.SelectNuclearType(std::nextafter(DomainMinimum[species],0.),pulls,false,0.)==-1,
                  "winner immediately below own minimum is rejected");
    DomainRequire(helper.SelectNuclearType(DomainMinimum[species],pulls,false,0.)==species,
                  "minimum boundary is inclusive");
    DomainRequire(helper.SelectNuclearType(std::nextafter(DomainMinimum[species],infinity),pulls,false,0.)==species,
                  "winner immediately above own minimum is accepted without low-R TOF");
  }
  const double competing[]={.2,999.,999.,.1};
  DomainRequire(helper.SelectNuclearType(.7,competing,false,0.)==-1,
                "cutoff-rejected He4 winner never promotes an otherwise allowed deuteron");
  const double tied[]={999.,999.,.1,.1};
  DomainRequire(helper.SelectNuclearType(.8,tied,false,0.)==2,"TPC tie retains lower species index");
}
void CheckHigh(FemtoLambdaLegacy& helper,bool extend) {
  const double infinity=std::numeric_limits<double>::infinity();
  const double nan=std::numeric_limits<double>::quiet_NaN();
  for(int species=0;species<4;++species) {
    double pulls[]={999.,999.,999.,999.};pulls[species]=.1;
    const double lo=2*species+1,hi=lo+1,inside=(lo+hi)/2;
    DomainRequire(helper.SelectNuclearType(std::nextafter(1.8,0.),pulls,false,nan)==species,
                  "all species remain TPC-only immediately below high-R threshold");
    DomainRequire(helper.SelectNuclearType(1.8,pulls,false,inside)==-1,"exact threshold requires TOF for every species");
    DomainRequire(helper.SelectNuclearType(1.8,pulls,true,inside)==species,"valid TOF selects every species at threshold");
    const double rejected[]={lo,hi,lo-.01,hi+.01,nan,infinity,-infinity};
    for(unsigned i=0;i<sizeof(rejected)/sizeof(rejected[0]);++i)
      DomainRequire(helper.SelectNuclearType(2.1,pulls,true,rejected[i])==-1,"strict TOF boundaries/nonfinite values rejected");
    DomainRequire(helper.SelectNuclearType(2.1,pulls,true,std::nextafter(lo,hi))==species,"TOF just inside lower boundary accepted");
    DomainRequire(helper.SelectNuclearType(2.1,pulls,true,std::nextafter(hi,lo))==species,"TOF just inside upper boundary accepted");
    DomainRequire(helper.SelectNuclearType(2.5,pulls,true,inside)==(species==3 && extend?3:-1),
                  "exact max rigidity is bypassed only by explicitly extended He4");
    DomainRequire(helper.SelectNuclearType(2.7,pulls,true,inside)==(species==3 && extend?3:-1),
                  "above max rigidity only explicitly extended He4 survives");
    if(species<3)pulls[3]=.2;else pulls[0]=.2;
    const double alternate=species<3?7.5:1.5;
    DomainRequire(helper.SelectNuclearType(2.1,pulls,true,alternate)==-1,
                  "TOF-rejected winner never promotes a hypothesis with matching TOF");
  }
  const double pulls[]={.1,999.,999.,999.};
  const double badR[]={0.,-1.,nan,infinity,-infinity};
  for(unsigned i=0;i<sizeof(badR)/sizeof(badR[0]);++i)
    DomainRequire(helper.SelectNuclearType(badR[i],pulls,true,1.5)==-1,"invalid rigidity fails closed");
  DomainRequire(!helper.PassNuclearHighRigidityTof(-1,2.,true,1.5),"invalid species index fails closed");
  DomainRequire(helper.PassHe4HighRigidityTof(2.,true,7.5)==helper.PassNuclearHighRigidityTof(3,2.,true,7.5),
                "He4 QA helper uses same all-species TOF policy");
}
void CheckLegacy(FemtoLambdaLegacy& helper) {
  const double rs[]={.3,.7,1.8,2.,2.4999,2.5,2.7};
  const double ms[]={-999.,1.5,7.5};
  const double inputs[][4]={{.1,.1,999.,999.},{999.,999.,.2,.1},{2.,2.,2.,2.},
                             {-1.99,999.,999.,999.},{999.,999.,999.,.1}};
  for(unsigned r=0;r<sizeof(rs)/sizeof(rs[0]);++r)
    for(unsigned p=0;p<sizeof(inputs)/sizeof(inputs[0]);++p)
      for(unsigned m=0;m<sizeof(ms)/sizeof(ms[0]);++m)for(int tof=0;tof<2;++tof) {
        int expected=-1;double closest=999.;
        if(rs[r]<2.5)for(int s=0;s<4;++s)
          if(std::fabs(inputs[p][s])<2. && std::fabs(inputs[p][s])<closest) {
            expected=s;closest=std::fabs(inputs[p][s]);
          }
        DomainRequire(helper.SelectNuclearType(rs[r],inputs[p],tof!=0,ms[m])==expected,
                      "absent/disabled blocks preserve legacy finite-input selection");
      }
}
void CheckVeto(FemtoLambdaLegacy& helper,const std::string& output) {
  const double infinity=std::numeric_limits<double>::infinity();
  const double minimum[]={.95,1.15,1.35,.95},maximum[]={1.05,1.25,1.45,1.05};
  for(int species=0;species<4;++species) {
    double pulls[]={999.,999.,999.,999.};pulls[species]=.1;
    const double mass2=2*species+1.5;
    DomainRequire(helper.SelectNuclearType(std::nextafter(minimum[species],0.),pulls,true,mass2)==species,
                  "rigidity immediately below veto minimum is allowed");
    DomainRequire(helper.SelectNuclearType(minimum[species],pulls,true,mass2)==-1,"veto lower boundary is inclusive");
    DomainRequire(helper.SelectNuclearType(std::nextafter(minimum[species],infinity),pulls,true,mass2)==-1,
                  "rigidity immediately inside veto lower boundary is rejected");
    DomainRequire(helper.SelectNuclearType(std::nextafter(maximum[species],0.),pulls,true,mass2)==-1,
                  "rigidity immediately inside veto upper boundary is rejected");
    DomainRequire(helper.SelectNuclearType(maximum[species],pulls,true,mass2)==species,"veto upper boundary is exclusive");
    DomainRequire(helper.SelectNuclearType(std::nextafter(maximum[species],infinity),pulls,true,mass2)==species,
                  "rigidity immediately above veto maximum is allowed");
  }
  const double deuteron[]={.1,999.,999.,999.};
  DomainRequire(helper.SelectNuclearType(1.1,deuteron,false,0.)==0,"gap between multiple veto intervals is allowed");
  DomainRequire(helper.SelectNuclearType(1.75,deuteron,false,0.)==-1,"second veto minimum is inclusive");
  DomainRequire(helper.SelectNuclearType(std::nextafter(2.,0.),deuteron,true,1.5)==-1,
                "unassessed region remains vetoed immediately below R=2");
  DomainRequire(helper.SelectNuclearType(2.,deuteron,false,1.5)==-1,"R=2 resumes TOF requirement rather than TPC-only acceptance");
  DomainRequire(helper.SelectNuclearType(2.,deuteron,true,1.5)==0,"R=2 is outside veto and accepted with valid TOF");
  DomainRequire(helper.SelectNuclearType(std::nextafter(2.,infinity),deuteron,true,1.5)==0,
                "rigidity immediately above R=2 is accepted with valid TOF");
  const double alternateHe4[]={.1,999.,999.,.2};
  DomainRequire(helper.SelectNuclearType(1.8,alternateHe4,true,7.5)==-1,
                "vetoed deuteron identity never promotes another species");
  const double alternateTriton[]={999.,.2,999.,.1};
  DomainRequire(helper.SelectNuclearType(1.,alternateTriton,false,0.)==-1,
                "vetoed He4 identity never promotes an allowed triton");
  {TFile file(output.c_str(),"CREATE");DomainRequire(!file.IsZombie() && helper.Write(&file),"write veto-domain fixture ROOT");}
  {TFile file(output.c_str(),"READ");TNamed* domain=dynamic_cast<TNamed*>(file.Get("FemtoNuclearPidDomain"));
   const std::string text=domain?domain->GetTitle():"";
   DomainRequire(domain && text.find("speciesPOverQVeto.d: [[")!=std::string::npos &&
                 text.find("[1.75, 2)")!=std::string::npos && text.find("min <= R < max")!=std::string::npos &&
                 text.find("not a statement of measured low purity")!=std::string::npos,
                 "veto intervals, half-open boundaries and non-physics reason are persisted");}
}
struct DomainPicoFixture {
  TClonesArray events,tracks,tofs;
  TClonesArray* arrays[StPicoArrays::NAllPicoArrays];StPicoDst dst;
  DomainPicoFixture():events("StPicoEvent"),tracks("StPicoTrack"),tofs("StPicoBTofPidTraits") {
    for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)arrays[i]=0;
    arrays[StPicoArrays::Event]=&events;arrays[StPicoArrays::Track]=&tracks;
    arrays[StPicoArrays::BTofPidTraits]=&tofs;
    new(events[0])StPicoEvent();StPicoDst::set(arrays);
  }
  ~DomainPicoFixture(){StPicoDst::unset();}
  void Add(FemtoLambdaLegacy& helper,double rigidity,bool tof) {
    const int index=tracks.GetEntriesFast();StPicoTrack* t=new(tracks[index])StPicoTrack();
    t->setId(index+20);t->setNHitsFit(25);t->setNHitsDedx(20);
    t->setPrimaryMomentum(rigidity,0.,0.);t->setGlobalMomentum(rigidity,0.,0.);
    // Recover the configured He4 center without assuming a particular calibration.
    const double p=t->pMom().Mag();
    const double pull1=helper.NuclearDedxPull(3,p,1.),pull2=helper.NuclearDedxPull(3,p,2.);
    DomainRequire(std::isfinite(pull1) && pull2>pull1,"valid fixture He4 calibration");
    const double mean=1.-pull1/(pull2-pull1);t->setDedx(mean*1.e-6);
    if(tof) {
      const int ti=tofs.GetEntriesFast();StPicoBTofPidTraits* b=new(tofs[ti])StPicoBTofPidTraits();
      b->setTrackIndex(index);b->setBTofMatchFlag(1);
      b->setBeta(static_cast<float>(1./std::sqrt(1.+7.5/(rigidity*rigidity))));
      t->setBTofPidTraitsIndex(ti);
    }
  }
};
void CheckCollection(FemtoLambdaLegacy& helper,const std::string& output,bool veto=false) {
  DomainPicoFixture fixture;
  fixture.Add(helper,.7,false);fixture.Add(helper,1.,false);fixture.Add(helper,2.2,false);
  fixture.Add(helper,2.2,true);fixture.Add(helper,2.6,true);
  FemtoCandidateStore store;helper.CollectNuclei(&fixture.dst,0,8,store);
  const unsigned selected=veto?2:3;
  DomainRequire(store["he4"].size()==selected,"CF candidates obey minimum/veto/TOF domain gates and explicit He4 extension");
  for(size_t i=0;i<store["he4"].size();++i) {
    const FemtoCandidate& c=store["he4"][i];
    DomainRequire(std::fabs(c.P4().P()-2.*fixture.dst.track(c.trk.trackIndex)->pMom().Mag())<1.e-6,
                  "He4 physical momentum receives Z=2 exactly once");
  }
  DomainRequire(helper.NuclearHistograms()->Get("hDedxP")->GetEntries()==5,"raw QA remains uncut by final domain gates");
  DomainRequire(helper.NuclearHistograms()->Get("hDedxP_4He_selected")->GetEntries()==selected,"selected He4 QA matches final CF candidates");
  DomainRequire(helper.RootHistograms()->Get("hNucleus_PtVsYLab_he4")->GetEntries()==selected,"acceptance matches final CF candidates");
  TH1* stages=helper.NuclearHistograms()->Get("hHe4PidStage");
  const int expected[]={5,2,3,2,2,int(selected)};
  for(int i=0;i<6;++i)DomainRequire(stages->GetBinContent(stages->FindBin(i))==expected[i],"pre-domain and final He4 stage semantics retained");
  {TFile file(output.c_str(),"CREATE");DomainRequire(!file.IsZombie() && helper.Write(&file),"write new-domain fixture ROOT");}
  {TFile file(output.c_str(),"READ");TNamed* domain=dynamic_cast<TNamed*>(file.Get("FemtoNuclearPidDomain"));
   DomainRequire(domain && std::string(domain->GetTitle()).find("never reassign")!=std::string::npos,"domain and no-reassignment provenance persisted");}
}
}

int test_femto_nuclear_pid_domain(const char* source,const char* outputDirectory) {
  try {
    DomainRequire(source && *source && outputDirectory && *outputDirectory,"explicit source and new output directory required");
    std::string output(outputDirectory);
    // Fixture mainconfs contain explicit nuclear-file references. Make those
    // references absolute rather than rebasing a relative output path twice.
    if(output[0]!='/')output=std::string(gSystem->WorkingDirectory())+"/"+output;
    DomainRequire(gSystem->AccessPathName(output.c_str()),"test output directory must not already exist");
    DomainRequire(gSystem->mkdir(output.c_str(),kTRUE)==0,"create new test output directory");
    const YAML::Node absent;
    {const std::string path=Fixture(source,output+"/legacy",absent,absent);
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"legacy fixture Init");CheckLegacy(helper);}
    {const std::string path=Fixture(source,output+"/disabled",absent,HighBlock(false));
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"explicitly disabled block Init");CheckLegacy(helper);}
    {const std::string path=Fixture(source,output+"/minimum",MinimumBlock(),absent);
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"minimum-only fixture Init");CheckMinimum(helper);}
    {const std::string path=Fixture(source,output+"/all_tof",MinimumBlock(),HighBlock());
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"all-species TOF fixture Init");
     CheckMinimum(helper);CheckHigh(helper,true);CheckCollection(helper,output+"/selected.root");}
    {const std::string path=Fixture(source,output+"/bounded",MinimumBlock(),HighBlock(true,false));
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"non-extended He4 fixture Init");CheckHigh(helper,false);}
    {const std::string path=Fixture(source,output+"/empty_veto",absent,absent,false,false,EmptyVetoBlock());
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"explicit empty veto lists Init");CheckLegacy(helper);}
    {YAML::Node high=HighBlock();high["thresholdPOverQ"]=2.;
     const std::string path=Fixture(source,output+"/veto",MinimumBlock(),high,false,false,VetoBlock());
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"species rigidity veto Init");
     CheckVeto(helper,output+"/veto_metadata.root");CheckCollection(helper,output+"/veto_selected.root",true);}
    {YAML::Node veto=EmptyVetoBlock();AddVeto(veto,0,.95,1.05);AddVeto(veto,0,1.05,1.15);
     const std::string path=Fixture(source,output+"/adjacent_veto",absent,absent,false,false,veto);
     FemtoLambdaLegacy helper;DomainRequire(helper.Init(path.c_str(),std::cerr),"adjacent half-open veto intervals are valid");
     const double pulls[]={.1,999.,999.,999.};
     DomainRequire(helper.SelectNuclearType(1.05,pulls,false,0.)==-1,"shared endpoint belongs to the following veto interval");
     DomainRequire(helper.SelectNuclearType(1.15,pulls,false,0.)==0,"final adjacent upper endpoint is not vetoed");}
    for(int variant=0;variant<25;++variant) {
      YAML::Node minimum=MinimumBlock(),high=HighBlock();
      if(variant==0)minimum.remove("d");
      if(variant==1)minimum["deuteron"]=.5;
      if(variant==2)minimum["he4"]=-.1;
      if(variant==3)minimum["t"]=std::numeric_limits<double>::infinity();
      if(variant==4)minimum["he3"]=std::numeric_limits<double>::quiet_NaN();
      if(variant==5)minimum="bad map";
      if(variant==6)high.remove("enabled");
      if(variant==7)high.remove("thresholdPOverQ");
      if(variant==8)high.remove("extendHe4BeyondMaxPOverQ");
      if(variant==9)high.remove("species");
      if(variant==10)high["thresholdPOverQ"]=0.;
      if(variant==11)high["thresholdPOverQ"]=std::numeric_limits<double>::quiet_NaN();
      if(variant==12)high["thresholdPOverq"]=1.8;
      if(variant==13)high["species"].remove("he4");
      if(variant==14)high["species"]["d"].remove("mass2Min");
      if(variant==15)high["species"]["t"]["mass2Max"]=3.;
      if(variant==16)high["species"]["he3"]["mass2Min"]=-1.;
      if(variant==17)high["species"]["he4"]["mass2Max"]=std::numeric_limits<double>::infinity();
      if(variant==18)high["species"]["he4"]["mass2Min"]=std::numeric_limits<double>::quiet_NaN();
      if(variant==19)high["species"]["he4"]["mass2min"]=7.;
      if(variant==20)high["species"]["deuteron"]=high["species"]["d"];
      if(variant==21)high="bad map";
      if(variant==22)high["species"]["he3"]="bad map";
      std::ostringstream suffix;suffix<<variant;
      const std::string path=Fixture(source,output+"/invalid_"+suffix.str(),minimum,high,variant==23,variant==24);
      FemtoLambdaLegacy helper;std::ostringstream diagnostics;
      DomainRequire(!helper.Init(path.c_str(),diagnostics),"malformed/missing/nonfinite/conflicting domain config fails closed");
      DomainRequire(diagnostics.str().find("Init failed")!=std::string::npos,"invalid domain config has explicit failure diagnostic");
    }
    for(int variant=0;variant<29;++variant) {
      YAML::Node veto=VetoBlock();
      if(variant<4)veto.remove(DomainSpecies[variant]);
      if(variant==4)veto["deuteron"]=YAML::Node(YAML::NodeType::Sequence);
      if(variant==5)veto="bad map";
      if(variant==6)veto["d"]=YAML::Node(YAML::NodeType::Map);
      if(variant==7)veto["d"][0]="bad interval map";
      if(variant==8)veto["d"][0].remove("min");
      if(variant==9)veto["d"][0].remove("max");
      if(variant==10)veto["d"][0]["minimum"]=.95;
      if(variant==11)veto["d"][0]["min"]=0.;
      if(variant==12)veto["d"][0]["min"]=-.1;
      if(variant==13)veto["d"][0]["max"]=.95;
      if(variant==14)veto["d"][0]["max"]=.9;
      if(variant==15)veto["d"][0]["min"]=std::numeric_limits<double>::quiet_NaN();
      if(variant==16)veto["d"][0]["max"]=std::numeric_limits<double>::infinity();
      if(variant==17)veto["d"][0]["min"]=std::numeric_limits<double>::infinity();
      if(variant==18)veto["d"][0]["max"]=std::numeric_limits<double>::quiet_NaN();
      if(variant==19)veto["d"][1]["min"]=1.;
      if(variant==20){veto["d"][1]["min"]=.6;veto["d"][1]["max"]=.7;}
      if(variant==21)veto=YAML::Load("d: []\nd: []\nt: []\nhe3: []\nhe4: []\n");
      if(variant==22)veto["d"]=YAML::Load("- {min: 0.95, min: 0.96, max: 1.05}\n");
      if(variant==23)veto["d"]=YAML::Load("- {min: 0.95, max: 1.05, max: 1.06}\n");
      if(variant==24)veto["d"]=YAML::Node(YAML::NodeType::Null);
      if(variant==25){veto["d"][1]["min"]=.95;veto["d"][1]["max"]=1.05;}
      if(variant==26)veto["d"][0]["max"]=-1.;
      if(variant==27)veto=YAML::Node(YAML::NodeType::Sequence);
      std::ostringstream suffix;suffix<<variant;
      const std::string path=Fixture(source,output+"/invalid_veto_"+suffix.str(),absent,absent,false,false,veto);
      if(variant==28) {
        const YAML::Node main=YAML::LoadFile(path);
        std::ofstream append(main["nuclearid"].as<std::string>().c_str(),std::ios::app);
        append<<"\nspeciesPOverQVeto: {d: [], t: [], he3: [], he4: []}\n";append.close();
        DomainRequire(append.good(),"write duplicate top-level veto block fixture");
      }
      FemtoLambdaLegacy helper;std::ostringstream diagnostics;
      DomainRequire(!helper.Init(path.c_str(),diagnostics),"invalid/missing/overlapping/unsorted/duplicate veto config fails closed");
      DomainRequire(diagnostics.str().find("Init failed")!=std::string::npos,"invalid veto config has explicit failure diagnostic");
    }
    std::cout<<"PASS nuclear PID domain: "<<DomainChecks<<" checks; legacy parity, winner cutoffs/veto intervals, all-species TOF, no reassignment, collection and provenance"<<std::endl;
    return 0;
  }catch(const std::exception& e) {
    std::cerr<<"FAIL nuclear PID domain after "<<DomainChecks<<" checks: "<<e.what()<<std::endl;return 1;
  }
}
#endif
