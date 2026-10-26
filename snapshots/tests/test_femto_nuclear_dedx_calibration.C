// ROOT5 ACLiC fixture: numerical values here are deliberately artificial test data.
int test_femto_nuclear_dedx_calibration(const char*,const char*,const char*);
#ifndef __CINT__
#include "StMaker/common/FemtoNuclearDedxCalibration.h"
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
#include "TFile.h"
#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
unsigned CalibrationChecks=0;
std::string TestDirectory;
void RequireCalibration(bool condition,const char* message) {
  ++CalibrationChecks;if(!condition)throw std::runtime_error(message);
}
void Near(double got,double expected,const char* message) {
  RequireCalibration(std::fabs(got-expected)<=1.e-11*std::max(1.,std::fabs(expected)),message);
}
std::string SaveYaml(const YAML::Node& node,const std::string& name) {
  const std::string path=TestDirectory+"/"+name+".yaml";
  std::ofstream file(path.c_str());file<<node;file.close();
  RequireCalibration(file.good(),"write isolated YAML fixture");return path;
}
YAML::Node Copy(const YAML::Node& node) { return YAML::Load(YAML::Dump(node)); }
YAML::Node Fixture() {
  YAML::Node doc;
  doc["schemaVersion"]=1;doc["profile"]="calibration-unit-test";
  doc["rigidityUnit"]="GeV/c";doc["dedxUnit"]="keV/cm";
  doc["interpolation"]="linear";doc["outsideRangePolicy"]="legacy";
  doc["edgePolicy"]="clampWithinDomain";doc["maxInterpolationGap"]=.015;
  YAML::Node source;
  source["collision"]="auau3p85";source["year"]=2021;source["production"]="P24iy";
  source["sourceRoot"]="synthetic-fixture.root";source["sourceRootFingerprint"]["fixture"]=true;
  const char* paths[]={"fitSummary","fitSource","settings"};
  const char* hashes[]={"fitSummarySha256","fitSourceSha256","settingsSha256"};
  for(int i=0;i<3;++i) { source[paths[i]]="synthetic";source[hashes[i]]=std::string(64,'0'); }
  doc["source"]=source;
  const char* species[]={"he3","he4"};
  const double positions[]={.415,.425,.455};
  for(int s=0;s<2;++s) {
    YAML::Node curve;curve["domain"].push_back(.41);curve["domain"].push_back(.46);
    curve["supportBinHalfWidth"]=.005;
    for(int k=0;k<3;++k) {
      YAML::Node knot;knot["rigidity"]=positions[k];
      knot["mean"]=60.+20.*s-2.*k;knot["sigma"]=2.+s+.5*k;
      knot["provenance"]["fixture"]=true;curve["knots"].push_back(knot);
    }
    doc["species"][species[s]]=curve;
  }
  return doc;
}
void LegacyParity(const FemtoNuclearDedxCalibration& provider) {
  const double r[]={.2,.399,.4,.405,.41,.415,.5,.7,1.3,2.,2.5,3.,5.,5.01,6.};
  for(int s=0;s<4;++s)for(unsigned j=0;j<sizeof(r)/sizeof(r[0]);++j) {
    const NuclearIdDeDxVsMom::SpeciesIndex species=static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s);
    Near(provider.Mean(s,r[j]),NuclearIdDeDxVsMom::GetMean(species,r[j]),"legacy mean unchanged");
    Near(provider.Sigma(s,r[j]),NuclearIdDeDxVsMom::GetSigma(species,r[j]),"legacy width unchanged");
    Near(provider.Pull(s,r[j],35.),NuclearIdDeDxVsMom::GetNSigma(species,r[j],35.),"legacy positive finite pull unchanged");
    RequireCalibration(!provider.UsesCalibration(s,r[j]),"legacy provider not flagged calibrated");
  }
}
void CheckInterpolation(FemtoNuclearDedxCalibration& provider) {
  for(int s=2;s<4;++s) {
    const double offset=20.*(s-2);
    Near(provider.Mean(s,.415),60.+offset,"first knot exact");
    Near(provider.Mean(s,.42),59.+offset,"adjacent mean linearly interpolated");
    Near(provider.Sigma(s,.42),2.25+s-2,"adjacent sigma linearly interpolated");
    Near(provider.Mean(s,.41),60.+offset,"lower supported domain edge clamp");
    Near(provider.Mean(s,.46),56.+offset,"upper supported domain edge clamp");
    Near(provider.Mean(s,.43),58.+offset,"gap left support edge clamp");
    Near(provider.Mean(s,.45),56.+offset,"gap right support edge clamp");
    RequireCalibration(provider.UsesCalibration(s,.43) && provider.UsesCalibration(s,.45),"supported gap half-bins calibrated");
    const double fallback[]={.409,.461,.4300000001,.435,.4499999999};
    for(unsigned i=0;i<sizeof(fallback)/sizeof(fallback[0]);++i) {
      const NuclearIdDeDxVsMom::SpeciesIndex sp=static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s);
      RequireCalibration(!provider.UsesCalibration(s,fallback[i]),"unsupported gap/outside falls back");
      Near(provider.Mean(s,fallback[i]),NuclearIdDeDxVsMom::GetMean(sp,fallback[i]),"declared fallback mean exact");
      Near(provider.Sigma(s,fallback[i]),NuclearIdDeDxVsMom::GetSigma(sp,fallback[i]),"declared fallback sigma exact");
    }
    Near(provider.Pull(s,.42,provider.Mean(s,.42)+1.5*provider.Sigma(s,.42)),1.5,"pull uses mean and sigma together");
  }
  for(int s=0;s<2;++s)for(int i=0;i<10;++i) {
    const double r=.3+.4*i;
    Near(provider.Mean(s,r),NuclearIdDeDxVsMom::GetMean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),r),"d/t never changed");
    RequireCalibration(!provider.UsesCalibration(s,r),"d/t calibration flag remains legacy");
  }
  const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
  RequireCalibration(provider.Pull(3,nan,80.)==999. && provider.Pull(3,inf,80.)==999.,"nonfinite rigidity cannot pass");
  RequireCalibration(provider.Pull(3,.42,nan)==999. && provider.Pull(3,.42,inf)==999.,"nonfinite dEdx cannot pass");
  RequireCalibration(provider.Pull(3,.42,0.)==999. && provider.Pull(3,-1.,80.)==999.,"nonpositive measurement cannot pass");
  RequireCalibration(provider.Pull(-1,.42,80.)==999. && provider.Pull(4,.42,80.)==999.,"invalid species rejected");
}
void CheckBadFiles(const YAML::Node& original) {
  for(int variant=0;variant<30;++variant) {
    YAML::Node doc=Copy(original);
    YAML::Node curve=doc["species"]["he4"];
    if(variant==0)doc["schemaVersion"]=2;
    if(variant==1)doc["unknownPolicy"]=true;
    if(variant==2)doc["dedxUnit"]="GeV";
    if(variant==3)doc["interpolation"]="spline";
    if(variant==4)doc["outsideRangePolicy"]="extrapolate";
    if(variant==5)doc["edgePolicy"]="extrapolate";
    if(variant==6)doc["maxInterpolationGap"]=0.;
    if(variant==7)doc["maxInterpolationGap"]=std::numeric_limits<double>::quiet_NaN();
    if(variant==8)doc["source"]["collision"]="auau13p5";
    if(variant==9)doc["source"]["year"]=2019;
    if(variant==10)doc["source"]["production"]="P23";
    if(variant==11)doc["source"].remove("fitSummarySha256");
    if(variant==12)doc["source"]["fitSourceSha256"]="bad";
    if(variant==13)doc["source"]["sourceRootFingerprint"]="bad";
    if(variant==14)doc["species"].remove("he3");
    if(variant==15)doc["species"]["triton"]=curve;
    if(variant==16)curve["domain"]=.4;
    if(variant==17)curve["domain"][1]=.4;
    if(variant==18)curve["domain"][0]=.40;
    if(variant==19)curve["supportBinHalfWidth"]=0.;
    if(variant==20)curve["supportBinHalfWidth"]=.02;
    if(variant==21)curve["knots"]="bad";
    if(variant==22)curve["knots"][1]["rigidity"]=.415;
    if(variant==23)curve["knots"][1]["mean"]=-1.;
    if(variant==24)curve["knots"][1]["sigma"]=0.;
    if(variant==25)curve["knots"][1]["sigma"]=std::numeric_limits<double>::infinity();
    if(variant==26)curve["knots"][1]["rigidity"]=.47;
    if(variant==27)curve["knots"][1]["provenance"]="bad";
    if(variant==28)curve["knots"][1]["sigm"]="typo";
    if(variant==29)doc["source"]["settings"]="";
    std::ostringstream name;name<<"bad-"<<variant;
    FemtoNuclearDedxCalibration provider;std::ostringstream diagnostic;
    RequireCalibration(!provider.Load(SaveYaml(doc,name.str()),diagnostic),"bad calibration rejected loudly");
    RequireCalibration(diagnostic.str().find("Load failed")!=std::string::npos,"bad calibration diagnostic");
  }
  FemtoNuclearDedxCalibration provider;std::ostringstream diagnostic;
  RequireCalibration(!provider.Load(TestDirectory+"/absent.yaml",diagnostic),"missing calibration rejected");
}
bool ExpectedCustom(const YAML::Node& curve,double r,double maximumGap,double& mean,double& sigma) {
  const double tol=1.e-12,min=curve["domain"][0].as<double>(),max=curve["domain"][1].as<double>();
  if(r<min-tol||r>max+tol)return false;
  const YAML::Node k=curve["knots"];
  if(r<=k[0]["rigidity"].as<double>()) { mean=k[0]["mean"].as<double>();sigma=k[0]["sigma"].as<double>();return true; }
  const unsigned last=k.size()-1;
  if(r>=k[last]["rigidity"].as<double>()) { mean=k[last]["mean"].as<double>();sigma=k[last]["sigma"].as<double>();return true; }
  for(unsigned i=1;i<k.size();++i) {
    const double left=k[i-1]["rigidity"].as<double>(),right=k[i]["rigidity"].as<double>();
    if(r>right)continue;
    unsigned endpoint=0;bool clamp=false;
    if(right-left>maximumGap+tol) {
      const double half=curve["supportBinHalfWidth"].as<double>();
      if(r<=left+half+tol) { endpoint=i-1;clamp=true; }
      else if(r>=right-half-tol) { endpoint=i;clamp=true; }
      else return false;
    }
    if(clamp) { mean=k[endpoint]["mean"].as<double>();sigma=k[endpoint]["sigma"].as<double>();return true; }
    const double f=(r-left)/(right-left);
    mean=(1.-f)*k[i-1]["mean"].as<double>()+f*k[i]["mean"].as<double>();
    sigma=(1.-f)*k[i-1]["sigma"].as<double>()+f*k[i]["sigma"].as<double>();return true;
  }
  return false;
}
void CheckExport(const char* file,FemtoNuclearDedxCalibration& provider) {
  YAML::Node doc=YAML::LoadFile(file);
  std::ofstream audit((TestDirectory+"/runtime_curve_audit.tsv").c_str());
  audit.precision(17);audit<<"species\trigidity\tmean\tsigma\tcustom\n";
  const char* names[]={"he3","he4"};
  for(int s=2;s<4;++s) {
    const YAML::Node curve=doc["species"][names[s-2]];
    RequireCalibration(curve["knots"].size()==205,"chosen205 export preserved");
    for(unsigned i=0;i<curve["knots"].size();++i) {
      const YAML::Node k=curve["knots"][i];const double r=k["rigidity"].as<double>();
      Near(provider.Mean(s,r),k["mean"].as<double>(),"export knot mean exact");
      Near(provider.Sigma(s,r),k["sigma"].as<double>(),"export knot width exact");
      audit<<s<<'\t'<<r<<'\t'<<provider.Mean(s,r)<<'\t'<<provider.Sigma(s,r)<<'\t'<<provider.UsesCalibration(s,r)<<'\n';
    }
    for(unsigned i=0;i<1250;++i) {
      const double r=.395+.0037*i;double mean=0,sigma=0;
      const bool expected=ExpectedCustom(curve,r,doc["maxInterpolationGap"].as<double>(),mean,sigma);
      RequireCalibration(provider.UsesCalibration(s,r)==expected,"export/runtime support domains and gaps agree");
      if(!expected) {
        mean=NuclearIdDeDxVsMom::GetMean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),r);
        sigma=NuclearIdDeDxVsMom::GetSigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),r);
      }
      Near(provider.Mean(s,r),mean,"export/runtime interpolated mean agree");
      Near(provider.Sigma(s,r),sigma,"export/runtime interpolated width agree");
      audit<<s<<'\t'<<r<<'\t'<<provider.Mean(s,r)<<'\t'<<provider.Sigma(s,r)<<'\t'<<provider.UsesCalibration(s,r)<<'\n';
    }
    const double half=curve["supportBinHalfWidth"].as<double>();
    for(unsigned i=0;i<curve["knots"].size();++i) {
      const double r=curve["knots"][i]["rigidity"].as<double>();
      const double offsets[]={-half-1.e-10,-half,-5.e-13,5.e-13,half,half+1.e-10};
      for(unsigned j=0;j<sizeof(offsets)/sizeof(offsets[0]);++j) {
        const double point=r+offsets[j];
        audit<<s<<'\t'<<point<<'\t'<<provider.Mean(s,point)<<'\t'<<provider.Sigma(s,point)<<'\t'<<provider.UsesCalibration(s,point)<<'\n';
      }
    }
  }
  audit.close();RequireCalibration(audit.good(),"write runtime/export cross-language audit");
}
std::string ChangedMain(const char* source,const YAML::Node& block,const char* wrongTag=0) {
  std::string mainPath=source;
  if(mainPath[0]!='/')mainPath=std::string(gSystem->WorkingDirectory())+"/"+mainPath;
  const size_t pos=mainPath.rfind("/config/");
  RequireCalibration(pos!=std::string::npos,"mainconf fixture under config/");
  const std::string base=mainPath.substr(0,pos)+"/config/";
  YAML::Node main=YAML::LoadFile(mainPath);
  for(YAML::const_iterator i=main.begin();i!=main.end();++i) {
    const std::string key=i->first.as<std::string>(),value=i->second.as<std::string>();
    if(!value.empty()&&value[0]!='/')main[key]=base+value;
  }
  YAML::Node nuclear=YAML::LoadFile(main["nuclearid"].as<std::string>());nuclear["dedxCalibration"]=block;
  static unsigned id=0;std::ostringstream prefix;prefix<<"config/main-"<<id++;
  main["nuclearid"]=SaveYaml(nuclear,prefix.str()+"-nuclear");
  if(wrongTag) {
    YAML::Node analysis=YAML::LoadFile(main["analysis"].as<std::string>());
    analysis["starTag"][wrongTag]=std::string(wrongTag)=="year"?"year2019":
      (std::string(wrongTag)=="collision"?"auau13p5":"P23");
    main["analysis"]=SaveYaml(analysis,prefix.str()+"-analysis");
  }
  return SaveYaml(main,prefix.str());
}
void CheckInit(const char* legacyconf,const char* calibratedconf,const char* calibration) {
  FemtoLambdaLegacy old;
  RequireCalibration(old.Init(legacyconf,std::cerr),"legacy profile Init");
  RequireCalibration(old.NuclearDedxCalibrationMetadata().find("mode: legacy")!=std::string::npos,"legacy metadata explicit");
  Near(old.NuclearDedxPull(3,.6,50.),NuclearIdDeDxVsMom::GetNSigma(NuclearIdDeDxVsMom::kHe4,.6,50.),"old helper legacy pull exact");
  YAML::Node block;block["enabled"]=true;block["calibrationFile"]=std::string(calibration);
  for(int i=0;i<7;++i) {
    YAML::Node bad=Copy(block);const char* wrong=0;
    if(i==0)bad="invalid";
    if(i==1)bad["fileTypo"]=true;
    if(i==2)bad.remove("calibrationFile");
    if(i==3)bad["calibrationFile"]="/unavailable-calibration.yaml";
    if(i==4)wrong="collision";
    if(i==5)wrong="year";
    if(i==6)wrong="productionTag";
    const std::string main=ChangedMain(calibratedconf,bad,wrong);
    FemtoLambdaLegacy helper;std::ostringstream diagnostic;
    RequireCalibration(!helper.Init(main.c_str(),diagnostic),"malformed calib config or dataset mismatch fails Init");
  }
  block["enabled"]=false;block.remove("calibrationFile");
  const std::string disabled=ChangedMain(calibratedconf,block);
  FemtoLambdaLegacy helper;
  RequireCalibration(helper.Init(disabled.c_str(),std::cerr),"disabled calibration Init");
  Near(helper.NuclearDedxPull(3,.6,50.),old.NuclearDedxPull(3,.6,50.),"disabled calibration reverts exactly");
  // Keep every unrelated reference valid, so these fixtures specifically test
  // name-versus-bachelor validation, not an incidental missing calibration.
  block["enabled"]=true;block["calibrationFile"]=std::string(calibration);
  for(int variant=0;variant<2;++variant) {
    const std::string initial=ChangedMain(calibratedconf,block);
    YAML::Node main=YAML::LoadFile(initial);
    YAML::Node maker=YAML::LoadFile(main["maker"].as<std::string>());
    if(variant==0) {
      maker["channel_0_partB"]="he3";maker["channel_4_partB"]="he4";
    } else {
      maker["channel_1_partB"]="he3";
    }
    std::ostringstream label;label<<"config/channel-bachelor-mismatch-"<<variant;
    main["maker"]=SaveYaml(maker,label.str()+"-maker");
    const std::string path=SaveYaml(main,label.str());
    FemtoLambdaLegacy rejected;std::ostringstream diagnostic;
    RequireCalibration(!rejected.Init(path.c_str(),diagnostic),"base/signal channel bachelor mismatch rejected");
    RequireCalibration(diagnostic.str().find("partB")!=std::string::npos,
                       "channel mismatch rejection specifically reports partB, not missing files");
  }
}
struct CalibrationPicoFixture {
  TClonesArray events,tracks,tofs;TClonesArray* arrays[StPicoArrays::NAllPicoArrays];StPicoDst dst;
  CalibrationPicoFixture():events("StPicoEvent"),tracks("StPicoTrack"),tofs("StPicoBTofPidTraits") {
    for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)arrays[i]=0;
    arrays[StPicoArrays::Event]=&events;arrays[StPicoArrays::Track]=&tracks;arrays[StPicoArrays::BTofPidTraits]=&tofs;
    new(events[0])StPicoEvent();StPicoDst::set(arrays);
  }
  ~CalibrationPicoFixture() { StPicoDst::unset(); }
  void Add(FemtoNuclearDedxCalibration& provider,int species,double r,bool tof=false) {
    const int index=tracks.GetEntriesFast();StPicoTrack* t=new(tracks[index])StPicoTrack();
    t->setId(index+100);t->setNHitsFit(25);t->setNHitsDedx(20);
    t->setPrimaryMomentum(r,0.,0.);t->setGlobalMomentum(r,0.,0.);
    t->setDedx(1.e-6*provider.Mean(species,t->pMom().Mag()));
    if(tof) {
      const int ti=tofs.GetEntriesFast();StPicoBTofPidTraits* b=new(tofs[ti])StPicoBTofPidTraits();
      b->setTrackIndex(index);b->setBTofMatchFlag(1);b->setBeta(1./std::sqrt(1.+3.5/(r*r)));
      t->setBTofPidTraitsIndex(ti);
    }
  }
};
void CheckCollect(const char* calibratedconf,FemtoNuclearDedxCalibration& provider) {
  FemtoLambdaLegacy helper;
  RequireCalibration(helper.Init(calibratedconf,std::cerr),"production calibrated dual-target Init");
  RequireCalibration(helper.IsPairNucleus(2)&&helper.IsPairNucleus(3),"calibrated He3/He4 dual targets enabled");
  CalibrationPicoFixture fixture;
  fixture.Add(provider,3,.6);fixture.Add(provider,2,.6);
  fixture.Add(provider,0,1.);fixture.Add(provider,1,1.);
  fixture.Add(provider,3,3.,true);fixture.Add(provider,2,3.,true);
  FemtoCandidateStore store;helper.CollectNuclei(&fixture.dst,17,8,store);
  RequireCalibration(store["he4"].size()==2,"calibrated lowR and explicit fallback highR He4 selected");
  RequireCalibration(store["he3"].size()==1,"calibrated He3 accepted but legacy highR cap maintained");
  RequireCalibration(store["deuteron"].size()==1&&store["triton"].size()==1,"d/t selection retained");
  const char* names[]={"deuteron","triton","he3","he4"};
  for(int s=0;s<4;++s)for(unsigned i=0;i<store[names[s]].size();++i) {
    const FemtoCandidate& c=store[names[s]][i];const StPicoTrack* t=fixture.dst.track(c.trk.trackIndex);
    const double r=t->pMom().Mag(),dedx=t->dEdx();
    const double pulls[]={c.trk.nSigmaDeuteron,c.trk.nSigmaTriton,c.trk.nSigmaHe3,c.trk.nSigmaHe4};
    for(int k=0;k<4;++k)Near(pulls[k],static_cast<float>(provider.Pull(k,r,dedx)),"all four stored pulls share calibration provider after declared Float_t storage");
    Near(c.P4().P(),r*(s>=2?2.:1.),"nuclear momentum Z correction exactly once");
  }
  RequireCalibration(helper.RootHistograms()->Get("hNucleus_PtVsYLab_he4")->GetEntries()==2,"He4 acceptance once per selected track");
  RequireCalibration(helper.RootHistograms()->Get("hNucleus_PtVsYLab_he3")->GetEntries()==1,"He3 acceptance once per selected track");
  TH2* control=dynamic_cast<TH2*>(helper.NuclearHistograms()->Get("hNSigmaHe4VsP_TOFControl"));
  RequireCalibration(control&&control->GetEntries()==6,"TOF control remains preTPC population");
  const StPicoTrack* t=fixture.dst.track(0);
  RequireCalibration(control->GetBinContent(control->FindBin(t->pMom().Mag(),provider.Pull(3,t->pMom().Mag(),t->dEdx())))>=1.,
                     "TOF control pull uses new calibration");
  const std::string root=TestDirectory+"/calibration_metadata.root";
  { TFile out(root.c_str(),"CREATE");RequireCalibration(!out.IsZombie()&&helper.Write(&out),"write calibrated metadata ROOT"); }
  { TFile input(root.c_str(),"READ");TNamed* metadata=dynamic_cast<TNamed*>(input.Get("FemtoNuclearDedxCalibration"));
    RequireCalibration(metadata&&std::string(metadata->GetTitle())==helper.NuclearDedxCalibrationMetadata(),"exact consumed calibration YAML saved in ROOT"); }
}
}
int test_femto_nuclear_dedx_calibration(const char* legacyconf,const char* calibratedconf,const char* calibration) {
  try {
    char temporary[]="/tmp/star-femto-calibration-test-XXXXXX";const char* directory=mkdtemp(temporary);
    RequireCalibration(directory!=0,"unique fixture directory");TestDirectory=directory;
    RequireCalibration(gSystem->mkdir((TestDirectory+"/config").c_str(),kTRUE)==0,"fixture config root");
    { FemtoNuclearDedxCalibration legacy;LegacyParity(legacy); }
    const YAML::Node fixture=Fixture();FemtoNuclearDedxCalibration provider;
    RequireCalibration(provider.Load(SaveYaml(fixture,"valid"),std::cerr),"valid fixture Load");
    CheckInterpolation(provider);CheckBadFiles(fixture);
    RequireCalibration(provider.Load(calibration,std::cerr),"production export Load");
    CheckExport(calibration,provider);CheckInit(legacyconf,calibratedconf,calibration);CheckCollect(calibratedconf,provider);
    std::cout<<"PASS helium dE/dx calibration: "<<CalibrationChecks<<" checks; chosen205, interpolation/gaps/edges, legacy fallback, strict config/dataset, all four pulls, QA, Z2 and ROOT metadata; fixture="<<TestDirectory<<std::endl;
    return 0;
  }catch(const std::exception& error) {
    std::cerr<<"FAIL helium calibration after "<<CalibrationChecks<<" checks: "<<error.what()<<std::endl;return 1;
  }
}
#endif
