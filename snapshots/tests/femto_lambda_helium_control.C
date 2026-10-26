int femto_lambda_helium_control(const char* mainconf,const char* rollback);
#ifndef __CINT__
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "HistManager.h"
#include "YamlParser.h"
#include "TH1.h"
#include "TH2.h"
#include "TMemFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace {
void RequireHelium(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
double HeliumKstarOracle(const FemtoCandidate& a,const FemtoCandidate& b){
  TLorentzVector boosted=a.P4();boosted.Boost(-(a.P4()+b.P4()).BoostVector());return boosted.Vect().Mag();
}
void PartitionHelium(const TH1& all,const TH1& low,const TH1& high){
  const int cells=(all.GetNbinsX()+2)*(all.GetDimension()>1?all.GetNbinsY()+2:1);
  RequireHelium(all.GetEntries()==low.GetEntries()+high.GetEntries(),"half entry closure");
  for(int bin=0;bin<cells;++bin){
    RequireHelium(all.GetBinContent(bin)==low.GetBinContent(bin)+high.GetBinContent(bin),"half bin closure including flows");
    const double variance=std::pow(all.GetBinError(bin),2),halves=std::pow(low.GetBinError(bin),2)+std::pow(high.GetBinError(bin),2);
    RequireHelium(std::fabs(variance-halves)<1.e-9*std::max(1.,variance),"half variance closure");
  }
}
double HeliumScalar(const std::map<std::string,std::string>& values,const char* key){
  std::map<std::string,std::string>::const_iterator i=values.find(key);RequireHelium(i!=values.end(),"missing test config key");
  std::istringstream parser(i->second);double value=0.;parser>>value;RequireHelium(bool(parser),"invalid test config scalar");return value;
}
}

int femto_lambda_helium_control(const char* mainconf,const char* rollback){
  try {
    std::string mainPath=mainconf;if(mainPath[0]!='/')mainPath=std::string(gSystem->WorkingDirectory())+"/"+mainPath;
    std::map<std::string,std::string> main,maker;
    RequireHelium(YamlParser::ParseFile(mainPath.c_str(),main),"mainconf parse");
    const std::string configBase=mainPath.substr(0,mainPath.rfind("/config/"))+"/config/";
    RequireHelium(YamlParser::ParseFile((configBase+main["maker"]).c_str(),maker),"maker parse");
    RequireHelium(maker["lambdaHe3ControlEnabled"]=="true","test requires explicit control profile");
    const double mean=HeliumScalar(maker,"lambdaSignalMean"),width=HeliumScalar(maker,"lambdaSignalSigma")*HeliumScalar(maker,"lambdaSignalNSigma");
    const double lambdaMass=HeliumScalar(maker,"lambdaPairMass");
    FemtoCandidate l,n3,n4,mapped;l.speciesKey="lambda";l.source=kFemtoCandResonance;l.eventIndex=10;
    l.reso.dau1EventIndex=10;l.reso.dau2EventIndex=10;l.reso.dau1Index=2;l.reso.dau2Index=3;
    TLorentzVector lambda;lambda.SetXYZM(.32,.14,.8,lambdaMass);l.SetP4(lambda);l.reso.invMass=mean;
    {
      FemtoLambdaLegacy old;RequireHelium(old.Init(rollback,std::cerr),"rollback Init");
      n3.speciesKey="he3";n3.source=kFemtoCandTrack;n3.eventIndex=10;n3.trk.trackIndex=9;n3.SetP4(old.NuclearP4(TVector3(.13,-.09,.43),2));
      RequireHelium(!old.IsPairNucleus(2)&&old.IsPairNucleus(3),"legacy remains a single He4 target");
      RequireHelium(!old.FillPair(l,n3,HeliumKstarOracle(l,n3),.2,8,false),"legacy rejects additional He3 pairs");
    }
    FemtoLambdaLegacy helper;RequireHelium(helper.Init(mainconf,std::cerr),"control Init");
    RequireHelium(helper.IsPairNucleus(2)&&helper.IsPairNucleus(3)&&!helper.IsPairNucleus(0),"dual targets only");
    n3.charge=2;n3.SetP4(helper.NuclearP4(TVector3(.13,-.09,.43),2));
    const double he3Mass=helper.NuclearP4(TVector3(),2).M(),he4Mass=helper.NuclearP4(TVector3(),3).M();
    mapped=n3;TLorentzVector mappedP4;mappedP4.SetVectM(n3.P4().Vect(),he4Mass);mapped.SetP4(mappedP4);
    n4=mapped;n4.speciesKey="he4";
    const double ks3=HeliumKstarOracle(l,n3),ks4=HeliumKstarOracle(l,mapped),qlab=(l.P4().Vect()-n3.P4().Vect()).Mag();
    RequireHelium(n3.px==mapped.px&&n3.py==mapped.py&&n3.pz==mapped.pz&&n3.energy!=mapped.energy,"mass mapping changes only energy, never Z=2 momentum");
    TH1* physical=helper.RootHistograms()->Get("hKstarSE_lambda_he3");TH1* control=helper.RootHistograms()->Get("hKstarSE_lambda_he3_as_he4");
    RequireHelium(physical->FindBin(ks3)!=control->FindBin(ks4),"fixture distinguishes physical and mapped kstar bins");
    helper.FillNucleusAcceptance(n3,2);helper.FillNucleusAcceptance(n4,3);
    RequireHelium(helper.RootHistograms()->Get("hNucleus_PtVsYLab_he3")->GetEntries()==1&&helper.RootHistograms()->Get("hNucleus_PtVsYLab_he3_as_he4")->GetEntries()==1&&helper.RootHistograms()->Get("hNucleus_PtVsYLab_he4")->GetEntries()==1,"acceptance once per track");
    const float masses[]={float(mean-.5*width),float(mean+.5*width),float(mean-2*width),float(mean+2*width),float(mean+5*width)};
    for(int mixed=0;mixed<2;++mixed)for(int cent=0;cent<9;++cent)for(unsigned i=0;i<sizeof(masses)/sizeof(masses[0]);++i){
      l.reso.invMass=masses[i];n3.eventIndex=n4.eventIndex=10+mixed;
      n3.trk.trackIndex=n4.trk.trackIndex=mixed?2:9;
      RequireHelium(helper.FillPair(l,n3,ks3,qlab,cent,mixed!=0),"physical He3 pair and mapped fanout accepted");
      RequireHelium(helper.FillPair(l,n4,ks4,qlab,cent,mixed!=0),"He4 pair accepted");
    }
    const char* species[]={"he3","he3_as_he4","he4"};const char* legacy[]={"3He","3HeAs4He","4He"};
    for(int mode=0;mode<2;++mode)for(int target=0;target<3;++target){
      const double ks=target==0?ks3:ks4,mass=target==0?he3Mass:he4Mass;
      const std::string canonical=std::string("hKstar")+(mode?"ME":"SE")+"_lambda_"+species[target];
      const std::string stem=std::string(mode?"Mixed_":"")+legacy[target];
      TH1* all=helper.RootHistograms()->Get(canonical.c_str());RequireHelium(all->GetEntries()==45&&all->GetBinContent(all->FindBin(ks))==45,"full same-pair kstar counts and independent boost oracle");
      PartitionHelium(*helper.RootHistograms()->Get((canonical+"_signal").c_str()),*helper.RootHistograms()->Get((canonical+"_signalLow").c_str()),*helper.RootHistograms()->Get((canonical+"_signalHigh").c_str()));
      for(int cent=0;cent<9;++cent){
        const std::string tail="_CentBin"+std::to_string(cent);
        TH2* full=dynamic_cast<TH2*>(helper.NuclearHistograms()->Get(("hKstarMass_"+stem+tail).c_str()));RequireHelium(full&&full->GetEntries()==5,"full-mass pair counts");
        TH2* mt=dynamic_cast<TH2*>(helper.NuclearHistograms()->Get(("hKstarMt_"+stem+tail).c_str()));
        const double px=double(l.px)+n3.px,py=double(l.py)+n3.py,expectedMt=std::sqrt(.25*(px*px+py*py+(lambdaMass+mass)*(lambdaMass+mass)));
        RequireHelium(mt&&mt->GetEntries()==2&&mt->GetBinContent(mt->FindBin(ks,expectedMt))==2,"mass-hypothesis mT oracle");
        TH1D* projection=mt->ProjectionX("helium_mt_projection",0,mt->GetNbinsY()+1,"e");projection->SetDirectory(0);
        TH1* reference=helper.NuclearHistograms()->Get(("hKstar_"+stem+tail).c_str());
        for(int bin=0;bin<=reference->GetNbinsX()+1;++bin)RequireHelium(projection->GetBinContent(bin)==reference->GetBinContent(bin)&&std::fabs(std::pow(projection->GetBinError(bin),2)-std::pow(reference->GetBinError(bin),2))<1.e-9,"mT projection closure including flow");
        delete projection;
        PartitionHelium(*mt,*helper.NuclearHistograms()->Get(("hKstarMt_"+stem+"_signalLow"+tail).c_str()),*helper.NuclearHistograms()->Get(("hKstarMt_"+stem+"_signalHigh"+tail).c_str()));
        RequireHelium(helper.NuclearHistograms()->Get(("hKstar_"+stem+"_SBNeg"+tail).c_str())->GetEntries()==1&&helper.NuclearHistograms()->Get(("hKstar_"+stem+"_SBPos"+tail).c_str())->GetEntries()==1,"both sideband populations retained");
      }
    }
    const double before=physical->GetEntries();n3.eventIndex=10;n3.trk.trackIndex=2;
    RequireHelium(!helper.FillPair(l,n3,ks3,qlab,8,false),"shared-track veto");
    n3.trk.trackIndex=9;RequireHelium(!helper.FillPair(l,n3,ks3,qlab,8,true),"ME same-event veto");
    RequireHelium(!helper.FillPair(l,n3,std::numeric_limits<double>::quiet_NaN(),qlab,8,false),"nonfinite pair veto");
    RequireHelium(physical->GetEntries()==before&&control->GetEntries()==before,"invalid pairs cannot partially fill the mapped control");
    TMemFile file("helium_control_test","RECREATE");RequireHelium(helper.Write(&file),"output write");
    TNamed* definition=dynamic_cast<TNamed*>(file.Get("FemtoNuclearDedxCalibration"));RequireHelium(definition&&definition->GetTitle()==helper.NuclearDedxCalibrationMetadata(),"exact PID metadata persisted");
    RequireHelium(file.Get("FemtoLambdaHe3ControlDefinition")&&file.Get("FemtoLambdaPairMtDefinition_he3")&&file.Get("FemtoLambdaPairMtDefinition_he3_as_he4"),"physical/mapped mass metadata present");
    for(int mode=0;mode<2;++mode)for(int target=0;target<3;++target){
      const std::string key=std::string("hKstarMass_")+(mode?"Mixed_":"")+legacy[target]+"_CentBin8";
      RequireHelium(file.Get((std::string(mode?"mix/":"true/")+key).c_str())!=0&&file.Get((std::string(mode?"true/":"mix/")+key).c_str())==0,"SE/ME directories preserved");
    }
    std::cout<<"PASS helium_control: same pair counts, independent kstar/mT mass mapping, SE/ME/SB/halves/cent/acceptance and metadata"<<std::endl;return 0;
  }catch(const std::exception& error){std::cerr<<"FAIL helium_control: "<<error.what()<<std::endl;return 1;}
}
#endif
