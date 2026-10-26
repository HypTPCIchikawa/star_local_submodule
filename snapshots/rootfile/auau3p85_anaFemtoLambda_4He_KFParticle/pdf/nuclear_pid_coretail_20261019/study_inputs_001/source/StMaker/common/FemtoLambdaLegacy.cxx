#include "FemtoLambdaLegacy.h"
#include "FemtoLambdaProvider.h"
#include "HistManager.h"
#include "NuclearIdDeDxVsMom.h"
#include "FemtoNuclearDedxCalibration.h"
#include "FemtoNuclearPairPidSparse.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "yaml-cpp/yaml.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TVector2.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace {
const char* const kSpecies[] = {"deuteron", "triton", "he3", "he4"};
const char* const kLegacy[] = {"d", "t", "3He", "4He"};
const char* const kMassKey[] = {"Deuteron", "Triton", "He3", "He4"};
const char* const kDomainSpecies[] = {"d", "t", "he3", "he4"};
std::string Index(int i) { std::ostringstream s; s << i; return s.str(); }
template<class T> T Required(const YAML::Node& n, const char* key) {
  if (!n[key]) throw std::runtime_error(std::string("Missing required key: ") + key);
  return n[key].as<T>();
}
double Number(const YAML::Node& n, const char* key) {
  const double v = Required<double>(n, key);
  if (!std::isfinite(v)) throw std::runtime_error(std::string("Nonfinite key: ") + key);
  return v;
}
void KnownKeys(const YAML::Node& n, const char* const* keys, unsigned count,
               const char* label) {
  if (!n.IsMap()) throw std::runtime_error(std::string(label)+" must be a map");
  std::vector<std::string> seen;
  for (YAML::const_iterator i=n.begin(); i!=n.end(); ++i) {
    const std::string key=i->first.as<std::string>();
    if(std::find(seen.begin(),seen.end(),key)!=seen.end())
      throw std::runtime_error(std::string("Duplicate ")+label+" key: "+key);
    seen.push_back(key);
    bool known=false;
    for(unsigned j=0;j<count;++j) known=known || key==keys[j];
    if(!known) throw std::runtime_error(std::string("Unknown ")+label+" key: "+key);
  }
}
std::string Resolve(const std::string& base, const YAML::Node& main, const char* key) {
  const std::string ref = Required<std::string>(main, key);
  if (ref.empty()) throw std::runtime_error(std::string("Empty reference: ") + key);
  return ref[0] == '/' ? ref : base + ref;
}
void KeysAndRequireLoaded(const std::string& path, HistManager& hist, std::vector<std::string>& keys) {
  YAML::Node doc = YAML::LoadFile(path);
  if (!doc["histograms"] || !doc["histograms"].IsMap()) throw std::runtime_error("Invalid histogram map: " + path);
  keys.clear();
  for (YAML::const_iterator i = doc["histograms"].begin(); i != doc["histograms"].end(); ++i) {
    const std::string name = i->first.as<std::string>();
    if (!hist.HasHistogram(name.c_str())) throw std::runtime_error("Histogram was not created: " + name);
    keys.push_back(name);
  }
}
void Need(HistManager& h, const std::string& key) {
  if (!h.HasHistogram(key.c_str())) throw std::runtime_error("Required histogram missing: " + key);
}
bool SameAxis(const TAxis& a, const TAxis& b) {
  if (a.GetNbins() != b.GetNbins()) return false;
  for (int i=1;i<=a.GetNbins()+1;++i)
    if (a.GetBinLowEdge(i) != b.GetBinLowEdge(i)) return false;
  return true;
}
void NeedPairMt(HistManager& hist, const std::string& key,
                const TH1& kstar, const TH2*& reference) {
  TH2* h = dynamic_cast<TH2*>(hist.Get(key.c_str()));
  if (!h || h->GetDimension()!=2 || !SameAxis(*h->GetXaxis(),*kstar.GetXaxis()))
    throw std::runtime_error("Missing/invalid kstar-mT histogram or mismatched kstar axis: "+key);
  if (reference && !SameAxis(*h->GetYaxis(),*reference->GetYaxis()))
    throw std::runtime_error("Inconsistent pair mT axes: "+key);
  reference=h;
}
double CurvatureAngle(double charge, double field, double radius, double pt) {
  const double a = 0.3 * charge * field * radius / (2.0 * pt);
  return std::asin(std::max(-1.0, std::min(1.0, a)));
}
double LambdaPairKstar(const TLorentzVector& a,const TLorentzVector& b) {
  // Same PRF q/2 convention as StFemtoMaker::ComputeKStar.
  TLorentzVector difference=a-b,pair=a+b;
  difference.Boost(-pair.BoostVector());
  return .5*difference.Vect().Mag();
}
}

FemtoLambdaLegacy::FemtoLambdaLegacy()
    : mRoot(0), mNuclear(0), mNuclearPidPairQa(0), mNuclearPidPairSparse(0),
      mDedxCalibration(new FemtoNuclearDedxCalibration),
      mSpeciesIndex(-1), mHe3ControlEnabled(false), mMean(0), mWindow(0),
      mOuterFactor(0), mFieldTesla(0), mRadiusMeters(0),
      mDaughterProtonMass(0), mDaughterPionMass(0), mLambdaPairMass(0), mMinNHitsDedx(0),
      mMinPt(0), mNSigmaFill(0), mNSigmaExclude(0), mMaxNSigma(0), mM2SigmaCut(0),
      mMinPTofQa(0), mMinM2(0), mMaxM2(0), mMaxRigidity(0), mM2Selection(false),
      mHe4HighRigidityTof(false), mHe4TofThreshold(0), mHe4TofMass2Min(0), mHe4TofMass2Max(0),
      mSpeciesMinRigidityEnabled(false), mHighRigidityTof(false),
      mExtendHe4BeyondMaxRigidity(false), mHighTofThreshold(0),
      mSpeciesRigidityVetoConfigured(false) {
  for(int s=0;s<4;++s) mMinRigidity[s]=mHighTofMass2Min[s]=mHighTofMass2Max[s]=0;
}
FemtoLambdaLegacy::~FemtoLambdaLegacy() {
  delete mRoot; delete mNuclear; delete mNuclearPidPairQa;
  delete mNuclearPidPairSparse; delete mDedxCalibration;
}

bool FemtoLambdaLegacy::Init(const char* mainconf, std::ostream& diagnostics) {
  try {
    if (!mainconf || !*mainconf) throw std::runtime_error("Explicit mainconf required");
    std::string mainPath(mainconf);
    if (mainPath[0] != '/') mainPath = std::string(gSystem->WorkingDirectory()) + "/" + mainPath;
    const size_t pos = mainPath.rfind("/config/");
    if (pos == std::string::npos) throw std::runtime_error("mainconf must reside below config/ (reference root)");
    const std::string base = mainPath.substr(0, pos) + "/config/";
    const YAML::Node main = YAML::LoadFile(mainPath);
    const std::string makerPath = Resolve(base, main, "maker");
    const std::string nuclearPath = Resolve(base, main, "nuclearid");
    const YAML::Node maker = YAML::LoadFile(makerPath), nuclear = YAML::LoadFile(nuclearPath);
    if (Required<std::string>(maker,"nuclearSelectionProfile") != "legacy_nuclearid")
      throw std::runtime_error("Only explicit nuclearSelectionProfile: legacy_nuclearid is implemented");
    std::string list = Required<std::string>(maker, "speciesKeys");
    std::replace(list.begin(), list.end(), ',', ' ');
    std::istringstream names(list);
    std::string name;
    mHe3ControlEnabled = maker["lambdaHe3ControlEnabled"] ? maker["lambdaHe3ControlEnabled"].as<bool>() : false;
    mPairSpeciesIndices.clear();
    while (names >> name) {
      if (name == "lambda") continue;
      int type=-1;
      for (int s = 0; s < 4; ++s) if (name == kSpecies[s]) type=s;
      if(type<0 || std::find(mPairSpeciesIndices.begin(),mPairSpeciesIndices.end(),type)!=mPairSpeciesIndices.end())
        throw std::runtime_error("Unknown/duplicate Lambda nuclear species");
      mPairSpeciesIndices.push_back(type);
    }
    if(mHe3ControlEnabled) {
      if(mPairSpeciesIndices.size()!=2 || !IsPairNucleus(2) || !IsPairNucleus(3))
        throw std::runtime_error("lambdaHe3ControlEnabled requires exactly he4 and he3 bachelors");
      mSpeciesIndex=3;
    } else {
      if(mPairSpeciesIndices.size()!=1)throw std::runtime_error("Lambda legacy profile requires exactly one nuclear species");
      mSpeciesIndex=mPairSpeciesIndices[0];
    }
    mSpecies=kSpecies[mSpeciesIndex];mLegacySpecies=kLegacy[mSpeciesIndex];
    if (Required<std::string>(maker,"lambdaPairMassMode") != "fixed")
      throw std::runtime_error("Lambda pair mT requires the fixed-mass pair convention");
    mLambdaPairMass = Number(maker,"lambdaPairMass");
    if (mLambdaPairMass <= 0) throw std::runtime_error("Lambda pair mass must be positive");
    mMean = Number(maker,"lambdaSignalMean");
    mWindow = Number(maker,"lambdaSignalSigma") * Number(maker,"lambdaSignalNSigma");
    mOuterFactor = Number(maker,"lambdaSidebandOuterFactor");
    mDaughterProtonMass = Number(maker,"lambdaDaughterProtonMass");
    mDaughterPionMass = Number(maker,"lambdaDaughterPionMass");
    if (mDaughterProtonMass <= 0 || mDaughterPionMass <= 0)
      throw std::runtime_error("Lambda daughter QA masses must be positive");
    const YAML::Node kf = YAML::LoadFile(Resolve(base, main, "kf"));
    const double lower[] = {Number(kf,"minMass"), mMean-mWindow, mMean-mOuterFactor*mWindow, mMean+mWindow};
    const double upper[] = {Number(kf,"maxMass"), mMean+mWindow, mMean-mWindow, mMean+mOuterFactor*mWindow};
    const char* regionSuffix[] = {"","_signal","_leftSB","_rightSB"};
    const int nChannels=Required<int>(maker,"nChannels");
    if(nChannels!=int(4*mPairSpeciesIndices.size()))throw std::runtime_error("Lambda channel count disagrees with enabled nuclei");
    for (int i=0;i<nChannels;++i) {
      const std::string prefix = "channel_" + Index(i) + "_";
      const std::string channel = Required<std::string>(maker,(prefix+"name").c_str());
      int region=-1;
      for(size_t t=0;t<mPairSpeciesIndices.size();++t)for(int r=0;r<4;++r)
        if(channel==std::string("lambda_")+kSpecies[mPairSpeciesIndices[t]]+regionSuffix[r]) region=r;
      if(region>=0 && channel!="lambda_"+Required<std::string>(maker,(prefix+"partB").c_str())+regionSuffix[region])
        throw std::runtime_error("Lambda channel name and partB nucleus disagree");
      if(region<0 || std::fabs(Number(maker,(prefix+"signalMin").c_str())-lower[region])>1.e-9 ||
         std::fabs(Number(maker,(prefix+"signalMax").c_str())-upper[region])>1.e-9)
        throw std::runtime_error("Channel mass windows disagree with the single Lambda window definition");
    }
    mFieldTesla = Number(maker,"mergingFieldTesla");
    mRadiusMeters = Number(maker,"mergingRadiusMeters");
    mMinNHitsDedx = Required<int>(nuclear,"nuclearMinNHitsDedx");
    mMinPt = Number(nuclear,"nuclearMinPt");
    mNSigmaFill = Number(nuclear,"nSigmaFill1D"); mNSigmaExclude = Number(nuclear,"nSigmaExclude");
    mMaxNSigma = Number(nuclear,"maxNSigmaNuclear"); mM2SigmaCut = Number(nuclear,"m2SigmaCut");
    mM2Selection = Required<bool>(nuclear,"m2_selection");
    mMinPTofQa = Number(nuclear,"minP_M2cut"); mMinM2 = Number(nuclear,"minM2_M2cut");
    mMaxM2 = Number(nuclear,"maxM2_M2cut"); mMaxRigidity = Number(nuclear,"maxPOverQ");
    // No block means the original shared tables. An explicit bad calibration
    // must fail Init, never quietly run with a different PID selection.
    mDedxCalibration->ResetLegacy();
    if(nuclear["dedxCalibration"]) {
      const YAML::Node calibration=nuclear["dedxCalibration"];
      if(!calibration.IsMap()) throw std::runtime_error("dedxCalibration must be a map");
      for(YAML::const_iterator i=calibration.begin();i!=calibration.end();++i) {
        const std::string key=i->first.as<std::string>();
        if(key!="enabled" && key!="calibrationFile")
          throw std::runtime_error("Unknown dedxCalibration key: "+key);
      }
      if(Required<bool>(calibration,"enabled")) {
        const YAML::Node analysis=YAML::LoadFile(Resolve(base,main,"analysis"));
        const YAML::Node tag=analysis["starTag"];
        if(!tag || !tag.IsMap() || Required<std::string>(tag,"collision")!="auau3p85" ||
           Required<std::string>(tag,"year")!="year2021" ||
           Required<std::string>(tag,"productionTag")!="P24iy")
          throw std::runtime_error("Explicit helium dE/dx calibration requires auau3p85/year2021/P24iy");
        const std::string calibrationPath=Resolve(base,calibration,"calibrationFile");
        if(!mDedxCalibration->Load(calibrationPath,diagnostics))
          throw std::runtime_error("Explicit helium dE/dx calibration could not be loaded");
      } else if(calibration["calibrationFile"]) {
        Required<std::string>(calibration,"calibrationFile");
      }
    }
    // An absent block deliberately preserves all historical mainconfs.
    mHe4HighRigidityTof = false;
    if (nuclear["he4HighRigidityTof"]) {
      const YAML::Node tof = nuclear["he4HighRigidityTof"];
      if (!tof.IsMap()) throw std::runtime_error("he4HighRigidityTof must be a map");
      for (YAML::const_iterator i=tof.begin(); i!=tof.end(); ++i) {
        const std::string key=i->first.as<std::string>();
        if (key!="enabled" && key!="thresholdPOverQ" && key!="mass2Min" && key!="mass2Max")
          throw std::runtime_error("Unknown he4HighRigidityTof key: "+key);
      }
      mHe4HighRigidityTof = Required<bool>(tof,"enabled");
      mHe4TofThreshold = Number(tof,"thresholdPOverQ");
      mHe4TofMass2Min = Number(tof,"mass2Min");
      mHe4TofMass2Max = Number(tof,"mass2Max");
      if (!(mHe4TofThreshold>0 && mHe4TofMass2Min>=0 && mHe4TofMass2Max>mHe4TofMass2Min))
        throw std::runtime_error("Invalid he4HighRigidityTof threshold / (m/q)^2 window");
      if (mHe4HighRigidityTof && (mSpeciesIndex!=3 || mM2Selection))
        throw std::runtime_error("he4HighRigidityTof requires the he4 target and m2_selection: false");
    }
    // These are final-identity gates, not changes to the TPC hypothesis list.
    // In particular, a rejected low-purity winner must not promote another nucleus.
    mSpeciesMinRigidityEnabled=false;
    for(int s=0;s<4;++s) mMinRigidity[s]=0;
    if(nuclear["speciesMinPOverQ"]) {
      const YAML::Node minimum=nuclear["speciesMinPOverQ"];
      KnownKeys(minimum,kDomainSpecies,4,"speciesMinPOverQ");
      for(int s=0;s<4;++s) {
        mMinRigidity[s]=Number(minimum,kDomainSpecies[s]);
        if(mMinRigidity[s]<0) throw std::runtime_error("speciesMinPOverQ values must be nonnegative");
      }
      mSpeciesMinRigidityEnabled=true;
    }
    mSpeciesRigidityVetoConfigured=false;
    for(int s=0;s<4;++s) mRigidityVeto[s].clear();
    unsigned vetoBlocks=0;
    for(YAML::const_iterator i=nuclear.begin();i!=nuclear.end();++i)
      if(i->first.as<std::string>()=="speciesPOverQVeto") ++vetoBlocks;
    if(vetoBlocks>1) throw std::runtime_error("Duplicate nuclear key: speciesPOverQVeto");
    if(nuclear["speciesPOverQVeto"]) {
      const YAML::Node veto=nuclear["speciesPOverQVeto"];
      KnownKeys(veto,kDomainSpecies,4,"speciesPOverQVeto");
      for(int s=0;s<4;++s) {
        const YAML::Node ranges=veto[kDomainSpecies[s]];
        if(!ranges || !ranges.IsSequence())
          throw std::runtime_error("speciesPOverQVeto requires an explicit sequence for every species");
        for(size_t i=0;i<ranges.size();++i) {
          const YAML::Node node=ranges[i];
          const char* keys[]={"min","max"};
          KnownKeys(node,keys,2,"speciesPOverQVeto interval");
          RigidityVetoRange range;range.minimum=Number(node,"min");range.maximum=Number(node,"max");
          if(!(range.minimum>0 && range.maximum>range.minimum))
            throw std::runtime_error("speciesPOverQVeto needs positive min < max");
          if(!mRigidityVeto[s].empty() && range.minimum<mRigidityVeto[s].back().maximum)
            throw std::runtime_error("speciesPOverQVeto intervals must be sorted and non-overlapping");
          mRigidityVeto[s].push_back(range);
        }
      }
      mSpeciesRigidityVetoConfigured=true;
    }
    mHighRigidityTof=false;mExtendHe4BeyondMaxRigidity=false;mHighTofThreshold=0;
    for(int s=0;s<4;++s) mHighTofMass2Min[s]=mHighTofMass2Max[s]=0;
    if(nuclear["highRigidityTof"]) {
      const YAML::Node tof=nuclear["highRigidityTof"];
      const char* keys[]={"enabled","thresholdPOverQ","extendHe4BeyondMaxPOverQ","species"};
      KnownKeys(tof,keys,4,"highRigidityTof");
      mHighRigidityTof=Required<bool>(tof,"enabled");
      mHighTofThreshold=Number(tof,"thresholdPOverQ");
      mExtendHe4BeyondMaxRigidity=Required<bool>(tof,"extendHe4BeyondMaxPOverQ");
      if(!(mHighTofThreshold>0)) throw std::runtime_error("highRigidityTof threshold must be positive");
      const YAML::Node species=tof["species"];
      KnownKeys(species,kDomainSpecies,4,"highRigidityTof species");
      for(int s=0;s<4;++s) {
        const YAML::Node window=species[kDomainSpecies[s]];
        const char* windowKeys[]={"mass2Min","mass2Max"};
        KnownKeys(window,windowKeys,2,"highRigidityTof window");
        mHighTofMass2Min[s]=Number(window,"mass2Min");
        mHighTofMass2Max[s]=Number(window,"mass2Max");
        if(!(mHighTofMass2Min[s]>=0 && mHighTofMass2Max[s]>mHighTofMass2Min[s]))
          throw std::runtime_error("Invalid highRigidityTof (m/q)^2 window");
      }
      if(mHighRigidityTof && (mHe4HighRigidityTof || mM2Selection))
        throw std::runtime_error("highRigidityTof cannot be combined with enabled he4HighRigidityTof or m2_selection");
    }
    mPidDomainMetadata.clear();
    if(nuclear["speciesMinPOverQ"] || nuclear["speciesPOverQVeto"] || nuclear["highRigidityTof"]) {
      std::ostringstream description;description.precision(17);
      description << "nuclearConfig: " << nuclearPath
                  << "\nrigidity: raw |p/q| [GeV/c]; mass2: (m/q)^2 [(GeV/c^2)^2]"
                  << "\nselectionOrder: legacy TPC best species, winner minimum rigidity, winner rigidity veto, winner high-rigidity TOF, maximum rigidity"
                  << "\nfailedWinnerPolicy: reject track; never reassign to another species"
                  << "\nminimumBoundary: R >= minimum (inclusive); TOF threshold: R >= threshold"
                  << "\nTOFWindowBoundary: mass2Min < mass2 < mass2Max (strict)"
                  << "\nmaxPOverQ: " << mMaxRigidity << "; R < maxPOverQ (strict when positive)"
                  << "\nspeciesMinPOverQEnabled: " << mSpeciesMinRigidityEnabled
                  << "\nspeciesPOverQVetoConfigured: " << mSpeciesRigidityVetoConfigured
                  << "\nvetoBoundary: min <= R < max; reason: configured veto (not a statement of measured low purity)"
                  << "\nhighRigidityTofEnabled: " << mHighRigidityTof
                  << "\nextendHe4BeyondMaxPOverQ: " << mExtendHe4BeyondMaxRigidity;
      for(int s=0;s<4;++s)
        description << "\n" << kDomainSpecies[s] << ": minimum=" << mMinRigidity[s]
                    << ", mass2Min=" << mHighTofMass2Min[s] << ", mass2Max=" << mHighTofMass2Max[s];
      for(int s=0;s<4;++s) {
        description << "\nspeciesPOverQVeto." << kDomainSpecies[s] << ": [";
        for(size_t i=0;i<mRigidityVeto[s].size();++i) {
          if(i) description << ", ";
          description << "[" << mRigidityVeto[s][i].minimum << ", "
                      << mRigidityVeto[s][i].maximum << ")";
        }
        description << "]";
      }
      description << "\nthresholdPOverQ: " << mHighTofThreshold
                  << "\nQA: raw hDedxP, per-hypothesis PID QA and TOF controls are before final domain gates; CF candidates and nucleus acceptance are after gates";
      mPidDomainMetadata=description.str();
    }
    // The original StNuclearIdMaker ALWAYS selects one best hypothesis.
    if (!Required<bool>(nuclear,"requireBestSpecies"))
      throw std::runtime_error("legacy_nuclearid requires requireBestSpecies: true");
    for (int s = 0; s < 4; ++s) {
      mNuclearMass[s] = Number(nuclear, (std::string("nuclearMass") + kMassKey[s]).c_str());
      if (mNuclearMass[s] <= 0) throw std::runtime_error("Nuclear masses must be positive");
      if (s < 3) {
        mTofMean[s] = Number(nuclear, (std::string("nuclearTofMean") + kMassKey[s]).c_str());
        mTofSigma[s] = Number(nuclear, (std::string("nuclearTofSigma") + kMassKey[s]).c_str());
        if (mTofSigma[s] <= 0) throw std::runtime_error("TOF sigmas must be positive");
      }
    }
    if (!(mWindow > 0 && mOuterFactor > 1 && mRadiusMeters > 0 && mMinNHitsDedx >= 0 &&
          mMinPt >= 0 && mMaxNSigma > 0 && mMaxM2 > mMinM2 && mM2SigmaCut > 0))
      throw std::runtime_error("Invalid legacy selection/window configuration");
    const std::string rootPath = Resolve(base, main, "femtoHist"), histPath = Resolve(base, main, "nuclearHist");
    delete mRoot; mRoot = new HistManager;
    delete mNuclear; mNuclear = new HistManager;
    if (!mRoot->LoadFromFile(rootPath.c_str()) || !mNuclear->LoadFromFile(histPath.c_str()))
      throw std::runtime_error("Cannot load required Lambda histograms");
    KeysAndRequireLoaded(rootPath, *mRoot, mRootKeys);
    KeysAndRequireLoaded(histPath, *mNuclear, mNuclearKeys);
    delete mNuclearPidPairQa; mNuclearPidPairQa=0;
    delete mNuclearPidPairSparse; mNuclearPidPairSparse=0;
    mNuclearPidPairQaKeys.clear();
    if(main["nuclearPidPairQaHist"]) {
      const std::string qaPath=Resolve(base,main,"nuclearPidPairQaHist");
      mNuclearPidPairQa=new HistManager;
      if(!mNuclearPidPairQa->LoadFromFile(qaPath.c_str()))
        throw std::runtime_error("Cannot load optional selected-pair nuclear PID QA: "+qaPath);
      KeysAndRequireLoaded(qaPath,*mNuclearPidPairQa,mNuclearPidPairQaKeys);
      mNuclearPidPairSparse=new FemtoNuclearPairPidSparse;
      if(!mNuclearPidPairSparse->Init(qaPath.c_str(),diagnostics))
        throw std::runtime_error("Cannot initialize joint nuclear PID pair sparse: "+qaPath);
      diagnostics << "[FemtoLambdaLegacy] selected/pre-PID nuclear PID QA=" << qaPath
                  << " (no change to nuclear selection or mixing)" << std::endl;
    }
    if (!ValidateHistograms(diagnostics)) return false;
    diagnostics << "[FemtoLambdaLegacy] mainconf(argument)=" << mainPath
                << "\n  nuclear=" << nuclearPath << "\n  rootHist=" << rootPath
                << "\n  nuclearHist=" << histPath << "\n  target=" << mSpecies
                << " he3Control=" << mHe3ControlEnabled
                << " signal=[" << mMean-mWindow << "," << mMean+mWindow
                << "] outerFactor=" << mOuterFactor << std::endl;
    if (mHe4HighRigidityTof)
      diagnostics << "  he4HighRigidityTof: R >= " << mHe4TofThreshold
                  << ", " << mHe4TofMass2Min << " < (m/q)^2 < " << mHe4TofMass2Max
                  << "; He4 bypasses maxPOverQ=" << mMaxRigidity
                  << "; TPC best-species first" << std::endl;
    if(!mPidDomainMetadata.empty()) diagnostics << mPidDomainMetadata << std::endl;
    return true;
  } catch (const std::exception& e) {
    diagnostics << "[FemtoLambdaLegacy] Init failed: " << e.what() << std::endl;
    return false;
  }
}

bool FemtoLambdaLegacy::ValidateHistograms(std::ostream& diagnostics) {
  try {
    const char* rootNames[] = {"hLambda_InvMass","hLambda_Pt","hLambda_Eta","hLambda_Phi",
      "hDCA12","hDCAV0","hCosPointing","hNSigmaProton","hNSigmaPion",
      "hLambda_InvMass_vs_Pt","hLambda_InvMass_vs_DecayLength","hLambda_InvMass_vs_Y",
      "hLambda_InvMass_vs_TransDecayLength","hDCAV0_vs_InvMass","hCosPointing_vs_InvMass",
      "hLambda_InvMass_vs_Cent9","hLambda_InvMass_vs_RefMultCorr","hVz","hRefMult",
      "hCentrality","hCentralityRaw","hCentrality16","hRefMultCorr","hRefMultWeight",
      "hRawMult","hRefMultVsNTOFMatch","hRefMultVsNTOFMatchAfter","hCentralityVsVz",
      "hRawMult_vs_RefMultCorr","hN","hRawMult_vs_Cent9","hRefMultCorr_vs_Cent9",
      "hNTracks_vs_Cent9","hTofMatchMult_vs_Cent9","hNProtonCand_vs_Cent9",
      "hNPionCand_vs_Cent9","hNLambdaPairs_vs_Cent9","hMixSamplerQA",
      "hLambdaKF_DecayLength","hLambdaKF_DecayLengthSignificance","hLambdaKF_Chi2Ndf","hLambdaKF_TopoChi2Ndf"};
    for (unsigned i=0;i<sizeof(rootNames)/sizeof(rootNames[0]);++i) Need(*mRoot,rootNames[i]);
    for (int c=0;c<9;++c) Need(*mRoot,"hLambda_InvMass_CentBin"+Index(c));
    Need(*mRoot,"hLambda_PtVsYLab_signal");
    Need(*mRoot,"hLambdaProton_PtVsYLab_signal");
    Need(*mRoot,"hLambdaPion_PtVsYLab_signal");
    for(size_t i=0;i<mPairSpeciesIndices.size();++i)
      Need(*mRoot,std::string("hNucleus_PtVsYLab_")+kSpecies[mPairSpeciesIndices[i]]);
    if(mHe3ControlEnabled)Need(*mRoot,"hNucleus_PtVsYLab_he3_as_he4");
    const char* common[]={"hDedxP","hDedxP_cut","hDedxP_e","hDedxP_pi","hDedxP_K","hDedxP_p","hDedxP_else","hM2P"};
    for (unsigned i=0;i<sizeof(common)/sizeof(common[0]);++i) Need(*mNuclear,common[i]);
    if (He4TofQaEnabled()) {
      const char* qa[]={"hDedxP_4He_TOFControl", "hNSigmaHe4VsP_TOFControl",
        "hM2P_4He_TPC", "hDedxP_4He_selected", "hM2P_4He_selected", "hHe4PidStage"};
      for (unsigned i=0;i<sizeof(qa)/sizeof(qa[0]);++i) Need(*mNuclear,qa[i]);
    }
    for (int s=0;s<4;++s) {
      Need(*mNuclear,std::string("hDedxP_")+kLegacy[s]);
      for (int q=0;q<3;++q) {
        const char* vars[]={"Eta","Y","Pt"};
        Need(*mNuclear,std::string("hPvs")+vars[q]+"_"+kLegacy[s]);
      }
      if (s<3) Need(*mNuclear,std::string("hDedxP_")+kLegacy[s]+"_m2");
    }
    const TH2* mtReference=0;
    const TH2* pidQaReference[2]={0,0};
    std::vector<std::string> expectedPidQaKeys;
    const size_t representations=mPairSpeciesIndices.size()+(mHe3ControlEnabled?1:0);
    for(size_t target=0;target<representations;++target) {
      const bool mapped=target>=mPairSpeciesIndices.size();
      const int type=mapped?3:mPairSpeciesIndices[target];
      const std::string species=mapped?"he3_as_he4":kSpecies[type];
      const std::string legacy=mapped?"3HeAs4He":kLegacy[type];
    for (int mixed=0;mixed<2;++mixed) {
      const std::string prefix=mixed?"Mixed_":"";
      const std::string mode=mixed?"ME":"SE";
      const char* suffixes[]={"","_SBPos","_SBNeg"};
      for(int r=0;r<3;++r) {
        const std::string stem=prefix+legacy+suffixes[r];
        Need(*mNuclear,"hKstar_"+stem); Need(*mNuclear,"hQlab_"+stem);
        for(int c=0;c<9;++c) Need(*mNuclear,"hKstar_"+stem+"_CentBin"+Index(c));
      }
      const char* mtSuffixes[]={"","_SBPos","_SBNeg","_signalLow","_signalHigh"};
      for(int r=0;r<5;++r) for(int c=0;c<9;++c) {
        const std::string cent="_CentBin"+Index(c);
        const std::string source="hKstar_"+prefix+legacy+(r<3?mtSuffixes[r]:"")+cent;
        NeedPairMt(*mNuclear,"hKstarMt_"+prefix+legacy+mtSuffixes[r]+cent,
                   *mNuclear->Get(source.c_str()),mtReference);
        if(mNuclearPidPairQa && type>=2) {
          const char* variables[]={"NSigma","POverQ"};
          for(int stage=0;stage<2;++stage) for(int v=0;v<2;++v) {
            const std::string key="hKstar"+std::string(variables[v])+(stage?"PrePid":"")+"_"+prefix+legacy+mtSuffixes[r]+cent;
            TH2* qa=dynamic_cast<TH2*>(mNuclearPidPairQa->Get(key.c_str()));
            if(!qa || qa->GetDimension()!=2 || !SameAxis(*qa->GetXaxis(),*mNuclear->Get(source.c_str())->GetXaxis()))
              throw std::runtime_error("Missing/invalid selected-pair PID QA or mismatched kstar axis: "+key);
            if(pidQaReference[v] && !SameAxis(*qa->GetYaxis(),*pidQaReference[v]->GetYaxis()))
              throw std::runtime_error("Inconsistent selected-pair PID QA axis: "+key);
            if(mNuclear->HasHistogram(key.c_str()) || mRoot->HasHistogram(key.c_str()))
              throw std::runtime_error("Selected-pair PID QA collides with an existing histogram: "+key);
            pidQaReference[v]=qa;
            expectedPidQaKeys.push_back(key);
          }
        }
      }
      for(int c=0;c<9;++c) Need(*mNuclear,"hKstarMass_"+prefix+legacy+"_CentBin"+Index(c));
      const char* channels[]={"","_signal","_leftSB","_rightSB","_signalLow","_signalHigh"};
      for(int r=0;r<6;++r) {
        Need(*mRoot,"hKstar"+mode+"_lambda_"+species+channels[r]);
        Need(*mRoot,"hKstar"+mode+"VsCent_lambda_"+species+channels[r]);
      }
    }
    Need(*mNuclear,"hDphiDeta_proton_"+legacy);
    Need(*mNuclear,"hDphiDeta_proton_"+legacy+"_EventField");
    }
    if(mNuclearPidPairQa) {
      std::sort(expectedPidQaKeys.begin(),expectedPidQaKeys.end());
      std::vector<std::string> actual=mNuclearPidPairQaKeys;
      std::sort(actual.begin(),actual.end());
      if(actual!=expectedPidQaKeys || actual.empty())
        throw std::runtime_error("Optional selected-pair PID QA keys must exactly match enabled helium representations");
      if(!mNuclearPidPairSparse || !mNuclearPidPairSparse->Histogram() ||
         !pidQaReference[0] || !pidQaReference[1] ||
         !SameAxis(*mNuclearPidPairSparse->Histogram()->GetAxis(0),*pidQaReference[0]->GetXaxis()) ||
         !SameAxis(*mNuclearPidPairSparse->Histogram()->GetAxis(1),*pidQaReference[1]->GetYaxis()) ||
         !SameAxis(*mNuclearPidPairSparse->Histogram()->GetAxis(2),*pidQaReference[0]->GetYaxis()))
        throw std::runtime_error("Joint and dense nuclear PID pair QA physical axes must match exactly");
    }
    return true;
  } catch (const std::exception& e) { diagnostics << e.what() << std::endl; return false; }
}

int FemtoLambdaLegacy::MassRegion(double mass) const {
  if (!std::isfinite(mass)) return 0;
  const double delta=mass-mMean;
  if (std::fabs(delta)<=mWindow) return 1;
  if (delta < -mWindow && delta >= -mOuterFactor*mWindow) return 2;
  if (delta > mWindow && delta <= mOuterFactor*mWindow) return 3;
  return 0;
}

int FemtoLambdaLegacy::SignalHalf(double mass) const {
  if (MassRegion(mass) != 1) return 0;
  return mass < mMean ? 1 : 2;
}

int FemtoLambdaLegacy::SelectNuclearType(double p, const double pulls[4], bool tof, double m2) const {
  if (!std::isfinite(p) || p<=0) return -1;
  const bool validTof=tof && std::isfinite(m2);
  const bool extendHe4=mHe4HighRigidityTof || (mHighRigidityTof && mExtendHe4BeyondMaxRigidity);
  if (!extendHe4 && mMaxRigidity > 0 && !(p < mMaxRigidity)) return -1;
  if (mM2Selection && !validTof) return -1;
  int best=-1; double minPull=999.;
  for(int s=0;s<4;++s) {
    bool pass=std::isfinite(pulls[s]) && std::fabs(pulls[s])<mMaxNSigma;
    if(mM2Selection) {
      const double lo=s<3?mTofMean[s]-mM2SigmaCut*mTofSigma[s]:mMinM2;
      const double hi=s<3?mTofMean[s]+mM2SigmaCut*mTofSigma[s]:mMaxM2;
      pass=pass && m2>lo && m2<hi;
    }
    if(pass && std::fabs(pulls[s])<minPull) { best=s; minPull=std::fabs(pulls[s]); }
  }
  if(best<0) return -1;
  if(mSpeciesMinRigidityEnabled && p<mMinRigidity[best]) return -1;
  for(size_t i=0;i<mRigidityVeto[best].size();++i)
    if(p>=mRigidityVeto[best][i].minimum && p<mRigidityVeto[best][i].maximum) return -1;
  if(mHighRigidityTof) {
    if(!PassNuclearHighRigidityTof(best,p,validTof,m2)) return -1;
    if(best==3 && mExtendHe4BeyondMaxRigidity) return best;
    if(mMaxRigidity>0 && !(p<mMaxRigidity)) return -1;
    return best;
  }
  // Determine the TPC identity before gating He4. A closer He3 hypothesis must
  // never become He4 merely because it failed an additional TOF requirement.
  if (best==3 && mHe4HighRigidityTof)
    return PassHe4HighRigidityTof(p,validTof,m2) ? best : -1;
  if (mMaxRigidity > 0 && !(p < mMaxRigidity)) return -1;
  return best;
}
bool FemtoLambdaLegacy::PassHe4HighRigidityTof(double p, bool tof, double m2) const {
  return PassNuclearHighRigidityTof(3,p,tof,m2);
}
bool FemtoLambdaLegacy::PassNuclearHighRigidityTof(int species,double p,bool tof,double m2) const {
  if(species<0 || species>=4) return false;
  if (!std::isfinite(p) || p<=0) return false;
  if(mHighRigidityTof) {
    if(p<mHighTofThreshold) return true;
    return tof && std::isfinite(m2) && m2>mHighTofMass2Min[species] && m2<mHighTofMass2Max[species];
  }
  if(species!=3) return true;
  if (!mHe4HighRigidityTof || p<mHe4TofThreshold) return true;
  return tof && std::isfinite(m2) && m2>mHe4TofMass2Min && m2<mHe4TofMass2Max;
}
bool FemtoLambdaLegacy::He4TofQaEnabled() const {
  return mSpeciesIndex==3 && (mHe4HighRigidityTof || mHighRigidityTof);
}
double FemtoLambdaLegacy::He4TofThreshold() const {
  return mHighRigidityTof?mHighTofThreshold:mHe4TofThreshold;
}
TLorentzVector FemtoLambdaLegacy::NuclearP4(const TVector3& raw, int type) const {
  TLorentzVector out;
  if(type>=0 && type<4) out.SetVectM(raw*(type>=2?2.:1.),mNuclearMass[type]);
  return out;
}

double FemtoLambdaLegacy::NuclearDedxPull(int species,double rigidity,double dedx) const {
  return mDedxCalibration->Pull(species,rigidity,dedx);
}
const std::string& FemtoLambdaLegacy::NuclearDedxCalibrationMetadata() const {
  return mDedxCalibration->Metadata();
}

bool FemtoLambdaLegacy::PassNuclearPidQaDomain(int type,double p,bool tof,double m2) const {
  if(type<0 || type>=4 || !std::isfinite(p) || p<=0) return false;
  const bool validTof=tof && std::isfinite(m2);
  if(mM2Selection) {
    const double lo=type<3?mTofMean[type]-mM2SigmaCut*mTofSigma[type]:mMinM2;
    const double hi=type<3?mTofMean[type]+mM2SigmaCut*mTofSigma[type]:mMaxM2;
    if(!validTof || !(m2>lo && m2<hi)) return false;
  }
  if(mSpeciesMinRigidityEnabled && p<mMinRigidity[type]) return false;
  for(size_t i=0;i<mRigidityVeto[type].size();++i)
    if(p>=mRigidityVeto[type][i].minimum && p<mRigidityVeto[type][i].maximum) return false;
  // Exactly the target's non-TPC final gates; no nSigma cut or winner decision.
  if(mHighRigidityTof) {
    if(!PassNuclearHighRigidityTof(type,p,validTof,m2)) return false;
    if(type==3 && mExtendHe4BeyondMaxRigidity) return true;
    return !(mMaxRigidity>0 && !(p<mMaxRigidity));
  }
  if(type==3 && mHe4HighRigidityTof)
    return PassHe4HighRigidityTof(p,validTof,m2);
  return !(mMaxRigidity>0 && !(p<mMaxRigidity));
}

void FemtoLambdaLegacy::CollectNuclei(StPicoDst* dst,int eventIndex,int cent9,
    FemtoCandidateStore& store,FemtoNuclearPidPairQaStore* pidQa) {
  if(!dst || !dst->event() || !mNuclear) return;
  int multiplicity[4]={0,0,0,0};
  for(unsigned it=0;it<dst->numberOfTracks();++it) {
    StPicoTrack* t=dst->track(it);
    if(!t || t->nHitsDedx()<mMinNHitsDedx || t->gPt()<mMinPt) continue;
    const TVector3 pvec=t->pMom();
    const double p=pvec.Mag(), dedx=t->dEdx();
    if(!std::isfinite(p) || p<=0 || !std::isfinite(dedx) || dedx<=0) continue;
    mNuclear->Fill("hDedxP",p,dedx);
    const double usual[]={t->nSigmaElectron(),t->nSigmaPion(),t->nSigmaKaon(),t->nSigmaProton()};
    const char* usualKey[]={"e","pi","K","p"};
    bool excluded=true;
    for(int s=0;s<4;++s) {
      if(std::fabs(usual[s])<mNSigmaFill) mNuclear->Fill((std::string("hDedxP_")+usualKey[s]).c_str(),p,dedx);
      excluded=excluded && std::fabs(usual[s])>mNSigmaExclude;
    }
    if(excluded) mNuclear->Fill("hDedxP_else",p,dedx);
    double pulls[4];
    for(int s=0;s<4;++s) {
      pulls[s]=NuclearDedxPull(s,p,dedx);
      if(std::fabs(pulls[s])<mMaxNSigma) {
        const TLorentzVector h=NuclearP4(pvec,s);
        mNuclear->Fill((std::string("hDedxP_")+kLegacy[s]).c_str(),p,dedx);
        mNuclear->Fill((std::string("hPvsEta_")+kLegacy[s]).c_str(),h.P(),pvec.Eta());
        mNuclear->Fill((std::string("hPvsY_")+kLegacy[s]).c_str(),h.P(),h.Rapidity());
        mNuclear->Fill((std::string("hPvsPt_")+kLegacy[s]).c_str(),h.P(),h.Pt());
      }
    }
    bool tof=false; double m2=-999.;
    if(t->isTofTrack() && t->bTofPidTraitsIndex()>=0 &&
       static_cast<unsigned>(t->bTofPidTraitsIndex())<dst->numberOfBTofPidTraits()) {
      StPicoBTofPidTraits* b=dst->btofPidTraits(t->bTofPidTraitsIndex());
      if(b && std::isfinite(b->btofBeta()) && b->btofBeta()>0) {
        const double beta=b->btofBeta(); m2=p*p*(1./(beta*beta)-1.);
        tof=std::isfinite(m2);
      }
    }
    if (He4TofQaEnabled()) {
      // A pre-TPC control sample can be fitted without truncating the He4 peak.
      // At low R it includes all quality tracks; at high R it includes the TOF gate.
      if (PassHe4HighRigidityTof(p,tof,m2)) {
        mNuclear->Fill("hDedxP_4He_TOFControl",p,dedx);
        mNuclear->Fill("hNSigmaHe4VsP_TOFControl",p,pulls[3]);
      }
      int tpcBest=-1; double closest=mMaxNSigma;
      for(int s=0;s<4;++s) if(std::isfinite(pulls[s]) && std::fabs(pulls[s])<closest) {
        tpcBest=s; closest=std::fabs(pulls[s]);
      }
      if(tpcBest==3) {
        mNuclear->Fill("hHe4PidStage",0.);
        if(p<He4TofThreshold()) mNuclear->Fill("hHe4PidStage",1.);
        else {
          mNuclear->Fill("hHe4PidStage",2.);
          if(tof) {
            mNuclear->Fill("hHe4PidStage",3.);
            mNuclear->Fill("hM2P_4He_TPC",p,m2);
          }
          if(PassHe4HighRigidityTof(p,tof,m2)) mNuclear->Fill("hHe4PidStage",4.);
        }
      }
    }
    const int best=SelectNuclearType(p,pulls,tof,m2);
    if(pidQa && mNuclearPidPairQa) for(size_t index=0;index<mPairSpeciesIndices.size();++index) {
      const int type=mPairSpeciesIndices[index];
      if(type<2 || !PassNuclearPidQaDomain(type,p,tof,m2)) continue;
      FemtoNuclearPidPairQaCandidate qa;
      FemtoCandidate& c=qa.candidate;
      c.eventIndex=eventIndex; c.source=kFemtoCandTrack;
      c.speciesKey=kSpecies[type]; c.charge=t->charge()*2;
      c.SetP4(NuclearP4(pvec,type)); c.trk.trackIndex=it; c.trk.nHitsFit=t->nHitsFit();
      c.trk.nSigmaDeuteron=pulls[0]; c.trk.nSigmaTriton=pulls[1];
      c.trk.nSigmaHe3=pulls[2]; c.trk.nSigmaHe4=pulls[3];
      c.trk.tofMatch=tof; c.trk.mass2=m2;
      qa.selected=best==type;
      (*pidQa)[c.speciesKey].push_back(qa);
    }
    if(best>=0) {
      ++multiplicity[best];
      FemtoCandidate c; c.eventIndex=eventIndex; c.source=kFemtoCandTrack;
      c.speciesKey=kSpecies[best]; c.charge=t->charge()*(best>=2?2:1);
      c.SetP4(NuclearP4(pvec,best)); c.trk.trackIndex=it; c.trk.nHitsFit=t->nHitsFit();
      c.trk.nSigmaDeuteron=pulls[0]; c.trk.nSigmaTriton=pulls[1];
      c.trk.nSigmaHe3=pulls[2]; c.trk.nSigmaHe4=pulls[3];
      c.trk.tofMatch=tof; c.trk.mass2=m2; store[c.speciesKey].push_back(c);
      if (best==3 && He4TofQaEnabled()) {
        mNuclear->Fill("hDedxP_4He_selected",p,dedx);
        if(tof) mNuclear->Fill("hM2P_4He_selected",p,m2);
        mNuclear->Fill("hHe4PidStage",5.);
      }
      // Final target nucleus once per track, independent of Lambda/pair multiplicity.
      // NuclearP4 already applied the physical mass and (for He) Z=2 exactly once.
      FillNucleusAcceptance(c,best);
      const std::string vzKey=std::string("hVz_")+kLegacy[best];
      if(mNuclear->HasHistogram(vzKey.c_str())) mNuclear->Fill(vzKey.c_str(),dst->event()->primaryVertex().Z());
    }
    // QA population is independent of final best-species / p-over-q selection.
    if(tof) {
      mNuclear->Fill("hM2P",p,m2);
      for(int s=0;s<3;++s) if(m2>mTofMean[s]-mM2SigmaCut*mTofSigma[s] && m2<mTofMean[s]+mM2SigmaCut*mTofSigma[s])
        mNuclear->Fill((std::string("hDedxP_")+kLegacy[s]+"_m2").c_str(),p,dedx);
    }
    if(!tof || !(m2>mMinM2 && m2<mMaxM2 && p>mMinPTofQa)) mNuclear->Fill("hDedxP_cut",p,dedx);
  }
  for(int s=0;s<4;++s) {
    const std::string key=std::string("hMult_")+kLegacy[s];
    if(mNuclear->HasHistogram(key.c_str())) mNuclear->Fill(key.c_str(),multiplicity[s]);
    const std::string ck=key+"_CentBin"+Index(cent9);
    if(cent9>=0 && cent9<9 && mNuclear->HasHistogram(ck.c_str())) mNuclear->Fill(ck.c_str(),multiplicity[s]);
  }
}

bool FemtoLambdaLegacy::IsPairNucleus(int type) const {
  return std::find(mPairSpeciesIndices.begin(),mPairSpeciesIndices.end(),type)!=mPairSpeciesIndices.end();
}

FemtoCandidate FemtoLambdaLegacy::MapHe3AsHe4(const FemtoCandidate& nucleus) const {
  FemtoCandidate mapped=nucleus;
  TLorentzVector p4;p4.SetVectM(nucleus.P4().Vect(),mNuclearMass[3]);
  mapped.SetP4(p4);mapped.speciesKey="he3_as_he4";
  // The stored physical momentum already has Z=2; do not multiply it again.
  return mapped;
}

void FemtoLambdaLegacy::FillNucleusAcceptance(const FemtoCandidate& nucleus,int type) {
  if(!mRoot||!IsPairNucleus(type))return;
  mRoot->Fill((std::string("hNucleus_PtVsYLab_")+kSpecies[type]).c_str(),nucleus.P4().Rapidity(),nucleus.P4().Pt());
  if(mHe3ControlEnabled&&type==2) {
    const FemtoCandidate mapped=MapHe3AsHe4(nucleus);
    mRoot->Fill("hNucleus_PtVsYLab_he3_as_he4",mapped.P4().Rapidity(),mapped.P4().Pt());
  }
}

void FemtoLambdaLegacy::FillLambda(const FemtoLambdaCandidate& c,int cent9,double refmult,const TVector3& pv) {
  if(!mRoot) return;
  const FemtoCandidate& f=c.candidate;
  const TVector3 displacement=TVector3(c.x,c.y,c.z)-pv;
  TLorentzVector measured; measured.SetVectM(f.P4().Vect(),c.mass);
  mRoot->Fill("hLambda_InvMass",c.mass); mRoot->Fill("hLambda_Pt",f.pt);
  mRoot->Fill("hLambda_Eta",f.eta); mRoot->Fill("hLambda_Phi",f.phi);
  mRoot->Fill("hDCA12",c.daughterDistance); mRoot->Fill("hDCAV0",c.distanceToPv);
  mRoot->Fill("hCosPointing",c.cosPointing);
  mRoot->Fill("hNSigmaProton",c.protonPidPull); mRoot->Fill("hNSigmaPion",c.pionPidPull);
  mRoot->Fill("hLambda_InvMass_vs_Pt",f.pt,c.mass);
  mRoot->Fill("hLambda_InvMass_vs_DecayLength",displacement.Mag(),c.mass);
  mRoot->Fill("hLambda_InvMass_vs_TransDecayLength",displacement.Perp(),c.mass);
  mRoot->Fill("hLambda_InvMass_vs_Y",measured.Rapidity(),c.mass);
  mRoot->Fill("hDCAV0_vs_InvMass",c.mass,c.distanceToPv);
  mRoot->Fill("hCosPointing_vs_InvMass",c.mass,c.cosPointing);
  if(cent9>=0 && cent9<9) {
    mRoot->Fill(("hLambda_InvMass_CentBin"+Index(cent9)).c_str(),c.mass);
    mRoot->Fill("hLambda_InvMass_vs_Cent9",cent9,c.mass);
    if(refmult>=0) mRoot->Fill("hLambda_InvMass_vs_RefMultCorr",refmult,c.mass);
  }
  mRoot->Fill("hLambdaKF_DecayLength",c.decayLength);
  mRoot->Fill("hLambdaKF_DecayLengthSignificance",c.decayLengthSignificance);
  mRoot->Fill("hLambdaKF_Chi2Ndf",c.chi2Ndf); mRoot->Fill("hLambdaKF_TopoChi2Ndf",c.topoChi2Ndf);
}

bool FemtoLambdaLegacy::FillLambdaAcceptance(const FemtoLambdaCandidate& c,StPicoDst* dst) {
  if (MassRegion(c.mass) != 1) return true;
  if (!mRoot || !dst || c.protonIndex < 0 || c.pionIndex < 0 ||
      static_cast<unsigned>(c.protonIndex) >= dst->numberOfTracks() ||
      static_cast<unsigned>(c.pionIndex) >= dst->numberOfTracks() ||
      !dst->track(c.protonIndex) || !dst->track(c.pionIndex)) {
    std::cerr << "[FemtoLambdaLegacy] Invalid daughter source for acceptance QA" << std::endl;
    return false;
  }
  TLorentzVector proton,pion;
  // Secondary daughters: global-track momentum, not PV-constrained primary pMom.
  proton.SetVectM(dst->track(c.protonIndex)->gMom(),mDaughterProtonMass);
  pion.SetVectM(dst->track(c.pionIndex)->gMom(),mDaughterPionMass);
  const TLorentzVector lambda=c.candidate.P4();
  if (!std::isfinite(lambda.Rapidity()) || !std::isfinite(lambda.Pt()) ||
      !std::isfinite(proton.Rapidity()) || !std::isfinite(proton.Pt()) ||
      !std::isfinite(pion.Rapidity()) || !std::isfinite(pion.Pt())) {
    std::cerr << "[FemtoLambdaLegacy] Nonfinite lab kinematics for acceptance QA" << std::endl;
    return false;
  }
  // No CM shift. One entry per Lambda candidate; daughters can be shared by
  // distinct Lambda candidates, but are never multiplied by nuclear pair count.
  mRoot->Fill("hLambda_PtVsYLab_signal",lambda.Rapidity(),lambda.Pt());
  mRoot->Fill("hLambdaProton_PtVsYLab_signal",proton.Rapidity(),proton.Pt());
  mRoot->Fill("hLambdaPion_PtVsYLab_signal",pion.Rapidity(),pion.Pt());
  return true;
}

bool FemtoLambdaLegacy::FillPair(const FemtoCandidate& l,const FemtoCandidate& n,
    double kstar,double qlab,int cent9,bool mixed,StPicoDst* dst) {
  int type=-1;
  for(size_t i=0;i<mPairSpeciesIndices.size();++i)
    if(n.speciesKey==kSpecies[mPairSpeciesIndices[i]])type=mPairSpeciesIndices[i];
  if(!mRoot || !mNuclear || type<0 || l.speciesKey!="lambda" ||
     !std::isfinite(kstar) || !std::isfinite(qlab) || !std::isfinite(l.reso.invMass) ||
     cent9<0 || cent9>=9 || FemtoCandidatesShareTrack(l,n)) return false;
  if((mixed && l.eventIndex==n.eventIndex) || (!mixed && l.eventIndex!=n.eventIndex)) return false;
  FemtoCandidate mapped;double mappedKstar=0.;
  if(mHe3ControlEnabled&&type==2) {
    mapped=MapHe3AsHe4(n);
    mappedKstar=LambdaPairKstar(l.P4(),mapped.P4());
    if(!std::isfinite(mappedKstar))return false;
  }
  if(!FillPairHistograms(l,n,kstar,qlab,cent9,mixed,dst,kSpecies[type],kLegacy[type],type,type))return false;
  if(mHe3ControlEnabled&&type==2)
    return FillPairHistograms(l,mapped,mappedKstar,qlab,cent9,mixed,dst,"he3_as_he4","3HeAs4He",3,2);
  return true;
}

bool FemtoLambdaLegacy::FillPairHistograms(const FemtoCandidate& l,const FemtoCandidate& n,
    double kstar,double qlab,int cent9,bool mixed,StPicoDst* dst,
    const std::string& species,const std::string& legacySpecies,int massIndex,int selectedNuclearType) {
  // User-defined femtoscopic mT: average fixed mass and HALF the vector
  // sum of transverse momenta; neither total pair Mt nor scalar Pt sum.
  // Nuclear momenta already include He's charge correction.
  const TLorentzVector pair=l.P4()+n.P4();
  const double massSum=mLambdaPairMass+mNuclearMass[massIndex];
  const double mt2=0.25*(pair.Px()*pair.Px()+pair.Py()*pair.Py()+massSum*massSum);
  if(!std::isfinite(mt2) || mt2<0) return false;
  const double pairMt=std::sqrt(mt2);
  const std::string mode=mixed?"ME":"SE", prefix=mixed?"Mixed_":"";
  const std::string channel="lambda_"+species;
  mRoot->Fill(("hKstar"+mode+"_"+channel).c_str(),kstar);
  mRoot->Fill(("hKstar"+mode+"VsCent_"+channel).c_str(),kstar,cent9);
  mNuclear->Fill(("hKstarMass_"+prefix+legacySpecies+"_CentBin"+Index(cent9)).c_str(),kstar,l.reso.invMass);
  const int region=MassRegion(l.reso.invMass);
  if(!region) return true; // Full-mass output is intentionally independent of narrow windows.
  const char* canonical[]={"","_signal","_leftSB","_rightSB"};
  const char* legacy[]={"","","_SBNeg","_SBPos"};
  const std::string stem=prefix+legacySpecies+legacy[region];
  mNuclear->Fill(("hKstar_"+stem).c_str(),kstar);
  mNuclear->Fill(("hQlab_"+stem).c_str(),qlab);
  mNuclear->Fill(("hKstar_"+stem+"_CentBin"+Index(cent9)).c_str(),kstar);
  mNuclear->Fill(("hKstarMt_"+stem+"_CentBin"+Index(cent9)).c_str(),kstar,pairMt);
  mRoot->Fill(("hKstar"+mode+"_"+channel+canonical[region]).c_str(),kstar);
  mRoot->Fill(("hKstar"+mode+"VsCent_"+channel+canonical[region]).c_str(),kstar,cent9);
  if (region==1) {
    // Subdivide the SAME accepted signal pair; never rerun/multiply mixing.
    // Midpoint belongs only to the high half, so low+high is exactly signal.
    const std::string half=SignalHalf(l.reso.invMass)==1?"_signalLow":"_signalHigh";
    mRoot->Fill(("hKstar"+mode+"_"+channel+half).c_str(),kstar);
    mRoot->Fill(("hKstar"+mode+"VsCent_"+channel+half).c_str(),kstar,cent9);
    mNuclear->Fill(("hKstarMt_"+prefix+legacySpecies+half+"_CentBin"+Index(cent9)).c_str(),kstar,pairMt);
  }
  // Diagnostic-only: never reject a previously accepted pair for an invalid
  // new diagnostic, and never choose a nuclear species from the mass mapping.
  FillNuclearPidPairQa(l,n,kstar,cent9,mixed,legacySpecies,selectedNuclearType,false);
  if(!mixed && region==1 && dst) FillMerging(l,n,dst,legacySpecies,massIndex);
  return true;
}

bool FemtoLambdaLegacy::FillPrePidNuclearPairQa(const FemtoCandidate& lambda,
    const FemtoCandidate& nucleus,double kstar,int cent9,bool mixed) {
  if(!mNuclearPidPairQa || lambda.speciesKey!="lambda" || !std::isfinite(kstar) ||
     !std::isfinite(lambda.reso.invMass) || cent9<0 || cent9>=9 ||
     FemtoCandidatesShareTrack(lambda,nucleus) ||
     (mixed && lambda.eventIndex==nucleus.eventIndex) ||
     (!mixed && lambda.eventIndex!=nucleus.eventIndex)) return false;
  int type=-1;
  for(size_t i=0;i<mPairSpeciesIndices.size();++i)
    if(nucleus.speciesKey==kSpecies[mPairSpeciesIndices[i]]) type=mPairSpeciesIndices[i];
  if(type<2) return false;
  FillNuclearPidPairQa(lambda,nucleus,kstar,cent9,mixed,kLegacy[type],type,true);
  if(mHe3ControlEnabled && type==2) {
    const FemtoCandidate mapped=MapHe3AsHe4(nucleus);
    const double mappedKstar=LambdaPairKstar(lambda.P4(),mapped.P4());
    if(std::isfinite(mappedKstar))
      FillNuclearPidPairQa(lambda,mapped,mappedKstar,cent9,mixed,"3HeAs4He",2,true);
  }
  return true;
}

void FemtoLambdaLegacy::FillNuclearPidPairQa(const FemtoCandidate& lambda,
    const FemtoCandidate& nucleus,double kstar,int cent9,bool mixed,
    const std::string& legacySpecies,int selectedNuclearType,bool prePid) {
  if(!mNuclearPidPairQa || selectedNuclearType<2 || selectedNuclearType>3) return;
  const int region=MassRegion(lambda.reso.invMass);
  if(!region) return;
  // P4 already contains the physical Z=2 momentum (also for He3AsHe4).
  // This is a Float-P4-derived QA estimate, not a newly reconstructed raw track.
  const double rigidity=nucleus.P4().P()/2.;
  const double pull=selectedNuclearType==2?nucleus.trk.nSigmaHe3:nucleus.trk.nSigmaHe4;
  // The joint tuple and both dense views share one QA-only validity gate.
  // Never propagate a bad diagnostic back to physical pair acceptance.
  if(!std::isfinite(kstar) || !std::isfinite(rigidity) || !std::isfinite(pull)) return;
  if(mNuclearPidPairSparse) {
    const int representation=legacySpecies=="3HeAs4He"?2:(selectedNuclearType==2?1:0);
    const int massClass=region==2?0:(region==3?3:SignalHalf(lambda.reso.invMass));
    mNuclearPidPairSparse->Fill(kstar,rigidity,pull,cent9,representation,mixed,massClass,prePid);
  }
  const std::string prefix=mixed?"Mixed_":"";
  const char* regionSuffix[]={"","","_SBNeg","_SBPos"};
  const std::string stems[]={prefix+legacySpecies+regionSuffix[region],
    prefix+legacySpecies+(SignalHalf(lambda.reso.invMass)==1?"_signalLow":"_signalHigh")};
  const int count=region==1?2:1;
  for(int i=0;i<count;++i) {
    const std::string suffix=stems[i]+"_CentBin"+Index(cent9);
    const std::string stage=prePid?"PrePid":"";
    if(std::isfinite(pull)) mNuclearPidPairQa->Fill(("hKstarNSigma"+stage+"_"+suffix).c_str(),kstar,pull);
    if(std::isfinite(rigidity)) mNuclearPidPairQa->Fill(("hKstarPOverQ"+stage+"_"+suffix).c_str(),kstar,rigidity);
  }
}

void FemtoLambdaLegacy::FillMerging(const FemtoCandidate& l,const FemtoCandidate& n,StPicoDst* dst,
    const std::string& legacySpecies,int massIndex) {
  const int ip=l.reso.dau1Index, in=n.trk.trackIndex;
  if(ip<0 || in<0 || ip>=int(dst->numberOfTracks()) || in>=int(dst->numberOfTracks())) return;
  StPicoTrack* p=dst->track(ip); StPicoTrack* a=dst->track(in);
  if(!p || !a || p->pPt()<=0 || a->pPt()<=0) return;
  const TVector3 pp=p->pMom(), pn=a->pMom();
  const double z=massIndex>=2?2.:1., deta=pp.Eta()-pn.Eta();
  // Legacy numerical convention retained separately: positive charges, raw rigidity,
  // fixed configured field. Only its original ID-as-index bug is corrected.
  const double old=pp.Phi()-pn.Phi()+CurvatureAngle(1.,mFieldTesla,mRadiusMeters,pp.Pt())-
    CurvatureAngle(z,mFieldTesla,mRadiusMeters,pn.Pt());
  mNuclear->Fill(("hDphiDeta_proton_"+legacySpecies).c_str(),TVector2::Phi_mpi_pi(old),deta);
  // Corrected QA uses event kG -> T and signed physical charge / physical pT.
  const double field=dst->event()->bField()*0.1;
  const double corrected=pp.Phi()-pn.Phi()+CurvatureAngle(p->charge(),field,mRadiusMeters,pp.Pt())-
    CurvatureAngle(a->charge()*z,field,mRadiusMeters,pn.Pt()*z);
  mNuclear->Fill(("hDphiDeta_proton_"+legacySpecies+"_EventField").c_str(),TVector2::Phi_mpi_pi(corrected),deta);
}

bool FemtoLambdaLegacy::Write(TDirectory* out) {
  if(!out || !mRoot || !mNuclear) return false;
  TDirectory* previous=gDirectory;
  out->cd();
  bool ok=true;
  std::ostringstream mtDescription;
  mtDescription.precision(17);
  mtDescription << "average_mass: sqrt(|pT_Lambda+pT_nucleus|^2/4+(m_Lambda+m_nucleus)^2/4); "
                << "fixed masses [GeV/c^2]: Lambda=" << mLambdaPairMass
                << ", nucleus=" << mNuclearMass[mSpeciesIndex]
                << "; X=kstar [GeV/c], Y=m_T [GeV/c^2]; not total pair transverse mass";
  TNamed mtDefinition("FemtoLambdaPairMtDefinition",mtDescription.str().c_str());
  if(mtDefinition.Write()<=0) ok=false;
  TNamed pidDefinition("FemtoNuclearDedxCalibration",NuclearDedxCalibrationMetadata().c_str());
  if(pidDefinition.Write()<=0)ok=false;
  if(!mPidDomainMetadata.empty()) {
    TNamed domainDefinition("FemtoNuclearPidDomain",mPidDomainMetadata.c_str());
    if(domainDefinition.Write()<=0)ok=false;
  }
  if(mNuclearPidPairQa) {
    const std::string description=
      "selected nuclear pairs: same event/track veto, Lambda mass regions, cent9 and SE/ME pairs as existing CF; "
      "X=kstar [GeV/c]; hKstarNSigma Y=selected source nuclear TPC pull; physical He3 and He3As4He both use nSigmaHe3, He4 uses nSigmaHe4; "
      "hKstarPOverQ Y=|physical Float-P4 momentum|/2 [GeV/c], physical Z=2 correction applied exactly once before storage; "
      "mapped He3As4He changes bachelor mass/energy and kstar only, not pull or rigidity; "
      "signal is +/-3sigma, signalLow and signalHigh partition the same signal pairs, SBPos/SBNeg retain existing windows; "
      "PrePid keys: same quality, target rigidity-domain and TOF gates, no TPC pull or best-species cut; "
      "pre-PID payload never admits events to the physical CF pool; QA ME uses only the same stored physical events, recent depth and both directions; "
      "histogram axes and under/overflow are stored; joint kstar-rigidity-nSigma in hNuclearPidPairJoint; no purity fit; "
      "unit raw-count TH2 fills use Poisson errors unless ROOT allocated Sumw2; sparse stores explicit Sumw2; "
      "a nonfinite QA tuple skips all new QA views only, not existing pair acceptance or histograms";
    TNamed definition("FemtoLambdaNuclearPairPidQaDefinition",description.c_str());
    if(definition.Write()<=0)ok=false;
    if(!mNuclearPidPairSparse || !mNuclearPidPairSparse->Write(out))ok=false;
  }
  if(mHe3ControlEnabled) {
    for(int representation=0;representation<2;++representation) {
      const int type=representation?3:2;
      const std::string species=representation?"he3_as_he4":"he3";
      std::ostringstream definition;definition.precision(17);
      definition<<"average_mass: sqrt(|pT_Lambda+pT_nucleus|^2/4+(m_Lambda+m_nucleus)^2/4); "
                <<"fixed masses [GeV/c^2]: Lambda="<<mLambdaPairMass<<", nucleus="<<mNuclearMass[type]
                <<"; namespace=lambda_"<<species<<"; physical He3 momentum unchanged (Z=2 already applied)";
      TNamed targetDefinition(("FemtoLambdaPairMtDefinition_"+species).c_str(),definition.str().c_str());
      if(targetDefinition.Write()<=0)ok=false;
    }
    std::ostringstream description;description.precision(17);
    description<<"selected physical He3: lambda_he3 / legacy 3He; mapped mass hypothesis: lambda_he3_as_he4 / legacy 3HeAs4He; "
               <<"same accepted SE/ME pairs, event provenance, Lambda mass windows and physical Z=2 momentum; "
               <<"recompute bachelor energy using He4 mass="<<mNuclearMass[3]<<" instead of He3 mass="<<mNuclearMass[2]
               <<"; no extra collection/mixing loop or pair cut; raw counts only, no inferred contamination coefficient";
    TNamed controlDefinition("FemtoLambdaHe3ControlDefinition",description.str().c_str());
    if(controlDefinition.Write()<=0)ok=false;
  }
  for(size_t i=0;i<mRootKeys.size();++i) if(mRoot->Get(mRootKeys[i].c_str())->Write()<=0) ok=false;
  TDirectory* se=out->GetDirectory("true"); if(!se) se=out->mkdir("true");
  TDirectory* me=out->GetDirectory("mix"); if(!me) me=out->mkdir("mix");
  if(!se || !me) ok=false;
  else for(size_t i=0;i<mNuclearKeys.size();++i) {
    const std::string& key=mNuclearKeys[i];
    (key.find("_Mixed_")!=std::string::npos?me:se)->cd();
    if(mNuclear->Get(key.c_str())->Write()<=0) ok=false;
  }
  if(se && me && mNuclearPidPairQa) for(size_t i=0;i<mNuclearPidPairQaKeys.size();++i) {
    const std::string& key=mNuclearPidPairQaKeys[i];
    (key.find("_Mixed_")!=std::string::npos?me:se)->cd();
    if(mNuclearPidPairQa->Get(key.c_str())->Write()<=0)ok=false;
  }
  if(previous) previous->cd();
  return ok;
}
