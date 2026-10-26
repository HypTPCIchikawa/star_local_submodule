#include "FemtoLambdaLegacy.h"
#include "FemtoLambdaProvider.h"
#include "HistManager.h"
#include "NuclearIdDeDxVsMom.h"
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
}

FemtoLambdaLegacy::FemtoLambdaLegacy()
    : mRoot(0), mNuclear(0), mSpeciesIndex(-1), mMean(0), mWindow(0),
      mOuterFactor(0), mFieldTesla(0), mRadiusMeters(0),
      mDaughterProtonMass(0), mDaughterPionMass(0), mLambdaPairMass(0), mMinNHitsDedx(0),
      mMinPt(0), mNSigmaFill(0), mNSigmaExclude(0), mMaxNSigma(0), mM2SigmaCut(0),
      mMinPTofQa(0), mMinM2(0), mMaxM2(0), mMaxRigidity(0), mM2Selection(false) {}
FemtoLambdaLegacy::~FemtoLambdaLegacy() { delete mRoot; delete mNuclear; }

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
    int nuclei = 0;
    while (names >> name) {
      if (name == "lambda") continue;
      for (int s = 0; s < 4; ++s) if (name == kSpecies[s]) {
        ++nuclei; mSpeciesIndex = s; mSpecies = name; mLegacySpecies = kLegacy[s];
      }
    }
    if (nuclei != 1) throw std::runtime_error("Lambda legacy profile requires exactly one nuclear species");
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
    for (int i=0;i<Required<int>(maker,"nChannels");++i) {
      const std::string prefix = "channel_" + Index(i) + "_";
      const std::string channel = Required<std::string>(maker,(prefix+"name").c_str());
      int region=-1;
      for(int r=0;r<4;++r) if(channel=="lambda_"+mSpecies+regionSuffix[r]) region=r;
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
    if (!ValidateHistograms(diagnostics)) return false;
    diagnostics << "[FemtoLambdaLegacy] mainconf(argument)=" << mainPath
                << "\n  nuclear=" << nuclearPath << "\n  rootHist=" << rootPath
                << "\n  nuclearHist=" << histPath << "\n  target=" << mSpecies
                << " signal=[" << mMean-mWindow << "," << mMean+mWindow
                << "] outerFactor=" << mOuterFactor << std::endl;
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
    Need(*mRoot,"hNucleus_PtVsYLab_"+mSpecies);
    const char* common[]={"hDedxP","hDedxP_cut","hDedxP_e","hDedxP_pi","hDedxP_K","hDedxP_p","hDedxP_else","hM2P"};
    for (unsigned i=0;i<sizeof(common)/sizeof(common[0]);++i) Need(*mNuclear,common[i]);
    for (int s=0;s<4;++s) {
      Need(*mNuclear,std::string("hDedxP_")+kLegacy[s]);
      for (int q=0;q<3;++q) {
        const char* vars[]={"Eta","Y","Pt"};
        Need(*mNuclear,std::string("hPvs")+vars[q]+"_"+kLegacy[s]);
      }
      if (s<3) Need(*mNuclear,std::string("hDedxP_")+kLegacy[s]+"_m2");
    }
    const TH2* mtReference=0;
    for (int mixed=0;mixed<2;++mixed) {
      const std::string prefix=mixed?"Mixed_":"";
      const std::string mode=mixed?"ME":"SE";
      const char* suffixes[]={"","_SBPos","_SBNeg"};
      for(int r=0;r<3;++r) {
        const std::string stem=prefix+mLegacySpecies+suffixes[r];
        Need(*mNuclear,"hKstar_"+stem); Need(*mNuclear,"hQlab_"+stem);
        for(int c=0;c<9;++c) Need(*mNuclear,"hKstar_"+stem+"_CentBin"+Index(c));
      }
      const char* mtSuffixes[]={"","_SBPos","_SBNeg","_signalLow","_signalHigh"};
      for(int r=0;r<5;++r) for(int c=0;c<9;++c) {
        const std::string cent="_CentBin"+Index(c);
        const std::string source="hKstar_"+prefix+mLegacySpecies+(r<3?mtSuffixes[r]:"")+cent;
        NeedPairMt(*mNuclear,"hKstarMt_"+prefix+mLegacySpecies+mtSuffixes[r]+cent,
                   *mNuclear->Get(source.c_str()),mtReference);
      }
      for(int c=0;c<9;++c) Need(*mNuclear,"hKstarMass_"+prefix+mLegacySpecies+"_CentBin"+Index(c));
      const char* channels[]={"","_signal","_leftSB","_rightSB","_signalLow","_signalHigh"};
      for(int r=0;r<6;++r) {
        Need(*mRoot,"hKstar"+mode+"_lambda_"+mSpecies+channels[r]);
        Need(*mRoot,"hKstar"+mode+"VsCent_lambda_"+mSpecies+channels[r]);
      }
    }
    Need(*mNuclear,"hDphiDeta_proton_"+mLegacySpecies);
    Need(*mNuclear,"hDphiDeta_proton_"+mLegacySpecies+"_EventField");
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
  if (mMaxRigidity > 0 && !(p < mMaxRigidity)) return -1;
  if (mM2Selection && !tof) return -1;
  int best=-1; double minPull=999.;
  for(int s=0;s<4;++s) {
    bool pass=std::fabs(pulls[s])<mMaxNSigma;
    if(mM2Selection) {
      const double lo=s<3?mTofMean[s]-mM2SigmaCut*mTofSigma[s]:mMinM2;
      const double hi=s<3?mTofMean[s]+mM2SigmaCut*mTofSigma[s]:mMaxM2;
      pass=pass && m2>lo && m2<hi;
    }
    if(pass && std::fabs(pulls[s])<minPull) { best=s; minPull=std::fabs(pulls[s]); }
  }
  return best;
}
TLorentzVector FemtoLambdaLegacy::NuclearP4(const TVector3& raw, int type) const {
  TLorentzVector out;
  if(type>=0 && type<4) out.SetVectM(raw*(type>=2?2.:1.),mNuclearMass[type]);
  return out;
}

void FemtoLambdaLegacy::CollectNuclei(StPicoDst* dst,int eventIndex,int cent9,FemtoCandidateStore& store) {
  if(!dst || !dst->event() || !mNuclear) return;
  int multiplicity[4]={0,0,0,0};
  for(unsigned it=0;it<dst->numberOfTracks();++it) {
    StPicoTrack* t=dst->track(it);
    if(!t || t->nHitsDedx()<mMinNHitsDedx || t->gPt()<mMinPt) continue;
    const TVector3 pvec=t->pMom();
    const double p=pvec.Mag(), dedx=t->dEdx();
    if(dedx<=0) continue;
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
      pulls[s]=NuclearIdDeDxVsMom::GetNSigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),p,dedx);
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
      if(b && b->btofBeta()>0) { const double beta=b->btofBeta(); tof=true; m2=p*p*(1./(beta*beta)-1.); }
    }
    const int best=SelectNuclearType(p,pulls,tof,m2);
    if(best>=0) {
      ++multiplicity[best];
      FemtoCandidate c; c.eventIndex=eventIndex; c.source=kFemtoCandTrack;
      c.speciesKey=kSpecies[best]; c.charge=t->charge()*(best>=2?2:1);
      c.SetP4(NuclearP4(pvec,best)); c.trk.trackIndex=it; c.trk.nHitsFit=t->nHitsFit();
      c.trk.nSigmaDeuteron=pulls[0]; c.trk.nSigmaTriton=pulls[1];
      c.trk.nSigmaHe3=pulls[2]; c.trk.nSigmaHe4=pulls[3];
      c.trk.tofMatch=tof; c.trk.mass2=m2; store[c.speciesKey].push_back(c);
      // Final target nucleus once per track, independent of Lambda/pair multiplicity.
      // NuclearP4 already applied the physical mass and (for He) Z=2 exactly once.
      if (best == mSpeciesIndex)
        mRoot->Fill(("hNucleus_PtVsYLab_"+mSpecies).c_str(),c.P4().Rapidity(),c.P4().Pt());
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
  if(!mRoot || !mNuclear || n.speciesKey!=mSpecies || l.speciesKey!="lambda" ||
     !std::isfinite(kstar) || !std::isfinite(qlab) || !std::isfinite(l.reso.invMass) ||
     cent9<0 || cent9>=9 || FemtoCandidatesShareTrack(l,n)) return false;
  if((mixed && l.eventIndex==n.eventIndex) || (!mixed && l.eventIndex!=n.eventIndex)) return false;
  // User-defined femtoscopic mT: average fixed mass and HALF the vector
  // sum of transverse momenta; neither total pair Mt nor scalar Pt sum.
  // Nuclear momenta already include He's charge correction.
  const TLorentzVector pair=l.P4()+n.P4();
  const double massSum=mLambdaPairMass+mNuclearMass[mSpeciesIndex];
  const double mt2=0.25*(pair.Px()*pair.Px()+pair.Py()*pair.Py()+massSum*massSum);
  if(!std::isfinite(mt2) || mt2<0) return false;
  const double pairMt=std::sqrt(mt2);
  const std::string mode=mixed?"ME":"SE", prefix=mixed?"Mixed_":"";
  const std::string channel="lambda_"+mSpecies;
  mRoot->Fill(("hKstar"+mode+"_"+channel).c_str(),kstar);
  mRoot->Fill(("hKstar"+mode+"VsCent_"+channel).c_str(),kstar,cent9);
  mNuclear->Fill(("hKstarMass_"+prefix+mLegacySpecies+"_CentBin"+Index(cent9)).c_str(),kstar,l.reso.invMass);
  const int region=MassRegion(l.reso.invMass);
  if(!region) return true; // Full-mass output is intentionally independent of narrow windows.
  const char* canonical[]={"","_signal","_leftSB","_rightSB"};
  const char* legacy[]={"","","_SBNeg","_SBPos"};
  const std::string stem=prefix+mLegacySpecies+legacy[region];
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
    mNuclear->Fill(("hKstarMt_"+prefix+mLegacySpecies+half+"_CentBin"+Index(cent9)).c_str(),kstar,pairMt);
  }
  if(!mixed && region==1 && dst) FillMerging(l,n,dst);
  return true;
}

void FemtoLambdaLegacy::FillMerging(const FemtoCandidate& l,const FemtoCandidate& n,StPicoDst* dst) {
  const int ip=l.reso.dau1Index, in=n.trk.trackIndex;
  if(ip<0 || in<0 || ip>=int(dst->numberOfTracks()) || in>=int(dst->numberOfTracks())) return;
  StPicoTrack* p=dst->track(ip); StPicoTrack* a=dst->track(in);
  if(!p || !a || p->pPt()<=0 || a->pPt()<=0) return;
  const TVector3 pp=p->pMom(), pn=a->pMom();
  const double z=mSpeciesIndex>=2?2.:1., deta=pp.Eta()-pn.Eta();
  // Legacy numerical convention retained separately: positive charges, raw rigidity,
  // fixed configured field. Only its original ID-as-index bug is corrected.
  const double old=pp.Phi()-pn.Phi()+CurvatureAngle(1.,mFieldTesla,mRadiusMeters,pp.Pt())-
    CurvatureAngle(z,mFieldTesla,mRadiusMeters,pn.Pt());
  mNuclear->Fill(("hDphiDeta_proton_"+mLegacySpecies).c_str(),TVector2::Phi_mpi_pi(old),deta);
  // Corrected QA uses event kG -> T and signed physical charge / physical pT.
  const double field=dst->event()->bField()*0.1;
  const double corrected=pp.Phi()-pn.Phi()+CurvatureAngle(p->charge(),field,mRadiusMeters,pp.Pt())-
    CurvatureAngle(a->charge()*z,field,mRadiusMeters,pn.Pt()*z);
  mNuclear->Fill(("hDphiDeta_proton_"+mLegacySpecies+"_EventField").c_str(),TVector2::Phi_mpi_pi(corrected),deta);
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
  for(size_t i=0;i<mRootKeys.size();++i) if(mRoot->Get(mRootKeys[i].c_str())->Write()<=0) ok=false;
  TDirectory* se=out->GetDirectory("true"); if(!se) se=out->mkdir("true");
  TDirectory* me=out->GetDirectory("mix"); if(!me) me=out->mkdir("mix");
  if(!se || !me) ok=false;
  else for(size_t i=0;i<mNuclearKeys.size();++i) {
    const std::string& key=mNuclearKeys[i];
    (key.find("_Mixed_")!=std::string::npos?me:se)->cd();
    if(mNuclear->Get(key.c_str())->Write()<=0) ok=false;
  }
  if(previous) previous->cd();
  return ok;
}
