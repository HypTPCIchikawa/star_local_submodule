// StFemtoMakerLambdaNuclear.cxx
// Single StMaker integrating:
//   1. Lambda (p+ pi-) V0 reconstruction via StLambdaV0Reconstruction (common/)
//   2. Nuclear candidate PID via StNuclearIdHelper (common/)
//   3. Same-event femto correlations (signal + sidebands)
//   4. Mixed-event femto correlations (event-mixing pool)
//
// Design follows StFemtoMaker.cxx.
// Histogram names are species-independent (no "_d", "_t", etc.) so that
// the same YAML hist config can be reused for d, t, 3He, 4He.

#include "StLambdaNuclearMaker.h"
#include "StLambdaV0Reconstruction.h"
#include "StNuclearIdHelper.h"
#include "ConfigManager.h"
#include "HistManager.h"
#include "CentralityHelper.h"
#include "cuts/EventCutConfig.h"
#include "cuts/LambdaCutConfig.h"
#include "cuts/NuclearIdCutConfig.h"
#include "cuts/CentralityCutConfig.h"
#include "cuts/MixingConfig.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>

namespace {
const Double_t kLambdaMass = 1.115683;  // PDG Lambda mass [GeV/c^2]
}  // namespace

//=============================================================================
// Constructor / Destructor
//=============================================================================

StFemtoMakerLambdaNuclear::StFemtoMakerLambdaNuclear(
    const char* name,
    StPicoDstMaker* picoMaker,
    const char* outName,
    NuclearSpecies targetSpecies)
    : StMaker(name),
      mPicoDstMaker(picoMaker),
      mPicoDst(0),
      mOutName(outName),
      mTargetSpecies(targetSpecies),
      mEventCounter(0),
      m_histManager(0),
      m_centrality(0),
      m_cent9(-1),
      m_cent16(-1),
      m_refMultCorr(-1.0),
      m_centWeight(1.0),
      m_centralityPercent(-1.0),
      m_bufferSize(150),
      m_maxMixEvents(100) {}

StFemtoMakerLambdaNuclear::~StFemtoMakerLambdaNuclear() {
  if (m_centrality) { delete m_centrality; m_centrality = 0; }
  if (m_histManager) { delete m_histManager; m_histManager = 0; }
}

//=============================================================================
// Factory functions
//=============================================================================

StFemtoMakerLambdaNuclear* createStFemtoMakerLambdaNuclear(
    const char* name, StPicoDstMaker* picoMaker, const char* outName,
    NuclearSpecies targetSpecies) {
  return new StFemtoMakerLambdaNuclear(name, picoMaker, outName, targetSpecies);
}

extern "C" void* createStFemtoMakerLambdaNuclearC(
    const char* name, void* picoMaker, const char* outName, int targetSpecies) {
  return (void*)createStFemtoMakerLambdaNuclear(
      name, (StPicoDstMaker*)picoMaker, outName,
      (NuclearSpecies)targetSpecies);
}

//=============================================================================
// Init
//=============================================================================

Int_t StFemtoMakerLambdaNuclear::Init() {
  // Load histogram config via the Maker name key (e.g. "femtoLambdaNuclear")
  std::string histPath = ConfigManager::GetInstance().GetHistConfigPath(GetName());
  if (histPath.empty()) {
    std::cerr << "[StFemtoMakerLambdaNuclear] GetHistConfigPath(\"" << GetName()
              << "\") returned empty; no histograms will be filled." << std::endl;
    m_histManager = 0;
  } else {
    m_histManager = new HistManager();
    if (!m_histManager->LoadFromFile(histPath.c_str())) {
      std::cerr << "[StFemtoMakerLambdaNuclear] Failed to load hist config from "
                << histPath << std::endl;
      delete m_histManager;
      m_histManager = 0;
    }
  }

  m_centrality = new CentralityHelper();
  if (!m_centrality->Init(ConfigManager::GetInstance().GetCentralityCuts())) {
    std::cerr << "[StFemtoMakerLambdaNuclear] CentralityHelper init failed" << std::endl;
  }

  // mixing buffer size from YAML
  const MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  if (mix.bufferSize > 0) m_bufferSize = mix.bufferSize;

  return kStOK;
}

//=============================================================================
// Clear (called each event by StChain)
//=============================================================================

void StFemtoMakerLambdaNuclear::Clear(Option_t* opt) {
  StMaker::Clear(opt);
  mLambdaCands.clear();
  mNuclearCands.clear();
}

//=============================================================================
// Event-level cuts (track multiplicity only; vertex cuts done in centrality)
//=============================================================================

Bool_t StFemtoMakerLambdaNuclear::PassEventCuts(Int_t nTracks) {
  EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
  if (ev.maxNTr > 0 && nTracks > ev.maxNTr) return kFALSE;
  return kTRUE;
}

//=============================================================================
// Lambda daughter PID cuts
//=============================================================================

Bool_t StFemtoMakerLambdaNuclear::PassLambdaProtonCuts(
    StPicoTrack* trk, const TVector3& pVtx) {
  if (!trk || trk->charge() <= 0) return kFALSE;
  LambdaCutConfig& lam = ConfigManager::GetInstance().GetLambdaCuts();
  if (trk->nHitsFit() < lam.minNHitsFit) return kFALSE;
  if (trk->nHitsMax() > 0 &&
      (Double_t)trk->nHitsFit() / (Double_t)trk->nHitsMax() < lam.minNHitsRatio)
    return kFALSE;
  if (TMath::Abs(trk->nSigmaProton()) > lam.nSigmaProton) return kFALSE;
  Double_t dca = trk->gDCA(pVtx.X(), pVtx.Y(), pVtx.Z());
  if (dca < lam.minDCAProton) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoMakerLambdaNuclear::PassLambdaPionCuts(
    StPicoTrack* trk, const TVector3& pVtx) {
  if (!trk || trk->charge() >= 0) return kFALSE;
  LambdaCutConfig& lam = ConfigManager::GetInstance().GetLambdaCuts();
  if (trk->nHitsFit() < lam.minNHitsFit) return kFALSE;
  if (trk->nHitsMax() > 0 &&
      (Double_t)trk->nHitsFit() / (Double_t)trk->nHitsMax() < lam.minNHitsRatio)
    return kFALSE;
  if (TMath::Abs(trk->nSigmaPion()) > lam.nSigmaPion) return kFALSE;
  Double_t dca = trk->gDCA(pVtx.X(), pVtx.Y(), pVtx.Z());
  if (dca < lam.minDCAPion) return kFALSE;
  return kTRUE;
}

//=============================================================================
// Lambda V0 reconstruction (fills mLambdaCands)
//=============================================================================

void StFemtoMakerLambdaNuclear::ReconstructLambdas(
    const TVector3& pVtx, Double_t bField) {
  Int_t nTr = mPicoDst->numberOfTracks();
  for (Int_t ip = 0; ip < nTr; ip++) {
    StPicoTrack* p = mPicoDst->track(ip);
    if (!p || !PassLambdaProtonCuts(p, pVtx)) continue;

    for (Int_t ii = 0; ii < nTr; ii++) {
      if (ii == ip) continue;
      StPicoTrack* pi = mPicoDst->track(ii);
      if (!pi || !PassLambdaPionCuts(pi, pVtx)) continue;

      // Fill track states for V0 reconstruction
      LambdaV0TrackState stP, stPi;
      StLambdaV0Reconstruction::FillFromPico(stP,  p,  bField);
      StLambdaV0Reconstruction::FillFromPico(stPi, pi, bField);

      LambdaV0Candidate cand;
      if (!StLambdaV0Reconstruction::ReconstructV0(stP, stPi, pVtx, cand)) continue;

      LambdaCand lc;
      lc.momentum = cand.momLambda;
      lc.invMass  = cand.invMass;
      lc.protonId = cand.protonId;
      lc.pionId   = cand.pionId;
      mLambdaCands.push_back(lc);

      // Lambda QA histograms
      if (m_histManager) {
        m_histManager->Fill("hLambda_InvMass", cand.invMass);
        m_histManager->Fill("hLambda_Pt",      cand.momLambda.Pt());
        m_histManager->Fill("hLambda_Eta",     cand.momLambda.PseudoRapidity());
        m_histManager->Fill("hDCA12",          cand.dca12);
        m_histManager->Fill("hDCAV0",          cand.dcaV0);
        m_histManager->Fill("hCosPointing",    cand.cosPoint);
        m_histManager->Fill("hLambda_InvMass_vs_Pt",          cand.momLambda.Pt(), cand.invMass);
        m_histManager->Fill("hLambda_InvMass_vs_DecayLength", cand.decayLength,    cand.invMass);
      }
    }
  }
}

//=============================================================================
// Nuclear candidate identification (fills mNuclearCands)
//=============================================================================

void StFemtoMakerLambdaNuclear::IdentifyNuclei(const TVector3& /*pVtx*/) {
  NuclearIdCutConfig& nucCuts = ConfigManager::GetInstance().GetNuclearIdCuts();
  Int_t nTr = mPicoDst->numberOfTracks();

  const bool hasCentBinDedxP  = m_histManager && m_histManager->HasHistogram("hDedxP_CentBin0");
  const bool hasCentBinPvsEta = m_histManager && m_histManager->HasHistogram("hPvsEta_CentBin0");

  for (Int_t it = 0; it < nTr; it++) {
    StPicoTrack* trk = mPicoDst->track(it);
    if (!trk) continue;
    if (trk->nHitsDedx() < nucCuts.minNHitsDedxNuclear) continue;
    if (trk->gPt() < 0.1) continue;

    Double_t p = trk->pMom().Mag();
    Double_t dedx = trk->dEdx();
    if (dedx <= 0) continue;

    // Fill NuclearTrackState for helper
    NuclearTrackState state;
    StNuclearIdHelper::FillFromPico(state, trk, mPicoDst);

    // Inclusive QA
    if (m_histManager && m_histManager->HasHistogram("hDedxP")) {
      m_histManager->Fill("hDedxP", p, dedx);
      if (hasCentBinDedxP && m_cent9 >= 0 && m_cent9 <= 8) {
        m_histManager->Fill(Form("hDedxP_CentBin%d", m_cent9), p, dedx);
      }
      if (fabs(trk->nSigmaElectron()) < nucCuts.nSigmaFill1D) m_histManager->Fill("hDedxP_e", p, dedx);
      if (fabs(trk->nSigmaPion()) < nucCuts.nSigmaFill1D) m_histManager->Fill("hDedxP_pi", p, dedx);
      if (fabs(trk->nSigmaKaon()) < nucCuts.nSigmaFill1D) m_histManager->Fill("hDedxP_K", p, dedx);
      if (fabs(trk->nSigmaProton()) < nucCuts.nSigmaFill1D) m_histManager->Fill("hDedxP_p", p, dedx);

      if (fabs(trk->nSigmaElectron()) > nucCuts.nSigmaExclude &&
          fabs(trk->nSigmaPion()) > nucCuts.nSigmaExclude &&
          fabs(trk->nSigmaKaon()) > nucCuts.nSigmaExclude &&
          fabs(trk->nSigmaProton()) > nucCuts.nSigmaExclude) {
        m_histManager->Fill("hDedxP_else", p, dedx);
      }
    }

    if (state.tofMatch && state.beta > 0) {
      if (m_histManager && m_histManager->HasHistogram("hM2P")) {
        m_histManager->Fill("hM2P", p, state.mass2);
      }
      if (nucCuts.m2_selection && StNuclearIdHelper::PassTof(mTargetSpecies, state)) {
        if (m_histManager && m_histManager->HasHistogram("hDedxP_m2")) m_histManager->Fill("hDedxP_m2", p, dedx);
      }
    } else {
      if (m_histManager && m_histManager->HasHistogram("hDedxP_cut")) {
        m_histManager->Fill("hDedxP_cut", p, dedx);
      }
    }

    // Check only the target species
    Bool_t isTarget = kFALSE;
    switch (mTargetSpecies) {
      case kNucDeuteron: isTarget = StNuclearIdHelper::IsDeuteron(state); break;
      case kNucTriton:   isTarget = StNuclearIdHelper::IsTriton  (state); break;
      case kNucHe3:      isTarget = StNuclearIdHelper::IsHe3     (state); break;
      case kNucHe4:      isTarget = StNuclearIdHelper::IsHe4     (state); break;
      default:           isTarget = kFALSE;                                break;
    }
    if (!isTarget) continue;

    // Target Nucleus Accepted: Fill Kinematic QA
    TVector3 pMom = trk->pMom();
    TVector3 actMom = StNuclearIdHelper::EffectiveMomentum(pMom, mTargetSpecies);
    Double_t actP = actMom.Mag();

    TLorentzVector lv = StNuclearIdHelper::NuclearP4(pMom, mTargetSpecies);
    Double_t eta = pMom.PseudoRapidity(); // match original raw eta
    Double_t pT = actMom.Perp();
    Double_t y = lv.Rapidity();

    if (m_histManager && m_histManager->HasHistogram("hPvsEta")) {
      m_histManager->Fill("hPvsEta", actP, eta);
      m_histManager->Fill("hPvsY", actP, y);
      m_histManager->Fill("hPvsPt", actP, pT);
      if (hasCentBinPvsEta && m_cent9 >= 0 && m_cent9 <= 8) {
        m_histManager->Fill(Form("hPvsEta_CentBin%d", m_cent9), actP, eta);
        m_histManager->Fill(Form("hPvsY_CentBin%d", m_cent9), actP, y);
        m_histManager->Fill(Form("hPvsPt_CentBin%d", m_cent9), actP, pT);
      }
    }

    NuclearCand nc;
    nc.momentum = actMom;
    nc.trackId  = trk->id();
    mNuclearCands.push_back(nc);
  }
}

//=============================================================================
// k* and q_lab calculation helper
//=============================================================================

static void ComputeFemtoKinematics(const TVector3& pLam, const TVector3& pNuc,
                                   NuclearSpecies sp,
                                   Double_t& k_star, Double_t& q_lab) {
  Double_t nucMass = StNuclearIdHelper::SpeciesMass(sp);

  TLorentzVector lLam, lNuc;
  lLam.SetVectM(pLam, kLambdaMass);
  lNuc.SetVectM(pNuc, nucMass);

  // Lab-frame relative momentum
  q_lab = (pLam - pNuc).Mag();

  // Rest-frame relative momentum (k*)
  TLorentzVector pair = lLam + lNuc;
  TVector3 boostVec   = -pair.BoostVector();
  TLorentzVector lamRest = lLam;
  lamRest.Boost(boostVec);
  k_star = lamRest.P();
}

//=============================================================================
// Lambda signal / sideband classification
//=============================================================================

static void ClassifyLambda(Double_t invMass,
                           Bool_t& isSig, Bool_t& isSBPos, Bool_t& isSBNeg) {
  NuclearIdCutConfig& nucCuts = ConfigManager::GetInstance().GetNuclearIdCuts();
  Double_t delta  = invMass - nucCuts.MeanLambda;
  Double_t window = nucCuts.nsigmaLambda * nucCuts.sigmaLambda;

  isSig   = (TMath::Abs(delta) <= window);
  isSBPos = (!isSig && delta >  window  && delta <= 4.0 * window);
  isSBNeg = (!isSig && delta < -window  && delta >= -4.0 * window);
}

//=============================================================================
// Fill a femto pair into histograms (same-event or mixed-event)
//=============================================================================

void StFemtoMakerLambdaNuclear::FillFemtoPair(
    Double_t k_star, Double_t q_lab, Double_t invMass,
    Bool_t isSig, Bool_t isSBPos, Bool_t isSBNeg,
    Bool_t isMixed) {
  if (!m_histManager) return;

  // k* vs invMass (all pairs, for signal extraction QA)
  if (!isMixed && m_cent9 >= 0 && m_cent9 <= 8) {
    m_histManager->Fill(Form("hKstarMass_CentBin%d", m_cent9), k_star, invMass);
  }

  if (!isSig && !isSBPos && !isSBNeg) return;

  if (!isMixed) {
    // same-event
    if (isSig) {
      m_histManager->Fill("hKstar", k_star);
      m_histManager->Fill("hQlab",  q_lab);
      if (m_cent9 >= 0 && m_cent9 <= 8) {
        m_histManager->Fill(Form("hKstar_CentBin%d", m_cent9), k_star);
      }
    } else if (isSBPos) {
      m_histManager->Fill("hKstar_SBPos", k_star);
      m_histManager->Fill("hQlab_SBPos",  q_lab);
      if (m_cent9 >= 0 && m_cent9 <= 8) {
        m_histManager->Fill(Form("hKstar_SBPos_CentBin%d", m_cent9), k_star);
      }
    } else if (isSBNeg) {
      m_histManager->Fill("hKstar_SBNeg", k_star);
      m_histManager->Fill("hQlab_SBNeg",  q_lab);
      if (m_cent9 >= 0 && m_cent9 <= 8) {
        m_histManager->Fill(Form("hKstar_SBNeg_CentBin%d", m_cent9), k_star);
      }
    }
  } else {
    // mixed-event
    if (isSig) {
      if (m_histManager->HasHistogram("hKstar_Mixed")) {
        m_histManager->Fill("hKstar_Mixed", k_star);
        m_histManager->Fill("hQlab_Mixed",  q_lab);
      }
      if (m_cent9 >= 0 && m_cent9 <= 8) {
        TString hname = Form("hKstar_Mixed_CentBin%d", m_cent9);
        if (m_histManager->HasHistogram(hname.Data()))
          m_histManager->Fill(hname.Data(), k_star);
      }
    } else if (isSBPos) {
      if (m_histManager->HasHistogram("hKstar_Mixed_SBPos")) {
        m_histManager->Fill("hKstar_Mixed_SBPos", k_star);
        m_histManager->Fill("hQlab_Mixed_SBPos",  q_lab);
      }
    } else if (isSBNeg) {
      if (m_histManager->HasHistogram("hKstar_Mixed_SBNeg")) {
        m_histManager->Fill("hKstar_Mixed_SBNeg", k_star);
        m_histManager->Fill("hQlab_Mixed_SBNeg",  q_lab);
      }
    }
  }
}

//=============================================================================
// Same-event femto pairs
//=============================================================================

void StFemtoMakerLambdaNuclear::FillSameEventPairs() {
  for (size_t il = 0; il < mLambdaCands.size(); il++) {
    const LambdaCand& lc = mLambdaCands[il];
    Bool_t isSig, isSBPos, isSBNeg;
    ClassifyLambda(lc.invMass, isSig, isSBPos, isSBNeg);

    for (size_t in = 0; in < mNuclearCands.size(); in++) {
      const NuclearCand& nc = mNuclearCands[in];
      // Self-correlation cut: nucleus track cannot be a Lambda daughter
      if (nc.trackId == lc.protonId || nc.trackId == lc.pionId) continue;

      Double_t k_star, q_lab;
      ComputeFemtoKinematics(lc.momentum, nc.momentum, mTargetSpecies, k_star, q_lab);

      // k* vs invMass (for signal extraction, all pairs)
      if (m_histManager && m_cent9 >= 0 && m_cent9 <= 8) {
        m_histManager->Fill(Form("hKstarMass_CentBin%d", m_cent9), k_star, lc.invMass);
      }

      FillFemtoPair(k_star, q_lab, lc.invMass,
                    isSig, isSBPos, isSBNeg, /*isMixed=*/kFALSE);
    }
  }
}

//=============================================================================
// Mixed-event femto pairs
//=============================================================================

void StFemtoMakerLambdaNuclear::FillMixedEventPairs(Float_t vz) {
  Int_t mixBin = GetMixingBin(vz, m_cent9);
  std::deque<MixedEvent>& pool = m_mixingPool[mixBin];

  Int_t mixCount = 0;
  for (auto it = pool.rbegin(); it != pool.rend(); ++it) {
    if (mixCount >= m_maxMixEvents) break;
    const MixedEvent& mixEvt = *it;

    // current Lambda × pool Nucleus
    for (size_t il = 0; il < mLambdaCands.size(); il++) {
      const LambdaCand& lc = mLambdaCands[il];
      Bool_t isSig, isSBPos, isSBNeg;
      ClassifyLambda(lc.invMass, isSig, isSBPos, isSBNeg);

      for (size_t in = 0; in < mixEvt.nuclei.size(); in++) {
        Double_t k_star, q_lab;
        ComputeFemtoKinematics(lc.momentum, mixEvt.nuclei[in].momentum,
                               mTargetSpecies, k_star, q_lab);

        // mixed-event k* vs invMass
        if (m_histManager && m_cent9 >= 0 && m_cent9 <= 8) {
          TString hname = Form("hKstarMass_Mixed_CentBin%d", m_cent9);
          if (m_histManager->HasHistogram(hname.Data()))
            m_histManager->Fill(hname.Data(), k_star, lc.invMass);
        }

        FillFemtoPair(k_star, q_lab, lc.invMass,
                      isSig, isSBPos, isSBNeg, /*isMixed=*/kTRUE);
      }
    }

    // pool Lambda × current Nucleus
    for (size_t il = 0; il < mixEvt.lambdas.size(); il++) {
      const LambdaCand& lc = mixEvt.lambdas[il];
      Bool_t isSig, isSBPos, isSBNeg;
      ClassifyLambda(lc.invMass, isSig, isSBPos, isSBNeg);

      for (size_t in = 0; in < mNuclearCands.size(); in++) {
        Double_t k_star, q_lab;
        ComputeFemtoKinematics(lc.momentum, mNuclearCands[in].momentum,
                               mTargetSpecies, k_star, q_lab);
        FillFemtoPair(k_star, q_lab, lc.invMass,
                      isSig, isSBPos, isSBNeg, /*isMixed=*/kTRUE);
      }
    }

    mixCount++;
  }
}

//=============================================================================
// Store current event into the mixing pool
//=============================================================================

void StFemtoMakerLambdaNuclear::StoreEventForMixing(Float_t vz) {
  if (mLambdaCands.empty() && mNuclearCands.empty()) return;

  Int_t mixBin = GetMixingBin(vz, m_cent9);
  std::deque<MixedEvent>& pool = m_mixingPool[mixBin];

  MixedEvent evt;
  evt.lambdas = mLambdaCands;
  evt.nuclei  = mNuclearCands;
  pool.push_back(evt);

  if ((Int_t)pool.size() > m_bufferSize) pool.pop_front();
}

//=============================================================================
// Mixing bin index
//=============================================================================

Int_t StFemtoMakerLambdaNuclear::GetMixingBin(Float_t vz, Int_t cent9) const {
  const MixingConfig& mix  = ConfigManager::GetInstance().GetMixingConfig();
  const EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();

  Int_t vzBin   = 0;
  Int_t centBin = 0;

  if (mix.nVzBins > 0) {
    Double_t vzSpan = ev.maxVz - ev.minVz;
    if (vzSpan > 0) {
      vzBin = (Int_t)((vz - ev.minVz) / vzSpan * mix.nVzBins);
      if (vzBin < 0)            vzBin = 0;
      if (vzBin >= mix.nVzBins) vzBin = mix.nVzBins - 1;
    }
  }

  if (mix.nCentralityBins > 0 && cent9 >= 0) {
    centBin = cent9;
    if (centBin >= mix.nCentralityBins) centBin = mix.nCentralityBins - 1;
  }

  return vzBin + mix.nVzBins * centBin;
}

//=============================================================================
// Centrality event-level QA
//=============================================================================

void StFemtoMakerLambdaNuclear::FillCentralityEventQA(
    Int_t rawMult, Double_t refMultCorr,
    Int_t nTracks, Int_t nBTOFMatch,
    Int_t nLambda, Int_t nNuclear) {
  if (!m_histManager || m_cent9 < 0) return;
  const Double_t centX = (Double_t)m_cent9;
  m_histManager->Fill("hRawMult_vs_Cent9",    centX, (Double_t)rawMult);
  m_histManager->Fill("hNTracks_vs_Cent9",    centX, (Double_t)nTracks);
  m_histManager->Fill("hTofMatchMult_vs_Cent9", centX, (Double_t)nBTOFMatch);
  m_histManager->Fill("hNLambda_vs_Cent9",    centX, (Double_t)nLambda);
  m_histManager->Fill("hNNuclear_vs_Cent9",   centX, (Double_t)nNuclear);
  if (refMultCorr >= 0.0) {
    m_histManager->Fill("hRefMultCorr_vs_Cent9", centX, refMultCorr);
  }
}

//=============================================================================
// Make (main event loop body)
//=============================================================================

Int_t StFemtoMakerLambdaNuclear::Make() {
  if (!mPicoDstMaker) return kStWarn;
  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) return kStWarn;

  StPicoEvent* event = mPicoDst->event();
  if (!event) return kStWarn;

  mEventCounter++;

  TVector3 pVtx       = event->primaryVertex();
  Int_t    refMult    = event->refMult();
  const Int_t runId   = event->runId();
  const Int_t nBTOFMatch = event->nBTOFMatch();
  const Double_t vz   = pVtx.Z();

  CentralityCutConfig& centCfg = ConfigManager::GetInstance().GetCentralityCuts();
  Int_t rawMult = refMult;
  if (centCfg.enabled) {
    TString mode(centCfg.mode.c_str());
    mode.ToLower();
    if (mode == "fxtmult") rawMult = event->fxtMult();
  }

  // Reset per-event centrality
  m_cent9            = -1;
  m_cent16           = -1;
  m_refMultCorr      = -1.0;
  m_centWeight       = 1.0;
  m_centralityPercent= -1.0;

  if (m_histManager && m_histManager->HasHistogram("hRefMultVsNTOFMatch")) {
    m_histManager->Fill("hRefMultVsNTOFMatch", (Double_t)nBTOFMatch, (Double_t)rawMult);
  }

  CentralityRejectReason centReason = kCentralityOk;

  if (m_centrality && m_centrality->IsEnabled()) {
    if (!m_centrality->CheckBadRun(runId, centReason))               return kStOK;
    if (!m_centrality->CheckPileup(rawMult, nBTOFMatch, vz, centReason)) return kStOK;
    if (!m_centrality->ComputeBins(event, rawMult, vz,
                                   m_cent9, m_cent16, m_refMultCorr,
                                   m_centWeight, centReason))        return kStOK;

    if (m_histManager && centCfg.fillCentralityQA) {
      m_histManager->Fill("hRawMult",             (Double_t)rawMult);
      m_histManager->Fill("hCentralityRaw",       (Double_t)m_cent9);
      m_histManager->Fill("hRefMultCorr",         m_refMultCorr);
      m_histManager->Fill("hCentralityVsVz",      vz, (Double_t)m_cent9);
      m_histManager->Fill("hRefMultVsNTOFMatchAfter", (Double_t)nBTOFMatch, (Double_t)rawMult);
    }

    if (!m_centrality->AcceptCentBin(m_cent9, m_refMultCorr, centReason)) return kStOK;

    m_centralityPercent = CentralityHelper::Cent9ToPercentile(m_cent9);
    if (m_histManager) {
      const Double_t w = centCfg.useWeight ? m_centWeight : 1.0;
      TH1* hCent = m_histManager->Get("hCentrality");
      if (hCent) hCent->Fill((Double_t)m_cent9, w);
      TH1* hCent16 = m_histManager->Get("hCentrality16");
      if (hCent16) hCent16->Fill((Double_t)m_cent16, w);
    }
  }

  if (m_histManager) {
    m_histManager->Fill("hVz",      pVtx.Z());
    m_histManager->Fill("hRefMult", refMult);
  }

  Int_t nTr = mPicoDst->numberOfTracks();
  if (!PassEventCuts(nTr)) return kStOK;

  Double_t bField = event->bField();

  // 1. Reconstruct Lambda candidates
  ReconstructLambdas(pVtx, bField);

  // 2. Identify nuclear candidates (target species only)
  IdentifyNuclei(pVtx);

  // 3. Same-event femto correlations
  FillSameEventPairs();

  // 4. Mixed-event femto correlations
  FillMixedEventPairs((Float_t)vz);

  // 5. Store this event for future mixing
  StoreEventForMixing((Float_t)vz);

  // Centrality QA
  if (m_cent9 >= 0) {
    FillCentralityEventQA(rawMult, m_refMultCorr, nTr, nBTOFMatch,
                          (Int_t)mLambdaCands.size(),
                          (Int_t)mNuclearCands.size());
  }

  if (m_histManager && m_histManager->HasHistogram("hMult")) {
    m_histManager->Fill("hMult", (Double_t)mNuclearCands.size());
    if (m_cent9 >= 0 && m_cent9 <= 8) {
      m_histManager->Fill(Form("hMult_CentBin%d", m_cent9), (Double_t)mNuclearCands.size());
    }
  }

  if (m_histManager) m_histManager->Fill("hN", 0);
  return kStOK;
}

//=============================================================================
// Finish
//=============================================================================

Int_t StFemtoMakerLambdaNuclear::Finish() {
  if (mOutName != "" && m_histManager) {
    gSystem->mkdir(gSystem->DirName(mOutName.Data()), kTRUE);
    TFile* fout = new TFile(mOutName.Data(), "RECREATE");
    if (fout && !fout->IsZombie()) {
      // Write true-event (non-mixed) histograms to "true/" subdirectory.
      // HistManager::WriteContaining writes only histograms whose name
      // contains the substring; we reuse it with a prefix common to all
      // non-mixed hists by writing everything first in true/ and then
      // overwriting mixed in mix/.
      // Strategy: write ALL to true/, then write Mixed subset to mix/.
      TDirectory* dirTrue = fout->mkdir("true");
      dirTrue->cd();
      if (m_histManager) m_histManager->Write();

      // Write mixed-event histograms to "mix/" subdirectory
      TDirectory* dirMix = fout->mkdir("mix");
      dirMix->cd();
      if (m_histManager) m_histManager->WriteContaining("Mixed");

      fout->Close();
      delete fout;
      std::cout << "[StFemtoMakerLambdaNuclear] Saved histograms to "
                << mOutName << std::endl;
    } else {
      std::cerr << "[StFemtoMakerLambdaNuclear] ERROR: failed to open "
                << mOutName << std::endl;
    }
  }

  std::cout << "[StFemtoMakerLambdaNuclear] Finish(): processed "
            << mEventCounter << " events" << std::endl;

  if (m_centrality && m_centrality->IsEnabled()) {
    std::cout << "[StFemtoMakerLambdaNuclear] centrality summary:"
              << " ok="          << m_centrality->CountOk()
              << " badRun="      << m_centrality->CountBadRun()
              << " pileup="      << m_centrality->CountPileup()
              << " invalidCent=" << m_centrality->CountInvalidCent()
              << " binRejected=" << m_centrality->CountBinRejected() << std::endl;
    m_centrality->Finish();
  }

  return kStOK;
}

//=============================================================================
// WriteHistograms (called externally when outName was set to "" in constructor)
//=============================================================================

void StFemtoMakerLambdaNuclear::WriteHistograms() {
  if (m_histManager) m_histManager->Write();
}
