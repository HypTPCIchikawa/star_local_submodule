// StLambdaV0Reconstruction.cxx
// Lambda (p+ pi-) V0 reconstruction utilities shared across StMakers.
// Ported from StLambdaMaker (MakeHelix, MakeLambdaHelix) into a standalone
// namespace, following the pattern of StPhiKKReconstruction.cxx.

#include "StLambdaV0Reconstruction.h"
#include "ConfigManager.h"
#include "cuts/LambdaCutConfig.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"

#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <utility>

namespace {
const Double_t kProtonMass = 0.938272;
const Double_t kPionMass   = 0.139570;
const Double_t kLambdaMass = 1.115683;
}  // namespace

namespace StLambdaV0Reconstruction {

//-----------------------------------------------------------------------------
Double_t ProtonMass() { return kProtonMass; }
Double_t PionMass()   { return kPionMass;   }
Double_t LambdaMass() { return kLambdaMass; }

//-----------------------------------------------------------------------------
StPhysicalHelixD BuildHelix(const LambdaV0TrackState& trk) {
  StThreeVectorF gmomSt(trk.momentumX, trk.momentumY, trk.momentumZ);
  StThreeVectorF orgSt (trk.originX,   trk.originY,   trk.originZ);
  return StPhysicalHelixD(gmomSt, orgSt,
                          trk.BField * units::kilogauss,
                          static_cast<float>(trk.charge));
}

//-----------------------------------------------------------------------------
void FillFromPico(LambdaV0TrackState& state, StPicoTrack* trk, Double_t bField) {
  if (!trk) return;
  state.momentumX = trk->gMom().X();
  state.momentumY = trk->gMom().Y();
  state.momentumZ = trk->gMom().Z();
  state.originX   = trk->origin().X();
  state.originY   = trk->origin().Y();
  state.originZ   = trk->origin().Z();
  state.BField     = static_cast<Float_t>(bField);
  state.charge     = static_cast<Short_t>(trk->charge());
  state.trackIndex = trk->id();
}

//-----------------------------------------------------------------------------
Bool_t ReconstructV0(const LambdaV0TrackState& proton,
                     const LambdaV0TrackState& pion,
                     const TVector3& primaryVertex,
                     LambdaV0Candidate& cand) {
  LambdaCutConfig& lam = ConfigManager::GetInstance().GetLambdaCuts();

  StPhysicalHelixD hp  = BuildHelix(proton);
  StPhysicalHelixD hpi = BuildHelix(pion);

  std::pair<Double_t, Double_t> s = hp.pathLengths(hpi);
  if (TMath::Abs(s.first)  > lam.maxPathLength) return kFALSE;
  if (TMath::Abs(s.second) > lam.maxPathLength) return kFALSE;

  StThreeVectorD dcaA = hp.at(s.first);
  StThreeVectorD dcaB = hpi.at(s.second);
  StThreeVectorD v0_((dcaA.x() + dcaB.x()) * 0.5,
                     (dcaA.y() + dcaB.y()) * 0.5,
                     (dcaA.z() + dcaB.z()) * 0.5);
  Double_t dca12 = (dcaA - dcaB).mag();

  if (dca12 < 0 || dca12 > lam.maxDaughterDCA) return kFALSE;

  StThreeVectorD pp  = hp.momentumAt (s.first,  proton.BField * units::kilogauss);
  StThreeVectorD ppi = hpi.momentumAt(s.second, pion.BField   * units::kilogauss);

  TVector3 momP (pp.x(),  pp.y(),  pp.z());
  TVector3 momPi(ppi.x(), ppi.y(), ppi.z());
  TVector3 v0   (v0_.x(), v0_.y(), v0_.z());

  TVector3 pLam = momP + momPi;
  Double_t pLamMag = pLam.Mag();
  if (pLamMag < 1e-5) return kFALSE;

  // DCA of Lambda candidate to primary vertex
  TVector3 pLamUnit = pLam * (1.0 / pLamMag);
  TVector3 diff     = primaryVertex - v0;
  Double_t dcaV0    = (diff.Cross(pLamUnit)).Mag();
  if (dcaV0 > lam.maxDCAV0) return kFALSE;

  // Pointing angle
  TVector3 flight = v0 - primaryVertex;
  Double_t cosPoint = flight.Dot(pLam) / (flight.Mag() * pLam.Mag() + 1e-10);
  if (cosPoint < lam.minCosPointing) return kFALSE;

  // Invariant mass
  TLorentzVector lp, lpi;
  lp.SetVectM (momP,  kProtonMass);
  lpi.SetVectM(momPi, kPionMass);
  TLorentzVector lLambda = lp + lpi;

  // Fill output candidate
  cand.v0          = v0;
  cand.momLambda   = pLam;
  cand.momP        = momP;
  cand.momPi       = momPi;
  cand.dca12       = dca12;
  cand.dcaV0       = dcaV0;
  cand.cosPoint    = cosPoint;
  cand.invMass     = lLambda.M();
  cand.decayLength = flight.Mag();
  cand.protonId    = proton.trackIndex;
  cand.pionId      = pion.trackIndex;

  return kTRUE;
}

}  // namespace StLambdaV0Reconstruction
