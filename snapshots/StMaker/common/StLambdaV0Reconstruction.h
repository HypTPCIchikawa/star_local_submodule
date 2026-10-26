#ifndef ST_LAMBDA_V0_RECONSTRUCTION_H
#define ST_LAMBDA_V0_RECONSTRUCTION_H

// StLambdaV0Reconstruction.h
// Lambda (p+ pi-) V0 reconstruction utilities shared across StMakers.
// Pattern follows StPhiKKReconstruction.h.

#include "Rtypes.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "TVector3.h"
class StPicoTrack;

/** Shared track state for Lambda V0 reconstruction (helix + kinematics). */
struct LambdaV0TrackState {
  Float_t momentumX;
  Float_t momentumY;
  Float_t momentumZ;
  Float_t originX;
  Float_t originY;
  Float_t originZ;
  Float_t BField;   // [kG]
  Short_t charge;
  Int_t   trackIndex;  // track id for self-correlation cut

  LambdaV0TrackState()
      : momentumX(0.f), momentumY(0.f), momentumZ(0.f),
        originX(0.f),   originY(0.f),   originZ(0.f),
        BField(0.f), charge(0), trackIndex(-1) {}
};

/** Result of a single Lambda candidate. */
struct LambdaV0Candidate {
  TVector3 v0;        // decay vertex position
  TVector3 momLambda; // Lambda 3-momentum (momP + momPi)
  TVector3 momP;      // proton 3-momentum at V0
  TVector3 momPi;     // pion   3-momentum at V0
  Double_t dca12;     // DCA between daughters at V0
  Double_t dcaV0;     // DCA of Lambda to primary vertex
  Double_t cosPoint;  // cos(pointing angle)
  Double_t invMass;   // invariant mass (p + pi-)
  Double_t decayLength;
  Int_t    protonId;  // proton  track index
  Int_t    pionId;    // pion    track index
};

namespace StLambdaV0Reconstruction {

Double_t ProtonMass();
Double_t PionMass();
Double_t LambdaMass();

/** Build a helix from a LambdaV0TrackState. */
StPhysicalHelixD BuildHelix(const LambdaV0TrackState& trk);

/**
 * Fill a LambdaV0TrackState from a StPicoTrack.
 * @param state   Output track state.
 * @param trk     Input pico track.
 * @param bField  Magnetic field [kG].
 */
void FillFromPico(LambdaV0TrackState& state, StPicoTrack* trk, Double_t bField);

/**
 * Attempt to reconstruct a Lambda from one proton candidate and one pion candidate.
 * Geometric cuts (maxPathLength, maxDaughterDCA, maxDCAV0, minCosPointing) are
 * read from LambdaCutConfig via ConfigManager.
 *
 * @param proton      Proton track state.
 * @param pion        Pion track state.
 * @param primaryVertex  Primary vertex position.
 * @param cand        [out] Filled on success.
 * @return kTRUE if a candidate passes all geometric cuts.
 */
Bool_t ReconstructV0(const LambdaV0TrackState& proton,
                     const LambdaV0TrackState& pion,
                     const TVector3& primaryVertex,
                     LambdaV0Candidate& cand);

}  // namespace StLambdaV0Reconstruction

#endif  // ST_LAMBDA_V0_RECONSTRUCTION_H
