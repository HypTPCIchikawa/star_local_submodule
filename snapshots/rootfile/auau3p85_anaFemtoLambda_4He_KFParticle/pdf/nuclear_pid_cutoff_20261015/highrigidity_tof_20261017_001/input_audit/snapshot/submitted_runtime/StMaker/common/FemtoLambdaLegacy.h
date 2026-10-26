#ifndef FEMTO_LAMBDA_LEGACY_H
#define FEMTO_LAMBDA_LEGACY_H

#include "FemtoCandidate.h"
#include <iosfwd>
#include <string>
#include <vector>

class HistManager;
class StPicoDst;
class TDirectory;
class TVector3;
struct FemtoLambdaCandidate;

// Lambda-only compatibility support. It owns its histograms; no KF/SIMD linkage.
// StFemtoMaker owns event acceptance, centrality, pair calculation and mixing.
class FemtoLambdaLegacy {
 public:
  FemtoLambdaLegacy();
  ~FemtoLambdaLegacy();
  bool Init(const char* mainconf, std::ostream& diagnostics);
  HistManager* RootHistograms() const { return mRoot; }
  HistManager* NuclearHistograms() const { return mNuclear; }
  const std::string& SpeciesKey() const { return mSpecies; }
  void CollectNuclei(StPicoDst*, int eventIndex, int cent9, FemtoCandidateStore&);
  void FillLambda(const FemtoLambdaCandidate&, int cent9, double refmultcorr,
                  const TVector3& primaryVertex);
  // Signal-candidate acceptance, not weighted by the number of nuclear pairs.
  bool FillLambdaAcceptance(const FemtoLambdaCandidate&, StPicoDst*);
  bool FillPair(const FemtoCandidate& lambda, const FemtoCandidate& nucleus,
                double kstar, double qlab, int cent9, bool mixed,
                StPicoDst* currentDst = 0);
  bool Write(TDirectory* destination);

  // 0=outside windows, 1=signal, 2=left sideband, 3=right sideband.
  int MassRegion(double mass) const;
  // 0=not signal, 1=low half (mass<mean), 2=high half (mass>=mean).
  int SignalHalf(double mass) const;
  int SelectNuclearType(double rigidity, const double pulls[4],
                        bool validTof, double mass2) const;
  TLorentzVector NuclearP4(const TVector3& rawRigidity, int type) const;

 private:
  FemtoLambdaLegacy(const FemtoLambdaLegacy&);
  FemtoLambdaLegacy& operator=(const FemtoLambdaLegacy&);
  void FillMerging(const FemtoCandidate&, const FemtoCandidate&, StPicoDst*);
  bool ValidateHistograms(std::ostream&);
  HistManager* mRoot;
  HistManager* mNuclear;
  std::vector<std::string> mRootKeys, mNuclearKeys;
  std::string mSpecies, mLegacySpecies;
  int mSpeciesIndex;
  double mMean, mWindow, mOuterFactor, mFieldTesla, mRadiusMeters;
  double mDaughterProtonMass, mDaughterPionMass, mLambdaPairMass;
  int mMinNHitsDedx;
  double mMinPt, mNuclearMass[4], mTofMean[3], mTofSigma[3];
  double mNSigmaFill, mNSigmaExclude, mMaxNSigma, mM2SigmaCut;
  double mMinPTofQa, mMinM2, mMaxM2, mMaxRigidity;
  bool mM2Selection;
};

#endif
