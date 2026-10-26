#ifndef StLambdaNuclearMaker_h
#define StLambdaNuclearMaker_h

// StFemtoMakerLambdaNuclear.hh
// Single StMaker integrating Lambda reconstruction, nuclear candidate ID,
// same-event and mixed-event femtoscopy for the Lambda x Nuclear system.
//
// Design pattern follows StFemtoMaker.h.
// Lambda V0 geometry is delegated to StLambdaV0Reconstruction (common/).
// Nuclear PID is delegated to StNuclearIdHelper (common/).
//
// Usage (in analysis macro):
//   StFemtoMakerLambdaNuclear* maker =
//       new StFemtoMakerLambdaNuclear("femtoLambdaNuclear",
//                                     picoMaker, outputFile,
//                                     kNucDeuteron);
//   chain->AddMaker(picoMaker);
//   chain->AddMaker(maker);

#include "StMaker.h"
#include "StNuclearIdHelper.h"

#include <deque>
#include <map>
#include <vector>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class TString;
class HistManager;
class CentralityHelper;
class TVector3;

class StFemtoMakerLambdaNuclear : public StMaker {
 public:
  StFemtoMakerLambdaNuclear(const char* name,
                             StPicoDstMaker* picoMaker,
                             const char* outName,
                             NuclearSpecies targetSpecies = kNucDeuteron);
  virtual ~StFemtoMakerLambdaNuclear();

  virtual Int_t  Init();
  virtual Int_t  Make();
  virtual void   Clear(Option_t* opt = "");
  virtual Int_t  Finish();

  void WriteHistograms();

 private:
  // ---- internal candidate structs -----------------------------------------

  struct LambdaCand {
    TVector3 momentum;
    Double_t invMass;
    Int_t    protonId;  // for self-correlation cut
    Int_t    pionId;    // for self-correlation cut
  };

  struct NuclearCand {
    TVector3 momentum;
    Int_t    trackId;   // for self-correlation cut
  };

  struct MixedEvent {
    std::vector<LambdaCand>  lambdas;
    std::vector<NuclearCand> nuclei;
  };

  // ---- data members -------------------------------------------------------

  StPicoDstMaker*  mPicoDstMaker;
  StPicoDst*       mPicoDst;
  TString          mOutName;
  NuclearSpecies   mTargetSpecies;   // which nucleus to correlate with Lambda
  Int_t            mEventCounter;

  HistManager*      m_histManager;
  CentralityHelper* m_centrality;
  Int_t             m_cent9;
  Int_t             m_cent16;
  Double_t          m_refMultCorr;
  Double_t          m_centWeight;
  Double_t          m_centralityPercent;

  // per-event candidate lists (cleared in Clear())
  std::vector<LambdaCand>  mLambdaCands;
  std::vector<NuclearCand> mNuclearCands;

  // mixing pool: key = mixing bin index
  std::map<Int_t, std::deque<MixedEvent> > m_mixingPool;
  Int_t m_bufferSize;
  Int_t m_maxMixEvents;

  // ---- private helpers ----------------------------------------------------

  Bool_t PassEventCuts(Int_t nTracks);
  Bool_t PassLambdaProtonCuts(StPicoTrack* trk, const TVector3& pVtx);
  Bool_t PassLambdaPionCuts  (StPicoTrack* trk, const TVector3& pVtx);

  void ReconstructLambdas(const TVector3& pVtx, Double_t bField);
  void IdentifyNuclei    (const TVector3& pVtx);

  void FillSameEventPairs ();
  void FillMixedEventPairs(Float_t vz);
  void StoreEventForMixing(Float_t vz);

  Int_t    GetMixingBin(Float_t vz, Int_t cent9) const;

  // Compute k* (rest-frame relative momentum) and q_lab, then fill histograms.
  void FillFemtoPair(Double_t k_star, Double_t q_lab,
                     Double_t invMass,
                     Bool_t isSig, Bool_t isSBPos, Bool_t isSBNeg,
                     Bool_t isMixed);

  // Centrality event-level QA histograms
  void FillCentralityEventQA(Int_t rawMult, Double_t refMultCorr,
                              Int_t nTracks, Int_t nBTOFMatch,
                              Int_t nLambda, Int_t nNuclear);
};

// Factory functions for CINT / shared-library loading
StFemtoMakerLambdaNuclear* createStFemtoMakerLambdaNuclear(
    const char* name, StPicoDstMaker* picoMaker, const char* outName,
    NuclearSpecies targetSpecies = kNucDeuteron);
extern "C" void* createStFemtoMakerLambdaNuclearC(
    const char* name, void* picoMaker, const char* outName, int targetSpecies);

#endif  // StFemtoMakerLambdaNuclear_hh
