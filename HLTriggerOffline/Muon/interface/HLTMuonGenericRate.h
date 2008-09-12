#ifndef HLTriggerOffline_Muon_HLTMuonGenericRate_H
#define HLTriggerOffline_Muon_HLTMuonGenericRate_H

/** \class HLTMuonGenericRate
 *  Get L1/HLT efficiency/rate plots
 *
 *  \author  M. Vander Donckt, J. Klukas  (copied from J. Alcaraz)
 */

// Base Class Headers

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <vector>
#include "TFile.h"
#include "TNtuple.h"


class HLTMuonGenericRate {

public:
  /// Constructor
  HLTMuonGenericRate(const edm::ParameterSet& pset, int triggerIndex);

  /// Destructor
  virtual ~HLTMuonGenericRate();

  // Operations

  void analyze          (const edm::Event & iEvent);
  void BookHistograms   ();
  void WriteHistograms  ();
  void SetCurrentFolder (TString folder );
  MonitorElement* BookIt(TString name, TString title, 
			 int Nbins, float Min, float Max);

private:

  // Input from cfg file
  edm::InputTag theL1CollectionLabel;
  edm::InputTag theGenLabel;
  edm::InputTag theRecoLabel;
  std::vector <edm::InputTag> theHLTCollectionLabels;

  double theL1ReferenceThreshold;
  double theHLTReferenceThreshold;
  double theLuminosity;
  double thePtMin;
  double thePtMax;
  double theMinPtCut;
  double theMaxEtaCut;

  std::vector<double> theNSigmas;
  unsigned int theNumberOfObjects;
  unsigned int theNbins;
  int thisEventWeight;
  int motherParticleId;
  bool useMuonFromGenerator;
  bool useMuonFromReco;

  // Struct for matching
  struct MatchStruct {
    const reco::GenParticle *genCand;
    const reco::Track*       recCand;
    const l1extra::L1MuonParticle* l1Cand;
    std::vector<const reco::RecoChargedCandidate*> hltCands;
  };
  std::vector<MatchStruct> genMatches;
  std::vector<MatchStruct> recMatches;
  int nL1Orphans ;
  int nHltOrphans;
  
  const reco::Candidate* findMother(const reco::Candidate*);

  // Histograms
  DQMStore* dbe_;

  std::vector <MonitorElement*> hPtPassGen ;
  std::vector <MonitorElement*> hEtaPassGen;
  std::vector <MonitorElement*> hPhiPassGen;
  std::vector <MonitorElement*> hPtPassRec ;
  std::vector <MonitorElement*> hEtaPassRec;
  std::vector <MonitorElement*> hPhiPassRec;

  MonitorElement *NumberOfEvents  ;
  MonitorElement *NumberOfL1Events;
  int theNumberOfEvents ;
  int theNumberOfL1Events;
  std::string theRootFileName;

  int findGenMatch( double eta, double phi, double maxDeltaR );
  int findRecMatch( double eta, double phi, double maxdeltaR );


  // ntuple
  bool    makeNtuple;
  TNtuple *theNtuple;
  TFile   *theFile;
  float   ntParams[18];

};
#endif
