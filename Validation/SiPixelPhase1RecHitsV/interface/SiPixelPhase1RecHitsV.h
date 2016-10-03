#ifndef SiPixelPhase1RecHitsV_h 
#define SiPixelPhase1RecHitsV_h 
// -*- C++ -*-
// 
// Package:     SiPixelPhase1RecHitsV
// Class  :     SiPixelPhase1RecHitsV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "Validation/SiPixelPhase1CommonV/interface/SiPixelPhase1BaseV.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

class SiPixelPhase1RecHitsV : public SiPixelPhase1BaseV {
  enum {
    NRECHITS,
    CLUST_X,
    CLUST_Y,
    ERROR_X,
    ERROR_Y,
    POS
  };

  public:
  explicit SiPixelPhase1RecHitsV(const edm::ParameterSet& conf);
  void analyze(const edm::Event&, const edm::EventSetup&);

  private:
  TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  edm::EDGetTokenT<SiPixelRecHitCollection> srcToken_;
};

#endif
