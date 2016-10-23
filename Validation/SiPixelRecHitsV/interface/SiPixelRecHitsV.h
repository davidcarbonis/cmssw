#ifndef SiPixelRecHitsV_h 
#define SiPixelRecHitsV_h 
// -*- C++ -*-
// 
// Package:     SiPixelRecHitsV
// Class  :     SiPixelRecHitsV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

class SiPixelRecHitsV : public SiPixelPhase1Base {
  enum {
    IN_TIME_BUNCH,
    OUT_TIME_BUNCH,
    NSIMHITS,
    RECHIT_X,
    RECHIT_Y,
    RES_X,
    RES_Y,
    ERROR_X,
    ERROR_Y,
    PULL_X,
    PULL_Y,
  };

  public:
  explicit SiPixelRecHitsV(const edm::ParameterSet& conf);
  void analyze(const edm::Event&, const edm::EventSetup&);

  private:
  TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  edm::EDGetTokenT<SiPixelRecHitCollection> srcToken_;
};

#endif
