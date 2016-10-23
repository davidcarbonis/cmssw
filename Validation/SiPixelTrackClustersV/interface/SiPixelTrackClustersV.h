#ifndef SiPixelTrackClustersV_h 
#define SiPixelTrackClustersV_h 
// -*- C++ -*-
// 
// Package:     SiPixelTrackClustersV
// Class  :     SiPixelTrackClustersV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

class SiPixelTrackClustersV : public SiPixelPhase1Base {
  enum {
    CHARGE,
    SIZE_X,
    SIZE_Y,
  };

  public:
  explicit SiPixelTrackClustersV(const edm::ParameterSet& conf);
  void analyze(const edm::Event&, const edm::EventSetup&);

  private:
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > clustersToken_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> trackAssociationToken_;
};

#endif
