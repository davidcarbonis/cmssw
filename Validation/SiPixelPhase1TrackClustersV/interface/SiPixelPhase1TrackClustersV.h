#ifndef SiPixelPhase1TrackClustersV_h 
#define SiPixelPhase1TrackClustersV_h 
// -*- C++ -*-
// 
// Package:     SiPixelPhase1TrackClustersV
// Class  :     SiPixelPhase1TrackClustersV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "Validation/SiPixelPhase1CommonV/interface/SiPixelPhase1BaseV.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

class SiPixelPhase1TrackClustersV : public SiPixelPhase1BaseV {
  enum {
    ONTRACK_CHARGE,
    ONTRACK_SIZE,
    ONTRACK_NCLUSTERS,
    ONTRACK_POSITION_B,
    ONTRACK_POSITION_F,

    OFFTRACK_CHARGE,
    OFFTRACK_SIZE,
    OFFTRACK_NCLUSTERS,
    OFFTRACK_POSITION_B,
    OFFTRACK_POSITION_F,

    NTRACKS,
    NTRACKS_VOLUME
  };

  public:
  explicit SiPixelPhase1TrackClustersV(const edm::ParameterSet& conf);
  void analyze(const edm::Event&, const edm::EventSetup&);

  private:
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > clustersToken_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> trackAssociationToken_;
};

#endif
