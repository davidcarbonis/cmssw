#ifndef SiPixelPhase1ClustersV_h 
#define SiPixelPhase1ClustersV_h 
// -*- C++ -*-
// 
// Package:     SiPixelPhase1ClustersV
// Class  :     SiPixelPhase1ClustersV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "Validation/SiPixelPhase1CommonV/interface/SiPixelPhase1BaseV.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

class SiPixelPhase1ClustersV : public SiPixelPhase1BaseV {
  enum {
    CHARGE,
    SIZE,
    NCLUSTERS,
    EVENTRATE,
    POSITION_B,
    POSITION_F,
    POSITION_XZ,
    POSITION_YZ,
    SIZE_VS_ETA
  };

  public:
  explicit SiPixelPhase1ClustersV(const edm::ParameterSet& conf);
  void analyze(const edm::Event&, const edm::EventSetup&);

  private:
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > srcToken_;
};

#endif
