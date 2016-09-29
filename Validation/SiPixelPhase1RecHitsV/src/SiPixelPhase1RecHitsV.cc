// -*- C++ -*-
//
// Package:     SiPixelPhase1RecHitsV
// Class:       SiPixelPhase1RecHitsV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "Validation/SiPixelPhase1RecHitsV/interface/SiPixelPhase1RecHitsV.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"


SiPixelPhase1RecHitsV::SiPixelPhase1RecHitsV(const edm::ParameterSet& iConfig) :
  SiPixelPhase1BaseV(iConfig) 
{
  srcToken_ = consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("src"));
}

void SiPixelPhase1RecHitsV::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<SiPixelRecHitCollection> input;
  iEvent.getByToken(srcToken_, input);
  if (!input.isValid()) return;

  SiPixelRecHitCollection::const_iterator it;
  for (it = input->begin(); it != input->end(); ++it) {
    auto id = DetId(it->detId());

    for(SiPixelRecHit const& rechit : *it) {
      SiPixelRecHit::ClusterRef const& clust = rechit.cluster();
      int sizeX = (*clust).sizeX();
      int sizeY = (*clust).sizeY();

      LocalPoint lp = rechit.localPosition();
      float rechit_x = lp.x();
      float rechit_y = lp.y();
      
      LocalError lerr = rechit.localPositionError();
      float lerr_x = sqrt(lerr.xx());
      float lerr_y = sqrt(lerr.yy());

      histo[NRECHITS].fill(id, &iEvent);

      histo[CLUST_X].fill(sizeX, id, &iEvent);
      histo[CLUST_Y].fill(sizeY, id, &iEvent);

      histo[ERROR_X].fill(lerr_x, id, &iEvent);
      histo[ERROR_Y].fill(lerr_y, id, &iEvent);

      histo[POS].fill(rechit_x, rechit_y, id, &iEvent);
    }
  }

  histo[NRECHITS].executePerEventHarvesting();
}

DEFINE_FWK_MODULE(SiPixelPhase1RecHitsV);

