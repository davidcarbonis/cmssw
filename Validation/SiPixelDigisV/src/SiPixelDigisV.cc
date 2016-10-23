// -*- C++ -*-
//
// Package:    SiPixelDigisV
// Class:      SiPixelDigisV
//

// Original Author: Marcel Schneider

#include "Validation/SiPixelDigisV/interface/SiPixelDigisV.h"
// Additional Authors: Alexander Morton - modifying code for validation use

// C++ stuff
#include <iostream>

// CMSSW stuff
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// DQM Stuff
#include "DQMServices/Core/interface/MonitorElement.h"

SiPixelDigisV::SiPixelDigisV(const edm::ParameterSet& iConfig) :
  SiPixelPhase1Base(iConfig)
{
  srcToken_ = consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("src"));
} 

void SiPixelDigisV::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::DetSetVector<PixelDigi>> input;
  iEvent.getByToken(srcToken_, input);
  if (!input.isValid()) return; 

  edm::DetSetVector<PixelDigi>::const_iterator it;
  for (it = input->begin(); it != input->end(); ++it) {
    for(PixelDigi const& digi : *it) {
      histo[ADC].fill((double) digi.adc(), DetId(it->detId()), &iEvent);
      histo[ROW].fill((double) digi.row(), DetId(it->detId()), &iEvent);
      histo[COLUMN].fill((double) digi.column(), DetId(it->detId()), &iEvent);
      histo[NDIGIS].fill(DetId(it->detId()), &iEvent); // count
    }
    histo[DEBUG].fill(geometryInterface.extract(geometryInterface.intern("PXLadder"), DetId(it->detId())), DetId(it->detId()));
  }
  histo[NDIGIS].executePerEventHarvesting();
}

DEFINE_FWK_MODULE(SiPixelDigisV);

