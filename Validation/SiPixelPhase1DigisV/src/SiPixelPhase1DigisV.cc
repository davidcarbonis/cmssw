// -*- C++ -*-
//
// Package:    SiPixelPhase1DigisV
// Class:      SiPixelPhase1DigisV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "Validation/SiPixelPhase1DigisV/interface/SiPixelPhase1DigisV.h"

// C++ stuff
#include <iostream>

// CMSSW stuff
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// DQM Stuff
#include "DQMServices/Core/interface/MonitorElement.h"

SiPixelPhase1DigisV::SiPixelPhase1DigisV(const edm::ParameterSet& iConfig) :
  SiPixelPhase1BaseV(iConfig)
{
  srcToken_ = consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("src"));
} 

void SiPixelPhase1DigisV::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::DetSetVector<PixelDigi>> input;
  iEvent.getByToken(srcToken_, input);
  if (!input.isValid()) return; 

  edm::DetSetVector<PixelDigi>::const_iterator it;
  for (it = input->begin(); it != input->end(); ++it) {
    for(PixelDigi const& digi : *it) {
      histo[ADC].fill((double) digi.adc(), DetId(it->detId()), &iEvent);
      histo[MAP].fill(DetId(it->detId()), &iEvent, digi.column(), digi.row()); 
      histo[NDIGIS    ].fill(DetId(it->detId()), &iEvent); // count
      histo[NDIGIS_FED].fill(DetId(it->detId()), &iEvent); 
    }
    histo[DEBUG].fill(geometryInterface.extract(geometryInterface.intern("PXLadder"), DetId(it->detId())), DetId(it->detId()));
  }
  histo[EVENT].fill(DetId(0), &iEvent);
  histo[NDIGIS    ].executePerEventHarvesting();
  histo[NDIGIS_FED].executePerEventHarvesting(); 
}

DEFINE_FWK_MODULE(SiPixelPhase1DigisV);

