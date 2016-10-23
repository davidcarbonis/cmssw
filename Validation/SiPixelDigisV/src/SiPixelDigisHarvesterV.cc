// -*- C++ -*-
//
// Package:     SiPixelDigisHarvesterV
// Class:       SiPixelDigisHarvesterV
//

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

#include "Validation/SiPixelDigisV/interface/SiPixelDigisV.h"
#include "FWCore/Framework/interface/MakerMacros.h"

SiPixelDigisHarvesterV::SiPixelDigisHarvesterV(const edm::ParameterSet& iConfig) :
  SiPixelPhase1Harvester(iConfig) 
{}

DEFINE_FWK_MODULE(SiPixelDigisHarvesterV);
