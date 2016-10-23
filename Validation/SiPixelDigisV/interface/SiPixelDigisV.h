#ifndef SiPixelDigisV_h // Can we use #pragma once?
#define SiPixelDigisV_h
// -*- C++ -*-
//
// Package:     SiPixelDigisV
// Class  :     SiPixelDigisV
// 

// Original Author: Marcel Schneider
// Additional Authors: Alexander Morton - modifying code for validation use

// Input data stuff
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

// PixelDQM Framework
#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"

class SiPixelDigisV : public SiPixelPhase1Base {
  // List of quantities to be plotted. 
  enum {
    ADC, // digi ADC readouts
    NDIGIS, // number of digis per event and module
    ROW, // number of digis per row
    COLUMN, // number of digis per column
    DEBUG, // geometry debugging

    MAX_HIST // a sentinel that gives the number of quantities (not a plot).
  };
  public:
  explicit SiPixelDigisV(const edm::ParameterSet& conf);

  void analyze(const edm::Event&, const edm::EventSetup&) ;

  private:
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> srcToken_;

};

class SiPixelDigisHarvesterV : public SiPixelPhase1Harvester {
  enum {
    ADC, // digi ADC readouts
    NDIGIS, // number of digis per event and module
    ROW, // number of digis per row
    COLUMN, // number of digis per column
    DEBUG, // geometry debugging

    MAX_HIST
  };
  public:
  explicit SiPixelDigisHarvesterV(const edm::ParameterSet& conf);

};

#endif
