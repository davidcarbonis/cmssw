import FWCore.ParameterSet.Config as cms

PilotBladeOccupancyFilter = cms.EDFilter("PilotBladeOccupancyFilter",
        digis = cms.InputTag("siPixelDigis"), 
        clusters = cms.InputTag("siPixelClusters"), 
        minclusters = cms.int32(20),
        mindigis = cms.int32(50)
)
