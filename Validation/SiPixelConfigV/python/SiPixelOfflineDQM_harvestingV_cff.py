import FWCore.ParameterSet.Config as cms

from  Validation.SiPixelConfigV.SiPixelOfflineDQM_sourceV_cff import *

siPixelOfflineDQM_harvestingV = cms.Sequence(SiPixelDigisHarvesterV
                                                + SiPixelRecHitsHarvesterV
                                                + SiPixelHitsHarvesterV
                                                + SiPixelRecHitsHarvesterV
                                                + SiPixelTrackClustersHarvesterV
                                                )
