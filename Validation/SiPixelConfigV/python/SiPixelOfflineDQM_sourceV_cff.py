import FWCore.ParameterSet.Config as cms

# Pixel Digi Monitoring
from Validation.SiPixelDigisV.SiPixelDigisV_cfi import *
# Hits
from Validation.SiPixelHitsV.SiPixelHitsV_cfi import *
# RecHit (clusters)
from Validation.SiPixelRecHitsV.SiPixelRecHitsV_cfi import *
# Clusters ontrack/offtrack (also general tracks)
from Validation.SiPixelTrackClustersV.SiPixelTrackClustersV_cfi import *

PerModule.enabled = False

siPixelOfflineDQM_sourceV = cms.Sequence(SiPixelDigisAnalyzerV
                                            + SiPixelHitsAnalyzerV
                                            + SiPixelRecHitsAnalyzerV
                                            + SiPixelTrackClustersAnalyzerV
                                            )
