import FWCore.ParameterSet.Config as cms
from Validation.SiPixelPhase1CommonV.HistogramManager_cfi import *

SiPixelPhase1RecHitsNRecHits = DefaultHisto.clone(
  name = "rechits",
  title = "RecHits",
  range_min = 0, range_max = 10, range_nbins = 10,
  xlabel = "number of rechits",
  dimensions = 0,
  topFolderName = "PixelPhase1V/RecHits",
  specs = cms.VPSet(
#    StandardSpecificationTrend_Num,
#    StandardSpecification2DProfile_Num
    Specification(PerLadder).groupBy("PXBarrel|PXForward/Shell|HalfCylinder/PXLayer|PXDisk/PXRing|/PXLadder|PXBlade") # per-ladder and profiles
                            .save(),
    Specification(PerLayer1D).groupBy(parent("PXBarrel|PXForward/Shell|HalfCylinder/PXLayer|PXDisk/PXRing|/PXLadder|PXBlade")) # per-layer
                             .save(),
    Specification(PerModule).groupBy("PXBarrel|PXForward/PXLayer|PXDisk/DetId").save()
  )
)

SiPixelPhase1RecHitsClustX = DefaultHisto.clone(
  name = "rechitsize_x",
  title = "X size of RecHit clusters",
  range_min = 0, range_max = 10, range_nbins = 10,
  xlabel = "RecHit size X dimension",
  dimensions = 1,
  topFolderName = "PixelPhase1V/RecHits",
  specs = cms.VPSet(
#    StandardSpecification2DProfile
    StandardSpecifications1D
  )
)

SiPixelPhase1RecHitsClustY = SiPixelPhase1RecHitsClustX.clone(
  name = "rechitsize_y",
  title = "Y size of RecHit clusters",
  xlabel = "RecHit size Y dimension"
)

SiPixelPhase1RecHitsErrorX = DefaultHisto.clone(
  name = "rechiterror_x",
  title = "RecHit Error in X-direction",
  range_min = 0, range_max = 0.02, range_nbins = 100,
  xlabel = "X error",
  dimensions = 1,
  topFolderName = "PixelPhase1V/RecHits",
  specs = cms.VPSet(
    StandardSpecifications1D
#    StandardSpecification2DProfile
  )
)

SiPixelPhase1RecHitsErrorY = SiPixelPhase1RecHitsErrorX.clone(
  name = "rechiterror_y",
  title = "RecHit Error in Y-direction",
  xlabel = "Y error"
)

SiPixelPhase1RecHitsPosition = DefaultHisto.clone(
  enabled = False,
  name = "rechit_pos",
  title = "Position of RecHits on Module",
  range_min   = -1, range_max   = 1, range_nbins   = 100,
  range_y_min = -4, range_y_max = 4, range_y_nbins = 100,
  xlabel = "x offset",
  ylabel = "y offset",
  dimensions = 2,
  topFolderName = "PixelPhase1V/RecHits",
  specs = cms.VPSet(
    Specification(PerModule).groupBy("PXBarrel|PXForward/PXLayer|PXDisk/DetId").save(),
  )
)

SiPixelPhase1RecHitsConf = cms.VPSet(
  SiPixelPhase1RecHitsNRecHits,
  SiPixelPhase1RecHitsClustX,
  SiPixelPhase1RecHitsClustY,
  SiPixelPhase1RecHitsErrorX,
  SiPixelPhase1RecHitsErrorY,
  SiPixelPhase1RecHitsPosition,
)

SiPixelPhase1RecHitsAnalyzerV = cms.EDAnalyzer("SiPixelPhase1RecHitsV",
        src = cms.InputTag("siPixelRecHits"),
        histograms = SiPixelPhase1RecHitsConf,
        geometry = SiPixelPhase1Geometry
)

SiPixelPhase1RecHitsHarvesterV = cms.EDAnalyzer("SiPixelPhase1HarvesterV",
        histograms = SiPixelPhase1RecHitsConf,
        geometry = SiPixelPhase1Geometry
)
