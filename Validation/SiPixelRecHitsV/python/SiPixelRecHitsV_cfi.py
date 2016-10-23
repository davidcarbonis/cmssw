import FWCore.ParameterSet.Config as cms
from DQM.SiPixelPhase1Common.HistogramManager_cfi import *

SiPixelRecHitsInTimeEvents = DefaultHisto.clone(
  name = "in_time_bunch",
  title = "Events (in-time bunch)",
  range_min = 0, range_max = 10, range_nbins = 10,
  xlabel = "number of in-time rechits events",
  dimensions = 1,
  topFolderName = "PixelV/RecHits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward").save(),
  )
)

SiPixelRecHitsOutTimeEvents = DefaultHisto.clone(
  name = "out_time_bunch",
  title = "Events (out-time bunch)",
  range_min = 0, range_max = 10, range_nbins = 10,
  xlabel = "number of out-time rechit events",
  dimensions = 1,
  topFolderName = "PixelV/RecHits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward").save(),
  )
)


SiPixelRecHitsNSimHits = DefaultHisto.clone(
  name = "nsimhits",
  title = "SimHits",
  range_min = 0, range_max = 100, range_nbins = 100,
  xlabel = "sim hit event number in event",
  dimensions = 1,
  topFolderName = "PixelV/RecHits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk").save(),
  )
)

SiPixelRecHitsPosX = DefaultHisto.clone(
  name = "rechit_x",
  title = "X position of RecHits",
  range_min = -2., range_max = 2., range_nbins = 80,
  xlabel = "RecHit position X dimension",
  dimensions = 1,
  topFolderName = "PixelV/RecHits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward").save(),
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelRecHitsPosY = SiPixelRecHitsPosX.clone(
  name = "rechit_y",
  title = "Y position of RecHits",
  xlabel = "RecHit position Y dimension",
  range_min = -4., range_max = 4., range_nbins = 80,
)

SiPixelRecHitsResX = DefaultHisto.clone(
  name = "res_x",
  title = "X resolution of RecHits",
  range_min = -200., range_max = 200., range_nbins = 200,
  xlabel = "RecHit resolution X dimension",
  dimensions = 1,
  topFolderName = "PixelV/RecHits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward").save(),
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelRecHitsResY = SiPixelRecHitsResX.clone(
  name = "res_y",
  title = "Y resolution of RecHits",
  xlabel = "RecHit resolution Y dimension"
)

SiPixelRecHitsErrorX = DefaultHisto.clone(
  name = "rechiterror_x",
  title = "RecHit Error in X-direction",
  range_min = 0, range_max = 0.02, range_nbins = 100,
  xlabel = "X error",
  dimensions = 1,
  topFolderName = "PixelV/RecHits",
  specs = cms.VPSet(
    Specification().groupBy("").save(),
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk").save(),
  )
)

SiPixelRecHitsErrorY = SiPixelRecHitsErrorX.clone(
  name = "rechiterror_y",
  title = "RecHit Error in Y-direction",
  xlabel = "Y error"
)

SiPixelRecHitsPullX = DefaultHisto.clone(
  name = "pull_x",
  title = "RecHit Pull in X-direction",
  range_min = -10., range_max = 10., range_nbins = 100,
  xlabel = "X Pull",
  dimensions = 1,
  topFolderName = "PixelV/RecHits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk").save(),
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelRecHitsPullY = SiPixelRecHitsPullX.clone(
  name = "pull_y",
  title = "RecHit Pull in Y-direction",
  xlabel = "Y Pull"
)

SiPixelRecHitsConf = cms.VPSet(
  SiPixelRecHitsInTimeEvents,
  SiPixelRecHitsOutTimeEvents,
  SiPixelRecHitsNSimHits,
  SiPixelRecHitsPosX,
  SiPixelRecHitsPosY,
  SiPixelRecHitsResX,
  SiPixelRecHitsResY,
  SiPixelRecHitsErrorX,
  SiPixelRecHitsErrorY,
  SiPixelRecHitsPullX,
  SiPixelRecHitsPullY,
)

SiPixelRecHitsAnalyzerV = cms.EDAnalyzer("SiPixelRecHitsV",
        src = cms.InputTag("siPixelRecHits"),
        # Track assoc. parameters
        associatePixel = cms.bool(True),
        ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof', 
            'g4SimHitsTrackerHitsPixelBarrelHighTof', 
            'g4SimHitsTrackerHitsPixelEndcapLowTof', 
            'g4SimHitsTrackerHitsPixelEndcapHighTof'),
        associateStrip = cms.bool(False),
        associateRecoTracks = cms.bool(False),
        pixelSimLinkSrc = cms.InputTag("simSiPixelDigis"),
        stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
        histograms = SiPixelRecHitsConf,
        geometry = SiPixelPhase1Geometry
)

SiPixelRecHitsHarvesterV = cms.EDAnalyzer("SiPixelPhase1Harvester",
        histograms = SiPixelRecHitsConf,
        geometry = SiPixelPhase1Geometry
)
