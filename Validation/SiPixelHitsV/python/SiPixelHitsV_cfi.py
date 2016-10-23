import FWCore.ParameterSet.Config as cms
from DQM.SiPixelPhase1Common.HistogramManager_cfi import *

SiPixelHitsEnergyLoss = DefaultHisto.clone(
  name = "eloss",
  title = "Energy loss",
  range_min = 0, range_max = 0.001, range_nbins = 10000,
  xlabel = "Energy Loss",
  dimensions = 1,
  topFolderName = "PixelV/Hits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelHitsEntryExitX = DefaultHisto.clone(
  name = "entry_exit_x",
  title = "Entryx-Exitx",
  range_min = -0.03, range_max = 0.03, range_nbins = 10000,
  xlabel = "",
  dimensions = 1,
  topFolderName = "PixelV/Hits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelHitsEntryExitY = SiPixelHitsEntryExitX.clone(
  name = "entry_exit_y",
  title = "Entryy-Exity",
  xlabel = "",
  range_min = -0.03, range_max = 0.03, range_nbins = 10000,
)

SiPixelHitsEntryExitZ = SiPixelHitsEntryExitX.clone(
  name = "entry_exit_z",
  title = "Entryz-Exitz",
  xlabel = "",
  range_min = 0.0, range_max = 0.05, range_nbins = 10000,
)

SiPixelHitsPosX = DefaultHisto.clone(
  name = "local_x",
  title = "X position of Hits",
  range_min = -3.5, range_max = 3.5, range_nbins = 10000,
  xlabel = "Hit position X dimension",
  dimensions = 1,
  topFolderName = "PixelV/Hits",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelHitsPosY = SiPixelHitsPosX.clone(
  name = "local_y",
  title = "Y position of Hits",
  xlabel = "Hit position Y dimension",
  range_min = -3.5, range_max = 3.5, range_nbins = 10000,
)

SiPixelHitsConf = cms.VPSet(
  SiPixelHitsEnergyLoss,
  SiPixelHitsEntryExitX,
  SiPixelHitsEntryExitY,
  SiPixelHitsEntryExitZ,
  SiPixelHitsPosX,
  SiPixelHitsPosY,
)

SiPixelHitsAnalyzerV = cms.EDAnalyzer("SiPixelHitsV",
        pixBarrelLowSrc = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
        pixBarrelHighSrc = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),
        pixForwardLowSrc = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
        pixForwardHighSrc = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"),
        # Track assoc. parameters
        histograms = SiPixelHitsConf,
        geometry = SiPixelPhase1Geometry
)

SiPixelHitsHarvesterV = cms.EDAnalyzer("SiPixelPhase1Harvester",
        histograms = SiPixelHitsConf,
        geometry = SiPixelPhase1Geometry
)
