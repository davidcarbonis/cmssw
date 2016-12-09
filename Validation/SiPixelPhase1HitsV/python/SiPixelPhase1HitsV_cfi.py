import FWCore.ParameterSet.Config as cms
from DQM.SiPixelPhase1Common.HistogramManager_cfi import *

SiPixelPhase1HitsEnergyLoss = DefaultHisto.clone(
  name = "eloss",
  title = "Energy loss",
  range_min = 0, range_max = 0.001, range_nbins = 10000,
  xlabel = "Energy Loss",
  dimensions = 1,
  topFolderName = "PixelPhase1V/Hits",
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer/P1PXModuleName").save(),
    Specification().groupBy("PXForward/PXDisk/P1PXModuleName").save(),
  )
)

SiPixelPhase1HitsEntryExitX = DefaultHisto.clone(
  name = "entry_exit_x",
  title = "Entryx-Exitx",
  range_min = -0.03, range_max = 0.03, range_nbins = 10000,
  xlabel = "",
  dimensions = 1,
  topFolderName = "PixelPhase1V/Hits",
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer/P1PXModuleName").save(),
    Specification().groupBy("PXForward/PXDisk/P1PXModuleName").save(),
  )
)

SiPixelPhase1HitsEntryExitY = SiPixelPhase1HitsEntryExitX.clone(
  name = "entry_exit_y",
  title = "Entryy-Exity",
  xlabel = "",
  range_min = -0.03, range_max = 0.03, range_nbins = 10000,
)

SiPixelPhase1HitsEntryExitZ = SiPixelPhase1HitsEntryExitX.clone(
  name = "entry_exit_z",
  title = "Entryz-Exitz",
  xlabel = "",
  range_min = 0.0, range_max = 0.05, range_nbins = 10000,
)

SiPixelPhase1HitsPosX = DefaultHisto.clone(
  name = "local_x",
  title = "X position of Hits",
  range_min = -3.5, range_max = 3.5, range_nbins = 10000,
  xlabel = "Hit position X dimension",
  dimensions = 1,
  topFolderName = "PixelPhase1V/Hits",
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer/P1PXModuleName").save(),
    Specification().groupBy("PXForward/PXDisk/P1PXModuleName").save(),
  )
)

SiPixelPhase1HitsPosY = SiPixelPhase1HitsPosX.clone(
  name = "local_y",
  title = "Y position of Hits",
  xlabel = "Hit position Y dimension",
  range_min = -3.5, range_max = 3.5, range_nbins = 10000,
)

SiPixelPhase1HitsPosZ = SiPixelPhase1HitsPosX.clone(
  name = "local_z",
  title = "Z position of Hits",
  xlabel = "Hit position Z dimension",
  range_min = -0.05, range_max = 0.05, range_nbins = 500,
)

SiPixelPhase1HitsPosPhi = SiPixelPhase1HitsPosX.clone(
  name = "local_phi",
  title = "Phi position of Hits",
  xlabel = "Hit position phi dimension",
  range_min = -3.5, range_max = 3.5, range_nbins = 10000,
)

SiPixelPhase1HitsPosEta = SiPixelPhase1HitsPosX.clone(
  name = "local_eta",
  title = "Eta position of Hits",
  xlabel = "Hit position Eta dimension",
  range_min = -0.1, range_max = 0.1, range_nbins = 1000,
)

SiPixelPhase1HitsConf = cms.VPSet(
  SiPixelPhase1HitsEnergyLoss,
  SiPixelPhase1HitsEntryExitX,
  SiPixelPhase1HitsEntryExitY,
  SiPixelPhase1HitsEntryExitZ,
  SiPixelPhase1HitsPosX,
  SiPixelPhase1HitsPosY,
  SiPixelPhase1HitsPosZ,
  SiPixelPhase1HitsPosPhi,
  SiPixelPhase1HitsPosEta,
)

SiPixelPhase1HitsAnalyzerV = cms.EDAnalyzer("SiPixelPhase1HitsV",
        pixBarrelLowSrc = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
        pixBarrelHighSrc = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),
        pixForwardLowSrc = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
        pixForwardHighSrc = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"),
        # Track assoc. parameters
        histograms = SiPixelPhase1HitsConf,
        geometry = SiPixelPhase1Geometry
)

SiPixelPhase1HitsHarvesterV = cms.EDAnalyzer("SiPixelPhase1Harvester",
        histograms = SiPixelPhase1HitsConf,
        geometry = SiPixelPhase1Geometry
)
