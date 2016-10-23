import FWCore.ParameterSet.Config as cms
from DQM.SiPixelPhase1Common.HistogramManager_cfi import *

SiPixelTrackClustersCharge = DefaultHisto.clone(
  name = "charge",
  title = "Corrected Cluster Charge",
  range_min = 0, range_max = 100, range_nbins = 200,
  xlabel = "Charge size (in ke)",
  topFolderName = "PixelV/Clusters",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelTrackClustersSizeX = DefaultHisto.clone(
  name = "size_x",
  title = "Cluster Size X",
  range_min = 0, range_max = 30, range_nbins = 30,
  xlabel = "Cluster size (in pixels)",
  topFolderName = "PixelV/Clusters",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelTrackClustersSizeY = DefaultHisto.clone(
  name = "size_y",
  title = "Cluster Size Y",
  range_min = 0, range_max = 30, range_nbins = 30,
  xlabel = "Cluster size (in pixels)",
  topFolderName = "PixelV/Clusters",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save(),
  )
)

SiPixelTrackClustersConf = cms.VPSet(
  SiPixelTrackClustersCharge,
  SiPixelTrackClustersSizeX,
  SiPixelTrackClustersSizeY
)


SiPixelTrackClustersAnalyzerV = cms.EDAnalyzer("SiPixelTrackClustersV",
        clusters = cms.InputTag("siPixelClusters"),
        trajectories = cms.InputTag("generalTracks"),
        histograms = SiPixelTrackClustersConf,
        geometry = SiPixelPhase1Geometry
)

SiPixelTrackClustersHarvesterV = cms.EDAnalyzer("SiPixelPhase1Harvester",
        histograms = SiPixelTrackClustersConf,
        geometry = SiPixelPhase1Geometry
)
