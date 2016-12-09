import FWCore.ParameterSet.Config as cms
from DQM.SiPixelPhase1Common.HistogramManager_cfi import *

SiPixelPhase1TrackClustersOnTrackCharge = DefaultHistoTrack.clone(
  name = "charge",
  title = "Corrected Cluster Charge",
  range_min = 0, range_max = 200, range_nbins = 200,
  xlabel = "Charge size (in ke)",
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").saveAll(),
    Specification().groupBy("PXForward/PXDisk").saveAll(),
    StandardSpecification2DProfile
  )
)

SiPixelPhase1TrackClustersOnTrackSize = DefaultHistoTrack.clone(
  name = "size",
  title = "Total Cluster Size",
  range_min = 0, range_max = 30, range_nbins = 30,
  xlabel = "Cluster size (in pixels)",
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").saveAll(),
    Specification().groupBy("PXForward/PXDisk").saveAll(),
  )
)

SiPixelPhase1TrackClustersOnTrackNClusters = DefaultHistoTrack.clone(
  name = "clusters",
  title = "Clusters",
  range_min = 0, range_max = 10, range_nbins = 10,
  xlabel = "Number of Clusters",
  dimensions = 0,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer" + "/DetId/Event") 
                   .reduce("COUNT") 
                   .groupBy("PXBarrel/PXLayer")
                   .saveAll(),
    Specification().groupBy("PXForward/PXDisk" + "/DetId/Event") 
                   .reduce("COUNT") 
                   .groupBy("PXForward/PXDisk")
                   .saveAll(),
  )
)

SiPixelPhase1TrackClustersOnTrackPositionB = DefaultHistoTrack.clone(
  name = "clusterposition_zphi",
  title = "Cluster Positions",
  range_min   =  -60, range_max   =  60, range_nbins   = 600,
  range_y_min = -3.2, range_y_max = 3.2, range_y_nbins = 200,
  xlabel = "Global Z (cm)", ylabel = "Global \phi",
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("").save(),
  )
)

SiPixelPhase1TrackClustersOnTrackPositionF = DefaultHistoTrack.clone(
  name = "clusterposition_xy",
  title = "Cluster Positions",
  xlabel = "Global X (cm)", ylabel = "Global Y (cm)",
  range_min   = -20, range_max   = 20, range_nbins   = 200,
  range_y_min = -20, range_y_max = 20, range_y_nbins = 200,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXForward/PXDisk").save(),
  )
)

SiPixelPhase1TrackClustersOffTrackCharge = \
  SiPixelPhase1TrackClustersOnTrackCharge.clone(topFolderName = "PixelPhase1/OffTrack", 
  enabled = False,
  title = "Cluster Charge")
SiPixelPhase1TrackClustersOffTrackSize = \
  SiPixelPhase1TrackClustersOnTrackSize.clone(topFolderName = "PixelPhase1/OffTrack",
  enabled = False)

SiPixelPhase1TrackClustersOffTrackNClusters = \
  SiPixelPhase1TrackClustersOnTrackNClusters.clone(topFolderName = "PixelPhase1/OffTrack",
  enabled = False)

SiPixelPhase1TrackClustersOffTrackPositionB = \
  SiPixelPhase1TrackClustersOnTrackPositionB.clone(topFolderName = "PixelPhase1/OffTrack",
  enabled = False)

SiPixelPhase1TrackClustersOffTrackPositionF = \
  SiPixelPhase1TrackClustersOnTrackPositionF.clone(topFolderName = "PixelPhase1/OffTrack",
  enabled = False)

SiPixelPhase1TrackClustersNTracks = DefaultHistoTrack.clone(
  name = "ntracks",
  title = "Number of Tracks",
  xlabel = "All - Pixel - BPIX - FPIX",
  range_min = 1, range_max = 5, range_nbins = 4,
  dimensions = 1,
  topFolderName = "PixelPhase1/Tracks",
  specs = VPSet(
    Specification().groupBy("").save()
  )
)

SiPixelPhase1TrackClustersNTracksInVolume = DefaultHistoTrack.clone(
  name = "ntracksinpixvolume",
  title = "Number of Tracks in Pixel fiducial Volume",
  xlabel = "with hits - without hits",
  range_min = 0, range_max = 2, range_nbins = 2,
  dimensions = 1,
  topFolderName = "PixelPhase1/Tracks",
  specs = VPSet(
    Specification().groupBy("").save()
  )

)

SiPixelPhase1TrackClustersProb = DefaultHistoTrack.clone(
  name = "clusterprob",
  title = "Cluster Probability",
  xlabel = "log_10(Pr)",
  range_min = -10, range_max = 0, range_nbins = 200,
  dimensions = 1,
  specs = VPSet(
    StandardSpecifications1D
  )
)


SiPixelPhase1TrackClustersConf = cms.VPSet(
  SiPixelPhase1TrackClustersOnTrackCharge,
  SiPixelPhase1TrackClustersOnTrackSize,
  SiPixelPhase1TrackClustersOnTrackNClusters,
  SiPixelPhase1TrackClustersOnTrackPositionB,
  SiPixelPhase1TrackClustersOnTrackPositionF,

  SiPixelPhase1TrackClustersOffTrackCharge,
  SiPixelPhase1TrackClustersOffTrackSize,
  SiPixelPhase1TrackClustersOffTrackNClusters,
  SiPixelPhase1TrackClustersOffTrackPositionB,
  SiPixelPhase1TrackClustersOffTrackPositionF,

  SiPixelPhase1TrackClustersNTracks,
  SiPixelPhase1TrackClustersNTracksInVolume,
  SiPixelPhase1TrackClustersProb,
)


SiPixelPhase1TrackClustersAnalyzer = cms.EDAnalyzer("SiPixelPhase1TrackClusters",
        clusters = cms.InputTag("siPixelClusters"),
        trajectories = cms.InputTag("generalTracks"),
        histograms = SiPixelPhase1TrackClustersConf,
        geometry = SiPixelPhase1Geometry
)

SiPixelPhase1TrackClustersHarvester = cms.EDAnalyzer("SiPixelPhase1Harvester",
        histograms = SiPixelPhase1TrackClustersConf,
        geometry = SiPixelPhase1Geometry
)

