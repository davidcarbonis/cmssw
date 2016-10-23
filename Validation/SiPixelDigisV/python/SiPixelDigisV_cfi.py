import FWCore.ParameterSet.Config as cms

# this might also go into te Common config,as we do not reference it
from DQM.SiPixelPhase1Common.HistogramManager_cfi import *

SiPixelDigisADC = DefaultHisto.clone(
  name = "adc",
  title = "Digi ADC values",
  xlabel = "ADC counts",
  range_min = 0,
  range_max = 300,
  range_nbins = 300,
  topFolderName = "PixelV/Digis",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save()
  )
)

SiPixelDigisNdigis = DefaultHisto.clone(
  name = "digis", # 'Count of' added automatically
  title = "Digis",
  xlabel = "Number of Digis",
  range_min = 0,
  range_max = 30,
  range_nbins = 30,
  dimensions = 0, # this is a count
  topFolderName = "PixelV/Digis",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save()
  )
)

SiPixelDigisRows = DefaultHisto.clone(
  name = "row",
  title = "Digi Rows",
  xlabel = "Row",
  range_min = 0,
  range_max = 200,
  range_nbins = 200,
  topFolderName = "PixelV/Digis",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save()
  )
)

SiPixelDigisColumns = DefaultHisto.clone(
  name = "column",
  title = "Digi Columns",
  xlabel = "Column",
  range_min = 0,
  range_max = 300,
  range_nbins = 300,
  topFolderName = "PixelV/Digis",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/PXLayer|PXDisk/PXBModule|PXFModule").save()
  )
)

SiPixelDigisDebug = DefaultHisto.clone(
  enabled = False,
  name = "debug",
  xlabel = "Ladder #",
  range_min = 1,
  range_max = 64,
  range_nbins = 64,
  topFolderName = "PixelV/Debug",
  specs = cms.VPSet(
    Specification().groupBy("PXBarrel|PXForward/Shell|HalfCylinder/PXLayer|PXDisk/PXRing|/PXLadder|PXBlade") 
                   .save()
                   .reduce("MEAN")
                   .groupBy(parent("PXBarrel|PXForward/Shell|HalfCylinder/PXLayer|PXDisk/PXRing|/PXLadder|PXBlade"), "EXTEND_X")
                   .saveAll(),
  )
)

# This has to match the order of the names in the C++ enum.
SiPixelDigisConf = cms.VPSet(
  SiPixelDigisADC,
  SiPixelDigisNdigis,
  SiPixelDigisRows,
  SiPixelDigisColumns,
  SiPixelDigisDebug
)

SiPixelDigisAnalyzerV = cms.EDAnalyzer("SiPixelDigisV",
        src = cms.InputTag("simSiPixelDigis"), 
        histograms = SiPixelDigisConf,
        geometry = SiPixelPhase1Geometry
)

SiPixelDigisHarvesterV = cms.EDAnalyzer("SiPixelDigisHarvesterV",
        histograms = SiPixelDigisConf,
        geometry = SiPixelPhase1Geometry
)
