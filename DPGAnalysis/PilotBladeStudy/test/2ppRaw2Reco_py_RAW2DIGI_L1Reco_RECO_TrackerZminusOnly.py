# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 -s RAW2DIGI,RECO:reconstruction_trackingOnly,DQM:@trackingOnlyDQM -n 100 --eventcontent RECO,DQM --datatier RECO,DQMIO --conditions auto:run2_data --filein file:JetHT_Run273447_raw.root --data --no_exec --python_filename=runFirst.py --scenario pp --era Run2_25ns
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO2',eras.Run2_25ns)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('L1Trigger.Configuration.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
##process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

PB = cms.bool(True)

if PB:
  process.load('DPGAnalysis.PilotBladeStudy.PilotBladeSetup_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

'''
process.Timing = cms.Service("Timing"
    ,summaryOnly = cms.untracked.bool(True)
)
'''

process.options = cms.untracked.PSet(
)


'''
#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
'''

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/data/vami/projects/pilotBlade/pp2016Processing/RAW/0AE3F2F3-F908-E611-B3A2-02163E014126.root'),
    ##fileNames = cms.untracked.vstring('file:/data/vami/projects/pilotBlade/pp2016Processing/FEVT/273730/00821BA0-381F-E611-BBC5-02163E0142C0.root'),
    secondaryFileNames = cms.untracked.vstring()
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('2ppRaw2Reco.py nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('2ppRaw2Reco_py_RAW2DIGI_L1Reco_RECO.root'),
    #outputCommands = process.RECOEventContent.outputCommands,
     outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_*lumi*_*_*',
    'keep *_*eamSpot*_*_*',
    'keep *_*onditions*_*_*',
    'keep *_*ixel*_*_*',
    'keep *_*PB*_*_*',
    'keep *_*racks*_*_*',
    'keep *_*iStripClusters*_*_*',
    'keep *awDataError*_*_*_*',
    'keep *_*rig*_*_*',
    'keep *rig*_*_*_*'
    ),
    splitLevel = cms.untracked.int32(0)
)


# Additional output definition

# kill event without hits in Pilot Blade
if PB :
  process.bysipixelclustmulteventfilter = cms.EDFilter('BySiPixelClusterMultiplicityEventFilter',
                                               multiplicityConfig = cms.PSet(
                                                                  collectionName = cms.InputTag("PBClusters"),
                                                                  moduleThreshold = cms.untracked.int32(-1),
                                                                  useQuality = cms.untracked.bool(False),
                                                                  qualityLabel = cms.untracked.string("")
                                                                  ),
                                               cut = cms.string("mult > 0")
                                               )

  process.PBClusterFilter_step = cms.Path(process.bysipixelclustmulteventfilter)
  process.RECOoutput.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("PBClusterFilter_step"))

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 
#  '80X_dataRun2_Prompt_v8_CustomTracker2016B_v0'
 '80X_dataRun2_Prompt_v8'
  , '')


##----------------------%<------------------
process.GlobalTag.toGet = cms.VPSet(
         cms.PSet(record = cms.string('TrackerAlignmentRcd'),
                  tag =  cms.string('Alignments'),
                  connect = cms.string('sqlite_file:./tracker_alignment_80X_dataRun2_Prompt_v8.db')
                  ),
         cms.PSet(record = cms.string('TrackerAlignmentErrorExtendedRcd'),
                  tag =  cms.string('AlignmentErrorsExtended'),
                  connect = cms.string('sqlite_file:./tracker_alignment_80X_dataRun2_Prompt_v8.db')
                  ),
         cms.PSet(record = cms.string('TrackerSurfaceDeformationRcd'),
                  tag =  cms.string('AlignmentSurfaceDeformations'),
                  connect = cms.string('sqlite_file:./tracker_alignment_80X_dataRun2_Prompt_v8.db')
                  )
         )
##--------------------->%---------------------

# --------------------- Reconstruction --------------------
#Standard Pixel Digis
process.siPixelDigis.UseQualityInfo = cms.bool(True)

if PB :
  #Pilot Blade Digis
  process.PBDigis = cms.EDProducer("SiPixelRawToDigi",
    InputLabel = cms.InputTag("rawDataCollector"),
    CablingMapLabel =  cms.string("pilotBlade"),
    UsePhase1 = cms.bool(False),
    UsePilotBlade = cms.bool(True),
    UseQualityInfo = cms.bool(False),
    IncludeErrors = cms.bool(True),
    UserErrorList = cms.vint32(25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)
  )
  #Pilot Blade Clusters
  from CondTools.SiPixel.SiPixelGainCalibrationService_cfi import *
  process.PBClusters = cms.EDProducer("SiPixelClusterProducer",
    SiPixelGainCalibrationServiceParameters,
    src = cms.InputTag("PBDigis"),
    ChannelThreshold = cms.int32(1000),
    MissCalibrate = cms.untracked.bool(False),
    SplitClusters = cms.bool(False),
    VCaltoElectronGain = cms.int32(65),
    VCaltoElectronOffset = cms.int32(-414),                          
    payloadType = cms.string('Offline'),
    SeedThreshold = cms.int32(1000),
    ClusterThreshold = cms.double(4000.0),
    maxNumberOfClusters = cms.int32(-1),
  )

  #Pilot Blade RecHits
  process.PBRecHits = cms.EDProducer("SiPixelRecHitConverter",
    src = cms.InputTag("PBClusters"),
    CPE = cms.string('PixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
  )
# endif
#-------------------------------------------------------

# Path and EndPath definitions
process.conditionsInEdm_step = cms.Path(process.conditionsInEdm)

# special reco step
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import siPixelClusters as _siPixelClusters
process.siPixelClusters = _siPixelClusters

from RecoLocalTracker.SiStripClusterizer.SiStripClustersFromRaw_cfi import SiStripClustersFromRawFacility as _siStripClustersFromRawFacility
process.SiStripClustersFromRawFacility = _siStripClustersFromRawFacility
process.SiStripClustersFromRawFacility.onDemand=True

process.MeasurementTrackerEvent.stripClusterProducer = "SiStripClustersFromRawFacility"


process.initialStepSeedLayers.layerList = [
    'BPix1+BPix2+FPix1_neg', 
    'BPix1+FPix1_neg+FPix2_neg'
]

process.InitialStepPreSplitting = cms.Sequence(
    process.siPixelDigis+
    process.siPixelClusters +
    process.siPixelRecHits +
    process.SiStripClustersFromRawFacility +
    process.MeasurementTrackerEvent +
    process.siPixelClusterShapeCache
)

from RecoTracker.FinalTrackSelectors.TrackCollectionMerger_cfi import *

import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
process.generalTracks =  TrackCollectionMerger.clone()
process.generalTracks.trackProducers = ['initialStepTracks',
                                     ]
process.generalTracks.inputClassifiers =["initialStep",
                                      ]

process.iterTracking = cms.Sequence(process.InitialStepPreSplitting*
                                    process.InitialStep*
                                    process.generalTracks
)

process.trk = cms.Path(
process.offlineBeamSpot*process.iterTracking
)

# endjob
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.RECOoutput)

if PB :
  process.PBDigi_step = cms.Path(process.PBDigis)
  process.pilotBladeReco_step = cms.Path(
    process.PBClusters*
    process.PBRecHits
  )

# Schedule definition
process.schedule = cms.Schedule(
    process.trk,
    process.conditionsInEdm_step,
    
    #PilotBlade part
    process.PBDigi_step,
    process.pilotBladeReco_step,
    process.PBClusterFilter_step,
  
    #end part 
    process.endjob_step,
    process.output_step
)


