# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: 2CosmicsRaw2Reco.py --data --conditions 80X_dataRun2_Prompt_v6 --era Run2_25ns --step RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO -n 100 --no_exec --filein root://xrootd.unl.edu//store/express/Run2016A/ExpressPhysics/FEVT/Express-v1/000/271/056/00000/0AE3F2F3-F908-E611-B3A2-02163E014126
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
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('L1Trigger.Configuration.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

PB = cms.bool(True)

if PB:
  process.load('DPGAnalysis.PilotBladeStudy.PilotBladeSetup_cfi')
#endif

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

# Input source 
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/data/vami/projects/pilotBlade/pp2016Processing/RAW/0AE3F2F3-F908-E611-B3A2-02163E014126.root'),
##    fileNames = cms.untracked.vstring('file:/data/vami/projects/pilotBlade/pp2016Processing/FEVT/273730/00821BA0-381F-E611-BBC5-02163E0142C0.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

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

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8', '')

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
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.conditionsInEdm_step = cms.Path(process.conditionsInEdm)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.RECOoutput)

if PB :
  process.PBDigi_step = cms.Path(process.PBDigis)
  process.pilotBladeReco_step = cms.Path(
    process.PBClusters*
    process.PBRecHits
  )
# endif

# Schedule definition
process.schedule = cms.Schedule(
  #standard part
  process.raw2digi_step,
  process.L1Reco_step,
  process.conditionsInEdm_step,
  process.reconstruction_step,
  
  #PilotBlade part
  process.PBDigi_step,
  process.pilotBladeReco_step,
  
  #end part
  process.endjob_step,
  process.output_step
)
