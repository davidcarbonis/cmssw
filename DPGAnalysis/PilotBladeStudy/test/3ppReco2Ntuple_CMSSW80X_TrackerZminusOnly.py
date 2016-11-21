import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('Analysis',eras.Run2_25ns)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#for the Pilot Blade reconstruction include the line below
process.load('DPGAnalysis.PilotBladeStudy.PilotBladeSetup_cfi')

# ------------------- Input ----------------------------
nEvents =  'All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:./2ppRaw2Reco_py_RAW2DIGI_L1Reco_RECO.root'),
    #eventsToProcess = cms.untracked.VEventRange('269794:286802','269807:198295'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
)
#-------------------------------------------------------

#----------------Production Info------------------------
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Ntuplizer.py nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)
#-------------------------------------------------------

## ------------------- Output ---------------------------
# we do not save the output from the standard sequance 
# only the output of the NTuplizer (PilotBladeStudy)
#-------------------------------------------------------

# --------------------- GlobalTag ----------------------
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

#---------------------- Refitter -----------------------
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.Refitter = process.TrackRefitterP5.clone()
process.Refitter.src = 'generalTracks'
process.Refitter.TrajectoryInEvent = True
process.MeasurementTrackerEvent.stripClusterProducer = "SiStripClustersFromRawFacility"

#-------------------------------------------------------

#------------------- PilotBladeStudy -------------------
process.PilotBladeStudy = cms.EDAnalyzer("PilotBladeStudy",
  trajectoryInput = cms.string('Refitter'),
  fileName = cms.string('nTuplePilotBlade_'+nEvents+'.root'),
  usePixelCPE= cms.bool(True),
  cosmicsCase = cms.bool(False),
  triggerNames=cms.vstring(
        "HLT_ZeroBias",
        "HLT_Physics",
        "HLT_Random",
        "HLT_MinBias"
        ),
  PositionCorrections = cms.untracked.VPSet(
       cms.PSet(
            id = cms.uint32(344132868),
            dx = cms.double(-0.1476),
            dy = cms.double(0.2455)
        ),
       cms.PSet(
            id = cms.uint32(344134148),
            dx = cms.double(-0.01271),
            dy = cms.double(0.3227)
        ),
       cms.PSet(
            id = cms.uint32(344131076),
            dx = cms.double(-0.2564),
            dy = cms.double(0.0032)
        ),
       cms.PSet(
            id = cms.uint32(344132100),
            dx = cms.double(-0.2183),
            dy = cms.double(-0.1538)
        ),
       cms.PSet(
            id = cms.uint32(344130820),
            dx = cms.double(0.1980),
            dy = cms.double(-0.086)
        ),
  ),
  FiducialRegions = cms.untracked.VPSet(
       cms.PSet(
            id = cms.uint32(344132868),
            marginX = cms.double(0.1),
            marginY = cms.double(0.1),
            rocX = cms.vint32(-1,-1,1,1),
            rocY = cms.vint32(1,2,1,2)
        )
  )                                       
)
#-------------------------------------------------------

#--------------------- EndPath -------------------------
process.endjob_step = cms.EndPath(process.endOfProcess)

process.Refitter_step = cms.Path(
  process.MeasurementTrackerEvent*
  process.Refitter
)

process.PilotBladeStudy_step = cms.Path(process.PilotBladeStudy)

#--------------------- Schedule ------------------------
process.schedule = cms.Schedule(
  process.Refitter_step,
  process.PilotBladeStudy_step,	
)
#-------------------------------------------------------
