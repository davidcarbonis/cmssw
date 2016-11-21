import FWCore.ParameterSet.Config as cms

from Configuration.Geometry.GeometrySimDB_cff import*
from Configuration.Geometry.GeometryRecoDB_cff import*
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import*

import CalibTracker.Configuration.Common.PoolDBESSource_cfi
# ----------------------- CONDITIONS -------------------

# --------------------- SiPixelQuality -----------------
SiPixelQualityDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  #connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelQuality_PilotBlade4.db'), #local DB
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), 
 #Fake Pilot Blade out, Pilot Blade in --> PilotBlade3.db
 #Fake Pilot Blade out, Pilot Blade out --> PilotBlade4.db
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelQualityFromDbRcd'),
      #tag = cms.string('SiPixelQuality_PilotBlade') #tagname for the local DB
      tag = cms.string('SiPixelQuality_PilotBlade_v1'),  
    )
  )
)
es_prefer_Quality = cms.ESPrefer("PoolDBESSource","SiPixelQualityDBReader")
#-------------------------------------------------------

# --------------------- CablingMap ---------------------
CablingMapDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  #connect = cms.string('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelCabling_PilotBlade_data.db'), #local DB
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), 
 #SiPixelCabling_Ph0andPilotBlade.db --> PB + BPix + FPix
 #SiPixelCabling_PilotBlade.db --> just PB, no BPix, no FPix
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelFedCablingMapRcd'),
      label = cms.untracked.string('pilotBlade'), 
      #tag = cms.string('SiPixelFedCablingMap_data') #tagname for the local DB
      tag = cms.string('SiPixelCabling_PilotBlade_data_v1'),  
    )
  )
)
es_prefer_CablingReader = cms.ESPrefer("PoolDBESSource","CablingMapDBReader")
#-------------------------------------------------------

# --------------------- Gain Calib DB ------------------
GainDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  #connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/GainDB_PilotBlade.db'),
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), #from Prep
  #connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/GainDB_PilotBlade_v2.db'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelGainCalibrationOfflineRcd'),
      #tag = cms.string('GainCalib_TEST_offline') #tagname for the local DB
      tag = cms.string('GainCalib_PilotBlade_offline_v1')
    )
  )
)
es_prefer_Gain = cms.ESPrefer("PoolDBESSource","GainDBReader")
#-------------------------------------------------------

# --------------------- Gain Calib HLT DB ------------------
GainDBReaderHLT = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), #from Prep
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelGainCalibrationForHLTRcd'),
      tag = cms.string('GainCalib_PilotBlade_hlt_v1')
    )
  )
)
es_prefer_GainHLT = cms.ESPrefer("PoolDBESSource","GainDBReaderHLT")
#-------------------------------------------------------

# ----------------------- GenError ---------------------
GenErrReader = cms.ESSource("PoolDBESSource",
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  #connect = cms.string('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelGenErrors38T3_PilotBlade.db'), #local DB
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), #from Prep
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelGenErrorDBObjectRcd'),
      #tag = cms.string('SiPixelGenErrorDBObject38T3') #tagname for the local DB
      tag = cms.string('SiPixelGenErrors38T_PilotBlade_v1'), 
    )
  )
)
es_prefer_GenErr = cms.ESPrefer("PoolDBESSource","GenErrReader")
#-------------------------------------------------------

# --------------------- LorentzAngle -------------------

LorentzAngleDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  #connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelLorentzAngle_PilotBlade3.db'), #local DB
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), #from Prep
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelLorentzAngleRcd'),
      label = cms.untracked.string(''),
      #tag = cms.string('SiPixelLorentzAngle_v02_mc') #tagname for the local DB
      tag = cms.string('SiPixelLorentzAngle_PilotBlade_v1'), 
    )
  )
)
es_prefer_LA = cms.ESPrefer("PoolDBESSource","LorentzAngleDBReader")
#-------------------------------------------------------

# ---------------------- PB Geometry -------------------
trackerGeometryDB.applyAlignment = cms.bool(True)
XMLFromDBSource.label=''
PoolDBESSourceGeometry = cms.ESSource("PoolDBESSource",
  #CondDBSetup,
  timetype = cms.string('runnumber'),
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'),
  #connect = cms.string('sqlite_file:/data/vami/projects/pilotBlade/1ConditionDBs/PilotGeometry.db') #local file
  #PilotGeometry.db --> with PB and Fake PB
  #PilotGeometry0.db --> with PB only
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('GeometryFileRcd'),
      #tag = cms.string('XMLFILE_Geometry_74YV2_Extended2015_mc') #for the local file
      tag = cms.string('SiPixelPilotGeometry_v1')
    ),
    cms.PSet(
      record = cms.string('IdealGeometryRecord'),
      tag = cms.string('TKRECO_Geometry_forPilotBlade_v1')
    ),
    cms.PSet(
      record = cms.string('PGeometricDetExtraRcd'),
      tag = cms.string('TKExtra_Geometry_forPilotBlade_v1')
    ),
    cms.PSet(
      record = cms.string('PTrackerParametersRcd'),
      tag = cms.string('TKParameters_Geometry_forPilotBlade_v1')
    ),
    #cms.PSet(
    #  record = cms.string('PEcalBarrelRcd'),
    #  tag = cms.string('EBRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('PEcalEndcapRcd'),
      #tag = cms.string('EERECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('PEcalPreshowerRcd'),
      #tag = cms.string('EPRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('PHcalRcd'),
      #tag = cms.string('HCALRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('PCaloTowerRcd'),
      #tag = cms.string('CTRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('PZdcRcd'),
      #tag = cms.string('ZDCRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('PCastorRcd'),
      #tag = cms.string('CASTORRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('CSCRecoGeometryRcd'),
      #tag = cms.string('CSCRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('CSCRecoDigiParametersRcd'),
      #tag = cms.string('CSCRECODIGI_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('DTRecoGeometryRcd'),
      #tag = cms.string('DTRECO_Geometry_forPilotBlade_v1')
    #),
    #cms.PSet(
      #record = cms.string('RPCRecoGeometryRcd'),
      #tag = cms.string('RPCRECO_Geometry_forPilotBlade_v1')
    #)
  )
)
es_prefer_geometry = cms.ESPrefer( "PoolDBESSource", "PoolDBESSourceGeometry" )
#-------------------------------------------------------

'''
# forWidth LA
# we dont use it
LAWidthReader = cms.ESSource("PoolDBESSource",
   BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    DBParameters = cms.PSet(
  messageLevel = cms.untracked.int32(0),
  authenticationPath = cms.untracked.string('')
   ),
   #connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelLorentzAngle_PilotBlade4.db'),
   #connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelLorentzAngle_PilotBlade5.db'),
   connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), 
   toGet = cms.VPSet(
   cms.PSet(
    record = cms.string('SiPixelLorentzAngleRcd'),
    label = cms.untracked.string('forWidth'),
    tag = cms.string('SiPixelLorentzAngle_PilotBlade_forWidth_v1')
    #tag = cms.string('SiPixelLorentzAngle_forWidth_v1_mc[cms_orcon_prod/CMS_COND_31X_PIXEL]') # local DB
 ),
   ),
)
es_prefer_LAWidth = cms.ESPrefer("PoolDBESSource","LAWidthReader")

# fromAlignment LA
# we dont use it
LAAlignmentReader = cms.ESSource("PoolDBESSource",
   BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    DBParameters = cms.PSet(
  messageLevel = cms.untracked.int32(0),
  authenticationPath = cms.untracked.string('')
   ),
   #connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelLorentzAngle_PilotBlade4.db'),
   connect = cms.string ('sqlite_file:/data/vami/projects/PilotBlade/1ConditionDBs/SiPixelLorentzAngle_PilotBlade6.db'),
    toGet = cms.VPSet(
   cms.PSet(
    record = cms.string('SiPixelLorentzAngleRcd'),
    label = cms.untracked.string('fromAlignment'),
    #tag = cms.string('SiPixelLorentzAngle_forWidth_v1_mc')
    tag = cms.string('SiPixelLorentzAngle_fromAlignment_v1_mc[cms_orcon_prod/CMS_COND_31X_PIXEL]')
 ),
   ),
)
es_prefer_LAAlignment = cms.ESPrefer("PoolDBESSource","LAAlignmentReader")
#-------------------------------------------------------
'''
