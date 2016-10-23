import FWCore.ParameterSet.Config as cms

#   module trackerHitsValid = TrackerHitProducer
trackerHitsValid = cms.EDAnalyzer("TrackerHitAnalyzer",
    G4TrkSrc = cms.InputTag("g4SimHits"),
    SiTIDLowSrc = cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"),
    Name = cms.untracked.string('TrackerHitAnalyzer'),
    Verbosity = cms.untracked.bool(False),
    runStandalone = cms.bool(False),
    outputFile =cms.untracked.string(''),

    Label = cms.string('TrkHits'),
    ProvenanceLookup = cms.PSet(
        PrintProvenanceInfo = cms.untracked.bool(False),
        GetAllProvenances = cms.untracked.bool(False)
    ),
    SiTOBLowSrc = cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"),
    SiTIBHighSrc = cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"),
    SiTIBLowSrc = cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"),
    SiTOBHighSrc = cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"),
    SiTECHighSrc = cms.InputTag("g4SimHits","TrackerHitsTECHighTof"),
    SiTECLowSrc = cms.InputTag("g4SimHits","TrackerHitsTECLowTof"),
    SiTIDHighSrc = cms.InputTag("g4SimHits","TrackerHitsTIDHighTof")
)


