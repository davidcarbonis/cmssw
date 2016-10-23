import FWCore.ParameterSet.Config as cms

from Validation.TrackerDigis.stripDigisValidation_cfi import *
trackerDigisValidation = cms.Sequence(stripDigisValid)


