import FWCore.ParameterSet.Config as cms

muonGlobalTrackProducer = cms.EDProducer("MuonGlobalTrackProducer",
    muons = cms.InputTag("muons"),
)
