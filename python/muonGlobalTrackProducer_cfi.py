import FWCore.ParameterSet.Config as cms

muonGlobalTrackProducer = cms.EDProducer("MuonGlobalTrackProducer",
    muons = cms.InputTag("muons"),
    displacedMuons = cms.InputTag("displacedMuons"),
    mode = cms.string("priority"),
    trackPriority = cms.vstring("globalTrack", "innerTrack", "outerTrack"),
)
