import FWCore.ParameterSet.Config as cms

fastSimMuonEnhancedTracks = cms.EDProducer("FastSimMuonEnhancedTracksProducer",
    generalTracks = cms.InputTag("generalTracks"),
    dsaMuonTracks = cms.InputTag("standAloneMuons"),
    dgMuonTracks = cms.InputTag("globalMuons"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
)
