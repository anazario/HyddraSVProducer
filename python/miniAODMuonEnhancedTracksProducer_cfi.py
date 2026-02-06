import FWCore.ParameterSet.Config as cms

miniAODMuonEnhancedTracks = cms.EDProducer("MiniAODMuonEnhancedTracksProducer",
    generalTracks = cms.InputTag("miniAODTrackProducer", "merged"),
    dsaMuonTracks = cms.InputTag("displacedStandAloneMuons"),
    dgMuonTracks = cms.InputTag("displacedGlobalMuons"),
    displacedTracks = cms.InputTag("displacedTracks"),
    pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
