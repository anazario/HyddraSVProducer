import FWCore.ParameterSet.Config as cms

miniAODTrackProducer = cms.EDProducer("MiniAODTrackProducer",
    pfCandidates = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    eleLostTracks = cms.InputTag("lostTracks", "eleTracks"),
)
