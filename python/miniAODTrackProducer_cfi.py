import FWCore.ParameterSet.Config as cms

miniAODTrackProducer = cms.EDProducer("MiniAODTrackProducer",
    pfCandidates = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    eleLostTracks = cms.InputTag("lostTracks", "eleTracks"),
    pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    # Cut configuration (disabled by default)
    applyCuts = cms.bool(False),
    minPt = cms.double(1.0),              # pT > 1 GeV
    minAbsSip2D = cms.double(4.0),        # |sip2D| >= 4 (select displaced tracks)
    maxNormalizedChi2 = cms.double(5.0),  # chi2/ndof < 5
)
