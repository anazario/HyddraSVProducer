import FWCore.ParameterSet.Config as cms

hyddraEXOAnalyzer = cms.EDAnalyzer("HyddraSVsEXOAnalyzer",
    inclusiveVertices = cms.InputTag("hyddraEXO", "inclusiveVertices"),
    isolatedVertices  = cms.InputTag("hyddraEXO", "isolatedVertices"),
    pvCollection      = cms.InputTag("offlinePrimaryVertices"),
    tracks            = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    genParticles      = cms.InputTag("genParticles"),
    hasGenInfo        = cms.bool(True),
    motherPdgId       = cms.int32(54),    # dark photon proxy; change per signal
    genDRCut          = cms.double(0.05), # gold/bronze track-to-gen matching threshold
    passSelDRCut      = cms.double(0.02), # passSelection acceptance threshold
)
