import FWCore.ParameterSet.Config as cms

hyddraEXOAnalyzer = cms.EDAnalyzer("HyddraSVsEXOAnalyzer",
    # Base collection: "seeds", "inclusive" (default), or "isolated"
    outputCollection  = cms.string("inclusive"),
    # Vertex collections
    seedVertices      = cms.InputTag("hyddraEXO", "seedVertices"),
    inclusiveVertices = cms.InputTag("hyddraEXO", "inclusiveVertices"),
    isolatedVertices  = cms.InputTag("hyddraEXO", "isolatedVertices"),
    # Flag vectors
    disambiguationFlags = cms.InputTag("hyddraEXO", "disambiguationFlags"),
    seedIsolationFlags  = cms.InputTag("hyddraEXO", "seedIsolationFlags"),
    isolationFlags      = cms.InputTag("hyddraEXO", "isolationFlags"),
    # Other
    pvCollection      = cms.InputTag("offlinePrimaryVertices"),
    tracks            = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    genParticles      = cms.InputTag("genParticles"),
    MET               = cms.InputTag("slimmedMETs"),
    hasGenInfo        = cms.bool(True),
    motherPdgId       = cms.int32(54),    # dark photon proxy; change per signal
    genDRCut          = cms.double(0.05), # gold/bronze track-to-gen matching threshold
    passSelDRCut      = cms.double(0.02), # passSelection acceptance threshold
)
