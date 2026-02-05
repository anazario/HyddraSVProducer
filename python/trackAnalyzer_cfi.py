import FWCore.ParameterSet.Config as cms

# =============================================================================
# Track Collection Configuration (shared with hyddraSVAnalyzer_cfi)
# =============================================================================
TRACK_COLLECTION_CONFIG = {
    'general': cms.InputTag("generalTracks"),
    'generalFiltered': cms.InputTag("filteredTrackProducer", "filteredTracks"),
    'selected': cms.InputTag("muonEnhancedTracks", "selectedTracks"),
    'muon': cms.InputTag("displacedGlobalMuons"),
    'muonGlobal': cms.InputTag("muonGlobalTrackProducer", "globalTracks"),
    'sip2D': cms.InputTag("muonEnhancedTracks", "sip2DTracks"),
    'sip2DMuonEnhanced': cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    'muonEnhanced': cms.InputTag("muonEnhancedTracks", "muonEnhancedTracks"),
}

# =============================================================================
# Default analyzer
# =============================================================================
trackAnalyzer = cms.EDAnalyzer("TrackAnalyzer",
    hasGenInfo = cms.bool(True),
    genMatchDeltaRCut = cms.double(0.02),
    tracks = TRACK_COLLECTION_CONFIG['sip2DMuonEnhanced'],
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
)

# =============================================================================
# Alternative configuration for data (no gen info)
# =============================================================================
trackAnalyzerData = trackAnalyzer.clone(
    hasGenInfo = cms.bool(False),
)
