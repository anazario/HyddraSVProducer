import FWCore.ParameterSet.Config as cms

trackAnalyzer = cms.EDAnalyzer("TrackAnalyzer",
    hasGenInfo = cms.bool(True),
    genMatchDeltaRCut = cms.double(0.02),
    tracks = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
)

# Alternative configuration for data (no gen info)
trackAnalyzerData = trackAnalyzer.clone(
    hasGenInfo = cms.bool(False),
)

# Pre-configured analyzers for different track collections
# These make it easy to switch between collections for efficiency studies

# General tracks (no selection)
trackAnalyzerGeneral = trackAnalyzer.clone(
    tracks = cms.InputTag("generalTracks"),
)

# Selected tracks (pt > 1 GeV from KUCMSNtupleizer)
trackAnalyzerSelected = trackAnalyzer.clone(
    tracks = cms.InputTag("muonEnhancedTracks", "selectedTracks"),
)

# Muon tracks (displaced global muons)
trackAnalyzerMuon = trackAnalyzer.clone(
    tracks = cms.InputTag("displacedGlobalMuons"),
)

# Sip2D tracks (sip2D > 3)
trackAnalyzerSip2D = trackAnalyzer.clone(
    tracks = cms.InputTag("muonEnhancedTracks", "sip2DTracks"),
)

# Sip2D muon-enhanced tracks (default for HYDDRA)
trackAnalyzerSip2DMuonEnhanced = trackAnalyzer.clone(
    tracks = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
)

# Muon-enhanced tracks (before sip2D cut)
trackAnalyzerMuonEnhanced = trackAnalyzer.clone(
    tracks = cms.InputTag("muonEnhancedTracks", "muonEnhancedTracks"),
)
