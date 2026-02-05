import FWCore.ParameterSet.Config as cms

# Import the standard TrackAssociatorParameters and customize for AOD
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

# Customize for AOD (reduced rechits)
tkAssocParamBlock = TrackAssociatorParameterBlock.clone()
tkAssocParamBlock.TrackAssociatorParameters.useMuon = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useCalo = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useHO = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.usePreshower = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEE")
tkAssocParamBlock.TrackAssociatorParameters.EBRecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEB")
tkAssocParamBlock.TrackAssociatorParameters.HBHERecHitCollectionLabel = cms.InputTag("reducedHcalRecHits", "hbhereco")
tkAssocParamBlock.TrackAssociatorParameters.HORecHitCollectionLabel = cms.InputTag("reducedHcalRecHits", "horeco")

# =============================================================================
# Track Collection Configuration
# =============================================================================
# IMPORTANT: The producer (hyddraSVs) and analyzer (hyddraSVAnalyzer) MUST use
# the same track collection! The producer builds vertices from tracks, and the
# analyzer does gen-matching on the same tracks.

TRACK_COLLECTION_CONFIG = {
    'general': cms.InputTag("generalTracks"),
    'generalFiltered': cms.InputTag("filteredTrackProducer", "filteredTracks"),
    'selected': cms.InputTag("muonEnhancedTracks", "selectedTracks"),
    'muon': cms.InputTag("displacedGlobalMuons"),
    'muonGlobal': cms.InputTag("muonGlobalTrackProducer", "globalTracks"),
    'displacedMuonGlobal': cms.InputTag("muonGlobalTrackProducer", "displacedGlobalTracks"),
    'sip2D': cms.InputTag("muonEnhancedTracks", "sip2DTracks"),
    'sip2DMuonEnhanced': cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    'muonEnhanced': cms.InputTag("muonEnhancedTracks", "muonEnhancedTracks"),
}

def configureTrackCollection(process, trackCollection):
    """
    Configure both the SV producer and analyzer to use the same track collection.
    Call this after loading both modules.

    Usage:
        from KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi import configureTrackCollection
        configureTrackCollection(process, 'general')
    """
    if trackCollection not in TRACK_COLLECTION_CONFIG:
        raise ValueError(f"Unknown track collection: {trackCollection}. "
                        f"Available: {list(TRACK_COLLECTION_CONFIG.keys())}")

    inputTag = TRACK_COLLECTION_CONFIG[trackCollection]

    # Configure producer
    if hasattr(process, 'hyddraSVs'):
        process.hyddraSVs.tracks = inputTag

    # Configure analyzer
    if hasattr(process, 'hyddraSVAnalyzer'):
        process.hyddraSVAnalyzer.tracks = inputTag

    print(f"Configured track collection: {trackCollection} -> {inputTag}")

# =============================================================================
# Default analyzer (sip2DMuonEnhanced tracks)
# =============================================================================
hyddraSVAnalyzer = cms.EDAnalyzer("HyddraSVAnalyzer",
    hasGenInfo = cms.bool(True),
    leptonicVertices = cms.InputTag("hyddraSVs", "leptonicVertices"),
    hadronicVertices = cms.InputTag("hyddraSVs", "hadronicVertices"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    tracks = TRACK_COLLECTION_CONFIG['sip2DMuonEnhanced'],
    muonTracks = cms.InputTag("displacedGlobalMuons"),
    mergedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs"),
    genParticles = cms.InputTag("genParticles"),
    TrackAssociatorParameters = tkAssocParamBlock.TrackAssociatorParameters.clone(),
)

# =============================================================================
# Alternative configuration for data (no gen info)
# =============================================================================
hyddraSVAnalyzerData = hyddraSVAnalyzer.clone(
    hasGenInfo = cms.bool(False),
)

# =============================================================================
# Leptonic-only and Hadronic-only variants
# =============================================================================
hyddraSVAnalyzerLeptonic = hyddraSVAnalyzer.clone(
    hadronicVertices = cms.InputTag(""),  # Empty tag disables hadronic
)

hyddraSVAnalyzerHadronic = hyddraSVAnalyzer.clone(
    leptonicVertices = cms.InputTag(""),  # Empty tag disables leptonic
)
