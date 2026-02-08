import FWCore.ParameterSet.Config as cms

# =============================================================================
# MiniAOD Track Collection Configuration
# =============================================================================
# These collections are produced by MiniAODTrackProducer from packed candidates.
# Use with MiniAOD/NanoAOD input files.

MINIAOD_TRACK_COLLECTION_CONFIG = {
    # From packed candidates
    'pf': cms.InputTag("miniAODTrackProducer", "pfCandidateTracks"),
    'lost': cms.InputTag("miniAODTrackProducer", "lostTracks"),
    'eleLost': cms.InputTag("miniAODTrackProducer", "eleLostTracks"),
    'merged': cms.InputTag("miniAODTrackProducer", "merged"),
    'mergedWithEle': cms.InputTag("miniAODTrackProducer", "mergedWithEle"),
    # From muons (global tracks)
    'promptMuonExtracted': cms.InputTag("miniAODTrackProducer", "muonGlobalTracks"),
    'displacedMuonExtracted': cms.InputTag("miniAODTrackProducer", "displacedMuonGlobalTracks"),
    # From MiniAODMuonEnhancedTracksProducer (muon-enhanced with sip2D selection)
    'sip2DMuonEnhanced': cms.InputTag("miniAODMuonEnhancedTracks", "sip2DMuonEnhancedTracks"),
}

def configureMiniAODTrackCollection(process, trackCollection):
    """
    Configure both the SV producer and analyzer to use the same MiniAOD track collection.
    Call this after loading both modules.

    Usage:
        from KUCMSNtupleizer.HyddraSVProducer.miniAODTrackCollections_cfi import configureMiniAODTrackCollection
        configureMiniAODTrackCollection(process, 'merged')
    """
    if trackCollection not in MINIAOD_TRACK_COLLECTION_CONFIG:
        raise ValueError(f"Unknown MiniAOD track collection: {trackCollection}. "
                        f"Available: {list(MINIAOD_TRACK_COLLECTION_CONFIG.keys())}")

    inputTag = MINIAOD_TRACK_COLLECTION_CONFIG[trackCollection]

    # Configure producer
    if hasattr(process, 'hyddraSVs'):
        process.hyddraSVs.tracks = inputTag

    # Configure analyzer
    if hasattr(process, 'hyddraSVAnalyzer'):
        process.hyddraSVAnalyzer.tracks = inputTag

    print(f"Configured MiniAOD track collection: {trackCollection} -> {inputTag}")
