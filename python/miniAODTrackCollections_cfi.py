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
    'mergedAll': cms.InputTag("miniAODTrackProducer", "mergedAll"),
    # From muons (global tracks)
    'promptMuonExtracted': cms.InputTag("miniAODTrackProducer", "muonGlobalTracks"),
    'displacedMuonExtracted': cms.InputTag("miniAODTrackProducer", "displacedMuonGlobalTracks"),
    # PF+Lost (muon-ID'd tracks removed, duplicates of slimmedMuons removed) + all slimmedMuon
    # bestTracks, with sip2D cuts applied uniformly to the full merged collection
    'mergedMuonSip2D': cms.InputTag("miniAODTrackProducer", "mergedMuonSip2D"),
    # From muons using LLPNanoAOD PatMuonVertex selection
    # (isGlobalMuon -> combinedMuon, isStandAloneMuon -> standAloneMuon, else -> tunePBestTrack)
    'promptMuonLLPNano': cms.InputTag("miniAODTrackProducer", "muonLLPNanoTracks"),
    'displacedMuonLLPNano': cms.InputTag("miniAODTrackProducer", "displacedMuonLLPNanoTracks"),
    # From MiniAODMuonEnhancedTracksProducer (muon-enhanced with sip2D selection)
    'sip2DMuonEnhanced': cms.InputTag("miniAODMuonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    'sip2DMuonEnhancedWithEle': cms.InputTag("miniAODMuonEnhancedTracks", "sip2DMuonEnhancedTracksWithEle"),
    'sip2DSlimmedDisplacedMuonEnhancedWithEle': cms.InputTag("miniAODMuonEnhancedTracks", "sip2DSlimmedDisplacedMuonEnhancedTracksWithEle"),
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


def configureMiniAODDiagnosticTrackCollection(process, trackCollection):
    """
    Configure the MiniAOD diagnostic producer and analyzer to use the same track collection.

    Usage:
        from KUCMSNtupleizer.HyddraSVProducer.miniAODTrackCollections_cfi import configureMiniAODDiagnosticTrackCollection
        configureMiniAODDiagnosticTrackCollection(process, 'sip2DMuonEnhanced')
    """
    if trackCollection not in MINIAOD_TRACK_COLLECTION_CONFIG:
        raise ValueError(f"Unknown MiniAOD track collection: {trackCollection}. "
                        f"Available: {list(MINIAOD_TRACK_COLLECTION_CONFIG.keys())}")

    inputTag = MINIAOD_TRACK_COLLECTION_CONFIG[trackCollection]

    if hasattr(process, 'hyddraSVsDiag'):
        process.hyddraSVsDiag.tracks = inputTag

    if hasattr(process, 'hyddraSVsDiagAnalyzer'):
        process.hyddraSVsDiagAnalyzer.tracks = inputTag

    print(f"Configured MiniAOD diagnostic track collection: {trackCollection} -> {inputTag}")


def configureMiniAODHadronicDiagnosticTrackCollection(process, trackCollection):
    """
    Configure the MiniAOD hadronic diagnostic producer and analyzer to use the same
    track collection.  Mirrors configureMiniAODDiagnosticTrackCollection() but targets
    hyddraSVsHadronicDiag and hyddraSVsHadronicDiagAnalyzer.

    Usage:
        from KUCMSNtupleizer.HyddraSVProducer.miniAODTrackCollections_cfi import configureMiniAODHadronicDiagnosticTrackCollection
        configureMiniAODHadronicDiagnosticTrackCollection(process, 'sip2DMuonEnhanced')
    """
    if trackCollection not in MINIAOD_TRACK_COLLECTION_CONFIG:
        raise ValueError(f"Unknown MiniAOD track collection: {trackCollection}. "
                        f"Available: {list(MINIAOD_TRACK_COLLECTION_CONFIG.keys())}")

    inputTag = MINIAOD_TRACK_COLLECTION_CONFIG[trackCollection]

    if hasattr(process, 'hyddraSVsHadronicDiag'):
        process.hyddraSVsHadronicDiag.tracks = inputTag

    if hasattr(process, 'hyddraSVsHadronicDiagAnalyzer'):
        process.hyddraSVsHadronicDiagAnalyzer.tracks = inputTag

    print(f"Configured MiniAOD hadronic diagnostic track collection: {trackCollection} -> {inputTag}")
