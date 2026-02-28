import FWCore.ParameterSet.Config as cms

from KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi import TRACK_COLLECTION_CONFIG
from KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi import leptonicHYDDRA

# =============================================================================
# Default diagnostic analyzer — promptMuonBestTrack collection
# =============================================================================
hyddraSVsDiagAnalyzer = cms.EDAnalyzer("HyddraSVsDiagnosticAnalyzer",

    hasGenInfo = cms.bool(True),
    isFullAOD  = cms.bool(True),
    genMatchDeltaRCut = cms.double(0.02),

    # Intermediate stage collections from HyddraSVsDiagnosticProducer
    seedVertices         = cms.InputTag("hyddraSVsDiag", "leptonicSeeds"),
    mergedVertices       = cms.InputTag("hyddraSVsDiag", "leptonicMerged"),
    cleanedVertices      = cms.InputTag("hyddraSVsDiag", "leptonicCleaned"),
    disambiguatedVertices= cms.InputTag("hyddraSVsDiag", "leptonicDisambiguated"),
    filteredVertices     = cms.InputTag("hyddraSVsDiag", "leptonicFiltered"),

    pvCollection     = cms.InputTag("offlinePrimaryVertices"),
    tracks           = TRACK_COLLECTION_CONFIG['promptMuonBestTrack'],
    genParticles     = cms.InputTag("genParticles"),
    packedGenParticles = cms.InputTag(""),

    leptonic = leptonicHYDDRA,
)

# =============================================================================
# Data variant (no gen info)
# =============================================================================
hyddraSVsDiagAnalyzerData = hyddraSVsDiagAnalyzer.clone(
    hasGenInfo = cms.bool(False),
)


def configureDiagnosticTrackCollection(process, trackCollection):
    """
    Configure the diagnostic producer and analyzer to use the same track collection.

    Extends configureTrackCollection() from hyddraSVAnalyzer_cfi to also handle
    the diagnostic modules.  Can be called independently or together with the
    main configureTrackCollection().

    Available collections: see TRACK_COLLECTION_CONFIG in hyddraSVAnalyzer_cfi.py

    Usage:
        from KUCMSNtupleizer.HyddraSVProducer.hyddraSVsDiagnosticAnalyzer_cfi \\
            import configureDiagnosticTrackCollection
        configureDiagnosticTrackCollection(process, 'displacedMuonBestTrack')
    """
    if trackCollection not in TRACK_COLLECTION_CONFIG:
        raise ValueError(
            f"Unknown track collection: {trackCollection}. "
            f"Available: {list(TRACK_COLLECTION_CONFIG.keys())}"
        )

    inputTag = TRACK_COLLECTION_CONFIG[trackCollection]

    if hasattr(process, 'hyddraSVsDiag'):
        process.hyddraSVsDiag.tracks = inputTag

    if hasattr(process, 'hyddraSVsDiagAnalyzer'):
        process.hyddraSVsDiagAnalyzer.tracks = inputTag

    print(f"[DiagnosticConfig] Track collection: {trackCollection} -> {inputTag}")
