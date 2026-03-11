import FWCore.ParameterSet.Config as cms

from KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi import TRACK_COLLECTION_CONFIG
from KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi import hadronicHYDDRA

# =============================================================================
# Default hadronic diagnostic analyzer — promptMuonBestTrack collection
# =============================================================================
hyddraSVsHadronicDiagAnalyzer = cms.EDAnalyzer("HyddraSVsHadronicDiagnosticAnalyzer",

    hasGenInfo = cms.bool(True),
    isFullAOD  = cms.bool(True),
    genMatchDeltaRCut = cms.double(0.02),

    # Intermediate stage collections from HyddraSVsHadronicDiagnosticProducer
    seedVertices          = cms.InputTag("hyddraSVsHadronicDiag", "hadronicSeeds"),
    mergedVertices        = cms.InputTag("hyddraSVsHadronicDiag", "hadronicMerged"),
    cleanedVertices       = cms.InputTag("hyddraSVsHadronicDiag", "hadronicCleaned"),
    disambiguatedVertices = cms.InputTag("hyddraSVsHadronicDiag", "hadronicDisambiguated"),
    filteredVertices      = cms.InputTag("hyddraSVsHadronicDiag", "hadronicFiltered"),

    pvCollection       = cms.InputTag("offlinePrimaryVertices"),
    tracks             = TRACK_COLLECTION_CONFIG['promptMuonBestTrack'],
    genParticles       = cms.InputTag("genParticles"),
    packedGenParticles = cms.InputTag(""),

    hadronic = hadronicHYDDRA,
)

# =============================================================================
# Data variant (no gen info)
# =============================================================================
hyddraSVsHadronicDiagAnalyzerData = hyddraSVsHadronicDiagAnalyzer.clone(
    hasGenInfo = cms.bool(False),
)


def configureHadronicDiagnosticTrackCollection(process, trackCollection):
    """
    Configure the hadronic diagnostic producer and analyzer to use the same track collection.

    Mirrors configureDiagnosticTrackCollection() from hyddraSVsDiagnosticAnalyzer_cfi but
    targets the hadronic diagnostic modules (hyddraSVsHadronicDiag and
    hyddraSVsHadronicDiagAnalyzer).

    Available collections: see TRACK_COLLECTION_CONFIG in hyddraSVAnalyzer_cfi.py

    Usage:
        from KUCMSNtupleizer.HyddraSVProducer.hyddraSVsHadronicDiagnosticAnalyzer_cfi \\
            import configureHadronicDiagnosticTrackCollection
        configureHadronicDiagnosticTrackCollection(process, 'displacedMuonBestTrack')
    """
    if trackCollection not in TRACK_COLLECTION_CONFIG:
        raise ValueError(
            f"Unknown track collection: {trackCollection}. "
            f"Available: {list(TRACK_COLLECTION_CONFIG.keys())}"
        )

    inputTag = TRACK_COLLECTION_CONFIG[trackCollection]

    if hasattr(process, 'hyddraSVsHadronicDiag'):
        process.hyddraSVsHadronicDiag.tracks = inputTag

    if hasattr(process, 'hyddraSVsHadronicDiagAnalyzer'):
        process.hyddraSVsHadronicDiagAnalyzer.tracks = inputTag

    print(f"[HadronicDiagnosticConfig] Track collection: {trackCollection} -> {inputTag}")
