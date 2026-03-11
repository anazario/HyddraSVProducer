import FWCore.ParameterSet.Config as cms

from KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi import hadronicHYDDRA

# =============================================================================
# Default hadronic diagnostic producer — promptMuonBestTrack collection
# =============================================================================
hyddraSVsHadronicDiag = cms.EDProducer("HyddraSVsHadronicDiagnosticProducer",

    tracks       = cms.InputTag("muonBestTrackProducer", "globalTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),

    hadronic = hadronicHYDDRA,
)
