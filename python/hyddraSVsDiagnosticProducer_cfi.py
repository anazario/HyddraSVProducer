import FWCore.ParameterSet.Config as cms

from KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi import leptonicHYDDRA

# =============================================================================
# Default diagnostic producer — promptMuonBestTrack collection
# =============================================================================
hyddraSVsDiag = cms.EDProducer("HyddraSVsDiagnosticProducer",

    tracks       = cms.InputTag("muonBestTrackProducer", "globalTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),

    leptonic = leptonicHYDDRA,
)
