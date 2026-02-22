import FWCore.ParameterSet.Config as cms

# =============================================================================
# Default diagnostic producer — promptMuonBestTrack collection
# =============================================================================
hyddraSVsDiag = cms.EDProducer("HyddraSVsDiagnosticProducer",

    tracks       = cms.InputTag("muonBestTrackProducer", "globalTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),

    # Leptonic algorithm parameters — mirror defaults in HyddraSVsProducer
    leptonic = cms.PSet(
        seedCosThetaCut          = cms.double(0.75),
        minMass                  = cms.double(2.0),
        minPOverE                = cms.double(0.6),
        maxNormChi2              = cms.double(5.0),
        minDxySignificance       = cms.double(25.0),
        maxCompatibility         = cms.double(1.5),
        minCleanCosTheta         = cms.double(0.5),
        useDiagonalCut           = cms.bool(False),
        cleanCutSlope            = cms.double(0.0),
        minTrackCosTheta         = cms.double(0.5),
        maxTrackCosThetaCM_Limit = cms.double(0.95),
        maxTrackCosThetaCM_Slope = cms.double(1.8),
    ),
)
