import FWCore.ParameterSet.Config as cms

leptonicHYDDRA = cms.PSet(
    # Seeding
    seedCosThetaCut = cms.double(0.75),

    # Shared kinematic cuts
    minMass            = cms.double(2.0),
    minPOverE          = cms.double(0.6),
    maxNormChi2        = cms.double(5.0),
    minDxySignificance = cms.double(25.0),

    # Cleaning (track removal for vertices with > 2 tracks)
    maxCompatibility = cms.double(1.5),
    minCleanCosTheta = cms.double(0.5),

    # Post-cleaning kinematic cuts
    minPostCosTheta = cms.double(0.75),

    # Final filtering (post-disambiguation, 2-track vertices only)
    minTrackCosTheta         = cms.double(0.5),
    maxTrackCosThetaCM_Limit = cms.double(0.95),
    maxTrackCosThetaCM_Slope = cms.double(1.8),
)

hadronicHYDDRA = cms.PSet(
    # Seeding
    seedCosThetaCut = cms.double(0.0),

    # Shared kinematic cuts
    minMass            = cms.double(2.0),
    minPOverE          = cms.double(0.6),
    maxNormChi2        = cms.double(5.0),
    minDxySignificance = cms.double(40.0),

    # Hadronic vertex cuts
    minSize       = cms.int32(5),
    minCosTheta   = cms.double(0.0),
    maxDecayAngle = cms.double(0.9),
)

hyddraSVs = cms.EDProducer("HyddraSVsProducer",
    tracks       = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    leptonic     = leptonicHYDDRA,
    hadronic     = hadronicHYDDRA,
)
