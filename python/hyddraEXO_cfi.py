import FWCore.ParameterSet.Config as cms

# Universal leptonic HYDDRA PSet for EXONanoAOD.
# No signal-dependent cuts. All kinematic variables stored in the ntuple.
leptonicHYDDRAEXO = cms.PSet(
    # Seeding: only chi2 gate; cosTheta cut disabled for universality
    seedCosThetaCut  = cms.double(-1.0),
    maxNormChi2      = cms.double(5.0),

    # The following parameters are unused (doCleaning=False, doFiltering=False)
    # but required by LeptonicHYDDRA constructor
    minMass                      = cms.double(0.0),
    minPOverE                    = cms.double(0.0),
    minDxySignificance           = cms.double(0.0),
    maxCompatibility             = cms.double(999.0),
    minCleanCosTheta             = cms.double(-1.0),
    maxCleanCosTheta             = cms.double(1.0),
    invertCleanCosThetaCut       = cms.bool(False),
    useDiagonalCut               = cms.bool(False),
    cleanCutSlope                = cms.double(0.0),
    minTrackCosTheta             = cms.double(-1.0),
    maxTrackCosThetaCM_Limit     = cms.double(1.0),
    maxTrackCosThetaCM_Intercept = cms.double(999.0),
    trackCosThetaCM_Slope        = cms.double(0.0),
    requireChargeNeutrality      = cms.bool(False),
    minVtxCosTheta               = cms.double(-1.0),
    maxVtxCosTheta               = cms.double(1.0),
    useAbsVtxCosTheta            = cms.bool(False),
    maxVtxDecayAngle             = cms.double(1.0),
    useAbsVtxDecayAngle          = cms.bool(False),
    applyVtxDecayAngleCleaning   = cms.bool(False),
    applyVtxDecayAngleFiltering  = cms.bool(False),
    minMassFilter                = cms.double(0.0),
    minBetaFilter                = cms.double(0.0),

    # Pipeline: merging and disambiguation enabled; cleaning and filtering disabled
    doMerging        = cms.bool(True),
    doCleaning       = cms.bool(False),
    doDisambiguation = cms.bool(True),
    doFiltering      = cms.bool(False),
    useVertexSmoothing = cms.bool(False),
)

hyddraEXO = cms.EDProducer("HyddraSVsEXOProducer",
    tracks       = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    leptonic     = leptonicHYDDRAEXO,
)
