import FWCore.ParameterSet.Config as cms

# Track Associator Parameters for SC matching
TrackAssociatorParameterBlock = cms.PSet(
    TrackAssociatorParameters = cms.PSet(
        useEcal = cms.bool(True),
        useHcal = cms.bool(False),
        useHO = cms.bool(False),
        useCalo = cms.bool(False),
        useMuon = cms.bool(False),
        usePreshower = cms.bool(False),
        dREcal = cms.double(1.0),
        dRHcal = cms.double(1.0),
        dRMuon = cms.double(1.0),
        dRPreshowerPreselection = cms.double(0.2),
        dRMuonPreselection = cms.double(0.2),
        dREcalPreselection = cms.double(0.5),
        dRHcalPreselection = cms.double(0.5),
        accountForTrajectoryChangeCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit", "EcalRecHitsEE"),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
    )
)

hyddraSVAnalyzer = cms.EDAnalyzer("HyddraSVAnalyzer",
    hasGenInfo = cms.bool(True),
    leptonicVertices = cms.InputTag("hyddraSVs", "leptonicVertices"),
    hadronicVertices = cms.InputTag("hyddraSVs", "hadronicVertices"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    tracks = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    muonTracks = cms.InputTag("displacedGlobalMuons"),
    mergedSCs = cms.InputTag("mergedSuperClusters"),
    genParticles = cms.InputTag("genParticles"),
    TrackAssociatorParameters = TrackAssociatorParameterBlock.TrackAssociatorParameters.clone(),
)

# Alternative configuration for data (no gen info)
hyddraSVAnalyzerData = hyddraSVAnalyzer.clone(
    hasGenInfo = cms.bool(False),
)
