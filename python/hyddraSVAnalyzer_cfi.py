import FWCore.ParameterSet.Config as cms

# Import the standard TrackAssociatorParameters and customize for AOD
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

# Customize for AOD (reduced rechits)
tkAssocParamBlock = TrackAssociatorParameterBlock.clone()
tkAssocParamBlock.TrackAssociatorParameters.useMuon = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useCalo = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useHO = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.usePreshower = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEE")
tkAssocParamBlock.TrackAssociatorParameters.EBRecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEB")
tkAssocParamBlock.TrackAssociatorParameters.HBHERecHitCollectionLabel = cms.InputTag("reducedHcalRecHits", "hbhereco")
tkAssocParamBlock.TrackAssociatorParameters.HORecHitCollectionLabel = cms.InputTag("reducedHcalRecHits", "horeco")

hyddraSVAnalyzer = cms.EDAnalyzer("HyddraSVAnalyzer",
    hasGenInfo = cms.bool(True),
    leptonicVertices = cms.InputTag("hyddraSVs", "leptonicVertices"),
    hadronicVertices = cms.InputTag("hyddraSVs", "hadronicVertices"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    tracks = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    muonTracks = cms.InputTag("displacedGlobalMuons"),
    mergedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs"),
    genParticles = cms.InputTag("genParticles"),
    TrackAssociatorParameters = tkAssocParamBlock.TrackAssociatorParameters.clone(),
)

# Alternative configuration for data (no gen info)
hyddraSVAnalyzerData = hyddraSVAnalyzer.clone(
    hasGenInfo = cms.bool(False),
)
