import FWCore.ParameterSet.Config as cms

# Import the standard TrackAssociatorParameters
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

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
