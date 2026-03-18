import FWCore.ParameterSet.Config as cms

from KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi import hyddraSVAnalyzer, tkAssocParamBlock
from KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi import leptonicHYDDRA

# =============================================================================
# Default diagnostic producer — shares TrackAssociatorParameters with analyzer
# =============================================================================
hyddraSVsDiag = cms.EDProducer("HyddraSVsDiagnosticProducer",

    tracks       = cms.InputTag("muonBestTrackProducer", "globalTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    muonTracks   = cms.InputTag("muonEnhancedTracks", "muonGlobalTracks"),
    mergedSCs    = cms.InputTag("ecalTracks", "displacedElectronSCs"),

    TrackAssociatorParameters = tkAssocParamBlock.TrackAssociatorParameters.clone(),

    leptonic = leptonicHYDDRA,
)
