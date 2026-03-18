import FWCore.ParameterSet.Config as cms

from KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi import hyddraSVAnalyzer

# =============================================================================
# Default diagnostic producer — clone TrackAssociatorParameters from analyzer
# =============================================================================
hyddraSVsDiag = cms.EDProducer("HyddraSVsDiagnosticProducer",

    tracks       = hyddraSVAnalyzer.tracks.clone(),
    pvCollection = hyddraSVAnalyzer.pvCollection.clone(),
    muonTracks   = hyddraSVAnalyzer.muonTracks.clone(),
    mergedSCs    = hyddraSVAnalyzer.mergedSCs.clone(),

    TrackAssociatorParameters = hyddraSVAnalyzer.TrackAssociatorParameters.clone(),

    leptonic = hyddraSVAnalyzer.leptonic.clone(),
)
