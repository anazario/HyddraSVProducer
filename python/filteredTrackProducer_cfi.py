import FWCore.ParameterSet.Config as cms

filteredTrackProducer = cms.EDProducer("FilteredTrackProducer",
    tracks = cms.InputTag("generalTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    # Quality cuts (optimized from signal vs background comparison)
    minPt = cms.double(5.0),
    maxPtResolution = cms.double(0.08),
    maxNormalizedChi2 = cms.double(5.0),
    maxAbsSip2D = cms.double(5.0),
)
