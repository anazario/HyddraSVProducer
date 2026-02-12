import FWCore.ParameterSet.Config as cms

filteredTrackProducer = cms.EDProducer("FilteredTrackProducer",
    tracks = cms.InputTag("generalTracks"),
    pvCollection = cms.InputTag("offlinePrimaryVertices"),
    # Quality cuts (optimized from signal vs background comparison)
    minPt = cms.double(1.0),
    maxPtResolution = cms.double(0.5),
    maxNormalizedChi2 = cms.double(5.0),
    # sip2D cut: invertSip2DCut=False keeps |sip2D| < sip2DCut (prompt-like)
    #            invertSip2DCut=True keeps |sip2D| >= sip2DCut (displaced)
    sip2DCut = cms.double(4.0),
    invertSip2DCut = cms.bool(True),
)
