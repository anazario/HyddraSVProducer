import FWCore.ParameterSet.Config as cms

# =============================================================================
# LLPNanoSVAnalyzer: reads LLPNanoAOD FlatTables and writes a TTree with the
# same branch structure as HyddraSVAnalyzer, enabling direct efficiency
# comparison between HYDDRA and the LLPNanoAOD dimuon vertex tools.
# =============================================================================

llpNanoSVAnalyzer = cms.EDAnalyzer("LLPNanoSVAnalyzer",
    collection = cms.string("PatMuonVertex"),
    motherPdgId = cms.int32(54),
    dsaDeltaR = cms.double(0.3),
    dsaRelPtDiff = cms.double(0.5),
    genPartTable = cms.InputTag("genParticleTable"),
    muonTable = cms.InputTag("muonTable"),
    dsaMuonTable = cms.InputTag("dsaMuonTable"),
    vertexTable = cms.InputTag("patMuonVertexTable"),
    refittedTracksTable = cms.InputTag("patMuonVertexRefittedTracksTable"),
)

# =============================================================================
# PatDSAMuonVertex variant
# =============================================================================
llpNanoSVAnalyzerDSA = llpNanoSVAnalyzer.clone(
    collection = "PatDSAMuonVertex",
    vertexTable = cms.InputTag("patDSAMuonVertexTable"),
    refittedTracksTable = cms.InputTag("patDSAMuonVertexRefittedTracksTable"),
)
