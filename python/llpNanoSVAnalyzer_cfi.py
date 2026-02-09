import FWCore.ParameterSet.Config as cms

# =============================================================================
# LLPNanoSVAnalyzer: reads LLPNanoAOD flat ROOT files directly (Events TTree)
# and writes a TTree with the same branch structure as HyddraSVAnalyzer,
# enabling direct efficiency comparison between HYDDRA and the LLPNanoAOD
# dimuon vertex tools.
#
# Uses EmptySource because LLPNanoAOD files may not be PoolSource-compatible
# across CMSSW versions.  Input files are passed to the analyzer PSet.
# =============================================================================

llpNanoSVAnalyzer = cms.EDAnalyzer("LLPNanoSVAnalyzer",
    collection = cms.string("PatMuonVertex"),
    motherPdgId = cms.int32(54),
    dsaDeltaR = cms.double(0.3),
    dsaRelPtDiff = cms.double(0.5),
    inputFiles = cms.vstring(),
)

# =============================================================================
# PatDSAMuonVertex variant
# =============================================================================
llpNanoSVAnalyzerDSA = llpNanoSVAnalyzer.clone(
    collection = "PatDSAMuonVertex",
)
