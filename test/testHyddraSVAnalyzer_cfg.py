import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# Setup command line options
options = VarParsing.VarParsing('analysis')
options.register('hasGenInfo',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Process gen-level information (set False for data)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraSV_ntuple.root')
options.parseArguments()

process = cms.Process("HYDDRAANALYSIS")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source - modify this path to your input file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles if options.inputFiles else ['file:input.root']
    )
)

# Load geometry and magnetic field (required for track-SC matching and track building)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Set the GlobalTag - modify for your data/MC era
from Configuration.AlCa.GlobalTag import GlobalTag
# For Run3 MC:
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')
# For Run2 UL MC, use something like:
# process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')

# TFileService for output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

# ============================================================================
# MuonEnhancedTracks Producer (from KUCMSNtupleizer)
# ============================================================================
process.load("KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi")

# ============================================================================
# HYDDRA SV Producer (uses sip2DMuonEnhancedTracks by default)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi")

# ============================================================================
# HYDDRA SV Analyzer (uses sip2DMuonEnhancedTracks by default)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi")
process.hyddraSVAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)

# ============================================================================
# Path: Run producer sequence then analyzer
# ============================================================================
process.p = cms.Path(
    process.muonEnhancedTracks +  # Produces sip2DMuonEnhancedTracks (and others)
    process.hyddraSVs +           # Produces leptonic/hadronic vertices
    process.hyddraSVAnalyzer      # Writes TTree output
)

process.schedule = cms.Schedule(process.p)

# ============================================================================
# Optional: Customize input tags if your collections have different names
# ============================================================================
# process.muonEnhancedTracks.generalTracks = cms.InputTag("generalTracks")
# process.muonEnhancedTracks.dsaMuonTracks = cms.InputTag("displacedStandAloneMuons")
# process.muonEnhancedTracks.dgMuonTracks = cms.InputTag("displacedGlobalMuons")
# process.muonEnhancedTracks.displacedTracks = cms.InputTag("displacedTracks")
# process.muonEnhancedTracks.displacedMuons = cms.InputTag("displacedMuons")
# process.muonEnhancedTracks.pvCollection = cms.InputTag("offlinePrimaryVertices")

# process.hyddraSVAnalyzer.leptonicVertices = cms.InputTag("hyddraSVs", "leptonicVertices")
# process.hyddraSVAnalyzer.hadronicVertices = cms.InputTag("hyddraSVs", "hadronicVertices")
# process.hyddraSVAnalyzer.muonTracks = cms.InputTag("displacedGlobalMuons")
# process.hyddraSVAnalyzer.mergedSCs = cms.InputTag("mergedSuperClusters")
# process.hyddraSVAnalyzer.genParticles = cms.InputTag("genParticles")
