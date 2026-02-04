import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# Setup command line options
options = VarParsing.VarParsing('analysis')
options.register('hasGenInfo',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Process gen-level information (set False for data)")
options.register('processMode',
                 'both',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Processing mode: both (default), leptonic, or hadronic")
options.register('trackCollection',
                 'sip2DMuonEnhanced',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection: general, selected, muon, sip2D, sip2DMuonEnhanced (default), muonEnhanced")
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
if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos:
    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')

process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

# Load TransientTrackBuilder (required for MuonEnhancedTracksProducer and vertex fitting)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

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
# ECALTracks Producer (produces displacedElectronSCs for SC matching)
# ============================================================================
process.load("KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi")

# ============================================================================
# HYDDRA SV Producer (uses sip2DMuonEnhancedTracks by default)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi")

# ============================================================================
# HYDDRA SV Analyzer (uses sip2DMuonEnhancedTracks by default)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi")
process.hyddraSVAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)
process.hyddraSVAnalyzer.mergedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs")

# ============================================================================
# Configure track collection for BOTH producer and analyzer
# ============================================================================
# IMPORTANT: Producer and analyzer MUST use the same track collection!
from KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi import configureTrackCollection
configureTrackCollection(process, options.trackCollection)

# Configure processing mode: both (default), leptonic, or hadronic
if options.processMode == 'leptonic':
    # Only process leptonic vertices
    process.hyddraSVAnalyzer.hadronicVertices = cms.InputTag("")
elif options.processMode == 'hadronic':
    # Only process hadronic vertices
    process.hyddraSVAnalyzer.leptonicVertices = cms.InputTag("")
# else: process both (default)

# ============================================================================
# Path: Run producer sequence then analyzer
# ============================================================================
process.p = cms.Path(
    process.muonEnhancedTracks +  # Produces sip2DMuonEnhancedTracks (and others)
    process.ecalTracks +          # Produces displacedElectronSCs for SC matching
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

# ============================================================================
# Usage examples:
# ============================================================================
# Default (sip2DMuonEnhanced tracks for both producer and analyzer):
#   cmsRun testHyddraSVAnalyzer_cfg.py inputFiles=file:myinput.root
#
# General tracks (configures BOTH producer and analyzer):
#   cmsRun testHyddraSVAnalyzer_cfg.py inputFiles=file:myinput.root trackCollection=general
#
# Muon-enhanced tracks (before sip2D cut):
#   cmsRun testHyddraSVAnalyzer_cfg.py inputFiles=file:myinput.root trackCollection=muonEnhanced
#
# Sip2D tracks only:
#   cmsRun testHyddraSVAnalyzer_cfg.py inputFiles=file:myinput.root trackCollection=sip2D
#
# Leptonic only with general tracks:
#   cmsRun testHyddraSVAnalyzer_cfg.py inputFiles=file:myinput.root trackCollection=general processMode=leptonic
#
# For data (no gen info):
#   cmsRun testHyddraSVAnalyzer_cfg.py inputFiles=file:myinput.root hasGenInfo=False
#
# NOTE: The trackCollection option configures BOTH the SV producer (hyddraSVs)
# and the analyzer (hyddraSVAnalyzer) to use the same track collection. This is
# required for gen-matching to work correctly.
