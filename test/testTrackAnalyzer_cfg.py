import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# Setup command line options
options = VarParsing.VarParsing('analysis')
options.register('hasGenInfo',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Process gen-level information (set False for data)")
options.register('trackCollection',
                 'sip2DMuonEnhanced',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection to analyze: general, selected, muon, sip2D, sip2DMuonEnhanced, muonEnhanced")
options.register('genMatchDeltaRCut',
                 0.02,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Delta R cut for gen-matching (default: 0.02)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'track_ntuple.root')
options.parseArguments()

process = cms.Process("TRACKANALYSIS")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles if options.inputFiles else ['file:input.root']
    )
)

# Load geometry and magnetic field
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos:
    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')

# Load TransientTrackBuilder
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Set the GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

# TFileService for output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

# ============================================================================
# MuonEnhancedTracks Producer (needed for most track collections)
# ============================================================================
process.load("KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi")

# ============================================================================
# Track Analyzer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.trackAnalyzer_cfi")
from KUCMSNtupleizer.HyddraSVProducer.trackAnalyzer_cfi import TRACK_COLLECTION_CONFIG

# Configure based on command line options
process.trackAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)
process.trackAnalyzer.genMatchDeltaRCut = cms.double(options.genMatchDeltaRCut)

# Set track collection based on option
if options.trackCollection in TRACK_COLLECTION_CONFIG:
    process.trackAnalyzer.tracks = TRACK_COLLECTION_CONFIG[options.trackCollection]
    print(f"Using track collection: {options.trackCollection}")
else:
    print(f"Warning: Unknown track collection '{options.trackCollection}'. Using default.")

# ============================================================================
# Path
# ============================================================================
# Build the path based on track collection
if options.trackCollection in ['general', 'muon']:
    # These collections don't need muonEnhancedTracks producer
    process.p = cms.Path(process.trackAnalyzer)
else:
    # Most collections need the muonEnhancedTracks producer
    process.p = cms.Path(
        process.muonEnhancedTracks +
        process.trackAnalyzer
    )

process.schedule = cms.Schedule(process.p)

# ============================================================================
# Usage examples:
# ============================================================================
# Default (sip2DMuonEnhanced tracks):
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root
#
# General tracks:
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root trackCollection=general
#
# With looser gen-matching:
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root genMatchDeltaRCut=0.05
#
# For data (no gen info):
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root hasGenInfo=False
