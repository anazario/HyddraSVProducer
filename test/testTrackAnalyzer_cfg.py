import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

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
                 "Track collection to analyze: general, generalFiltered, selected, displacedGlobalMuon, globalMuon, displacedMuonGlobal, sip2D, sip2DMuonEnhanced, muonEnhanced")
options.register('genMatchDeltaRCut',
                 0.02,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Delta R cut for gen-matching (default: 0.02)")
options.register('inputFileList',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to text file containing list of input files (one per line)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'track_ntuple.root')
options.parseArguments()

# Build input file list
def getInputFiles():
    files = []
    # First check if a file list was provided
    if options.inputFileList:
        if os.path.exists(options.inputFileList):
            with open(options.inputFileList, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        files.append(line)
            print(f"Loaded {len(files)} files from {options.inputFileList}")
        else:
            raise FileNotFoundError(f"Input file list not found: {options.inputFileList}")
    # Then add any files from inputFiles option
    if options.inputFiles:
        files.extend(options.inputFiles)
    # Default fallback
    if not files:
        files = ['file:input.root']
    return files

inputFiles = getInputFiles()

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
    fileNames = cms.untracked.vstring(inputFiles)
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
# MuonGlobalTrackProducer (extracts global tracks from AOD muon collection)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.muonGlobalTrackProducer_cfi")

# ============================================================================
# FilteredTrackProducer (applies quality cuts to general tracks)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.filteredTrackProducer_cfi")

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
if options.trackCollection in ['general', 'displacedGlobalMuon']:
    # These collections don't need any producer
    process.p = cms.Path(process.trackAnalyzer)
elif options.trackCollection in ['globalMuon', 'displacedMuonGlobal']:
    # globalMuon and displacedMuonGlobal need the muonGlobalTrackProducer
    process.p = cms.Path(
        process.muonGlobalTrackProducer +
        process.trackAnalyzer
    )
elif options.trackCollection == 'generalFiltered':
    # generalFiltered needs the filteredTrackProducer
    process.p = cms.Path(
        process.filteredTrackProducer +
        process.trackAnalyzer
    )
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
# Single file:
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root
#
# Multiple files (comma-separated):
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:file1.root,file:file2.root,file:file3.root
#
# Multiple files (repeated option):
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:file1.root inputFiles=file:file2.root
#
# From a file list (text file with one path per line):
#   cmsRun testTrackAnalyzer_cfg.py inputFileList=myfiles.txt
#
# General tracks:
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root trackCollection=general
#
# With looser gen-matching:
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root genMatchDeltaRCut=0.05
#
# For data (no gen info):
#   cmsRun testTrackAnalyzer_cfg.py inputFiles=file:myinput.root hasGenInfo=False
