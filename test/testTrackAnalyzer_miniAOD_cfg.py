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
                 'merged',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection: pf, lost, eleLost, merged (default), mergedWithEle")
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
options.setDefault('outputFile', 'track_miniAOD_ntuple.root')
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
# MiniAOD Track Producer (extracts reco::Track from packed candidates)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.miniAODTrackProducer_cfi")

# ============================================================================
# Track Analyzer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.trackAnalyzer_cfi")

# Configure based on command line options
process.trackAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)
process.trackAnalyzer.genMatchDeltaRCut = cms.double(options.genMatchDeltaRCut)
# Override PV collection for MiniAOD
process.trackAnalyzer.pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
# Use pruned genParticles for MiniAOD
process.trackAnalyzer.genParticles = cms.InputTag("prunedGenParticles")

# Set track collection based on option
from KUCMSNtupleizer.HyddraSVProducer.miniAODTrackCollections_cfi import MINIAOD_TRACK_COLLECTION_CONFIG
if options.trackCollection in MINIAOD_TRACK_COLLECTION_CONFIG:
    process.trackAnalyzer.tracks = MINIAOD_TRACK_COLLECTION_CONFIG[options.trackCollection]
    print(f"Using MiniAOD track collection: {options.trackCollection}")
else:
    print(f"Warning: Unknown track collection '{options.trackCollection}'. Using default 'merged'.")
    process.trackAnalyzer.tracks = MINIAOD_TRACK_COLLECTION_CONFIG['merged']

# ============================================================================
# Path
# ============================================================================
process.p = cms.Path(
    process.miniAODTrackProducer +  # Produces track collections from packed candidates
    process.trackAnalyzer            # Writes TTree output
)

process.schedule = cms.Schedule(process.p)

# ============================================================================
# Usage examples:
# ============================================================================
# Single file:
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root
#
# Multiple files (comma-separated):
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:file1.root,file:file2.root
#
# From a file list (text file with one path per line):
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFileList=myfiles.txt
#
# Different track collections:
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=pf
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=lost
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=merged
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=mergedWithEle
#
# With looser gen-matching:
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root genMatchDeltaRCut=0.05
#
# For data (no gen info):
#   cmsRun testTrackAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root hasGenInfo=False
