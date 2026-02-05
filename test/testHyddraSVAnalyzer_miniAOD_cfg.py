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
options.register('processMode',
                 'both',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Processing mode: both (default), leptonic, or hadronic")
options.register('trackCollection',
                 'merged',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection: pf, lost, eleLost, merged (default), mergedWithEle")
options.register('inputFileList',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to text file containing list of input files (one per line)")
options.register('applyCuts',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Apply track quality cuts (default: False)")
options.register('minPt',
                 1.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Minimum track pT in GeV (default: 1.0)")
options.register('minAbsSip2D',
                 4.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Minimum |sip2D| for displaced track selection (default: 4.0)")
options.register('maxNormalizedChi2',
                 5.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Maximum normalized chi2 (default: 5.0)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraSV_miniAOD_ntuple.root')
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

process = cms.Process("HYDDRAANALYSIS")

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

# Load geometry and magnetic field (required for track-SC matching and track building)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos:
    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')

process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

# Load TransientTrackBuilder (required for vertex fitting)
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
# MiniAOD Track Producer (extracts reco::Track from packed candidates)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.miniAODTrackProducer_cfi")

# Configure track cuts
process.miniAODTrackProducer.applyCuts = cms.bool(options.applyCuts)
process.miniAODTrackProducer.minPt = cms.double(options.minPt)
process.miniAODTrackProducer.minAbsSip2D = cms.double(options.minAbsSip2D)
process.miniAODTrackProducer.maxNormalizedChi2 = cms.double(options.maxNormalizedChi2)

if options.applyCuts:
    print(f"Track cuts enabled: pT > {options.minPt} GeV, |sip2D| >= {options.minAbsSip2D}, chi2/ndof < {options.maxNormalizedChi2}")

# ============================================================================
# HYDDRA SV Producer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi")
# Override PV collection for MiniAOD
process.hyddraSVs.pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices")

# ============================================================================
# HYDDRA SV Analyzer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi")
process.hyddraSVAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)
# Override PV collection for MiniAOD
process.hyddraSVAnalyzer.pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
# Use pruned genParticles for MiniAOD
process.hyddraSVAnalyzer.genParticles = cms.InputTag("prunedGenParticles")
# Use slimmed muons for MiniAOD
process.hyddraSVAnalyzer.muonTracks = cms.InputTag("slimmedMuons")
# Disable SC matching for now (would need MiniAOD-compatible producer)
process.hyddraSVAnalyzer.mergedSCs = cms.InputTag("")

# ============================================================================
# Configure track collection for BOTH producer and analyzer
# ============================================================================
from KUCMSNtupleizer.HyddraSVProducer.miniAODTrackCollections_cfi import configureMiniAODTrackCollection
configureMiniAODTrackCollection(process, options.trackCollection)

# Configure processing mode: both (default), leptonic, or hadronic
if options.processMode == 'leptonic':
    # Only process leptonic vertices
    process.hyddraSVAnalyzer.hadronicVertices = cms.InputTag("")
elif options.processMode == 'hadronic':
    # Only process hadronic vertices
    process.hyddraSVAnalyzer.leptonicVertices = cms.InputTag("")
# else: process both (default)

# ============================================================================
# Path: Run MiniAOD track producer, then SV producer, then analyzer
# ============================================================================
process.p = cms.Path(
    process.miniAODTrackProducer +  # Produces track collections from packed candidates
    process.hyddraSVs +              # Produces leptonic/hadronic vertices
    process.hyddraSVAnalyzer         # Writes TTree output
)

process.schedule = cms.Schedule(process.p)

# ============================================================================
# Usage examples:
# ============================================================================
# Single file:
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root
#
# Multiple files (comma-separated):
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:file1.root,file:file2.root
#
# From a file list (text file with one path per line):
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFileList=myfiles.txt
#
# Different track collections:
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=pf
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=lost
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=merged
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=mergedWithEle
#
# With track cuts enabled (selects displaced tracks):
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root applyCuts=True
#
# With custom cut values:
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root applyCuts=True minPt=2.0 minAbsSip2D=3.0
#
# Leptonic only:
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root processMode=leptonic
#
# For data (no gen info):
#   cmsRun testHyddraSVAnalyzer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root hasGenInfo=False
