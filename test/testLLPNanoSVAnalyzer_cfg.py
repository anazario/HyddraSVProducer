import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

# Setup command line options
options = VarParsing.VarParsing('analysis')
options.register('collection',
                 'PatMuonVertex',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Vertex collection: PatMuonVertex (default) or PatDSAMuonVertex")
options.register('motherPdgId',
                 54,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Mother particle PDG ID for gen vertex identification (default: 54)")
options.register('inputFileList',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to text file containing list of input files (one per line)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'llpNanoSV_ntuple.root')
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

process = cms.Process("LLPNANOANALYSIS")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# EmptySource: the analyzer reads NanoAOD files directly via TFile::Open.
# Each CMSSW "event" processes one input file, giving per-file progress.
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(len(inputFiles))
)

process.source = cms.Source("EmptySource",
    numberEventsInLuminosityBlock = cms.untracked.uint32(max(len(inputFiles), 1)),
    firstRun = cms.untracked.uint32(1),
    firstLuminosityBlock = cms.untracked.uint32(1),
    firstEvent = cms.untracked.uint32(1),
)

# TFileService for output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

# ============================================================================
# LLPNanoSV Analyzer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.llpNanoSVAnalyzer_cfi")

# Configure based on collection choice
if options.collection == 'PatDSAMuonVertex':
    process.llpNanoSVAnalyzer = process.llpNanoSVAnalyzerDSA.clone()
# else: use the default PatMuonVertex configuration

# Override motherPdgId if specified
process.llpNanoSVAnalyzer.motherPdgId = cms.int32(options.motherPdgId)

# Pass input files to the analyzer (it reads NanoAOD TTrees directly)
process.llpNanoSVAnalyzer.inputFiles = cms.vstring(inputFiles)

# ============================================================================
# Path: just the analyzer (no producers needed for NanoAOD)
# ============================================================================
process.p = cms.Path(process.llpNanoSVAnalyzer)
process.schedule = cms.Schedule(process.p)

# ============================================================================
# Usage examples:
# ============================================================================
# PatMuonVertex (default):
#   cmsRun testLLPNanoSVAnalyzer_cfg.py inputFiles=file:output20.root
#
# PatDSAMuonVertex:
#   cmsRun testLLPNanoSVAnalyzer_cfg.py inputFiles=file:output20.root collection=PatDSAMuonVertex
#
# Remote file via xrootd:
#   cmsRun testLLPNanoSVAnalyzer_cfg.py inputFiles=root://xrootd.path/sample.root
#
# Multiple files from a file list:
#   cmsRun testLLPNanoSVAnalyzer_cfg.py inputFileList=myfiles.txt collection=PatMuonVertex
#
# Custom output file name:
#   cmsRun testLLPNanoSVAnalyzer_cfg.py inputFiles=file:output20.root outputFile=myoutput.root
#
# Compatible with parallelRun.sh:
#   ./scripts/parallelRun.sh files.txt -c testLLPNanoSVAnalyzer_cfg.py
