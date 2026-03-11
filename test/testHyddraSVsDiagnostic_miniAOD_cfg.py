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
                 "Track collection: pf, lost, eleLost, merged, mergedWithEle, mergedAll, "
                 "promptMuonExtracted, displacedMuonExtracted, "
                 "sip2DMuonEnhanced (default), sip2DMuonEnhancedWithEle, "
                 "sip2DSlimmedDisplacedMuonEnhancedWithEle")
options.register('processMode',
                 'leptonic',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Processing mode: leptonic (default), hadronic, or both")
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
options.setDefault('outputFile', 'hyddra_diagnostic_miniAOD.root')
options.parseArguments()

# Build input file list
def getInputFiles():
    files = []
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
    if options.inputFiles:
        files.extend(options.inputFiles)
    if not files:
        files = ['file:input.root']
    return files

inputFiles = getInputFiles()

process = cms.Process("HYDDRADIAG")

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

# Geometry, magnetic field, and conditions (required for vertex fitting)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos:
    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')

# TransientTrackBuilder (required for vertex fitting)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Set the GlobalTag — adjust for your data/MC era
from Configuration.AlCa.GlobalTag import GlobalTag
# For Run3 MC:
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')
# For Run2 UL MC:
# process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')

# TFileService for output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

# ============================================================================
# MiniAOD Track Producer (extracts reco::Track from packed candidates)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.miniAODTrackProducer_cfi")

process.miniAODTrackProducer.applyCuts = cms.bool(options.applyCuts)
process.miniAODTrackProducer.minPt = cms.double(options.minPt)
process.miniAODTrackProducer.minAbsSip2D = cms.double(options.minAbsSip2D)
process.miniAODTrackProducer.maxNormalizedChi2 = cms.double(options.maxNormalizedChi2)

if options.applyCuts:
    print(f"Track cuts enabled: pT > {options.minPt} GeV, |sip2D| >= {options.minAbsSip2D}, chi2/ndof < {options.maxNormalizedChi2}")

# ============================================================================
# MiniAOD Muon Enhanced Tracks Producer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.miniAODMuonEnhancedTracksProducer_cfi")

# ============================================================================
# Determine which modes to run
# ============================================================================
_run_leptonic = options.processMode in ('leptonic', 'both')
_run_hadronic = options.processMode in ('hadronic', 'both')
if not _run_leptonic and not _run_hadronic:
    raise ValueError(f"Unknown processMode '{options.processMode}'. "
                     "Valid options: leptonic, hadronic, both")

# ============================================================================
# Leptonic diagnostic producer and analyzer
# ============================================================================
if _run_leptonic:
    process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVsDiagnosticProducer_cfi")
    process.hyddraSVsDiag.pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices")

    process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVsDiagnosticAnalyzer_cfi")
    process.hyddraSVsDiagAnalyzer.hasGenInfo        = cms.bool(options.hasGenInfo)
    process.hyddraSVsDiagAnalyzer.isFullAOD         = cms.bool(False)
    process.hyddraSVsDiagAnalyzer.pvCollection      = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.hyddraSVsDiagAnalyzer.genParticles      = cms.InputTag("prunedGenParticles")
    process.hyddraSVsDiagAnalyzer.packedGenParticles= cms.InputTag("packedGenParticles")

    from KUCMSNtupleizer.HyddraSVProducer.miniAODTrackCollections_cfi import (
        configureMiniAODDiagnosticTrackCollection
    )
    configureMiniAODDiagnosticTrackCollection(process, options.trackCollection)

# ============================================================================
# Hadronic diagnostic producer and analyzer
# ============================================================================
if _run_hadronic:
    process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVsHadronicDiagnosticProducer_cfi")
    process.hyddraSVsHadronicDiag.pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices")

    process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVsHadronicDiagnosticAnalyzer_cfi")
    process.hyddraSVsHadronicDiagAnalyzer.hasGenInfo        = cms.bool(options.hasGenInfo)
    process.hyddraSVsHadronicDiagAnalyzer.isFullAOD         = cms.bool(False)
    process.hyddraSVsHadronicDiagAnalyzer.pvCollection      = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.hyddraSVsHadronicDiagAnalyzer.genParticles      = cms.InputTag("prunedGenParticles")
    process.hyddraSVsHadronicDiagAnalyzer.packedGenParticles= cms.InputTag("packedGenParticles")

    from KUCMSNtupleizer.HyddraSVProducer.miniAODTrackCollections_cfi import (
        configureMiniAODHadronicDiagnosticTrackCollection
    )
    configureMiniAODHadronicDiagnosticTrackCollection(process, options.trackCollection)

# ============================================================================
# Build combined diagnostic sequence
# ============================================================================
_lep_seq = cms.Sequence(
    process.hyddraSVsDiag + process.hyddraSVsDiagAnalyzer
) if _run_leptonic else None

_had_seq = cms.Sequence(
    process.hyddraSVsHadronicDiag + process.hyddraSVsHadronicDiagAnalyzer
) if _run_hadronic else None

def _diag_seq():
    seqs = [s for s in [_lep_seq, _had_seq] if s is not None]
    result = seqs[0]
    for s in seqs[1:]:
        result = result + s
    return result

# ============================================================================
# Path: MiniAOD track producers -> diagnostic producer(s) -> diagnostic analyzer(s)
# ============================================================================
process.p = cms.Path(
    process.miniAODTrackProducer +      # Produces track collections from packed candidates
    process.miniAODMuonEnhancedTracks + # Produces sip2DMuonEnhancedTracks and variants
    _diag_seq()
)

process.schedule = cms.Schedule(process.p)

# ============================================================================
# Usage examples:
# ============================================================================
# Default (leptonic, sip2DMuonEnhanced, gen-matched signal):
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root
#
# Hadronic mode only:
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root processMode=hadronic
#
# Both leptonic and hadronic:
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root processMode=both
#
# From a file list:
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFileList=myfiles.txt
#
# Custom output file:
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root outputFile=diag_output.root
#
# Different track collection:
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root trackCollection=merged
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root trackCollection=sip2DMuonEnhancedWithEle
#
# With track cuts enabled:
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root applyCuts=True
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root applyCuts=True minPt=2.0 minAbsSip2D=3.0
#
# Data (no gen info — disables genFunnel tree, stageCounts still filled):
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:data_miniAOD.root hasGenInfo=False
#
# Limit events for a quick check:
#   cmsRun testHyddraSVsDiagnostic_miniAOD_cfg.py inputFiles=file:miniAOD.root maxEvents=1000
#
# Then make diagnostic plots:
#   python3 scripts/hyddraDiagPlots/run.py --signal hyddra_diagnostic_miniAOD.root
#   python3 scripts/hyddraDiagPlots/run.py --signal signal_diag.root --background bkg_diag.root
