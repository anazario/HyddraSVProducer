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
options.register('hyddraPreset',
                 'default',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "HYDDRA leptonic algorithm preset: default, NonIso, TightIso")
options.register('trackCollection',
                 'promptMuonBestTrack',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection: promptMuonBestTrack (default), displacedMuonBestTrack, "
                 "promptMuonExtracted, displacedMuonExtracted, "
                 "promptMuonPriority, displacedMuonPriority, "
                 "promptMuonLLPNano, displacedMuonLLPNano, "
                 "general, generalFiltered, selected, sip2D, sip2DMuonEnhanced, muonEnhanced")
options.register('processMode',
                 'leptonic',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Processing mode: leptonic or hadronic (diagnostic is leptonic-only; "
                 "this argument is accepted for script compatibility)")
options.register('inputFileList',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to text file containing list of input files (one per line)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddra_diagnostic.root')
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
# Track producers (only the ones relevant to the diagnostic — no ecalTracks,
# no MuonEnhancedTracks for the muon-best-track path)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.muonGlobalTrackProducer_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.filteredTrackProducer_cfi")
process.load("KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi")

# ============================================================================
# Diagnostic producer — snapshots each leptonic stage
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVsDiagnosticProducer_cfi")

# Apply leptonic algorithm preset if requested
if options.hyddraPreset != 'default':
    from KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi import (
        leptonicHYDDRA_NonIso, leptonicHYDDRA_TightIso)
    _presets = {'NonIso': leptonicHYDDRA_NonIso, 'TightIso': leptonicHYDDRA_TightIso}
    if options.hyddraPreset not in _presets:
        raise ValueError(f"Unknown hyddraPreset '{options.hyddraPreset}'. "
                         f"Valid options: default, NonIso, TightIso")
    process.hyddraSVsDiag.leptonic = _presets[options.hyddraPreset]

# ============================================================================
# Diagnostic analyzer — per-gen-vertex funnel TTrees
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVsDiagnosticAnalyzer_cfi")
process.hyddraSVsDiagAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)

# ============================================================================
# Configure track collection for BOTH diagnostic producer and analyzer
# ============================================================================
from KUCMSNtupleizer.HyddraSVProducer.hyddraSVsDiagnosticAnalyzer_cfi import (
    configureDiagnosticTrackCollection
)
configureDiagnosticTrackCollection(process, options.trackCollection)

# ============================================================================
# Maps track collection names to (module name, mode, trackPriority list)
# ============================================================================
MUON_TRACK_COLLECTIONS = {
    'promptMuonExtracted':    ('muonGlobalTrackProducer', 'globalTrack', []),
    'displacedMuonExtracted': ('muonGlobalTrackProducer', 'globalTrack', []),
    'promptMuonBestTrack':    ('muonBestTrackProducer',   'bestTrack',   []),
    'displacedMuonBestTrack': ('muonBestTrackProducer',   'bestTrack',   []),
    'promptMuonPriority':     ('muonPriorityTrackProducer', 'priority',
                               ['globalTrack', 'innerTrack', 'outerTrack']),
    'displacedMuonPriority':  ('muonPriorityTrackProducer', 'priority',
                               ['globalTrack', 'innerTrack', 'outerTrack']),
    'promptMuonLLPNano':      ('muonLLPNanoTrackProducer',  'llpNano', [],
                               {'minMuonPt': cms.double(3.0)}),
    'displacedMuonLLPNano':   ('muonLLPNanoTrackProducer',  'llpNano', [],
                               {'minMuonPt': cms.double(3.0)}),
}

# ============================================================================
# Build the processing path based on track collection
# ============================================================================
if options.trackCollection in MUON_TRACK_COLLECTIONS:
    entry = MUON_TRACK_COLLECTIONS[options.trackCollection]
    moduleName, mode, priority = entry[0], entry[1], entry[2]
    extraKwargs = entry[3] if len(entry) > 3 else {}
    if not hasattr(process, moduleName):
        cloneArgs = dict(mode=cms.string(mode))
        if priority:
            cloneArgs['trackPriority'] = cms.vstring(*priority)
        cloneArgs.update(extraKwargs)
        setattr(process, moduleName,
                process.muonGlobalTrackProducer.clone(**cloneArgs))
    else:
        getattr(process, moduleName).mode = cms.string(mode)
        if priority:
            getattr(process, moduleName).trackPriority = cms.vstring(*priority)
        for k, v in extraKwargs.items():
            setattr(getattr(process, moduleName), k, v)
    process.p = cms.Path(
        getattr(process, moduleName) +   # Extracts/builds the muon track collection
        process.hyddraSVsDiag +          # Snapshots each leptonic stage
        process.hyddraSVsDiagAnalyzer    # Writes genFunnel / stageCounts / cleaningTracks
    )
elif options.trackCollection == 'generalFiltered':
    process.p = cms.Path(
        process.filteredTrackProducer +  # Applies quality cuts to general tracks
        process.hyddraSVsDiag +
        process.hyddraSVsDiagAnalyzer
    )
elif options.trackCollection in ['general', 'displacedGlobalMuon', 'displacedStandAloneMuon']:
    process.p = cms.Path(
        process.hyddraSVsDiag +
        process.hyddraSVsDiagAnalyzer
    )
else:
    # sip2D, sip2DMuonEnhanced, muonEnhanced, selected, etc.
    process.p = cms.Path(
        process.muonEnhancedTracks +     # Produces sip2DMuonEnhancedTracks (and variants)
        process.hyddraSVsDiag +
        process.hyddraSVsDiagAnalyzer
    )

process.schedule = cms.Schedule(process.p)

# ============================================================================
# Usage examples:
# ============================================================================
# Default (promptMuonBestTrack, gen-matched signal):
#   cmsRun testHyddraSVsDiagnostic_cfg.py inputFiles=file:signal.root
#
# From a file list:
#   cmsRun testHyddraSVsDiagnostic_cfg.py inputFileList=myfiles.txt
#
# Custom output file:
#   cmsRun testHyddraSVsDiagnostic_cfg.py inputFiles=file:signal.root outputFile=diag_output.root
#
# Displaced muon best tracks:
#   cmsRun testHyddraSVsDiagnostic_cfg.py inputFiles=file:signal.root trackCollection=displacedMuonBestTrack
#
# Data (no gen info — disables genFunnel tree, stageCounts still filled):
#   cmsRun testHyddraSVsDiagnostic_cfg.py inputFiles=file:data.root hasGenInfo=False
#
# Limit events for a quick check:
#   cmsRun testHyddraSVsDiagnostic_cfg.py inputFiles=file:signal.root maxEvents=1000
#
# Then make funnel plots:
#   python3 scripts/hyddraSVFunnelPlots.py --signal hyddra_diagnostic.root
#   python3 scripts/hyddraSVFunnelPlots.py --signal signal_diag.root --background bkg_diag.root
