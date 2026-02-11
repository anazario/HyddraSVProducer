import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

# ============================================================================
# FastSim AOD Configuration for HYDDRA SV Analyzer
# ============================================================================
# FastSim AOD does not contain several displaced reconstruction collections
# that are present in full sim AOD. This config remaps to available equivalents:
#
# Full Sim Collection          -> FastSim Equivalent
# -----------------------------------------------
# displacedGlobalMuons         -> globalMuons
# displacedStandAloneMuons     -> standAloneMuons
# displacedTracks              -> (not available)
# displacedMuons               -> (not available)
# particleFlowEGamma (SCs)     -> particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel
#                                  (barrel-only; endcap lost unless particleFlowEGamma exists)
# particleFlowSuperClusterOOT  -> same InputTag (exists in FastSim)
#
# Supported trackCollection options:
#   general, generalFiltered,
#   fastSimSip2DMuonEnhanced (default), fastSimMuonEnhanced, fastSimSip2D, fastSimSelected
# ============================================================================

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
                 'fastSimSip2DMuonEnhanced',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection: general, generalFiltered, "
                 "fastSimSip2DMuonEnhanced (default), fastSimMuonEnhanced, fastSimSip2D, fastSimSelected")
options.register('inputFileList',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to text file containing list of input files (one per line)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraSV_fastSim_ntuple.root')
options.parseArguments()

# Validate track collection choice for FastSim
FASTSIM_ALLOWED_COLLECTIONS = [
    'general', 'generalFiltered',
    'fastSimSip2DMuonEnhanced', 'fastSimMuonEnhanced', 'fastSimSip2D', 'fastSimSelected',
]
if options.trackCollection not in FASTSIM_ALLOWED_COLLECTIONS:
    raise ValueError(
        f"Track collection '{options.trackCollection}' is not available in FastSim AOD. "
        f"Allowed: {FASTSIM_ALLOWED_COLLECTIONS}. "
        f"Collections requiring full sim displaced reconstruction (sip2DMuonEnhanced, muonEnhanced, "
        f"displacedGlobalMuon, etc.) are not available in FastSim."
    )

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

# Load geometry and magnetic field
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos:
    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')

process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Set the GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

# TFileService for output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

# ============================================================================
# FastSim MuonEnhancedTracks Producer
# ============================================================================
# Same as MuonEnhancedTracksProducer but without displacedTracks/displacedMuons.
# Uses globalMuons and standAloneMuons as proxies.
process.load("KUCMSNtupleizer.HyddraSVProducer.fastSimMuonEnhancedTracksProducer_cfi")

# ============================================================================
# FilteredTrackProducer (applies quality cuts to general tracks)
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.filteredTrackProducer_cfi")

# ============================================================================
# ECALTracks Producer (produces displacedElectronSCs for SC matching)
# ============================================================================
process.load("KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi")

# FastSim SC remapping:
# In-time SCs: particleFlowEGamma may not exist in FastSim.
# Use particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel as fallback.
# NOTE: This is barrel-only. Endcap SCs from
# particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower
# are not included (the producer only takes one InputTag for in-time SCs).
# If particleFlowEGamma exists in your FastSim sample, comment out this override.
process.ecalTracks.superClusters = cms.InputTag(
    "particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"
)
# OOT SCs: same InputTag works in FastSim
# process.ecalTracks.ootSuperClusters is already correct (particleFlowSuperClusterOOTECAL:...Barrel)

# ============================================================================
# HYDDRA SV Producer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi")

# ============================================================================
# HYDDRA SV Analyzer
# ============================================================================
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi")
process.hyddraSVAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)
process.hyddraSVAnalyzer.mergedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs")

# FastSim muon track remapping:
# displacedGlobalMuons does not exist in FastSim -> use globalMuons
# (prompt global muon tracks; less efficient for displaced signatures)
process.hyddraSVAnalyzer.muonTracks = cms.InputTag("globalMuons")

# ============================================================================
# Configure track collection for BOTH producer and analyzer
# ============================================================================
from KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi import configureTrackCollection
configureTrackCollection(process, options.trackCollection)

# Configure processing mode
if options.processMode == 'leptonic':
    process.hyddraSVAnalyzer.hadronicVertices = cms.InputTag("")
elif options.processMode == 'hadronic':
    process.hyddraSVAnalyzer.leptonicVertices = cms.InputTag("")

# ============================================================================
# Path
# ============================================================================
# FastSim track collections that need the FastSimMuonEnhancedTracksProducer
FASTSIM_MUON_ENHANCED_COLLECTIONS = [
    'fastSimSip2DMuonEnhanced', 'fastSimMuonEnhanced', 'fastSimSip2D', 'fastSimSelected',
]

if options.trackCollection in FASTSIM_MUON_ENHANCED_COLLECTIONS:
    process.p = cms.Path(
        process.fastSimMuonEnhancedTracks +  # Produces muon-enhanced tracks from FastSim collections
        process.ecalTracks +                  # Produces displacedElectronSCs for SC matching
        process.hyddraSVs +                   # Produces leptonic/hadronic vertices
        process.hyddraSVAnalyzer              # Writes TTree output
    )
elif options.trackCollection == 'generalFiltered':
    process.p = cms.Path(
        process.filteredTrackProducer +
        process.ecalTracks +
        process.hyddraSVs +
        process.hyddraSVAnalyzer
    )
elif options.trackCollection == 'general':
    process.p = cms.Path(
        process.ecalTracks +
        process.hyddraSVs +
        process.hyddraSVAnalyzer
    )

process.schedule = cms.Schedule(process.p)

# ============================================================================
# Usage examples:
# ============================================================================
# Default (fastSimSip2DMuonEnhanced tracks — muon-enhanced with sip2D selection):
#   cmsRun testHyddraSVAnalyzer_fastSim_cfg.py inputFiles=file:fastsim_aod.root
#
# FastSim muon-enhanced tracks (before sip2D cut):
#   cmsRun testHyddraSVAnalyzer_fastSim_cfg.py inputFiles=file:fastsim_aod.root trackCollection=fastSimMuonEnhanced
#
# General tracks with filtering:
#   cmsRun testHyddraSVAnalyzer_fastSim_cfg.py inputFiles=file:fastsim_aod.root trackCollection=generalFiltered
#
# General tracks (no filtering):
#   cmsRun testHyddraSVAnalyzer_fastSim_cfg.py inputFiles=file:fastsim_aod.root trackCollection=general
#
# Leptonic only:
#   cmsRun testHyddraSVAnalyzer_fastSim_cfg.py inputFiles=file:fastsim_aod.root processMode=leptonic
#
# From a file list:
#   cmsRun testHyddraSVAnalyzer_fastSim_cfg.py inputFileList=fastsim_files.txt
#
# If particleFlowEGamma exists in your FastSim sample, comment out the
# ecalTracks.superClusters override above to use the full (barrel+endcap) collection.
