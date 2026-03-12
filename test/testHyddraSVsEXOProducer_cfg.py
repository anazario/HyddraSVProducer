import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

options = VarParsing.VarParsing('analysis')
options.register('hasGenInfo',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Process gen-level information (set False for data)")
options.register('motherPdgId',
                 54,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "PDG ID of the signal mother particle for gen matching")
options.register('trackCollection',
                 'sip2DMuonEnhanced',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection to use (same options as testHyddraSVAnalyzer_cfg.py)")
options.register('inputFileList',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to text file with input files (one per line)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraEXO_ntuple.root')
options.parseArguments()

def getInputFiles():
    files = []
    if options.inputFileList:
        if os.path.exists(options.inputFileList):
            with open(options.inputFileList) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        files.append(line)
        else:
            raise FileNotFoundError(f"Input file list not found: {options.inputFileList}")
    if options.inputFiles:
        files.extend(options.inputFiles)
    if not files:
        files = ['file:input.root']
    return files

process = cms.Process("HYDDRAEXO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(getInputFiles()))

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos:
    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile))

# ── Track producers ───────────────────────────────────────────────────────────
process.load("KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.muonGlobalTrackProducer_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.filteredTrackProducer_cfi")

# ── EXO producer + analyzer ───────────────────────────────────────────────────
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraEXO_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraEXOAnalyzer_cfi")

process.hyddraEXOAnalyzer.hasGenInfo  = cms.bool(options.hasGenInfo)
process.hyddraEXOAnalyzer.motherPdgId = cms.int32(options.motherPdgId)

# ── Configure track collection ────────────────────────────────────────────────
def configureTrackCollection(trackCol):
    """Point both producer and analyzer at the chosen track collection."""
    tags = {
        'general':                 cms.InputTag("generalTracks"),
        'generalFiltered':         cms.InputTag("filteredTrackProducer"),
        'sip2D':                   cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
        'sip2DMuonEnhanced':       cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
        'muonEnhanced':            cms.InputTag("muonEnhancedTracks", "muonEnhancedTracks"),
        'displacedGlobalMuon':     cms.InputTag("displacedGlobalMuons"),
        'displacedStandAloneMuon': cms.InputTag("displacedStandAloneMuons"),
    }
    if trackCol not in tags:
        raise ValueError(f"Unknown trackCollection '{trackCol}'")
    tag = tags[trackCol]
    process.hyddraEXO.tracks            = tag
    process.hyddraEXOAnalyzer.tracks    = tag

configureTrackCollection(options.trackCollection)

# ── Build path ────────────────────────────────────────────────────────────────
MUON_COLLECTIONS = {
    'displacedGlobalMuon', 'displacedStandAloneMuon', 'general'
}

if options.trackCollection == 'generalFiltered':
    process.p = cms.Path(
        process.filteredTrackProducer +
        process.hyddraEXO +
        process.hyddraEXOAnalyzer
    )
elif options.trackCollection in MUON_COLLECTIONS:
    process.p = cms.Path(
        process.hyddraEXO +
        process.hyddraEXOAnalyzer
    )
else:
    process.p = cms.Path(
        process.muonEnhancedTracks +
        process.hyddraEXO +
        process.hyddraEXOAnalyzer
    )

process.schedule = cms.Schedule(process.p)

# ── Usage ─────────────────────────────────────────────────────────────────────
# Single file:
#   cmsRun testHyddraSVsEXOProducer_cfg.py inputFiles=file:myinput.root
#
# From file list:
#   cmsRun testHyddraSVsEXOProducer_cfg.py inputFileList=myfiles.txt
#
# Without gen info (data):
#   cmsRun testHyddraSVsEXOProducer_cfg.py inputFiles=file:myinput.root hasGenInfo=False
#
# Custom signal mother PDG:
#   cmsRun testHyddraSVsEXOProducer_cfg.py inputFiles=file:myinput.root motherPdgId=9000006
