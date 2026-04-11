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
                 'mergedElectronTracks',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection to use. Electron options: gedElectronTracks, "
                 "lowPtElectronTracks, mergedElectronTracks (default). "
                 "Muon/packed options: sip2DMuonEnhanced, merged, mergedAll, pf, "
                 "lost, eleLost, promptMuonExtracted, displacedMuonExtracted.")
options.register('inputFileList',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Path to text file with input files (one per line)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraEXO_miniAOD_ntuple.root')
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
process.load("KUCMSNtupleizer.HyddraSVProducer.miniAODTrackProducer_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.miniAODElectronTrackProducer_cfi")

# ── EXO producer + analyzer ───────────────────────────────────────────────────
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraEXO_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraEXOAnalyzer_cfi")

# Override full-AOD defaults to MiniAOD collections
process.hyddraEXO.pvCollection          = cms.InputTag("offlineSlimmedPrimaryVertices")
process.hyddraEXOAnalyzer.pvCollection  = cms.InputTag("offlineSlimmedPrimaryVertices")
process.hyddraEXOAnalyzer.genParticles  = cms.InputTag("prunedGenParticles")
process.hyddraEXOAnalyzer.hasGenInfo    = cms.bool(options.hasGenInfo)
process.hyddraEXOAnalyzer.motherPdgId   = cms.int32(options.motherPdgId)

# ── Configure track collection ────────────────────────────────────────────────
ELECTRON_TRACK_COLLECTIONS = {
    'gedElectronTracks':    cms.InputTag("miniAODElectronTrackProducer", "gedElectronTracks"),
    'lowPtElectronTracks':  cms.InputTag("miniAODElectronTrackProducer", "lowPtElectronTracks"),
    'mergedElectronTracks': cms.InputTag("miniAODElectronTrackProducer", "mergedElectronTracks"),
}

MINIAOD_MUON_TRACK_COLLECTIONS = {
    'sip2DMuonEnhanced':       cms.InputTag("miniAODMuonEnhancedTracks", "sip2DMuonEnhancedTracks"),
    'merged':                  cms.InputTag("miniAODTrackProducer", "merged"),
    'mergedAll':               cms.InputTag("miniAODTrackProducer", "mergedAll"),
    'pf':                      cms.InputTag("miniAODTrackProducer", "pfCandidateTracks"),
    'lost':                    cms.InputTag("miniAODTrackProducer", "lostTracks"),
    'eleLost':                 cms.InputTag("miniAODTrackProducer", "eleLostTracks"),
    'promptMuonExtracted':     cms.InputTag("miniAODTrackProducer", "muonGlobalTracks"),
    'displacedMuonExtracted':  cms.InputTag("miniAODTrackProducer", "displacedMuonGlobalTracks"),
}

ALL_COLLECTIONS = {**ELECTRON_TRACK_COLLECTIONS, **MINIAOD_MUON_TRACK_COLLECTIONS}

if options.trackCollection not in ALL_COLLECTIONS:
    raise ValueError(f"Unknown trackCollection '{options.trackCollection}'. "
                     f"Available: {list(ALL_COLLECTIONS.keys())}")

tag = ALL_COLLECTIONS[options.trackCollection]
process.hyddraEXO.tracks         = tag
process.hyddraEXOAnalyzer.tracks = tag

# Enable per-SV GED/LowPt categorization for merged-collection runs
if options.trackCollection == 'mergedElectronTracks':
    process.hyddraEXOAnalyzer.gedTracks = cms.InputTag(
        "miniAODElectronTrackProducer", "gedElectronTracks")

# ── Build path ────────────────────────────────────────────────────────────────
# Only load the producers actually needed for the chosen track collection.
# miniAODMuonEnhancedTracks requires displacedGlobalMuons/displacedTracks which
# are full-AOD collections not present in MiniAOD, so it is excluded here.

NEEDS_PACKED_TRACKS = set(MINIAOD_MUON_TRACK_COLLECTIONS.keys())  # all need miniAODTrackProducer

if options.trackCollection in ELECTRON_TRACK_COLLECTIONS:
    process.p = cms.Path(
        process.miniAODElectronTrackProducer +
        process.hyddraEXO +
        process.hyddraEXOAnalyzer
    )
else:
    process.p = cms.Path(
        process.miniAODTrackProducer +
        process.miniAODElectronTrackProducer +
        process.hyddraEXO +
        process.hyddraEXOAnalyzer
    )

process.schedule = cms.Schedule(process.p)

# ── Notes ─────────────────────────────────────────────────────────────────────
# GED MVA ID and low-pT ID re-running:
#   By default the electron producer uses stored IDs. If your MiniAOD samples do
#   not have the correct IDs embedded (e.g. Run2 UL CMSSW 10_6_X), you can:
#     1. Disable the ID cuts while testing:
#          process.miniAODElectronTrackProducer.gedMVALabel   = cms.string("")
#          process.miniAODElectronTrackProducer.lowPtMVALabel = cms.string("")
#     2. For Run2 UL with correct IDs, add EGammaPostRecoTools before the path:
#          from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#          setupEgammaPostRecoSeq(process, runVID=True, runEnergyCorrections=False,
#                                 era='2018-UL')  # adjust era as needed
#          process.p.insert(0, process.egammaPostRecoSeq)
#
# Usage:
#   cmsRun testHyddraSVsEXOProducer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root
#   cmsRun testHyddraSVsEXOProducer_miniAOD_cfg.py inputFiles=file:myMiniAOD.root trackCollection=gedElectronTracks
#   cmsRun testHyddraSVsEXOProducer_miniAOD_cfg.py inputFileList=myfiles.txt hasGenInfo=True motherPdgId=9000006
