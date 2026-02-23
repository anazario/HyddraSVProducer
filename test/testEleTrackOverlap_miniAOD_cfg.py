import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('inputFileList', '',
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Path to text file with input files (one per line)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'eleTrackOverlap.root')
options.parseArguments()

import os
def getInputFiles():
    files = []
    if options.inputFileList and os.path.exists(options.inputFileList):
        with open(options.inputFileList) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    files.append(line)
    if options.inputFiles:
        files.extend(options.inputFiles)
    if not files:
        files = ['file:input_miniAOD.root']
    return files

process = cms.Process("ELEOVERLAPCHECK")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(getInputFiles())
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.eleTrackOverlapAnalyzer = cms.EDAnalyzer("EleTrackOverlapAnalyzer",
    pfCandidates  = cms.InputTag("packedPFCandidates"),
    lostTracks    = cms.InputTag("lostTracks"),
    eleLostTracks = cms.InputTag("lostTracks", "eleTracks"),
    drNoMatchThreshold = cms.double(0.4),
)

process.p = cms.Path(process.eleTrackOverlapAnalyzer)

# Usage:
#   cmsRun testEleTrackOverlap_miniAOD_cfg.py inputFiles=file:myMiniAOD.root
#   cmsRun testEleTrackOverlap_miniAOD_cfg.py inputFileList=myfiles.txt maxEvents=1000
