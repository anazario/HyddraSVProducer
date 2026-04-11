import FWCore.ParameterSet.Config as cms

miniAODElectronTrackProducer = cms.EDProducer("MiniAODElectronTrackProducer",
    # Input collections (MiniAOD defaults)
    gedElectrons   = cms.InputTag("slimmedElectrons"),
    lowPtElectrons = cms.InputTag("slimmedLowPtElectrons"),

    # Kinematic baseline (iDM analysis: pT > 1 GeV, |eta| < 2.4)
    minPt     = cms.double(1.0),
    maxAbsEta = cms.double(2.4),

    # HEM veto: off by default (2018 data artifact; not applicable to MC efficiency studies)
    applyHEMVeto = cms.bool(False),

    # GED electron MVA ID working point.
    # For Run2 UL MiniAOD (after EGammaPostRecoTools): "mvaEleID-Fall17-noIso-V2-wp90"
    # Set to "" to disable the ID cut while getting started, or for Run3.
    gedMVALabel = cms.string("mvaEleID-Fall17-noIso-V2-wp90"),

    # Low-pT electron MVA score key (stored via electronID(), not userFloat()).
    # Available keys on Run2 UL MiniAOD: "ID", "unbiased", "ptbiased". Set to "" to disable.
    lowPtMVALabel     = cms.string("ID"),
    lowPtMVAThreshold = cms.double(-0.25),

    # Cross-cleaning: reject low-pT electrons within deltaR < 0.05 of a GED electron
    crossCleanDR = cms.double(0.05),
)
