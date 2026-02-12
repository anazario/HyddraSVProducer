// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      TrackAnalyzer
//
// Description: Standalone EDAnalyzer for studying tracks with gen-matching
//              for single-track signal efficiency studies.
//
// Original Author:  Andres Abreu
//

#include <memory>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// Helper classes
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenTools.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"

// ROOT
#include "TTree.h"

class TrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit TrackAnalyzer(const edm::ParameterSet&);
  ~TrackAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  void clearBranches();
  reco::GenParticleCollection getStableChargedDaughtersFromPacked(
      const GenVertex& genVertex,
      const std::vector<pat::PackedGenParticle>& packedGenParticles) const;

  // Configuration
  bool hasGenInfo_;
  bool doChargedHadronMatching_;
  double genMatchDeltaRCut_;
  bool isFullAOD_;

  // Tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> packedGenToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;

  // Handles
  edm::Handle<reco::TrackCollection> tracksHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;
  edm::Handle<reco::GenParticleCollection> genHandle_;
  edm::Handle<std::vector<pat::PackedGenParticle>> packedGenHandle_;

  // TTree
  TTree* tree_;

  // Event info
  unsigned int run_;
  unsigned int lumi_;
  unsigned long long event_;

  // Track branches
  unsigned int nTracks_;
  std::vector<float> pt_;
  std::vector<float> eta_;
  std::vector<float> phi_;
  std::vector<float> p_;
  std::vector<float> px_;
  std::vector<float> py_;
  std::vector<float> pz_;
  std::vector<int> charge_;
  std::vector<float> vx_;
  std::vector<float> vy_;
  std::vector<float> vz_;
  std::vector<float> dxy_;
  std::vector<float> dxyError_;
  std::vector<float> dxySignificance_;
  std::vector<float> dz_;
  std::vector<float> dzError_;
  std::vector<float> ip2D_;
  std::vector<float> sip2D_;
  std::vector<float> ip3D_;
  std::vector<float> sip3D_;
  std::vector<float> chi2_;
  std::vector<float> ndof_;
  std::vector<float> normalizedChi2_;
  std::vector<float> ptError_;
  std::vector<float> etaError_;
  std::vector<float> phiError_;
  std::vector<float> qoverp_;
  std::vector<float> qoverpError_;
  std::vector<int> qualityMask_;
  std::vector<float> validHitFraction_;
  std::vector<int> nValidHits_;
  std::vector<int> nLostHits_;
  std::vector<int> nValidPixelHits_;
  std::vector<int> nValidStripHits_;
  std::vector<int> nValidPixelBarrelHits_;
  std::vector<int> nValidPixelEndcapHits_;
  std::vector<int> missingInnerHits_;
  std::vector<int> missingOuterHits_;
  std::vector<int> trackerLayersWithMeasurement_;
  std::vector<int> pixelLayersWithMeasurement_;
  std::vector<bool> isHighPurity_;

  // Gen-matching branches (per track)
  std::vector<bool> isSignal_;
  std::vector<bool> isSignalElectron_;
  std::vector<bool> isSignalMuon_;
  std::vector<float> genMatchDeltaR_;
  std::vector<float> genMatchPt_;
  std::vector<float> genMatchEta_;
  std::vector<float> genMatchPhi_;
  std::vector<float> genMatchP_;
  std::vector<int> genMatchPdgId_;
  std::vector<int> genMatchMotherPdgId_;
  std::vector<float> genMatchVx_;
  std::vector<float> genMatchVy_;
  std::vector<float> genMatchVz_;
  std::vector<float> genMatchDxy_;
  std::vector<float> genMatchRelPtDiff_;

  // Charged hadron matching (gated by doChargedHadronMatching_)
  GenVertices genVertices_;
  std::map<GenVertex, GenMatches> chargedMatches_;
  reco::TrackCollection signalTracks_;

  std::vector<bool> isSignalHadron_;

  std::vector<int> chargedHadron_genVertexIndex_;
  std::vector<float> chargedHadron_genDxy_;
  std::vector<float> chargedHadron_p_;
  std::vector<float> chargedHadron_pt_;
  std::vector<float> chargedHadron_eta_;
  std::vector<float> chargedHadron_phi_;
  std::vector<int> chargedHadron_charge_;
  std::vector<float> chargedHadron_deltaR_;
  std::vector<int> chargedHadron_trackIndex_;
  std::vector<float> chargedHadron_trackP_;
  std::vector<float> chargedHadron_trackPt_;
  std::vector<float> chargedHadron_trackEta_;
  std::vector<float> chargedHadron_trackPhi_;
  std::vector<int> chargedHadron_trackCharge_;
  std::vector<float> chargedHadron_trackNormChi2_;

  // Gen-level summary (all signal particles)
  unsigned int nSignalGen_;
  unsigned int nSignalElectronGen_;
  unsigned int nSignalMuonGen_;
  std::vector<float> signalGen_pt_;
  std::vector<float> signalGen_eta_;
  std::vector<float> signalGen_phi_;
  std::vector<float> signalGen_p_;
  std::vector<int> signalGen_pdgId_;
  std::vector<int> signalGen_motherPdgId_;
  std::vector<float> signalGen_vx_;
  std::vector<float> signalGen_vy_;
  std::vector<float> signalGen_vz_;
  std::vector<float> signalGen_dxy_;
  std::vector<bool> signalGen_isMatched_;
  std::vector<float> signalGen_matchedTrackPt_;
  std::vector<float> signalGen_matchedDeltaR_;
};

TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig) :
  hasGenInfo_(iConfig.getParameter<bool>("hasGenInfo")),
  doChargedHadronMatching_(iConfig.getParameter<bool>("doChargedHadronMatching")),
  genMatchDeltaRCut_(iConfig.getParameter<double>("genMatchDeltaRCut")),
  isFullAOD_(iConfig.getParameter<bool>("isFullAOD")),
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  transientTrackBuilder_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder")))
{
  usesResource("TFileService");

  if(hasGenInfo_) {
    genToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
    if(!isFullAOD_) {
      packedGenToken_ = consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("packedGenParticles"));
    }
  }
}

void TrackAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "Track Tree");

  // Event info
  tree_->Branch("run", &run_);
  tree_->Branch("lumi", &lumi_);
  tree_->Branch("event", &event_);

  // Track branches
  tree_->Branch("Track_nTotal", &nTracks_);
  tree_->Branch("Track_pt", &pt_);
  tree_->Branch("Track_eta", &eta_);
  tree_->Branch("Track_phi", &phi_);
  tree_->Branch("Track_p", &p_);
  tree_->Branch("Track_px", &px_);
  tree_->Branch("Track_py", &py_);
  tree_->Branch("Track_pz", &pz_);
  tree_->Branch("Track_charge", &charge_);
  tree_->Branch("Track_vx", &vx_);
  tree_->Branch("Track_vy", &vy_);
  tree_->Branch("Track_vz", &vz_);
  tree_->Branch("Track_dxy", &dxy_);
  tree_->Branch("Track_dxyError", &dxyError_);
  tree_->Branch("Track_dxySignificance", &dxySignificance_);
  tree_->Branch("Track_dz", &dz_);
  tree_->Branch("Track_dzError", &dzError_);
  tree_->Branch("Track_ip2D", &ip2D_);
  tree_->Branch("Track_sip2D", &sip2D_);
  tree_->Branch("Track_ip3D", &ip3D_);
  tree_->Branch("Track_sip3D", &sip3D_);
  tree_->Branch("Track_chi2", &chi2_);
  tree_->Branch("Track_ndof", &ndof_);
  tree_->Branch("Track_normalizedChi2", &normalizedChi2_);
  tree_->Branch("Track_ptError", &ptError_);
  tree_->Branch("Track_etaError", &etaError_);
  tree_->Branch("Track_phiError", &phiError_);
  tree_->Branch("Track_qoverp", &qoverp_);
  tree_->Branch("Track_qoverpError", &qoverpError_);
  tree_->Branch("Track_qualityMask", &qualityMask_);
  tree_->Branch("Track_validHitFraction", &validHitFraction_);
  tree_->Branch("Track_nValidHits", &nValidHits_);
  tree_->Branch("Track_nLostHits", &nLostHits_);
  tree_->Branch("Track_nValidPixelHits", &nValidPixelHits_);
  tree_->Branch("Track_nValidStripHits", &nValidStripHits_);
  tree_->Branch("Track_nValidPixelBarrelHits", &nValidPixelBarrelHits_);
  tree_->Branch("Track_nValidPixelEndcapHits", &nValidPixelEndcapHits_);
  tree_->Branch("Track_missingInnerHits", &missingInnerHits_);
  tree_->Branch("Track_missingOuterHits", &missingOuterHits_);
  tree_->Branch("Track_trackerLayersWithMeasurement", &trackerLayersWithMeasurement_);
  tree_->Branch("Track_pixelLayersWithMeasurement", &pixelLayersWithMeasurement_);
  tree_->Branch("Track_isHighPurity", &isHighPurity_);

  // Gen-matching branches
  if(hasGenInfo_) {
    tree_->Branch("Track_isSignal", &isSignal_);
    tree_->Branch("Track_isSignalElectron", &isSignalElectron_);
    tree_->Branch("Track_isSignalMuon", &isSignalMuon_);
    if(doChargedHadronMatching_) {
      tree_->Branch("Track_isSignalHadron", &isSignalHadron_);
    }
    tree_->Branch("Track_genMatchDeltaR", &genMatchDeltaR_);
    tree_->Branch("Track_genMatchPt", &genMatchPt_);
    tree_->Branch("Track_genMatchEta", &genMatchEta_);
    tree_->Branch("Track_genMatchPhi", &genMatchPhi_);
    tree_->Branch("Track_genMatchP", &genMatchP_);
    tree_->Branch("Track_genMatchPdgId", &genMatchPdgId_);
    tree_->Branch("Track_genMatchMotherPdgId", &genMatchMotherPdgId_);
    tree_->Branch("Track_genMatchVx", &genMatchVx_);
    tree_->Branch("Track_genMatchVy", &genMatchVy_);
    tree_->Branch("Track_genMatchVz", &genMatchVz_);
    tree_->Branch("Track_genMatchDxy", &genMatchDxy_);
    tree_->Branch("Track_genMatchRelPtDiff", &genMatchRelPtDiff_);

    // Gen-level summary
    tree_->Branch("SignalGen_nTotal", &nSignalGen_);
    tree_->Branch("SignalGen_nElectron", &nSignalElectronGen_);
    tree_->Branch("SignalGen_nMuon", &nSignalMuonGen_);
    tree_->Branch("SignalGen_pt", &signalGen_pt_);
    tree_->Branch("SignalGen_eta", &signalGen_eta_);
    tree_->Branch("SignalGen_phi", &signalGen_phi_);
    tree_->Branch("SignalGen_p", &signalGen_p_);
    tree_->Branch("SignalGen_pdgId", &signalGen_pdgId_);
    tree_->Branch("SignalGen_motherPdgId", &signalGen_motherPdgId_);
    tree_->Branch("SignalGen_vx", &signalGen_vx_);
    tree_->Branch("SignalGen_vy", &signalGen_vy_);
    tree_->Branch("SignalGen_vz", &signalGen_vz_);
    tree_->Branch("SignalGen_dxy", &signalGen_dxy_);
    tree_->Branch("SignalGen_isMatched", &signalGen_isMatched_);
    tree_->Branch("SignalGen_matchedTrackPt", &signalGen_matchedTrackPt_);
    tree_->Branch("SignalGen_matchedDeltaR", &signalGen_matchedDeltaR_);

    if(doChargedHadronMatching_) {
      tree_->Branch("ChargedHadron_genVertexIndex", &chargedHadron_genVertexIndex_);
      tree_->Branch("ChargedHadron_genDxy", &chargedHadron_genDxy_);
      tree_->Branch("ChargedHadron_p", &chargedHadron_p_);
      tree_->Branch("ChargedHadron_pt", &chargedHadron_pt_);
      tree_->Branch("ChargedHadron_eta", &chargedHadron_eta_);
      tree_->Branch("ChargedHadron_phi", &chargedHadron_phi_);
      tree_->Branch("ChargedHadron_charge", &chargedHadron_charge_);
      tree_->Branch("ChargedHadron_deltaR", &chargedHadron_deltaR_);
      tree_->Branch("ChargedHadron_trackIndex", &chargedHadron_trackIndex_);
      tree_->Branch("ChargedHadron_trackP", &chargedHadron_trackP_);
      tree_->Branch("ChargedHadron_trackPt", &chargedHadron_trackPt_);
      tree_->Branch("ChargedHadron_trackEta", &chargedHadron_trackEta_);
      tree_->Branch("ChargedHadron_trackPhi", &chargedHadron_trackPhi_);
      tree_->Branch("ChargedHadron_trackCharge", &chargedHadron_trackCharge_);
      tree_->Branch("ChargedHadron_trackNormChi2", &chargedHadron_trackNormChi2_);
    }
  }
}

void TrackAnalyzer::clearBranches() {
  run_ = 0;
  lumi_ = 0;
  event_ = 0;

  nTracks_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  p_.clear();
  px_.clear();
  py_.clear();
  pz_.clear();
  charge_.clear();
  vx_.clear();
  vy_.clear();
  vz_.clear();
  dxy_.clear();
  dxyError_.clear();
  dxySignificance_.clear();
  dz_.clear();
  dzError_.clear();
  ip2D_.clear();
  sip2D_.clear();
  ip3D_.clear();
  sip3D_.clear();
  chi2_.clear();
  ndof_.clear();
  normalizedChi2_.clear();
  ptError_.clear();
  etaError_.clear();
  phiError_.clear();
  qoverp_.clear();
  qoverpError_.clear();
  qualityMask_.clear();
  validHitFraction_.clear();
  nValidHits_.clear();
  nLostHits_.clear();
  nValidPixelHits_.clear();
  nValidStripHits_.clear();
  nValidPixelBarrelHits_.clear();
  nValidPixelEndcapHits_.clear();
  missingInnerHits_.clear();
  missingOuterHits_.clear();
  trackerLayersWithMeasurement_.clear();
  pixelLayersWithMeasurement_.clear();
  isHighPurity_.clear();

  if(hasGenInfo_) {
    isSignal_.clear();
    isSignalElectron_.clear();
    isSignalMuon_.clear();
    if(doChargedHadronMatching_) {
      isSignalHadron_.clear();
    }
    genMatchDeltaR_.clear();
    genMatchPt_.clear();
    genMatchEta_.clear();
    genMatchPhi_.clear();
    genMatchP_.clear();
    genMatchPdgId_.clear();
    genMatchMotherPdgId_.clear();
    genMatchVx_.clear();
    genMatchVy_.clear();
    genMatchVz_.clear();
    genMatchDxy_.clear();
    genMatchRelPtDiff_.clear();

    nSignalGen_ = 0;
    nSignalElectronGen_ = 0;
    nSignalMuonGen_ = 0;
    signalGen_pt_.clear();
    signalGen_eta_.clear();
    signalGen_phi_.clear();
    signalGen_p_.clear();
    signalGen_pdgId_.clear();
    signalGen_motherPdgId_.clear();
    signalGen_vx_.clear();
    signalGen_vy_.clear();
    signalGen_vz_.clear();
    signalGen_dxy_.clear();
    signalGen_isMatched_.clear();
    signalGen_matchedTrackPt_.clear();
    signalGen_matchedDeltaR_.clear();

    if(doChargedHadronMatching_) {
      genVertices_.clear();
      chargedMatches_.clear();
      signalTracks_.clear();
      chargedHadron_genVertexIndex_.clear();
      chargedHadron_genDxy_.clear();
      chargedHadron_p_.clear();
      chargedHadron_pt_.clear();
      chargedHadron_eta_.clear();
      chargedHadron_phi_.clear();
      chargedHadron_charge_.clear();
      chargedHadron_deltaR_.clear();
      chargedHadron_trackIndex_.clear();
      chargedHadron_trackP_.clear();
      chargedHadron_trackPt_.clear();
      chargedHadron_trackEta_.clear();
      chargedHadron_trackPhi_.clear();
      chargedHadron_trackCharge_.clear();
      chargedHadron_trackNormChi2_.clear();
    }
  }
}

void TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clearBranches();

  // Event info
  run_ = iEvent.id().run();
  lumi_ = iEvent.luminosityBlock();
  event_ = iEvent.id().event();

  // Get handles
  iEvent.getByToken(tracksToken_, tracksHandle_);
  iEvent.getByToken(pvToken_, pvHandle_);

  if(!tracksHandle_.isValid() || tracksHandle_->empty()) {
    tree_->Fill();
    return;
  }

  const reco::Vertex& pv = pvHandle_->at(0);
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);

  // Gen-matching setup
  typedef PairedObjectCollection<reco::TransientTrack, reco::GenParticle> TTGenMatches;
  TTGenMatches genMatches;
  reco::GenParticleCollection signalGenParticles;

  if(hasGenInfo_) {
    iEvent.getByToken(genToken_, genHandle_);
    if(!isFullAOD_)
      iEvent.getByToken(packedGenToken_, packedGenHandle_);

    // Build transient tracks for gen-matching
    std::vector<reco::TransientTrack> ttracks;
    for(const auto& track : *tracksHandle_) {
      ttracks.emplace_back(ttBuilder->build(track));
    }

    if(doChargedHadronMatching_) {
      // Full GenVertices procedure (following HadronSVs.hh)
      GenVertices allSignalSVs(*genHandle_);
      DeltaRGenMatchHungarian<reco::TransientTrack> assigner(ttracks, allSignalSVs.getAllGenParticles());
      genMatches = assigner.GetPairedObjects();

      genVertices_ = GenVertices(genMatches.ConvertFromTTracks(), genMatchDeltaRCut_);
      allSignalSVs += genVertices_;
      genVertices_ = allSignalSVs;

      // Per-vertex charged hadron matching and signal track collection
      for(const auto& genVertex : genVertices_) {
        reco::GenParticleCollection stableChargedDaughters = isFullAOD_ ?
            genVertex.getStableChargedDaughters(*genHandle_) :
            getStableChargedDaughtersFromPacked(genVertex, *packedGenHandle_);
        DeltaRGenMatchHungarian<reco::Track> chargedParticleAssigner(*tracksHandle_, stableChargedDaughters);
        chargedMatches_[genVertex] = chargedParticleAssigner.GetPairedObjects();

        if(!genVertex.hasTracks()) continue;
        for(const auto& pair : genVertex.genMatches())
          signalTracks_.emplace_back(pair.GetObjectA());
      }

      // Build signalGenParticles from leptonic gen vertices for the SignalGen summary
      for(const auto& gen : *genHandle_) {
        if(gen.status() == 1 && (isSignalGenMuon(gen) || isSignalGenElectron(gen))) {
          signalGenParticles.emplace_back(gen);
        }
      }
    } else {
      // Simple lepton-only matching (original behavior)
      for(const auto& gen : *genHandle_) {
        if(gen.status() == 1 && (isSignalGenMuon(gen) || isSignalGenElectron(gen))) {
          signalGenParticles.emplace_back(gen);
        }
      }

      DeltaRGenMatchHungarian<reco::TransientTrack> assigner(ttracks, signalGenParticles);
      genMatches = assigner.GetPairedObjects();
    }
  }

  nTracks_ = tracksHandle_->size();

  // Loop over tracks
  for(size_t i = 0; i < tracksHandle_->size(); ++i) {
    const reco::Track& track = tracksHandle_->at(i);
    const reco::TransientTrack ttrack = ttBuilder->build(track);

    // Calculate IP parameters
    GlobalVector direction(track.px(), track.py(), track.pz());
    std::pair<bool, Measurement1D> ip3DResult = IPTools::signedImpactParameter3D(ttrack, direction, pv);
    std::pair<bool, Measurement1D> ip2DResult = IPTools::signedTransverseImpactParameter(ttrack, direction, pv);

    // Fill track branches
    pt_.push_back(float(track.pt()));
    eta_.push_back(float(track.eta()));
    phi_.push_back(float(track.phi()));
    p_.push_back(float(track.p()));
    px_.push_back(float(track.px()));
    py_.push_back(float(track.py()));
    pz_.push_back(float(track.pz()));
    charge_.push_back(int(track.charge()));
    vx_.push_back(float(track.vx()));
    vy_.push_back(float(track.vy()));
    vz_.push_back(float(track.vz()));
    dxy_.push_back(float(track.dxy(pv.position())));
    dxyError_.push_back(float(track.dxyError()));
    dxySignificance_.push_back(float(fabs(track.dxy(pv.position())) / track.dxyError()));
    dz_.push_back(float(track.dz(pv.position())));
    dzError_.push_back(float(track.dzError()));
    ip2D_.push_back(float(ip2DResult.second.value()));
    sip2D_.push_back(float(ip2DResult.second.value() / ip2DResult.second.error()));
    ip3D_.push_back(float(ip3DResult.second.value()));
    sip3D_.push_back(float(ip3DResult.second.value() / ip3DResult.second.error()));
    chi2_.push_back(float(track.chi2()));
    ndof_.push_back(float(track.ndof()));
    normalizedChi2_.push_back(float(track.normalizedChi2()));
    ptError_.push_back(float(track.ptError()));
    etaError_.push_back(float(track.etaError()));
    phiError_.push_back(float(track.phiError()));
    qoverp_.push_back(float(track.qoverp()));
    qoverpError_.push_back(float(track.qoverpError()));
    qualityMask_.push_back(int(track.qualityMask()));
    validHitFraction_.push_back(float(track.validFraction()));
    nValidHits_.push_back(int(track.numberOfValidHits()));
    nLostHits_.push_back(int(track.numberOfLostHits()));
    nValidPixelHits_.push_back(int(track.hitPattern().numberOfValidPixelHits()));
    nValidStripHits_.push_back(int(track.hitPattern().numberOfValidStripHits()));
    nValidPixelBarrelHits_.push_back(int(track.hitPattern().numberOfValidPixelBarrelHits()));
    nValidPixelEndcapHits_.push_back(int(track.hitPattern().numberOfValidPixelEndcapHits()));
    missingInnerHits_.push_back(int(track.missingInnerHits()));
    missingOuterHits_.push_back(int(track.missingOuterHits()));
    trackerLayersWithMeasurement_.push_back(int(track.hitPattern().trackerLayersWithMeasurement()));
    pixelLayersWithMeasurement_.push_back(int(track.hitPattern().pixelLayersWithMeasurement()));
    isHighPurity_.push_back(track.quality(reco::TrackBase::highPurity));

    // Gen-matching for this track
    if(hasGenInfo_) {
      bool foundMatch = false;
      bool isElectron = false;
      bool isMuon = false;
      float matchDeltaR = -1.0;
      float matchPt = -1.0;
      float matchEta = -999.0;
      float matchPhi = -999.0;
      float matchP = -1.0;
      int matchPdgId = 0;
      int matchMotherPdgId = 0;
      float matchVx = -999.0;
      float matchVy = -999.0;
      float matchVz = -999.0;
      float matchDxy = -1.0;
      float matchRelPtDiff = -999.0;

      // Find the gen match for this track
      for(const auto& pair : genMatches) {
        if(TrackHelper::SameTrack(track, pair.GetObjectA().track())) {
          matchDeltaR = float(pair.GetDeltaR());
          const reco::GenParticle& gen = pair.GetObjectB();

          matchPt = float(gen.pt());
          matchEta = float(gen.eta());
          matchPhi = float(gen.phi());
          matchP = float(gen.p());
          matchPdgId = int(gen.pdgId());
          matchMotherPdgId = gen.mother() ? int(gen.mother()->pdgId()) : 0;
          matchVx = float(gen.vx());
          matchVy = float(gen.vy());
          matchVz = float(gen.vz());
          matchDxy = float(sqrt(gen.vx() * gen.vx() + gen.vy() * gen.vy()));
          matchRelPtDiff = float((gen.pt() - track.pt()) / gen.pt());

          // Apply deltaR cut to determine if it's a signal match
          if(pair.GetDeltaR() < genMatchDeltaRCut_) {
            foundMatch = true;
            isElectron = isSignalGenElectron(gen);
            isMuon = isSignalGenMuon(gen);
          }
          break;
        }
      }

      isSignal_.push_back(foundMatch);
      isSignalElectron_.push_back(isElectron);
      isSignalMuon_.push_back(isMuon);
      if(doChargedHadronMatching_) {
        bool isHadron = false;
        for(const auto& [genVertex, matches] : chargedMatches_) {
          for(const auto& chPair : matches) {
            if(TrackHelper::SameTrack(track, chPair.GetObjectA()) && chPair.GetDeltaR() < genMatchDeltaRCut_) {
              isHadron = true;
              break;
            }
          }
          if(isHadron) break;
        }
        isSignalHadron_.push_back(isHadron);
      }
      genMatchDeltaR_.push_back(matchDeltaR);
      genMatchPt_.push_back(matchPt);
      genMatchEta_.push_back(matchEta);
      genMatchPhi_.push_back(matchPhi);
      genMatchP_.push_back(matchP);
      genMatchPdgId_.push_back(matchPdgId);
      genMatchMotherPdgId_.push_back(matchMotherPdgId);
      genMatchVx_.push_back(matchVx);
      genMatchVy_.push_back(matchVy);
      genMatchVz_.push_back(matchVz);
      genMatchDxy_.push_back(matchDxy);
      genMatchRelPtDiff_.push_back(matchRelPtDiff);
    }
  }

  // Fill gen-level summary (from signal gen particle perspective)
  if(hasGenInfo_) {
    nSignalGen_ = signalGenParticles.size();
    unsigned int nElectron = 0;
    unsigned int nMuon = 0;

    for(const auto& gen : signalGenParticles) {
      bool isMatched = false;
      float matchedTrackPt = -1.0;
      float matchedDeltaR = -1.0;

      // Find if this gen particle was matched to a track
      for(const auto& pair : genMatches) {
        // Check if this gen particle matches (by comparing pointers isn't reliable, use kinematics)
        if(fabs(pair.GetObjectB().pt() - gen.pt()) < 1e-6 &&
           fabs(pair.GetObjectB().eta() - gen.eta()) < 1e-6) {
          matchedDeltaR = float(pair.GetDeltaR());
          if(pair.GetDeltaR() < genMatchDeltaRCut_) {
            isMatched = true;
            matchedTrackPt = float(pair.GetObjectA().track().pt());
          }
          break;
        }
      }

      float genDxy = sqrt(gen.vx() * gen.vx() + gen.vy() * gen.vy());

      signalGen_pt_.push_back(float(gen.pt()));
      signalGen_eta_.push_back(float(gen.eta()));
      signalGen_phi_.push_back(float(gen.phi()));
      signalGen_p_.push_back(float(gen.p()));
      signalGen_pdgId_.push_back(int(gen.pdgId()));
      signalGen_motherPdgId_.push_back(gen.mother() ? int(gen.mother()->pdgId()) : 0);
      signalGen_vx_.push_back(float(gen.vx()));
      signalGen_vy_.push_back(float(gen.vy()));
      signalGen_vz_.push_back(float(gen.vz()));
      signalGen_dxy_.push_back(float(genDxy));
      signalGen_isMatched_.push_back(isMatched);
      signalGen_matchedTrackPt_.push_back(matchedTrackPt);
      signalGen_matchedDeltaR_.push_back(matchedDeltaR);

      if(isSignalGenElectron(gen)) nElectron++;
      if(isSignalGenMuon(gen)) nMuon++;
    }

    nSignalElectronGen_ = nElectron;
    nSignalMuonGen_ = nMuon;

    // Fill ChargedHadron branches (per gen charged daughter)
    if(doChargedHadronMatching_) {
      int genVertexIndex = 0;
      for(const auto& genVertex : genVertices_) {
        if(chargedMatches_.find(genVertex) == chargedMatches_.end()) {
          genVertexIndex++;
          continue;
        }
        for(const auto& pair : chargedMatches_.at(genVertex)) {
          const reco::Track& matchedTrack = pair.GetObjectA();
          const reco::GenParticle& gen = pair.GetObjectB();

          chargedHadron_genVertexIndex_.push_back(genVertexIndex);
          chargedHadron_genDxy_.push_back(float(genVertex.dxy()));
          chargedHadron_p_.push_back(float(gen.p()));
          chargedHadron_pt_.push_back(float(gen.pt()));
          chargedHadron_eta_.push_back(float(gen.eta()));
          chargedHadron_phi_.push_back(float(gen.phi()));
          chargedHadron_charge_.push_back(int(gen.charge()));
          chargedHadron_deltaR_.push_back(float(pair.GetDeltaR()));
          chargedHadron_trackIndex_.push_back(int(TrackHelper::FindTrackIndex(matchedTrack, *tracksHandle_)));
          chargedHadron_trackP_.push_back(float(matchedTrack.p()));
          chargedHadron_trackPt_.push_back(float(matchedTrack.pt()));
          chargedHadron_trackEta_.push_back(float(matchedTrack.eta()));
          chargedHadron_trackPhi_.push_back(float(matchedTrack.phi()));
          chargedHadron_trackCharge_.push_back(int(matchedTrack.charge()));
          chargedHadron_trackNormChi2_.push_back(float(matchedTrack.normalizedChi2()));
        }
        genVertexIndex++;
      }
    }
  }

  tree_->Fill();
}

reco::GenParticleCollection TrackAnalyzer::getStableChargedDaughtersFromPacked(
    const GenVertex& genVertex,
    const std::vector<pat::PackedGenParticle>& packedGenParticles) const {

  const reco::Candidate* genZ = genVertex.genPair().first.mother();
  reco::GenParticleCollection result;

  if(!genZ) return result;

  for(const auto& packed : packedGenParticles) {
    if(packed.status() != 1 || packed.charge() == 0) continue;

    const reco::Candidate* mom = packed.mother(0);
    if(!mom) continue;

    bool isFromSameZ = false;
    const reco::Candidate* prev = nullptr;
    while(mom && mom != prev) {
      if(mom->pdgId() == 23) {
        if(mom == genZ) isFromSameZ = true;
        break;
      }
      prev = mom;
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }

    if(isFromSameZ) {
      result.emplace_back(reco::GenParticle(
          packed.charge(), packed.p4(), packed.vertex(),
          packed.pdgId(), packed.status(), true));
    }
  }
  return result;
}

void TrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<bool>("hasGenInfo", true);
  desc.add<bool>("doChargedHadronMatching", false);
  desc.add<double>("genMatchDeltaRCut", 0.02);
  desc.add<bool>("isFullAOD", true);
  desc.add<edm::InputTag>("tracks", edm::InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("packedGenParticles", edm::InputTag(""));

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TrackAnalyzer);
