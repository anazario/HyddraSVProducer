// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      HyddraSVAnalyzer
//
// Description: Standalone EDAnalyzer that reads HYDDRA leptonic and hadronic
//              displaced vertex collections and writes them to a flat TTree,
//              producing output identical to KUCMSHyddraSV.
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <limits>

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
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Helper classes from KUCMSNtupleizer
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedTrackSCPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchTracksToSC.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"

// Local helper classes
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

// ROOT
#include "TTree.h"

class HyddraSVAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit HyddraSVAnalyzer(const edm::ParameterSet&);
  ~HyddraSVAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  void clearBranches();
  void fillVertexBranches(const reco::Vertex& vertex, const reco::Vertex& pv, bool isLeptonic, int vtxIndex);

  // Gen matching helpers
  bool IsBronze(const reco::Vertex& vertex) const;
  bool IsSilver(const reco::Vertex& vertex) const;
  bool IsGold(const reco::Vertex& vertex) const;
  int FindGenVertexIndex(const reco::Vertex& vertex) const;
  int FindNearestGenVertexIndex(const reco::Vertex& vertex, double& distance) const;
  double matchRatio(const reco::Vertex& vertex, const GenVertex& genVertex) const;
  bool getSCMatch(const reco::Track& track, reco::SuperCluster& sc, double& deltaR) const;

  // Configuration
  bool hasGenInfo_;

  // Tokens
  edm::EDGetTokenT<reco::VertexCollection> leptonicVerticesToken_;
  edm::EDGetTokenT<reco::VertexCollection> hadronicVerticesToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<reco::TrackCollection> muonEnhancedTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> muonTracksToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> mergedSCsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;

  // Handles
  edm::Handle<reco::VertexCollection> leptonicVerticesHandle_;
  edm::Handle<reco::VertexCollection> hadronicVerticesHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;
  edm::Handle<reco::TrackCollection> muonEnhancedTracksHandle_;
  edm::Handle<reco::TrackCollection> muonTracksHandle_;
  edm::Handle<reco::SuperClusterCollection> mergedSCsHandle_;
  edm::Handle<reco::GenParticleCollection> genHandle_;

  // Track association
  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  // Track-SC matching (rebuilt each event)
  std::vector<MatchedTrackSCPair<reco::Track>> trackSCPairs_;
  reco::TrackCollection electronTracks_;

  // Gen matching state (rebuilt each event)
  GenVertices genVertices_;
  std::map<GenVertex, GenMatches> chargedMatches_;
  reco::TrackCollection signalTracks_;
  PairedObjectCollection<reco::Track, reco::GenParticle> allGenMatches_;  // All matches before deltaR cut

  // TTree
  TTree* tree_;

  // Branch variables - Reco vertex branches
  unsigned int nVertices_;
  std::vector<bool> isLeptonic_;
  std::vector<unsigned int> nTracks_;
  std::vector<float> x_;
  std::vector<float> y_;
  std::vector<float> z_;
  std::vector<float> mass_;
  std::vector<float> p_;
  std::vector<float> px_;
  std::vector<float> py_;
  std::vector<float> pz_;
  std::vector<float> pt_;
  std::vector<float> eta_;
  std::vector<float> phi_;
  std::vector<float> cosTheta_;
  std::vector<float> decayAngle_;
  std::vector<float> dxy_;
  std::vector<float> dxyError_;
  std::vector<float> chi2_;
  std::vector<float> normalizedChi2_;
  std::vector<float> ndof_;
  std::vector<int> sumCharge_;
  std::vector<float> cxx_;
  std::vector<float> cyy_;
  std::vector<float> czz_;
  std::vector<float> cxy_;
  std::vector<float> cxz_;
  std::vector<float> cyz_;
  std::vector<float> scMatchRatio_;
  std::vector<bool> passLooseMuonID_;
  std::vector<bool> passLooseElectronID_;
  std::vector<bool> passLooseID_;

  // Track branches
  std::vector<unsigned int> track_vertexIndex_;
  std::vector<unsigned int> track_trackIndex_;
  std::vector<float> track_trackCosTheta_;
  std::vector<float> track_trackCosThetaAtCM_;
  std::vector<float> track_SCDR_;
  std::vector<float> track_energySC_;
  std::vector<float> track_ratioPToEnergySC_;

  // Gen matching branches (reco vertex)
  std::vector<int> genVertexIndex_;
  std::vector<int> nearestGenVertexIndex_;
  std::vector<float> min3D_;
  std::vector<float> vtxMatchRatio_;
  std::vector<bool> isBronze_;
  std::vector<bool> isSilver_;
  std::vector<bool> isGold_;

  // Gen matching branches (tracks)
  std::vector<bool> track_isSignalTrack_;
  std::vector<bool> track_isSignalElectron_;
  std::vector<bool> track_isSignalMuon_;
  std::vector<float> track_genMatchDeltaR_;  // deltaR to matched gen particle

  // Gen matching diagnostic branches (all track-to-gen matches, before cuts)
  std::vector<float> genMatch_deltaR_;
  std::vector<float> genMatch_trackPt_;
  std::vector<float> genMatch_trackEta_;
  std::vector<float> genMatch_genPt_;
  std::vector<float> genMatch_genEta_;
  std::vector<int> genMatch_genPdgId_;
  std::vector<int> genMatch_genMotherPdgId_;

  // Gen vertex branches
  unsigned int genVertex_nTotal_;
  std::vector<unsigned int> genVertex_nTracks_;
  unsigned int genVertex_nElectron_;
  unsigned int genVertex_nMuon_;
  unsigned int genVertex_nHadronic_;
  std::vector<float> genVertex_mass_;
  std::vector<float> genVertex_x_;
  std::vector<float> genVertex_y_;
  std::vector<float> genVertex_z_;
  std::vector<float> genVertex_p_;
  std::vector<float> genVertex_px_;
  std::vector<float> genVertex_py_;
  std::vector<float> genVertex_pz_;
  std::vector<float> genVertex_pt_;
  std::vector<float> genVertex_eta_;
  std::vector<float> genVertex_phi_;
  std::vector<float> genVertex_dxy_;
  std::vector<bool> genVertex_isElectron_;
  std::vector<bool> genVertex_isMuon_;
  std::vector<bool> genVertex_isHadronic_;
  std::vector<bool> genVertex_passSelection_;
  std::vector<bool> genVertex_passSelectionAndCuts_;
};

HyddraSVAnalyzer::HyddraSVAnalyzer(const edm::ParameterSet& iConfig) :
  hasGenInfo_(iConfig.getParameter<bool>("hasGenInfo")),
  leptonicVerticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("leptonicVertices"))),
  hadronicVerticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("hadronicVertices"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  muonEnhancedTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  muonTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("muonTracks"))),
  mergedSCsToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("mergedSCs"))),
  transientTrackBuilder_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  magneticFieldToken_(esConsumes()),
  caloGeometryToken_(esConsumes())
{
  usesResource("TFileService");

  if(hasGenInfo_) {
    genToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  }

  // Initialize track associator parameters from config
  edm::ParameterSet trackAssocPSet = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  edm::ConsumesCollector iC = consumesCollector();
  trackAssocParameters_.loadParameters(trackAssocPSet, iC);
  trackAssociator_.useDefaultPropagator();
}

void HyddraSVAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "HyddraSV Tree");

  // Reco vertex branches
  tree_->Branch("HyddraSV_nVertices", &nVertices_);
  tree_->Branch("HyddraSV_isLeptonic", &isLeptonic_);
  tree_->Branch("HyddraSV_nTracks", &nTracks_);
  tree_->Branch("HyddraSV_x", &x_);
  tree_->Branch("HyddraSV_y", &y_);
  tree_->Branch("HyddraSV_z", &z_);
  tree_->Branch("HyddraSV_mass", &mass_);
  tree_->Branch("HyddraSV_p", &p_);
  tree_->Branch("HyddraSV_px", &px_);
  tree_->Branch("HyddraSV_py", &py_);
  tree_->Branch("HyddraSV_pz", &pz_);
  tree_->Branch("HyddraSV_pt", &pt_);
  tree_->Branch("HyddraSV_eta", &eta_);
  tree_->Branch("HyddraSV_phi", &phi_);
  tree_->Branch("HyddraSV_cosTheta", &cosTheta_);
  tree_->Branch("HyddraSV_decayAngle", &decayAngle_);
  tree_->Branch("HyddraSV_dxy", &dxy_);
  tree_->Branch("HyddraSV_dxyError", &dxyError_);
  tree_->Branch("HyddraSV_chi2", &chi2_);
  tree_->Branch("HyddraSV_normalizedChi2", &normalizedChi2_);
  tree_->Branch("HyddraSV_ndof", &ndof_);
  tree_->Branch("HyddraSV_sumCharge", &sumCharge_);
  tree_->Branch("HyddraSV_cxx", &cxx_);
  tree_->Branch("HyddraSV_cyy", &cyy_);
  tree_->Branch("HyddraSV_czz", &czz_);
  tree_->Branch("HyddraSV_cxy", &cxy_);
  tree_->Branch("HyddraSV_cxz", &cxz_);
  tree_->Branch("HyddraSV_cyz", &cyz_);
  tree_->Branch("HyddraSV_scMatchRatio", &scMatchRatio_);
  tree_->Branch("HyddraSV_passLooseMuonID", &passLooseMuonID_);
  tree_->Branch("HyddraSV_passLooseElectronID", &passLooseElectronID_);
  tree_->Branch("HyddraSV_passLooseID", &passLooseID_);

  // Track branches
  tree_->Branch("HyddraSVTrack_vertexIndex", &track_vertexIndex_);
  tree_->Branch("HyddraSVTrack_trackIndex", &track_trackIndex_);
  tree_->Branch("HyddraSVTrack_trackCosTheta", &track_trackCosTheta_);
  tree_->Branch("HyddraSVTrack_trackCosThetaAtCM", &track_trackCosThetaAtCM_);
  tree_->Branch("HyddraSVTrack_SCDR", &track_SCDR_);
  tree_->Branch("HyddraSVTrack_energySC", &track_energySC_);
  tree_->Branch("HyddraSVTrack_ratioPToEnergySC", &track_ratioPToEnergySC_);

  // Gen matching branches
  if(hasGenInfo_) {
    tree_->Branch("HyddraSV_genVertexIndex", &genVertexIndex_);
    tree_->Branch("HyddraSV_nearestGenVertexIndex", &nearestGenVertexIndex_);
    tree_->Branch("HyddraSV_min3D", &min3D_);
    tree_->Branch("HyddraSV_matchRatio", &vtxMatchRatio_);
    tree_->Branch("HyddraSV_isBronze", &isBronze_);
    tree_->Branch("HyddraSV_isSilver", &isSilver_);
    tree_->Branch("HyddraSV_isGold", &isGold_);

    tree_->Branch("HyddraSVTrack_isSignalTrack", &track_isSignalTrack_);
    tree_->Branch("HyddraSVTrack_isSignalElectron", &track_isSignalElectron_);
    tree_->Branch("HyddraSVTrack_isSignalMuon", &track_isSignalMuon_);
    tree_->Branch("HyddraSVTrack_genMatchDeltaR", &track_genMatchDeltaR_);

    // Gen matching diagnostics (all track-to-gen matches before deltaR cut)
    tree_->Branch("HyddraGenMatch_deltaR", &genMatch_deltaR_);
    tree_->Branch("HyddraGenMatch_trackPt", &genMatch_trackPt_);
    tree_->Branch("HyddraGenMatch_trackEta", &genMatch_trackEta_);
    tree_->Branch("HyddraGenMatch_genPt", &genMatch_genPt_);
    tree_->Branch("HyddraGenMatch_genEta", &genMatch_genEta_);
    tree_->Branch("HyddraGenMatch_genPdgId", &genMatch_genPdgId_);
    tree_->Branch("HyddraGenMatch_genMotherPdgId", &genMatch_genMotherPdgId_);

    tree_->Branch("HyddraGenVertex_nTotal", &genVertex_nTotal_);
    tree_->Branch("HyddraGenVertex_nTracks", &genVertex_nTracks_);
    tree_->Branch("HyddraGenVertex_nElectron", &genVertex_nElectron_);
    tree_->Branch("HyddraGenVertex_nMuon", &genVertex_nMuon_);
    tree_->Branch("HyddraGenVertex_nHadronic", &genVertex_nHadronic_);
    tree_->Branch("HyddraGenVertex_mass", &genVertex_mass_);
    tree_->Branch("HyddraGenVertex_x", &genVertex_x_);
    tree_->Branch("HyddraGenVertex_y", &genVertex_y_);
    tree_->Branch("HyddraGenVertex_z", &genVertex_z_);
    tree_->Branch("HyddraGenVertex_p", &genVertex_p_);
    tree_->Branch("HyddraGenVertex_px", &genVertex_px_);
    tree_->Branch("HyddraGenVertex_py", &genVertex_py_);
    tree_->Branch("HyddraGenVertex_pz", &genVertex_pz_);
    tree_->Branch("HyddraGenVertex_pt", &genVertex_pt_);
    tree_->Branch("HyddraGenVertex_eta", &genVertex_eta_);
    tree_->Branch("HyddraGenVertex_phi", &genVertex_phi_);
    tree_->Branch("HyddraGenVertex_dxy", &genVertex_dxy_);
    tree_->Branch("HyddraGenVertex_isElectron", &genVertex_isElectron_);
    tree_->Branch("HyddraGenVertex_isMuon", &genVertex_isMuon_);
    tree_->Branch("HyddraGenVertex_isHadronic", &genVertex_isHadronic_);
    tree_->Branch("HyddraGenVertex_passSelection", &genVertex_passSelection_);
    tree_->Branch("HyddraGenVertex_passSelectionAndCuts", &genVertex_passSelectionAndCuts_);
  }
}

void HyddraSVAnalyzer::clearBranches() {
  nVertices_ = 0;
  isLeptonic_.clear();
  nTracks_.clear();
  x_.clear();
  y_.clear();
  z_.clear();
  mass_.clear();
  p_.clear();
  px_.clear();
  py_.clear();
  pz_.clear();
  pt_.clear();
  eta_.clear();
  phi_.clear();
  cosTheta_.clear();
  decayAngle_.clear();
  dxy_.clear();
  dxyError_.clear();
  chi2_.clear();
  normalizedChi2_.clear();
  ndof_.clear();
  sumCharge_.clear();
  cxx_.clear();
  cyy_.clear();
  czz_.clear();
  cxy_.clear();
  cxz_.clear();
  cyz_.clear();
  scMatchRatio_.clear();
  passLooseMuonID_.clear();
  passLooseElectronID_.clear();
  passLooseID_.clear();

  track_vertexIndex_.clear();
  track_trackIndex_.clear();
  track_trackCosTheta_.clear();
  track_trackCosThetaAtCM_.clear();
  track_SCDR_.clear();
  track_energySC_.clear();
  track_ratioPToEnergySC_.clear();

  if(hasGenInfo_) {
    genVertexIndex_.clear();
    nearestGenVertexIndex_.clear();
    min3D_.clear();
    vtxMatchRatio_.clear();
    isBronze_.clear();
    isSilver_.clear();
    isGold_.clear();

    track_isSignalTrack_.clear();
    track_isSignalElectron_.clear();
    track_isSignalMuon_.clear();
    track_genMatchDeltaR_.clear();

    genMatch_deltaR_.clear();
    genMatch_trackPt_.clear();
    genMatch_trackEta_.clear();
    genMatch_genPt_.clear();
    genMatch_genEta_.clear();
    genMatch_genPdgId_.clear();
    genMatch_genMotherPdgId_.clear();

    genVertex_nTotal_ = 0;
    genVertex_nTracks_.clear();
    genVertex_nElectron_ = 0;
    genVertex_nMuon_ = 0;
    genVertex_nHadronic_ = 0;
    genVertex_mass_.clear();
    genVertex_x_.clear();
    genVertex_y_.clear();
    genVertex_z_.clear();
    genVertex_p_.clear();
    genVertex_px_.clear();
    genVertex_py_.clear();
    genVertex_pz_.clear();
    genVertex_pt_.clear();
    genVertex_eta_.clear();
    genVertex_phi_.clear();
    genVertex_dxy_.clear();
    genVertex_isElectron_.clear();
    genVertex_isMuon_.clear();
    genVertex_isHadronic_.clear();
    genVertex_passSelection_.clear();
    genVertex_passSelectionAndCuts_.clear();
  }
}

void HyddraSVAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clearBranches();

  // Get handles
  iEvent.getByToken(leptonicVerticesToken_, leptonicVerticesHandle_);
  iEvent.getByToken(hadronicVerticesToken_, hadronicVerticesHandle_);
  iEvent.getByToken(pvToken_, pvHandle_);
  iEvent.getByToken(muonEnhancedTracksToken_, muonEnhancedTracksHandle_);
  iEvent.getByToken(muonTracksToken_, muonTracksHandle_);
  iEvent.getByToken(mergedSCsToken_, mergedSCsHandle_);

  // Clear per-event state
  genVertices_.clear();
  chargedMatches_.clear();
  signalTracks_.clear();
  electronTracks_.clear();
  trackSCPairs_.clear();
  allGenMatches_.clear();

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  // Collect all tracks from both vertex collections for SC matching
  reco::TrackCollection allVertexTracks;
  if(leptonicVerticesHandle_.isValid()) {
    for(const auto& vertex : *leptonicVerticesHandle_) {
      for(auto it = vertex.tracks_begin(); it != vertex.tracks_end(); ++it) {
        allVertexTracks.emplace_back(**it);
      }
    }
  }
  if(hadronicVerticesHandle_.isValid()) {
    for(const auto& vertex : *hadronicVerticesHandle_) {
      for(auto it = vertex.tracks_begin(); it != vertex.tracks_end(); ++it) {
        allVertexTracks.emplace_back(**it);
      }
    }
  }

  // Build track-SC matching
  MatchTracksToSC<reco::Track> scAssigner(iEvent, iSetup, magfield, ecalGeometry, trackAssocParameters_, allVertexTracks, *mergedSCsHandle_);
  trackSCPairs_ = scAssigner.GetMatchedTrackSCPairs();

  // Build electron tracks collection
  for(const auto& pair : trackSCPairs_) {
    if(pair.GetDeltaR() < 0.04)
      electronTracks_.emplace_back(pair.GetTrack());
  }

  if(hasGenInfo_) {
    iEvent.getByToken(genToken_, genHandle_);

    // Build transient tracks from the muon-enhanced track collection
    std::vector<reco::TransientTrack> ttracks;
    for(const auto &track : *muonEnhancedTracksHandle_)
      ttracks.emplace_back(ttBuilder->build(track));

    // Build gen vertices using the same approach as KUCMSDisplacedVertex
    GenVertices allSignalSVs(*genHandle_);
    DeltaRGenMatchHungarian<reco::TransientTrack> assigner(ttracks, allSignalSVs.getAllGenParticles());

    // Do a separate matching with reco::Track for diagnostics (no deltaR cut applied here)
    std::vector<reco::Track> trackVec(muonEnhancedTracksHandle_->begin(), muonEnhancedTracksHandle_->end());
    DeltaRGenMatchHungarian<reco::Track> trackAssigner(trackVec, allSignalSVs.getAllGenParticles());
    allGenMatches_ = trackAssigner.GetPairedObjects();

    // Fill gen match diagnostic branches (all matches, no deltaR cut)
    for(const auto &pair : allGenMatches_) {
      genMatch_deltaR_.push_back(float(pair.GetDeltaR()));
      genMatch_trackPt_.push_back(float(pair.GetObjectA().pt()));
      genMatch_trackEta_.push_back(float(pair.GetObjectA().eta()));
      genMatch_genPt_.push_back(float(pair.GetObjectB().pt()));
      genMatch_genEta_.push_back(float(pair.GetObjectB().eta()));
      genMatch_genPdgId_.push_back(int(pair.GetObjectB().pdgId()));
      // Get mother pdgId
      int motherPdgId = 0;
      if(pair.GetObjectB().mother()) {
        motherPdgId = pair.GetObjectB().mother()->pdgId();
      }
      genMatch_genMotherPdgId_.push_back(motherPdgId);
    }

    genVertices_ = GenVertices(assigner.GetPairedObjects().ConvertFromTTracks(), 0.02);
    allSignalSVs += genVertices_;
    genVertices_ = allSignalSVs;

    // Build charged particle matches for matchRatio computation
    for(const auto &genVertex : genVertices_) {
      DeltaRGenMatchHungarian<reco::Track> chargedParticleAssigner(*muonEnhancedTracksHandle_, genVertex.getStableChargedDaughters(*genHandle_));
      chargedMatches_[genVertex] = chargedParticleAssigner.GetPairedObjects();

      if(!genVertex.hasTracks()) continue;

      for(const auto &pair : genVertex.genMatches())
        signalTracks_.emplace_back(pair.GetObjectA());
    }
  }

  // Process vertices
  const reco::Vertex& pv = pvHandle_->at(0);

  unsigned int nTotal = 0;
  int vtxIndex = 0;

  if(leptonicVerticesHandle_.isValid()) {
    for(const auto& vertex : *leptonicVerticesHandle_) {
      fillVertexBranches(vertex, pv, true, vtxIndex);
      nTotal++;
      vtxIndex++;
    }
  }

  if(hadronicVerticesHandle_.isValid()) {
    for(const auto& vertex : *hadronicVerticesHandle_) {
      fillVertexBranches(vertex, pv, false, vtxIndex);
      nTotal++;
      vtxIndex++;
    }
  }

  nVertices_ = nTotal;

  // Gen vertex branches
  if(hasGenInfo_) {
    GenVertices allSignalSVs(*genHandle_);
    allSignalSVs += genVertices_;

    int nSigElectrons(0), nSigMuons(0), nSigHadrons(0);
    genVertex_nTotal_ = unsigned(allSignalSVs.size());

    for(const auto &genVertex : allSignalSVs) {

      // Check if this gen vertex passes selection and cuts
      bool passSelectionAndCuts = false;
      if(genVertex.hasTracks()) {
        auto checkVertex = [&](const reco::Vertex& recoVtx) {
          for(const auto& pair : genVertex.genMatches()) {
            for(auto it = recoVtx.tracks_begin(); it != recoVtx.tracks_end(); ++it) {
              if(TrackHelper::SameTrack(pair.GetObjectA(), **it))
                return true;
            }
          }
          return false;
        };

        if(leptonicVerticesHandle_.isValid()) {
          for(const auto& vtx : *leptonicVerticesHandle_) {
            if(checkVertex(vtx)) { passSelectionAndCuts = true; break; }
          }
        }
        if(!passSelectionAndCuts && hadronicVerticesHandle_.isValid()) {
          for(const auto& vtx : *hadronicVerticesHandle_) {
            if(checkVertex(vtx)) { passSelectionAndCuts = true; break; }
          }
        }
      }

      genVertex_nTracks_.push_back(unsigned(genVertex.tracks().size()));
      genVertex_mass_.push_back(float(genVertex.mass()));
      genVertex_x_.push_back(float(genVertex.x()));
      genVertex_y_.push_back(float(genVertex.y()));
      genVertex_z_.push_back(float(genVertex.z()));
      genVertex_p_.push_back(float(genVertex.p()));
      genVertex_px_.push_back(float(genVertex.px()));
      genVertex_py_.push_back(float(genVertex.py()));
      genVertex_pz_.push_back(float(genVertex.pz()));
      genVertex_pt_.push_back(float(genVertex.pt()));
      genVertex_eta_.push_back(float(genVertex.eta()));
      genVertex_phi_.push_back(float(genVertex.phi()));
      genVertex_dxy_.push_back(float(genVertex.dxy()));
      genVertex_isElectron_.push_back(bool(genVertex.isGenElectron()));
      genVertex_isMuon_.push_back(bool(genVertex.isGenMuon()));
      genVertex_isHadronic_.push_back(bool(genVertex.isGenHadronic()));
      genVertex_passSelection_.push_back(bool(genVertex.hasTracks()));
      genVertex_passSelectionAndCuts_.push_back(passSelectionAndCuts);

      if(genVertex.isGenElectron()) nSigElectrons++;
      if(genVertex.isGenMuon()) nSigMuons++;
      if(genVertex.isGenHadronic()) nSigHadrons++;
    }

    genVertex_nElectron_ = unsigned(nSigElectrons);
    genVertex_nMuon_ = unsigned(nSigMuons);
    genVertex_nHadronic_ = unsigned(nSigHadrons);
  }

  tree_->Fill();
}

void HyddraSVAnalyzer::fillVertexBranches(const reco::Vertex& vertex, const reco::Vertex& pv, bool isLeptonic, int vtxIndex) {

  LorentzVec vertex4Vec(VertexHelper::GetVertex4Vector(vertex));

  const bool passLooseMuonID(vertex.tracksSize() == 2 && VertexHelper::CountInstances(vertex, *muonTracksHandle_) == 2);
  const bool passLooseElectronID(vertex.tracksSize() == 2 && VertexHelper::CountInstances(vertex, electronTracks_) == 2);
  const bool passLooseID(vertex.tracksSize() < 4 ? (passLooseMuonID || passLooseElectronID) : true);

  isLeptonic_.push_back(isLeptonic);
  nTracks_.push_back(unsigned(vertex.tracksSize()));
  x_.push_back(float(vertex.x()));
  y_.push_back(float(vertex.y()));
  z_.push_back(float(vertex.z()));
  mass_.push_back(float(vertex4Vec.M()));
  p_.push_back(float(vertex4Vec.P()));
  px_.push_back(float(vertex4Vec.px()));
  py_.push_back(float(vertex4Vec.py()));
  pz_.push_back(float(vertex4Vec.pz()));
  pt_.push_back(float(vertex4Vec.pt()));
  eta_.push_back(float(vertex4Vec.eta()));
  phi_.push_back(float(vertex4Vec.phi()));
  cosTheta_.push_back(float(VertexHelper::CalculateCosTheta(pv, vertex)));
  decayAngle_.push_back(float(VertexHelper::CalculateDecayAngle(vertex)));
  dxy_.push_back(float(VertexHelper::CalculateDxy(vertex, pv)));
  dxyError_.push_back(float(VertexHelper::CalculateDxyError(vertex, pv)));
  chi2_.push_back(float(vertex.chi2()));
  normalizedChi2_.push_back(float(vertex.normalizedChi2()));
  ndof_.push_back(float(vertex.ndof()));
  sumCharge_.push_back(int(VertexHelper::CalculateTotalCharge(vertex)));
  cxx_.push_back(float(vertex.error().At(0,0)));
  cyy_.push_back(float(vertex.error().At(1,1)));
  czz_.push_back(float(vertex.error().At(2,2)));
  cxy_.push_back(float(vertex.error().At(0,1)));
  cxz_.push_back(float(vertex.error().At(0,2)));
  cyz_.push_back(float(vertex.error().At(1,2)));
  scMatchRatio_.push_back(float(VertexHelper::CountInstances(vertex, electronTracks_) / float(vertex.tracksSize())));
  passLooseMuonID_.push_back(passLooseMuonID);
  passLooseElectronID_.push_back(passLooseElectronID);
  passLooseID_.push_back(passLooseID);

  // Gen matching
  if(hasGenInfo_) {
    double minDistance;
    int nearestIdx = FindNearestGenVertexIndex(vertex, minDistance);
    isBronze_.push_back(bool(IsBronze(vertex)));
    isSilver_.push_back(bool(IsSilver(vertex)));
    isGold_.push_back(bool(IsGold(vertex)));
    genVertexIndex_.push_back(int(FindGenVertexIndex(vertex)));
    nearestGenVertexIndex_.push_back(int(nearestIdx));
    min3D_.push_back(float(minDistance));
    vtxMatchRatio_.push_back(float(nearestIdx >= 0 ? matchRatio(vertex, genVertices_[nearestIdx]) : -1));
  }

  // Fill per-track branches
  for(auto it = vertex.tracks_begin(); it != vertex.tracks_end(); ++it) {
    const reco::Track track(**it);

    double deltaR(-1.);
    reco::SuperCluster sc;
    const bool isSCMatched(getSCMatch(track, sc, deltaR));

    track_vertexIndex_.push_back(unsigned(vtxIndex));
    track_trackIndex_.push_back(unsigned(TrackHelper::FindTrackIndex(track, *muonEnhancedTracksHandle_)));
    track_trackCosTheta_.push_back(float(TrackHelper::CalculateCosTheta(pv, vertex, track)));
    track_trackCosThetaAtCM_.push_back(float(VertexHelper::CalculateCMCosTheta(vertex, track)));
    track_SCDR_.push_back(float(deltaR));
    track_energySC_.push_back(float(isSCMatched ? sc.correctedEnergy() : -1.));
    track_ratioPToEnergySC_.push_back(float(isSCMatched ? track.p() / sc.correctedEnergy() : -1.));

    if(hasGenInfo_) {
      const bool isSignal(TrackHelper::FindTrackIndex(track, signalTracks_) >= 0);
      track_isSignalTrack_.push_back(isSignal);
      track_isSignalElectron_.push_back(isSignal ? genVertices_.getGenVertexFromTrack(track).isGenElectron() : false);
      track_isSignalMuon_.push_back(isSignal ? genVertices_.getGenVertexFromTrack(track).isGenMuon() : false);

      // Find deltaR for this track from the gen matching
      float matchDeltaR = -1.0;
      for(const auto &pair : allGenMatches_) {
        if(TrackHelper::SameTrack(track, pair.GetObjectA())) {
          matchDeltaR = float(pair.GetDeltaR());
          break;
        }
      }
      track_genMatchDeltaR_.push_back(matchDeltaR);
    }
  }
}

// ---- Gen matching methods --------------------------------------------------------

bool HyddraSVAnalyzer::IsGold(const reco::Vertex& vertex) const {
  for(const auto &genVertex : genVertices_) {
    if(genVertex.isGold(vertex)) return true;
  }
  return false;
}

bool HyddraSVAnalyzer::IsSilver(const reco::Vertex& vertex) const {
  for(const auto &genVertex : genVertices_) {
    if(genVertex.isSilver(vertex)) return true;
  }
  return false;
}

bool HyddraSVAnalyzer::IsBronze(const reco::Vertex& vertex) const {
  for(const auto &genVertex : genVertices_) {
    if(genVertex.isBronze(vertex)) return true;
  }
  return false;
}

int HyddraSVAnalyzer::FindGenVertexIndex(const reco::Vertex& vertex) const {
  int index(0);
  for(const auto &genVertex : genVertices_) {
    if(genVertex.hasTracks() && (genVertex.isBronze(vertex) || genVertex.isSilver(vertex) || genVertex.isGold(vertex)))
      return index;
    index++;
  }
  return -1;
}

int HyddraSVAnalyzer::FindNearestGenVertexIndex(const reco::Vertex& vertex, double& distance) const {
  int index(0), minIndex(-1);
  distance = std::numeric_limits<double>::max();
  for(const auto &genVertex : genVertices_) {
    if(!genVertex.isGenHadronic()) continue;

    const double tempMin(genVertex.distance3D(vertex));
    if(tempMin < distance) {
      distance = tempMin;
      minIndex = index;
    }
    index++;
  }

  if(minIndex < 0)
    distance = -1.;

  return minIndex;
}

double HyddraSVAnalyzer::matchRatio(const reco::Vertex& vertex, const GenVertex& genVertex) const {

  if(chargedMatches_.find(genVertex) == chargedMatches_.end())
    return -1.;

  int count(0);
  for(const auto &pair : chargedMatches_.at(genVertex)) {
    reco::TrackRef ref(TrackHelper::GetTrackRef(pair.GetObjectA(), muonEnhancedTracksHandle_));

    const double genPt(pair.GetObjectB().pt());
    const double trackPt(pair.GetObjectA().pt());
    const double relPtDiff((genPt - trackPt) / genPt);
    const double deltaR(pair.GetDeltaR());

    // Check if this track is in the reco vertex
    bool inVertex(false);
    for(auto it = vertex.tracks_begin(); it != vertex.tracks_end(); ++it) {
      if(TrackHelper::SameTrack(*ref, **it)) {
        inVertex = true;
        break;
      }
    }

    if(inVertex && ((deltaR < 0.01 && fabs(relPtDiff) < (-10 * deltaR + 0.15)) ||
                    (deltaR > 0.01 && deltaR < 0.1 && fabs(relPtDiff) < (-deltaR / 2 + 0.055)))) {
      count++;
    }
  }
  return double(count) / vertex.tracksSize();
}

bool HyddraSVAnalyzer::getSCMatch(const reco::Track& track, reco::SuperCluster& sc, double& deltaR) const {

  bool isMatched(false);
  for(const auto &pair : trackSCPairs_) {
    if(TrackHelper::SameTrack(track, pair.GetTrack())) {
      isMatched = true;
      sc = pair.GetSuperCluster();
      deltaR = pair.GetDeltaR();
      break;
    }
  }

  return isMatched;
}

void HyddraSVAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<bool>("hasGenInfo", true);
  desc.add<edm::InputTag>("leptonicVertices", edm::InputTag("hyddraSVs", "leptonicVertices"));
  desc.add<edm::InputTag>("hadronicVertices", edm::InputTag("hyddraSVs", "hadronicVertices"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"));
  desc.add<edm::InputTag>("muonTracks", edm::InputTag("displacedGlobalMuons"));
  desc.add<edm::InputTag>("mergedSCs", edm::InputTag("ecalTracks", "displacedElectronSCs"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));

  // Track associator parameters
  edm::ParameterSetDescription trackAssocDesc;
  trackAssocDesc.add<bool>("useEcal", true);
  trackAssocDesc.add<bool>("useHcal", false);
  trackAssocDesc.add<bool>("useHO", false);
  trackAssocDesc.add<bool>("useCalo", false);
  trackAssocDesc.add<bool>("useMuon", false);
  trackAssocDesc.add<bool>("usePreshower", false);
  trackAssocDesc.add<double>("dREcal", 1.0);
  trackAssocDesc.add<double>("dRHcal", 1.0);
  trackAssocDesc.add<double>("dRMuon", 1.0);
  trackAssocDesc.add<double>("dRPreshowerPreselection", 0.2);
  trackAssocDesc.add<double>("dRMuonPreselection", 0.2);
  trackAssocDesc.add<double>("dREcalPreselection", 0.5);
  trackAssocDesc.add<double>("dRHcalPreselection", 0.5);
  trackAssocDesc.add<bool>("accountForTrajectoryChangeCalo", false);
  trackAssocDesc.add<edm::InputTag>("EBRecHitCollectionLabel", edm::InputTag("ecalRecHit", "EcalRecHitsEB"));
  trackAssocDesc.add<edm::InputTag>("EERecHitCollectionLabel", edm::InputTag("ecalRecHit", "EcalRecHitsEE"));
  trackAssocDesc.add<edm::InputTag>("HBHERecHitCollectionLabel", edm::InputTag("hbhereco"));
  trackAssocDesc.add<edm::InputTag>("HORecHitCollectionLabel", edm::InputTag("horeco"));
  trackAssocDesc.setAllowAnything();
  desc.add<edm::ParameterSetDescription>("TrackAssociatorParameters", trackAssocDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVAnalyzer);
