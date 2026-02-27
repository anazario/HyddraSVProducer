// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      HyddraSVsDiagnosticAnalyzer
//
// Description: Reads the 5 intermediate leptonic vertex collections produced by
//              HyddraSVsDiagnosticProducer and gen-matches each one to trace
//              where signal di-lepton vertices are lost at each pipeline stage.
//              Writes three TTrees:
//                genFunnel    — one row per gen vertex per event
//                stageCounts  — one row per event (total + matched counts)
//                cleaningTracks — one row per event, per-track cleaning variables
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <limits>
#include <array>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// Gen matching and helpers from KUCMSNtupleizer
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

// Kalman fitter (for cleaningTracks compatibility computation)
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// ROOT
#include "TTree.h"

// ─── Stage index constants ────────────────────────────────────────────────────
static constexpr size_t kSeed    = 0;
static constexpr size_t kMerged  = 1;
static constexpr size_t kCleaned = 2;
static constexpr size_t kDisambig= 3;
static constexpr size_t kFiltered= 4;
static constexpr size_t kNStages = 5;

// ─── Helper: result of gen-matching against one stage collection ──────────────
struct StageMatch {
  bool isGold   = false;
  bool isSilver = false;
  bool isBronze = false;
  const reco::Vertex* bestVertex = nullptr;  // highest-quality match
};

// ─────────────────────────────────────────────────────────────────────────────
class HyddraSVsDiagnosticAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit HyddraSVsDiagnosticAnalyzer(const edm::ParameterSet&);
  ~HyddraSVsDiagnosticAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  void clearBranches();

  // Gen matching helpers
  StageMatch findStageMatch(const GenVertex& genVtx,
                            const reco::VertexCollection& coll) const;
  float computeMatchRatio(const reco::Vertex& vtx) const;
  float computeDxySignif(const reco::Vertex& vtx, const reco::Vertex& pv) const;
  void fillRecoPropsAt(size_t stageIdx, const reco::Vertex& vtx,
                       const reco::Vertex& pv);

  reco::GenParticleCollection getStableChargedDaughtersFromPacked(
      const GenVertex& genVertex,
      const std::vector<pat::PackedGenParticle>& packed) const;

  // ── Configuration ──────────────────────────────────────────────────────────
  bool hasGenInfo_;
  bool isFullAOD_;

  // ── Tokens ─────────────────────────────────────────────────────────────────
  std::array<edm::EDGetTokenT<reco::VertexCollection>, kNStages> stageTokens_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<reco::TrackCollection>  tracksToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> packedGenToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // ── Handles ────────────────────────────────────────────────────────────────
  std::array<edm::Handle<reco::VertexCollection>, kNStages> stageHandles_;
  edm::Handle<reco::VertexCollection> pvHandle_;
  edm::Handle<reco::TrackCollection>  tracksHandle_;
  edm::Handle<reco::GenParticleCollection> genHandle_;
  edm::Handle<std::vector<pat::PackedGenParticle>> packedGenHandle_;

  // ── Per-event gen matching state ────────────────────────────────────────────
  GenVertices       genVertices_;
  reco::TrackCollection signalTracks_;

  // ── TTrees ─────────────────────────────────────────────────────────────────
  TTree* genFunnelTree_;
  TTree* stageCountsTree_;
  TTree* cleaningTracksTree_;
  TTree* allStageVtxTree_;
  TTree* seedTracksTree_;

  // ═══════════════════════════════════════════════════════════════════════════
  // genFunnel branches (vectors — one element per gen vertex)
  // ═══════════════════════════════════════════════════════════════════════════
  unsigned int genFunnel_n_;

  // Gen-level truth
  std::vector<float> genFunnel_dxy_;
  std::vector<float> genFunnel_pt_;
  std::vector<float> genFunnel_eta_;
  std::vector<float> genFunnel_mass_;
  std::vector<float> genFunnel_phi_;
  std::vector<bool>  genFunnel_isMuon_;
  std::vector<bool>  genFunnel_isElectron_;
  std::vector<bool>  genFunnel_hasTracks_;

  // Match quality booleans per stage
  std::vector<bool> genFunnel_gold_seed_,    genFunnel_gold_merged_,
                    genFunnel_gold_cleaned_,  genFunnel_gold_disambig_,
                    genFunnel_gold_filtered_;
  std::vector<bool> genFunnel_silver_seed_,  genFunnel_silver_merged_,
                    genFunnel_silver_cleaned_,genFunnel_silver_disambig_,
                    genFunnel_silver_filtered_;
  std::vector<bool> genFunnel_bronze_seed_,  genFunnel_bronze_merged_,
                    genFunnel_bronze_cleaned_,genFunnel_bronze_disambig_,
                    genFunnel_bronze_filtered_;

  // Loss stage (int): -1=survived, 0=never found, 1-4=lost after that stage index
  std::vector<int>  genFunnel_lostGoldAtStage_;
  std::vector<int>  genFunnel_lostSilverAtStage_;

  // Reco properties at each stage (-1 if not found at that stage)
  std::vector<float> genFunnel_mass_seed_,      genFunnel_mass_merged_,
                     genFunnel_mass_cleaned_,    genFunnel_mass_disambig_,
                     genFunnel_mass_filtered_;
  std::vector<float> genFunnel_dxySignif_seed_,  genFunnel_dxySignif_merged_,
                     genFunnel_dxySignif_cleaned_,genFunnel_dxySignif_disambig_,
                     genFunnel_dxySignif_filtered_;
  std::vector<float> genFunnel_normChi2_seed_,   genFunnel_normChi2_merged_,
                     genFunnel_normChi2_cleaned_, genFunnel_normChi2_disambig_,
                     genFunnel_normChi2_filtered_;
  std::vector<int>   genFunnel_nTracks_seed_,    genFunnel_nTracks_merged_,
                     genFunnel_nTracks_cleaned_,  genFunnel_nTracks_disambig_,
                     genFunnel_nTracks_filtered_;
  std::vector<float> genFunnel_cosTheta_seed_,   genFunnel_cosTheta_merged_,
                     genFunnel_cosTheta_cleaned_, genFunnel_cosTheta_disambig_,
                     genFunnel_cosTheta_filtered_;
  std::vector<float> genFunnel_matchRatio_seed_,  genFunnel_matchRatio_merged_,
                     genFunnel_matchRatio_cleaned_,genFunnel_matchRatio_disambig_,
                     genFunnel_matchRatio_filtered_;
  std::vector<float> genFunnel_decayAngle_seed_,   genFunnel_decayAngle_merged_,
                     genFunnel_decayAngle_cleaned_,  genFunnel_decayAngle_disambig_,
                     genFunnel_decayAngle_filtered_;
  std::vector<float> genFunnel_pOverE_seed_,       genFunnel_pOverE_merged_,
                     genFunnel_pOverE_cleaned_,      genFunnel_pOverE_disambig_,
                     genFunnel_pOverE_filtered_;

  // ═══════════════════════════════════════════════════════════════════════════
  // stageCounts branches (scalars — one per event)
  // ═══════════════════════════════════════════════════════════════════════════
  unsigned int stage_n_seed_,     stage_n_merged_,  stage_n_cleaned_,
               stage_n_disambig_, stage_n_filtered_;
  unsigned int stage_nGold_seed_,     stage_nGold_merged_,   stage_nGold_cleaned_,
               stage_nGold_disambig_, stage_nGold_filtered_;
  unsigned int stage_nSilver_seed_,   stage_nSilver_merged_,  stage_nSilver_cleaned_,
               stage_nSilver_disambig_,stage_nSilver_filtered_;
  unsigned int stage_nBronze_seed_,   stage_nBronze_merged_,  stage_nBronze_cleaned_,
               stage_nBronze_disambig_,stage_nBronze_filtered_;

  // ═══════════════════════════════════════════════════════════════════════════
  // cleaningTracks branches (vectors — one element per track in multi-track
  //   post-merge vertices)
  // ═══════════════════════════════════════════════════════════════════════════
  std::vector<float> cleanTrack_compatibility_;
  std::vector<float> cleanTrack_cosTheta_;
  std::vector<bool>  cleanTrack_isSignal_;
  std::vector<bool>  cleanTrack_isRemoved_;
  std::vector<bool>  cleanTrack_vtxIsGold_;
  std::vector<bool>  cleanTrack_vtxIsSilver_;
  std::vector<unsigned int> cleanTrack_vtxNTracks_;

  // ═══════════════════════════════════════════════════════════════════════════
  // allStageVtx branches (scalar — one TTree::Fill per vertex)
  // ═══════════════════════════════════════════════════════════════════════════
  int   stageVtx_stageIdx_;
  float stageVtx_mass_;
  float stageVtx_cosTheta_;
  float stageVtx_minTrackCosTheta_;
  float stageVtx_decayAngle_;
  float stageVtx_pOverE_;
  float stageVtx_dxySignif_;
  int   stageVtx_charge_;
  bool  stageVtx_isGold_;
  bool  stageVtx_isSilver_;
  bool  stageVtx_isBronze_;
  std::vector<float> stageVtx_trackCosTheta_;
  std::vector<float> stageVtx_trackCosThetaCM_;

  // ═══════════════════════════════════════════════════════════════════════════
  // seedTracks branches (vectors — one element per track in 2-track OS seeds)
  // ═══════════════════════════════════════════════════════════════════════════
  std::vector<float> seedTrack_cosTheta_;
  std::vector<bool>  seedTrack_isSignal_;
  std::vector<bool>  seedTrack_vtxIsGold_;

  // ═══════════════════════════════════════════════════════════════════════════
  // leptonicConfig branches (scalars — single Fill in beginJob)
  // ═══════════════════════════════════════════════════════════════════════════
  TTree* leptonicConfigTree_;
  double cfg_maxCompatibility_;
  double cfg_minCleanCosTheta_;
  double cfg_maxCleanCosTheta_;
  bool   cfg_invertCleanCosThetaCut_;
  bool   cfg_useDiagonalCut_;
  double cfg_cleanCutSlope_;
  double cfg_minTrackCosTheta_;
  double cfg_maxTrackCosThetaCM_Limit_;
  double cfg_maxTrackCosThetaCM_Intercept_;
  double cfg_trackCosThetaCM_Slope_;
  bool   cfg_requireChargeNeutrality_;
  double cfg_minDxySignificance_;
  double cfg_minVtxCosTheta_;
  double cfg_maxVtxCosTheta_;
  bool   cfg_useAbsVtxCosTheta_;
  double cfg_maxVtxDecayAngle_;
  bool   cfg_useAbsVtxDecayAngle_;
  bool   cfg_applyVtxDecayAngleCleaning_;
  bool   cfg_applyVtxDecayAngleFiltering_;
  double cfg_minMassFilter_;
  double cfg_minBetaFilter_;

};

// ─────────────────────────────────────────────────────────────────────────────
HyddraSVsDiagnosticAnalyzer::HyddraSVsDiagnosticAnalyzer(const edm::ParameterSet& iConfig) :
  hasGenInfo_(iConfig.getParameter<bool>("hasGenInfo")),
  isFullAOD_( iConfig.getParameter<bool>("isFullAOD")),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder")))
{
  usesResource("TFileService");

  // 5 stage input collections
  static const std::array<std::string, kNStages> kTagNames = {{
    "seedVertices", "mergedVertices", "cleanedVertices",
    "disambiguatedVertices", "filteredVertices"
  }};
  for (size_t i = 0; i < kNStages; ++i) {
    stageTokens_[i] = consumes<reco::VertexCollection>(
        iConfig.getParameter<edm::InputTag>(kTagNames[i]));
  }

  if (hasGenInfo_) {
    genToken_ = consumes<reco::GenParticleCollection>(
        iConfig.getParameter<edm::InputTag>("genParticles"));
    if (!isFullAOD_) {
      packedGenToken_ = consumes<std::vector<pat::PackedGenParticle>>(
          iConfig.getParameter<edm::InputTag>("packedGenParticles"));
    }
  }

  // Leptonic algorithm parameters — persisted as metadata in the output file.
  // Read from an optional 'leptonic' sub-PSet; falls back to defaults otherwise.
  cfg_maxCompatibility_             = 1.5;
  cfg_minCleanCosTheta_             = 0.5;
  cfg_maxCleanCosTheta_             = 1.0;
  cfg_invertCleanCosThetaCut_       = false;
  cfg_useDiagonalCut_               = false;
  cfg_cleanCutSlope_                = 0.0;
  cfg_minTrackCosTheta_             = 0.5;
  cfg_maxTrackCosThetaCM_Limit_     = 0.95;
  cfg_maxTrackCosThetaCM_Intercept_ = 1.8;
  cfg_trackCosThetaCM_Slope_        = -1.0;
  cfg_requireChargeNeutrality_      = true;
  cfg_minDxySignificance_           = 25.0;
  cfg_minVtxCosTheta_               = -1.0;
  cfg_maxVtxCosTheta_               =  1.0;
  cfg_useAbsVtxCosTheta_            = false;
  cfg_maxVtxDecayAngle_             = 1.0;
  cfg_useAbsVtxDecayAngle_          = false;
  cfg_applyVtxDecayAngleCleaning_   = false;
  cfg_applyVtxDecayAngleFiltering_  = false;
  cfg_minMassFilter_                = 0.0;
  cfg_minBetaFilter_                = 0.0;
  if (iConfig.exists("leptonic")) {
    const edm::ParameterSet& lep = iConfig.getParameterSet("leptonic");
    cfg_maxCompatibility_             = lep.getParameter<double>("maxCompatibility");
    cfg_minCleanCosTheta_             = lep.getParameter<double>("minCleanCosTheta");
    cfg_maxCleanCosTheta_             = lep.getParameter<double>("maxCleanCosTheta");
    cfg_invertCleanCosThetaCut_       = lep.getParameter<bool>("invertCleanCosThetaCut");
    cfg_useDiagonalCut_               = lep.getParameter<bool>("useDiagonalCut");
    cfg_cleanCutSlope_                = lep.getParameter<double>("cleanCutSlope");
    cfg_minTrackCosTheta_             = lep.getParameter<double>("minTrackCosTheta");
    cfg_maxTrackCosThetaCM_Limit_     = lep.getParameter<double>("maxTrackCosThetaCM_Limit");
    cfg_maxTrackCosThetaCM_Intercept_ = lep.getParameter<double>("maxTrackCosThetaCM_Intercept");
    cfg_trackCosThetaCM_Slope_        = lep.getParameter<double>("trackCosThetaCM_Slope");
    cfg_requireChargeNeutrality_      = lep.getParameter<bool>("requireChargeNeutrality");
    cfg_minDxySignificance_           = lep.getParameter<double>("minDxySignificance");
    cfg_minVtxCosTheta_               = lep.getParameter<double>("minVtxCosTheta");
    cfg_maxVtxCosTheta_               = lep.getParameter<double>("maxVtxCosTheta");
    cfg_useAbsVtxCosTheta_            = lep.getParameter<bool>("useAbsVtxCosTheta");
    cfg_maxVtxDecayAngle_             = lep.getParameter<double>("maxVtxDecayAngle");
    cfg_useAbsVtxDecayAngle_          = lep.getParameter<bool>("useAbsVtxDecayAngle");
    cfg_applyVtxDecayAngleCleaning_   = lep.getParameter<bool>("applyVtxDecayAngleCleaning");
    cfg_applyVtxDecayAngleFiltering_  = lep.getParameter<bool>("applyVtxDecayAngleFiltering");
    cfg_minMassFilter_                = lep.getParameter<double>("minMassFilter");
    cfg_minBetaFilter_                = lep.getParameter<double>("minBetaFilter");
  }
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsDiagnosticAnalyzer::beginJob() {
  edm::Service<TFileService> fs;

  // ── genFunnel ──────────────────────────────────────────────────────────────
  genFunnelTree_ = fs->make<TTree>("genFunnel", "Per-gen-vertex pipeline funnel");

  genFunnelTree_->Branch("GenFunnel_n",           &genFunnel_n_);
  genFunnelTree_->Branch("GenFunnel_dxy",          &genFunnel_dxy_);
  genFunnelTree_->Branch("GenFunnel_pt",           &genFunnel_pt_);
  genFunnelTree_->Branch("GenFunnel_eta",          &genFunnel_eta_);
  genFunnelTree_->Branch("GenFunnel_mass",         &genFunnel_mass_);
  genFunnelTree_->Branch("GenFunnel_phi",          &genFunnel_phi_);
  genFunnelTree_->Branch("GenFunnel_isMuon",       &genFunnel_isMuon_);
  genFunnelTree_->Branch("GenFunnel_isElectron",   &genFunnel_isElectron_);
  genFunnelTree_->Branch("GenFunnel_hasTracks",    &genFunnel_hasTracks_);

  genFunnelTree_->Branch("GenFunnel_gold_seed",        &genFunnel_gold_seed_);
  genFunnelTree_->Branch("GenFunnel_gold_merged",      &genFunnel_gold_merged_);
  genFunnelTree_->Branch("GenFunnel_gold_cleaned",     &genFunnel_gold_cleaned_);
  genFunnelTree_->Branch("GenFunnel_gold_disambig",    &genFunnel_gold_disambig_);
  genFunnelTree_->Branch("GenFunnel_gold_filtered",    &genFunnel_gold_filtered_);

  genFunnelTree_->Branch("GenFunnel_silver_seed",      &genFunnel_silver_seed_);
  genFunnelTree_->Branch("GenFunnel_silver_merged",    &genFunnel_silver_merged_);
  genFunnelTree_->Branch("GenFunnel_silver_cleaned",   &genFunnel_silver_cleaned_);
  genFunnelTree_->Branch("GenFunnel_silver_disambig",  &genFunnel_silver_disambig_);
  genFunnelTree_->Branch("GenFunnel_silver_filtered",  &genFunnel_silver_filtered_);

  genFunnelTree_->Branch("GenFunnel_bronze_seed",      &genFunnel_bronze_seed_);
  genFunnelTree_->Branch("GenFunnel_bronze_merged",    &genFunnel_bronze_merged_);
  genFunnelTree_->Branch("GenFunnel_bronze_cleaned",   &genFunnel_bronze_cleaned_);
  genFunnelTree_->Branch("GenFunnel_bronze_disambig",  &genFunnel_bronze_disambig_);
  genFunnelTree_->Branch("GenFunnel_bronze_filtered",  &genFunnel_bronze_filtered_);

  genFunnelTree_->Branch("GenFunnel_lostGoldAtStage",   &genFunnel_lostGoldAtStage_);
  genFunnelTree_->Branch("GenFunnel_lostSilverAtStage", &genFunnel_lostSilverAtStage_);

  genFunnelTree_->Branch("GenFunnel_mass_seed",         &genFunnel_mass_seed_);
  genFunnelTree_->Branch("GenFunnel_mass_merged",       &genFunnel_mass_merged_);
  genFunnelTree_->Branch("GenFunnel_mass_cleaned",      &genFunnel_mass_cleaned_);
  genFunnelTree_->Branch("GenFunnel_mass_disambig",     &genFunnel_mass_disambig_);
  genFunnelTree_->Branch("GenFunnel_mass_filtered",     &genFunnel_mass_filtered_);

  genFunnelTree_->Branch("GenFunnel_dxySignif_seed",    &genFunnel_dxySignif_seed_);
  genFunnelTree_->Branch("GenFunnel_dxySignif_merged",  &genFunnel_dxySignif_merged_);
  genFunnelTree_->Branch("GenFunnel_dxySignif_cleaned", &genFunnel_dxySignif_cleaned_);
  genFunnelTree_->Branch("GenFunnel_dxySignif_disambig",&genFunnel_dxySignif_disambig_);
  genFunnelTree_->Branch("GenFunnel_dxySignif_filtered",&genFunnel_dxySignif_filtered_);

  genFunnelTree_->Branch("GenFunnel_normChi2_seed",     &genFunnel_normChi2_seed_);
  genFunnelTree_->Branch("GenFunnel_normChi2_merged",   &genFunnel_normChi2_merged_);
  genFunnelTree_->Branch("GenFunnel_normChi2_cleaned",  &genFunnel_normChi2_cleaned_);
  genFunnelTree_->Branch("GenFunnel_normChi2_disambig", &genFunnel_normChi2_disambig_);
  genFunnelTree_->Branch("GenFunnel_normChi2_filtered", &genFunnel_normChi2_filtered_);

  genFunnelTree_->Branch("GenFunnel_nTracks_seed",      &genFunnel_nTracks_seed_);
  genFunnelTree_->Branch("GenFunnel_nTracks_merged",    &genFunnel_nTracks_merged_);
  genFunnelTree_->Branch("GenFunnel_nTracks_cleaned",   &genFunnel_nTracks_cleaned_);
  genFunnelTree_->Branch("GenFunnel_nTracks_disambig",  &genFunnel_nTracks_disambig_);
  genFunnelTree_->Branch("GenFunnel_nTracks_filtered",  &genFunnel_nTracks_filtered_);

  genFunnelTree_->Branch("GenFunnel_cosTheta_seed",     &genFunnel_cosTheta_seed_);
  genFunnelTree_->Branch("GenFunnel_cosTheta_merged",   &genFunnel_cosTheta_merged_);
  genFunnelTree_->Branch("GenFunnel_cosTheta_cleaned",  &genFunnel_cosTheta_cleaned_);
  genFunnelTree_->Branch("GenFunnel_cosTheta_disambig", &genFunnel_cosTheta_disambig_);
  genFunnelTree_->Branch("GenFunnel_cosTheta_filtered", &genFunnel_cosTheta_filtered_);

  genFunnelTree_->Branch("GenFunnel_matchRatio_seed",    &genFunnel_matchRatio_seed_);
  genFunnelTree_->Branch("GenFunnel_matchRatio_merged",  &genFunnel_matchRatio_merged_);
  genFunnelTree_->Branch("GenFunnel_matchRatio_cleaned", &genFunnel_matchRatio_cleaned_);
  genFunnelTree_->Branch("GenFunnel_matchRatio_disambig",&genFunnel_matchRatio_disambig_);
  genFunnelTree_->Branch("GenFunnel_matchRatio_filtered",&genFunnel_matchRatio_filtered_);

  genFunnelTree_->Branch("GenFunnel_decayAngle_seed",     &genFunnel_decayAngle_seed_);
  genFunnelTree_->Branch("GenFunnel_decayAngle_merged",   &genFunnel_decayAngle_merged_);
  genFunnelTree_->Branch("GenFunnel_decayAngle_cleaned",  &genFunnel_decayAngle_cleaned_);
  genFunnelTree_->Branch("GenFunnel_decayAngle_disambig", &genFunnel_decayAngle_disambig_);
  genFunnelTree_->Branch("GenFunnel_decayAngle_filtered", &genFunnel_decayAngle_filtered_);

  genFunnelTree_->Branch("GenFunnel_pOverE_seed",         &genFunnel_pOverE_seed_);
  genFunnelTree_->Branch("GenFunnel_pOverE_merged",       &genFunnel_pOverE_merged_);
  genFunnelTree_->Branch("GenFunnel_pOverE_cleaned",      &genFunnel_pOverE_cleaned_);
  genFunnelTree_->Branch("GenFunnel_pOverE_disambig",     &genFunnel_pOverE_disambig_);
  genFunnelTree_->Branch("GenFunnel_pOverE_filtered",     &genFunnel_pOverE_filtered_);

  // ── stageCounts ────────────────────────────────────────────────────────────
  stageCountsTree_ = fs->make<TTree>("stageCounts", "Per-event vertex counts per stage");

  stageCountsTree_->Branch("Stage_n_seed",         &stage_n_seed_);
  stageCountsTree_->Branch("Stage_n_merged",        &stage_n_merged_);
  stageCountsTree_->Branch("Stage_n_cleaned",       &stage_n_cleaned_);
  stageCountsTree_->Branch("Stage_n_disambig",      &stage_n_disambig_);
  stageCountsTree_->Branch("Stage_n_filtered",      &stage_n_filtered_);

  stageCountsTree_->Branch("Stage_nGold_seed",      &stage_nGold_seed_);
  stageCountsTree_->Branch("Stage_nGold_merged",    &stage_nGold_merged_);
  stageCountsTree_->Branch("Stage_nGold_cleaned",   &stage_nGold_cleaned_);
  stageCountsTree_->Branch("Stage_nGold_disambig",  &stage_nGold_disambig_);
  stageCountsTree_->Branch("Stage_nGold_filtered",  &stage_nGold_filtered_);

  stageCountsTree_->Branch("Stage_nSilver_seed",    &stage_nSilver_seed_);
  stageCountsTree_->Branch("Stage_nSilver_merged",  &stage_nSilver_merged_);
  stageCountsTree_->Branch("Stage_nSilver_cleaned", &stage_nSilver_cleaned_);
  stageCountsTree_->Branch("Stage_nSilver_disambig",&stage_nSilver_disambig_);
  stageCountsTree_->Branch("Stage_nSilver_filtered",&stage_nSilver_filtered_);

  stageCountsTree_->Branch("Stage_nBronze_seed",    &stage_nBronze_seed_);
  stageCountsTree_->Branch("Stage_nBronze_merged",  &stage_nBronze_merged_);
  stageCountsTree_->Branch("Stage_nBronze_cleaned", &stage_nBronze_cleaned_);
  stageCountsTree_->Branch("Stage_nBronze_disambig",&stage_nBronze_disambig_);
  stageCountsTree_->Branch("Stage_nBronze_filtered",&stage_nBronze_filtered_);

  // ── cleaningTracks ─────────────────────────────────────────────────────────
  cleaningTracksTree_ = fs->make<TTree>("cleaningTracks",
      "Per-track cleaning variables for post-merge multi-track vertices");

  cleaningTracksTree_->Branch("CleanTrack_compatibility", &cleanTrack_compatibility_);
  cleaningTracksTree_->Branch("CleanTrack_cosTheta",      &cleanTrack_cosTheta_);
  cleaningTracksTree_->Branch("CleanTrack_isSignal",      &cleanTrack_isSignal_);
  cleaningTracksTree_->Branch("CleanTrack_isRemoved",     &cleanTrack_isRemoved_);
  cleaningTracksTree_->Branch("CleanTrack_vtxIsGold",     &cleanTrack_vtxIsGold_);
  cleaningTracksTree_->Branch("CleanTrack_vtxIsSilver",   &cleanTrack_vtxIsSilver_);
  cleaningTracksTree_->Branch("CleanTrack_vtxNTracks",    &cleanTrack_vtxNTracks_);

  // ── allStageVtx ────────────────────────────────────────────────────────────
  allStageVtxTree_ = fs->make<TTree>("allStageVtx",
      "All vertices at each stage with match quality flags");
  allStageVtxTree_->Branch("StageVtx_stageIdx",        &stageVtx_stageIdx_);
  allStageVtxTree_->Branch("StageVtx_mass",             &stageVtx_mass_);
  allStageVtxTree_->Branch("StageVtx_cosTheta",         &stageVtx_cosTheta_);
  allStageVtxTree_->Branch("StageVtx_minTrackCosTheta", &stageVtx_minTrackCosTheta_);
  allStageVtxTree_->Branch("StageVtx_decayAngle",       &stageVtx_decayAngle_);
  allStageVtxTree_->Branch("StageVtx_pOverE",           &stageVtx_pOverE_);
  allStageVtxTree_->Branch("StageVtx_dxySignif",        &stageVtx_dxySignif_);
  allStageVtxTree_->Branch("StageVtx_charge",           &stageVtx_charge_);
  allStageVtxTree_->Branch("StageVtx_isGold",           &stageVtx_isGold_);
  allStageVtxTree_->Branch("StageVtx_isSilver",         &stageVtx_isSilver_);
  allStageVtxTree_->Branch("StageVtx_isBronze",         &stageVtx_isBronze_);
  allStageVtxTree_->Branch("StageVtx_trackCosTheta",    &stageVtx_trackCosTheta_);
  allStageVtxTree_->Branch("StageVtx_trackCosThetaCM",  &stageVtx_trackCosThetaCM_);

  // ── seedTracks ─────────────────────────────────────────────────────────────
  seedTracksTree_ = fs->make<TTree>("seedTracks",
      "Per-track info for tracks in 2-track OS seeded vertices");
  seedTracksTree_->Branch("SeedTrack_cosTheta",  &seedTrack_cosTheta_);
  seedTracksTree_->Branch("SeedTrack_isSignal",  &seedTrack_isSignal_);
  seedTracksTree_->Branch("SeedTrack_vtxIsGold", &seedTrack_vtxIsGold_);

  // ── leptonicConfig ─────────────────────────────────────────────────────────
  leptonicConfigTree_ = fs->make<TTree>("leptonicConfig",
      "Leptonic HYDDRA configuration parameters (single row)");
  leptonicConfigTree_->Branch("maxCompatibility",          &cfg_maxCompatibility_);
  leptonicConfigTree_->Branch("minCleanCosTheta",          &cfg_minCleanCosTheta_);
  leptonicConfigTree_->Branch("maxCleanCosTheta",          &cfg_maxCleanCosTheta_);
  leptonicConfigTree_->Branch("invertCleanCosThetaCut",    &cfg_invertCleanCosThetaCut_);
  leptonicConfigTree_->Branch("useDiagonalCut",            &cfg_useDiagonalCut_);
  leptonicConfigTree_->Branch("cleanCutSlope",             &cfg_cleanCutSlope_);
  leptonicConfigTree_->Branch("minTrackCosTheta",          &cfg_minTrackCosTheta_);
  leptonicConfigTree_->Branch("maxTrackCosThetaCM_Limit",  &cfg_maxTrackCosThetaCM_Limit_);
  leptonicConfigTree_->Branch("maxTrackCosThetaCM_Intercept", &cfg_maxTrackCosThetaCM_Intercept_);
  leptonicConfigTree_->Branch("trackCosThetaCM_Slope",     &cfg_trackCosThetaCM_Slope_);
  leptonicConfigTree_->Branch("requireChargeNeutrality",   &cfg_requireChargeNeutrality_);
  leptonicConfigTree_->Branch("minDxySignificance",        &cfg_minDxySignificance_);
  leptonicConfigTree_->Branch("minVtxCosTheta",            &cfg_minVtxCosTheta_);
  leptonicConfigTree_->Branch("maxVtxCosTheta",            &cfg_maxVtxCosTheta_);
  leptonicConfigTree_->Branch("useAbsVtxCosTheta",         &cfg_useAbsVtxCosTheta_);
  leptonicConfigTree_->Branch("maxVtxDecayAngle",          &cfg_maxVtxDecayAngle_);
  leptonicConfigTree_->Branch("useAbsVtxDecayAngle",       &cfg_useAbsVtxDecayAngle_);
  leptonicConfigTree_->Branch("applyVtxDecayAngleCleaning",  &cfg_applyVtxDecayAngleCleaning_);
  leptonicConfigTree_->Branch("applyVtxDecayAngleFiltering", &cfg_applyVtxDecayAngleFiltering_);
  leptonicConfigTree_->Branch("minMassFilter",               &cfg_minMassFilter_);
  leptonicConfigTree_->Branch("minBetaFilter",               &cfg_minBetaFilter_);
  leptonicConfigTree_->Fill();  // values set in constructor; written once here

}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsDiagnosticAnalyzer::clearBranches() {
  genFunnel_n_ = 0;
  genFunnel_dxy_.clear();   genFunnel_pt_.clear();  genFunnel_eta_.clear();
  genFunnel_mass_.clear();  genFunnel_phi_.clear();
  genFunnel_isMuon_.clear();genFunnel_isElectron_.clear();genFunnel_hasTracks_.clear();

  genFunnel_gold_seed_.clear();   genFunnel_gold_merged_.clear();
  genFunnel_gold_cleaned_.clear();genFunnel_gold_disambig_.clear();
  genFunnel_gold_filtered_.clear();
  genFunnel_silver_seed_.clear(); genFunnel_silver_merged_.clear();
  genFunnel_silver_cleaned_.clear();genFunnel_silver_disambig_.clear();
  genFunnel_silver_filtered_.clear();
  genFunnel_bronze_seed_.clear(); genFunnel_bronze_merged_.clear();
  genFunnel_bronze_cleaned_.clear();genFunnel_bronze_disambig_.clear();
  genFunnel_bronze_filtered_.clear();

  genFunnel_lostGoldAtStage_.clear();
  genFunnel_lostSilverAtStage_.clear();

  genFunnel_mass_seed_.clear();     genFunnel_mass_merged_.clear();
  genFunnel_mass_cleaned_.clear();  genFunnel_mass_disambig_.clear();
  genFunnel_mass_filtered_.clear();
  genFunnel_dxySignif_seed_.clear();  genFunnel_dxySignif_merged_.clear();
  genFunnel_dxySignif_cleaned_.clear();genFunnel_dxySignif_disambig_.clear();
  genFunnel_dxySignif_filtered_.clear();
  genFunnel_normChi2_seed_.clear();  genFunnel_normChi2_merged_.clear();
  genFunnel_normChi2_cleaned_.clear();genFunnel_normChi2_disambig_.clear();
  genFunnel_normChi2_filtered_.clear();
  genFunnel_nTracks_seed_.clear();   genFunnel_nTracks_merged_.clear();
  genFunnel_nTracks_cleaned_.clear();genFunnel_nTracks_disambig_.clear();
  genFunnel_nTracks_filtered_.clear();
  genFunnel_cosTheta_seed_.clear();  genFunnel_cosTheta_merged_.clear();
  genFunnel_cosTheta_cleaned_.clear();genFunnel_cosTheta_disambig_.clear();
  genFunnel_cosTheta_filtered_.clear();
  genFunnel_matchRatio_seed_.clear();  genFunnel_matchRatio_merged_.clear();
  genFunnel_matchRatio_cleaned_.clear();genFunnel_matchRatio_disambig_.clear();
  genFunnel_matchRatio_filtered_.clear();
  genFunnel_decayAngle_seed_.clear();  genFunnel_decayAngle_merged_.clear();
  genFunnel_decayAngle_cleaned_.clear();genFunnel_decayAngle_disambig_.clear();
  genFunnel_decayAngle_filtered_.clear();
  genFunnel_pOverE_seed_.clear();      genFunnel_pOverE_merged_.clear();
  genFunnel_pOverE_cleaned_.clear();   genFunnel_pOverE_disambig_.clear();
  genFunnel_pOverE_filtered_.clear();

  stage_n_seed_ = stage_n_merged_ = stage_n_cleaned_ =
  stage_n_disambig_ = stage_n_filtered_ = 0;
  stage_nGold_seed_ = stage_nGold_merged_ = stage_nGold_cleaned_ =
  stage_nGold_disambig_ = stage_nGold_filtered_ = 0;
  stage_nSilver_seed_ = stage_nSilver_merged_ = stage_nSilver_cleaned_ =
  stage_nSilver_disambig_ = stage_nSilver_filtered_ = 0;
  stage_nBronze_seed_ = stage_nBronze_merged_ = stage_nBronze_cleaned_ =
  stage_nBronze_disambig_ = stage_nBronze_filtered_ = 0;

  cleanTrack_compatibility_.clear();cleanTrack_cosTheta_.clear();
  cleanTrack_isSignal_.clear();     cleanTrack_isRemoved_.clear();
  cleanTrack_vtxIsGold_.clear();    cleanTrack_vtxIsSilver_.clear();
  cleanTrack_vtxNTracks_.clear();

  seedTrack_cosTheta_.clear();
  seedTrack_isSignal_.clear();
  seedTrack_vtxIsGold_.clear();
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: find the highest-quality match for a gen vertex in a stage collection
StageMatch HyddraSVsDiagnosticAnalyzer::findStageMatch(
    const GenVertex& genVtx,
    const reco::VertexCollection& coll) const {

  StageMatch result;
  int bestQuality = 0;

  for (const auto& vtx : coll) {
    int quality = 0;
    if      (genVtx.isGold(vtx))   quality = 3;
    else if (genVtx.isSilver(vtx)) quality = 2;
    else if (genVtx.isBronze(vtx)) quality = 1;

    if (quality == 3) result.isGold   = true;
    if (quality == 2) result.isSilver = true;
    if (quality == 1) result.isBronze = true;

    if (quality > bestQuality) {
      bestQuality = quality;
      result.bestVertex = &vtx;
    }
  }
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
float HyddraSVsDiagnosticAnalyzer::computeMatchRatio(const reco::Vertex& vtx) const {
  if (vtx.tracksSize() == 0) return 0.f;
  int nSignal = 0;
  for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
    reco::Track t(**it);
    if (TrackHelper::FindTrackIndex(t, signalTracks_) >= 0) nSignal++;
  }
  return float(nSignal) / float(vtx.tracksSize());
}

// ─────────────────────────────────────────────────────────────────────────────
float HyddraSVsDiagnosticAnalyzer::computeDxySignif(
    const reco::Vertex& vtx, const reco::Vertex& pv) const {
  float dxy    = float(VertexHelper::CalculateDxy(vtx, pv));
  float dxyErr = float(VertexHelper::CalculateDxyError(vtx, pv));
  return (dxyErr > 0) ? dxy / dxyErr : -1.f;
}

// ─────────────────────────────────────────────────────────────────────────────
// Fill reco property branches at a given stage index from the matched vertex
void HyddraSVsDiagnosticAnalyzer::fillRecoPropsAt(
    size_t s, const reco::Vertex& vtx, const reco::Vertex& pv) {

  // Indexed by stage: kSeed=0, kMerged=1, kCleaned=2, kDisambig=3, kFiltered=4
  std::array<std::vector<float>*, kNStages> massVecs = {{
      &genFunnel_mass_seed_,  &genFunnel_mass_merged_,
      &genFunnel_mass_cleaned_,&genFunnel_mass_disambig_,&genFunnel_mass_filtered_}};
  std::array<std::vector<float>*, kNStages> dxsVecs = {{
      &genFunnel_dxySignif_seed_,  &genFunnel_dxySignif_merged_,
      &genFunnel_dxySignif_cleaned_,&genFunnel_dxySignif_disambig_,&genFunnel_dxySignif_filtered_}};
  std::array<std::vector<float>*, kNStages> chi2Vecs = {{
      &genFunnel_normChi2_seed_,  &genFunnel_normChi2_merged_,
      &genFunnel_normChi2_cleaned_,&genFunnel_normChi2_disambig_,&genFunnel_normChi2_filtered_}};
  std::array<std::vector<int>*, kNStages> ntrkVecs = {{
      &genFunnel_nTracks_seed_,  &genFunnel_nTracks_merged_,
      &genFunnel_nTracks_cleaned_,&genFunnel_nTracks_disambig_,&genFunnel_nTracks_filtered_}};
  std::array<std::vector<float>*, kNStages> ctVecs = {{
      &genFunnel_cosTheta_seed_,  &genFunnel_cosTheta_merged_,
      &genFunnel_cosTheta_cleaned_,&genFunnel_cosTheta_disambig_,&genFunnel_cosTheta_filtered_}};
  std::array<std::vector<float>*, kNStages> mrVecs = {{
      &genFunnel_matchRatio_seed_,  &genFunnel_matchRatio_merged_,
      &genFunnel_matchRatio_cleaned_,&genFunnel_matchRatio_disambig_,&genFunnel_matchRatio_filtered_}};
  std::array<std::vector<float>*, kNStages> daVecs = {{
      &genFunnel_decayAngle_seed_,  &genFunnel_decayAngle_merged_,
      &genFunnel_decayAngle_cleaned_,&genFunnel_decayAngle_disambig_,&genFunnel_decayAngle_filtered_}};
  std::array<std::vector<float>*, kNStages> poeVecs = {{
      &genFunnel_pOverE_seed_,  &genFunnel_pOverE_merged_,
      &genFunnel_pOverE_cleaned_,&genFunnel_pOverE_disambig_,&genFunnel_pOverE_filtered_}};

  auto p4 = VertexHelper::GetVertex4Vector(vtx);
  double p  = p4.P();
  double m  = p4.M();
  double E  = std::sqrt(p*p + m*m);
  float  pOverE = (E > 1e-6) ? float(p / E) : -1.f;

  massVecs[s]->push_back(float(m));
  dxsVecs [s]->push_back(computeDxySignif(vtx, pv));
  chi2Vecs[s]->push_back(float(vtx.normalizedChi2()));
  ntrkVecs[s]->push_back(int(vtx.tracksSize()));
  ctVecs  [s]->push_back(float(VertexHelper::CalculateCosTheta(pv, vtx)));
  mrVecs  [s]->push_back(computeMatchRatio(vtx));
  daVecs  [s]->push_back(float(VertexHelper::CalculateDecayAngle(vtx)));
  poeVecs [s]->push_back(pOverE);
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsDiagnosticAnalyzer::analyze(const edm::Event& iEvent,
                                          const edm::EventSetup& iSetup) {
  clearBranches();

  // Get stage vertex collection handles
  for (size_t i = 0; i < kNStages; ++i)
    iEvent.getByToken(stageTokens_[i], stageHandles_[i]);

  iEvent.getByToken(pvToken_, pvHandle_);
  iEvent.getByToken(tracksToken_, tracksHandle_);

  const reco::Vertex& pv = pvHandle_->at(0);
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);

  // ── Fill stageCounts totals ────────────────────────────────────────────────
  auto safeSize = [&](size_t i) {
    return stageHandles_[i].isValid() ? unsigned(stageHandles_[i]->size()) : 0u;
  };
  stage_n_seed_     = safeSize(kSeed);
  stage_n_merged_   = safeSize(kMerged);
  stage_n_cleaned_  = safeSize(kCleaned);
  stage_n_disambig_ = safeSize(kDisambig);
  stage_n_filtered_ = safeSize(kFiltered);

  // ── Build gen vertices and signal tracks ───────────────────────────────────
  genVertices_.clear();
  signalTracks_.clear();

  if (hasGenInfo_) {
    iEvent.getByToken(genToken_, genHandle_);
    if (!isFullAOD_)
      iEvent.getByToken(packedGenToken_, packedGenHandle_);

    // Build transient tracks for Hungarian matching
    std::vector<reco::TransientTrack> ttracks;
    for (const auto& track : *tracksHandle_)
      ttracks.emplace_back(ttBuilder->build(track));

    GenVertices allSignalSVs(*genHandle_);
    DeltaRGenMatchHungarian<reco::TransientTrack> assigner(
        ttracks, allSignalSVs.getAllGenParticles());

    genVertices_ = GenVertices(assigner.GetPairedObjects().ConvertFromTTracks(), 0.02);
    allSignalSVs += genVertices_;
    genVertices_ = allSignalSVs;

    for (const auto& gv : genVertices_) {
      if (!gv.hasTracks()) continue;
      for (const auto& pair : gv.genMatches())
        signalTracks_.emplace_back(pair.GetObjectA());
    }
  }

  // ── Fill stageCounts gold/silver/bronze per vertex ─────────────────────────
  if (hasGenInfo_) {
    // Helper: count matches across a stage collection
    auto countMatches = [&](size_t si,
                            unsigned& nG, unsigned& nS, unsigned& nB) {
      if (!stageHandles_[si].isValid()) return;
      for (const auto& vtx : *stageHandles_[si]) {
        bool anyGold = false, anySilver = false, anyBronze = false;
        for (const auto& gv : genVertices_) {
          if (gv.isGold(vtx))        { anyGold   = true; break; }
          if (gv.isSilver(vtx))      anySilver = true;
          else if (gv.isBronze(vtx)) anyBronze = true;
        }
        if (anyGold)        nG++;
        else if (anySilver) nS++;
        else if (anyBronze) nB++;
      }
    };

    countMatches(kSeed,    stage_nGold_seed_,     stage_nSilver_seed_,    stage_nBronze_seed_);
    countMatches(kMerged,  stage_nGold_merged_,   stage_nSilver_merged_,  stage_nBronze_merged_);
    countMatches(kCleaned, stage_nGold_cleaned_,  stage_nSilver_cleaned_, stage_nBronze_cleaned_);
    countMatches(kDisambig,stage_nGold_disambig_,  stage_nSilver_disambig_,stage_nBronze_disambig_);
    countMatches(kFiltered,stage_nGold_filtered_,  stage_nSilver_filtered_,stage_nBronze_filtered_);
  }

  // ── genFunnel: one entry per gen vertex ────────────────────────────────────
  if (hasGenInfo_) {
    for (const auto& gv : genVertices_) {

      // Gen-level truth
      genFunnel_dxy_       .push_back(float(gv.dxy()));
      genFunnel_pt_        .push_back(float(gv.pt()));
      genFunnel_eta_       .push_back(float(gv.eta()));
      genFunnel_mass_      .push_back(float(gv.mass()));
      genFunnel_phi_       .push_back(float(gv.phi()));
      genFunnel_isMuon_    .push_back(bool(gv.isGenMuon()));
      genFunnel_isElectron_.push_back(bool(gv.isGenElectron()));
      genFunnel_hasTracks_ .push_back(bool(gv.hasTracks()));

      // Match at each stage
      std::array<StageMatch, kNStages> m;
      for (size_t si = 0; si < kNStages; ++si) {
        if (stageHandles_[si].isValid())
          m[si] = findStageMatch(gv, *stageHandles_[si]);
      }

      // Push booleans
      genFunnel_gold_seed_    .push_back(m[kSeed].isGold);
      genFunnel_gold_merged_  .push_back(m[kMerged].isGold);
      genFunnel_gold_cleaned_ .push_back(m[kCleaned].isGold);
      genFunnel_gold_disambig_.push_back(m[kDisambig].isGold);
      genFunnel_gold_filtered_.push_back(m[kFiltered].isGold);

      genFunnel_silver_seed_    .push_back(m[kSeed].isSilver);
      genFunnel_silver_merged_  .push_back(m[kMerged].isSilver);
      genFunnel_silver_cleaned_ .push_back(m[kCleaned].isSilver);
      genFunnel_silver_disambig_.push_back(m[kDisambig].isSilver);
      genFunnel_silver_filtered_.push_back(m[kFiltered].isSilver);

      genFunnel_bronze_seed_    .push_back(m[kSeed].isBronze);
      genFunnel_bronze_merged_  .push_back(m[kMerged].isBronze);
      genFunnel_bronze_cleaned_ .push_back(m[kCleaned].isBronze);
      genFunnel_bronze_disambig_.push_back(m[kDisambig].isBronze);
      genFunnel_bronze_filtered_.push_back(m[kFiltered].isBronze);

      // Loss stage encoding:
      //  -1 = survived (gold at filter)
      //   0 = never gold at any stage
      //   1 = was gold at seed, lost gold after seeding
      //   2 = was gold at merge, lost gold after merging  (etc.)
      auto computeLostStage = [&](bool useGold) -> int {
        auto hasMatch = [&](size_t si) {
          return useGold ? m[si].isGold : (m[si].isGold || m[si].isSilver);
        };
        if (hasMatch(kFiltered)) return -1;
        for (int si = int(kNStages) - 2; si >= 0; --si) {
          if (hasMatch(size_t(si))) return si + 1;
        }
        return 0;
      };
      genFunnel_lostGoldAtStage_  .push_back(computeLostStage(true));
      genFunnel_lostSilverAtStage_.push_back(computeLostStage(false));

      // Reco properties per stage (-1 sentinel if no match at that stage)
      // Uses the same indexed-pointer pattern as fillRecoPropsAt().
      std::array<std::vector<float>*, kNStages> massVecs = {{
          &genFunnel_mass_seed_,  &genFunnel_mass_merged_,
          &genFunnel_mass_cleaned_,&genFunnel_mass_disambig_,&genFunnel_mass_filtered_}};
      std::array<std::vector<float>*, kNStages> dxsVecs = {{
          &genFunnel_dxySignif_seed_,  &genFunnel_dxySignif_merged_,
          &genFunnel_dxySignif_cleaned_,&genFunnel_dxySignif_disambig_,&genFunnel_dxySignif_filtered_}};
      std::array<std::vector<float>*, kNStages> chi2Vecs = {{
          &genFunnel_normChi2_seed_,  &genFunnel_normChi2_merged_,
          &genFunnel_normChi2_cleaned_,&genFunnel_normChi2_disambig_,&genFunnel_normChi2_filtered_}};
      std::array<std::vector<int>*, kNStages> ntrkVecs = {{
          &genFunnel_nTracks_seed_,  &genFunnel_nTracks_merged_,
          &genFunnel_nTracks_cleaned_,&genFunnel_nTracks_disambig_,&genFunnel_nTracks_filtered_}};
      std::array<std::vector<float>*, kNStages> ctVecs = {{
          &genFunnel_cosTheta_seed_,  &genFunnel_cosTheta_merged_,
          &genFunnel_cosTheta_cleaned_,&genFunnel_cosTheta_disambig_,&genFunnel_cosTheta_filtered_}};
      std::array<std::vector<float>*, kNStages> mrVecs = {{
          &genFunnel_matchRatio_seed_,  &genFunnel_matchRatio_merged_,
          &genFunnel_matchRatio_cleaned_,&genFunnel_matchRatio_disambig_,&genFunnel_matchRatio_filtered_}};
      std::array<std::vector<float>*, kNStages> daVecs = {{
          &genFunnel_decayAngle_seed_,  &genFunnel_decayAngle_merged_,
          &genFunnel_decayAngle_cleaned_,&genFunnel_decayAngle_disambig_,&genFunnel_decayAngle_filtered_}};
      std::array<std::vector<float>*, kNStages> poeVecs = {{
          &genFunnel_pOverE_seed_,  &genFunnel_pOverE_merged_,
          &genFunnel_pOverE_cleaned_,&genFunnel_pOverE_disambig_,&genFunnel_pOverE_filtered_}};

      for (size_t si = 0; si < kNStages; ++si) {
        if (m[si].bestVertex) {
          fillRecoPropsAt(si, *m[si].bestVertex, pv);
        } else {
          massVecs[si]->push_back(-1.f);
          dxsVecs [si]->push_back(-1.f);
          chi2Vecs[si]->push_back(-1.f);
          ntrkVecs[si]->push_back(-1);
          ctVecs  [si]->push_back(-1.f);
          mrVecs  [si]->push_back(-1.f);
          daVecs  [si]->push_back(-1.f);
          poeVecs [si]->push_back(-1.f);
        }
      }
    } // end gen vertex loop

    genFunnel_n_ = unsigned(genVertices_.size());
  }

  // ── cleaningTracks: tracks in post-merge vertices with >2 tracks ──────────
  if (hasGenInfo_ && stageHandles_[kMerged].isValid() && stageHandles_[kCleaned].isValid()) {

    // Build set of all track refs present in the cleaned collection
    std::set<reco::TrackRef> cleanedTrackRefs;
    for (const auto& vtx : *stageHandles_[kCleaned]) {
      for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
        cleanedTrackRefs.insert(it->castTo<reco::TrackRef>());
      }
    }

    for (const auto& vtx : *stageHandles_[kMerged]) {
      if (vtx.tracksSize() <= 2) continue;  // 2-track vertices skip track removal

      // Determine gen match quality for this vertex
      bool vtxIsGold = false, vtxIsSilver = false;
      for (const auto& gv : genVertices_) {
        if (gv.isGold(vtx))        { vtxIsGold = true; break; }
        if (gv.isSilver(vtx))      vtxIsSilver = true;
      }

      // Collect track refs and build transient tracks for Kalman compatibility
      std::vector<reco::TrackRef> trackRefs;
      std::vector<reco::TransientTrack> ttracks;
      for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
        trackRefs.push_back(it->castTo<reco::TrackRef>());
        ttracks.push_back(ttBuilder->build(*trackRefs.back()));
      }

      // Refit vertex with Kalman fitter for compatibility estimation
      KalmanVertexFitter fitter;
      TransientVertex tv;
      bool validFit = (ttracks.size() >= 2);
      if (validFit) tv = fitter.vertex(ttracks);
      validFit = validFit && tv.isValid();

      const KalmanVertexTrackCompatibilityEstimator<5> estimator;

      for (size_t ti = 0; ti < trackRefs.size(); ++ti) {
        const reco::TrackRef& tref = trackRefs[ti];
        bool isSignal = (TrackHelper::FindTrackIndex(*tref, signalTracks_) >= 0);
        bool isRemoved = (cleanedTrackRefs.find(tref) == cleanedTrackRefs.end());

        float compat = -1.f;
        if (validFit) {
          auto res = estimator.estimate(tv, ttracks[ti]);
          if (res.first) compat = float(std::sqrt(res.second));
        }

        float cosTheta = float(TrackHelper::CalculateCosTheta(pv, vtx, *tref));

        cleanTrack_compatibility_.push_back(compat);
        cleanTrack_cosTheta_     .push_back(cosTheta);
        cleanTrack_isSignal_     .push_back(isSignal);
        cleanTrack_isRemoved_    .push_back(isRemoved);
        cleanTrack_vtxIsGold_    .push_back(vtxIsGold);
        cleanTrack_vtxIsSilver_  .push_back(vtxIsSilver);
        cleanTrack_vtxNTracks_   .push_back(unsigned(vtx.tracksSize()));
      }
    }
  }

  // ── allStageVtx: one row per vertex per stage ─────────────────────────────
  for (size_t si = 0; si < kNStages; ++si) {
    if (!stageHandles_[si].isValid()) continue;
    stageVtx_stageIdx_ = int(si);
    for (const auto& vtx : *stageHandles_[si]) {
      stageVtx_cosTheta_   = float(VertexHelper::CalculateCosTheta(pv, vtx));
      stageVtx_decayAngle_ = float(VertexHelper::CalculateDecayAngle(vtx));
      {
        stageVtx_trackCosTheta_.clear();
        stageVtx_trackCosThetaCM_.clear();
        for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
          reco::TrackRef tref = it->castTo<reco::TrackRef>();
          stageVtx_trackCosTheta_.push_back(
              float(TrackHelper::CalculateCosTheta(pv, vtx, *tref)));
          stageVtx_trackCosThetaCM_.push_back(
              float(VertexHelper::CalculateCMCosTheta(vtx, *tref)));
        }
        stageVtx_minTrackCosTheta_ = stageVtx_trackCosTheta_.empty() ? -1.0f
            : *std::min_element(stageVtx_trackCosTheta_.begin(), stageVtx_trackCosTheta_.end());
      }
      {
        auto p4  = VertexHelper::GetVertex4Vector(vtx);
        double p = p4.P(), mv = p4.M();
        double E = std::sqrt(p*p + mv*mv);
        stageVtx_mass_   = float(mv);
        stageVtx_pOverE_ = (E > 1e-6) ? float(p / E) : -1.f;
      }
      stageVtx_dxySignif_ = computeDxySignif(vtx, pv);
      stageVtx_charge_    = VertexHelper::CalculateTotalCharge(vtx);
      stageVtx_isGold_    = false;
      stageVtx_isSilver_  = false;
      stageVtx_isBronze_  = false;
      if (hasGenInfo_) {
        for (const auto& gv : genVertices_) {
          if      (gv.isGold(vtx))   { stageVtx_isGold_   = true; break; }
          else if (gv.isSilver(vtx))   stageVtx_isSilver_ = true;
          else if (gv.isBronze(vtx))   stageVtx_isBronze_ = true;
        }
      }
      allStageVtxTree_->Fill();
    }
  }

  // ── seedTracks: per-track info in 2-track OS seeded vertices ──────────────
  if (stageHandles_[kSeed].isValid()) {
    for (const auto& vtx : *stageHandles_[kSeed]) {
      if (vtx.tracksSize() != 2) continue;
      auto it0 = vtx.tracks_begin();
      reco::TrackRef tr0 = it0->castTo<reco::TrackRef>();
      auto it1 = std::next(it0);
      reco::TrackRef tr1 = it1->castTo<reco::TrackRef>();
      if (tr0->charge() * tr1->charge() > 0) continue;  // same sign: skip

      bool vtxIsGold = false;
      if (hasGenInfo_) {
        for (const auto& gv : genVertices_) {
          if (gv.isGold(vtx)) { vtxIsGold = true; break; }
        }
      }

      for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
        reco::TrackRef tref = it->castTo<reco::TrackRef>();
        seedTrack_cosTheta_ .push_back(float(TrackHelper::CalculateCosTheta(pv, vtx, *tref)));
        seedTrack_isSignal_ .push_back(TrackHelper::FindTrackIndex(*tref, signalTracks_) >= 0);
        seedTrack_vtxIsGold_.push_back(vtxIsGold);
      }
    }
  }

  // ── Fill all trees ─────────────────────────────────────────────────────────
  genFunnelTree_->Fill();
  stageCountsTree_->Fill();
  cleaningTracksTree_->Fill();
  seedTracksTree_->Fill();
}

// ─────────────────────────────────────────────────────────────────────────────
reco::GenParticleCollection HyddraSVsDiagnosticAnalyzer::getStableChargedDaughtersFromPacked(
    const GenVertex& genVertex,
    const std::vector<pat::PackedGenParticle>& packed) const {

  const reco::Candidate* genZ = genVertex.genPair().first.mother();
  reco::GenParticleCollection result;
  if (!genZ) return result;

  for (const auto& p : packed) {
    if (p.status() != 1 || p.charge() == 0) continue;
    const reco::Candidate* mom = p.mother(0);
    if (!mom) continue;

    bool isFromSameZ = false;
    const reco::Candidate* prev = nullptr;
    while (mom && mom != prev) {
      if (mom->pdgId() == 23) {
        if (mom == genZ) isFromSameZ = true;
        break;
      }
      prev = mom;
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }
    if (isFromSameZ) {
      result.emplace_back(reco::GenParticle(
          p.charge(), p.p4(), p.vertex(), p.pdgId(), p.status(), true));
    }
  }
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsDiagnosticAnalyzer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.add<bool>("hasGenInfo", true);
  desc.add<bool>("isFullAOD",  true);
  desc.add<edm::InputTag>("seedVertices",          edm::InputTag("hyddraSVsDiag", "leptonicSeeds"));
  desc.add<edm::InputTag>("mergedVertices",         edm::InputTag("hyddraSVsDiag", "leptonicMerged"));
  desc.add<edm::InputTag>("cleanedVertices",        edm::InputTag("hyddraSVsDiag", "leptonicCleaned"));
  desc.add<edm::InputTag>("disambiguatedVertices",  edm::InputTag("hyddraSVsDiag", "leptonicDisambiguated"));
  desc.add<edm::InputTag>("filteredVertices",       edm::InputTag("hyddraSVsDiag", "leptonicFiltered"));
  desc.add<edm::InputTag>("pvCollection",           edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("tracks",                 edm::InputTag("muonBestTrackProducer", "globalTracks"));
  desc.add<edm::InputTag>("genParticles",           edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("packedGenParticles",     edm::InputTag(""));

  // Optional leptonic algorithm config — when provided, parameters are stored
  // as metadata in the output ROOT file for use by the diagnostic plot scripts.
  edm::ParameterSetDescription lepDesc;
  lepDesc.add<double>("seedCosThetaCut",              0.75);
  lepDesc.add<double>("minMass",                      2.0);
  lepDesc.add<double>("minPOverE",                    0.6);
  lepDesc.add<double>("maxNormChi2",                  5.0);
  lepDesc.add<double>("minDxySignificance",           25.0);
  lepDesc.add<double>("maxCompatibility",             1.5);
  lepDesc.add<double>("minCleanCosTheta",             0.5);
  lepDesc.add<double>("maxCleanCosTheta",             1.0);
  lepDesc.add<bool>  ("invertCleanCosThetaCut",       false);
  lepDesc.add<bool>  ("useDiagonalCut",               false);
  lepDesc.add<double>("cleanCutSlope",                0.0);
  lepDesc.add<double>("minTrackCosTheta",             0.5);
  lepDesc.add<double>("maxTrackCosThetaCM_Limit",     0.95);
  lepDesc.add<double>("maxTrackCosThetaCM_Intercept", 1.8);
  lepDesc.add<double>("trackCosThetaCM_Slope",        -1.0);
  lepDesc.add<bool>  ("requireChargeNeutrality",      true);
  lepDesc.add<double>("minVtxCosTheta",               -1.0);
  lepDesc.add<double>("maxVtxCosTheta",                1.0);
  lepDesc.add<bool>  ("useAbsVtxCosTheta",            false);
  lepDesc.add<double>("maxVtxDecayAngle",              1.0);
  lepDesc.add<bool>  ("useAbsVtxDecayAngle",          false);
  lepDesc.add<bool>  ("applyVtxDecayAngleCleaning",   false);
  lepDesc.add<bool>  ("applyVtxDecayAngleFiltering",  false);
  lepDesc.add<double>("minMassFilter",                 0.0);
  lepDesc.add<double>("minBetaFilter",                 0.0);
  lepDesc.add<bool>  ("doMerging",                    true);
  lepDesc.add<bool>  ("doCleaning",                   true);
  lepDesc.add<bool>  ("doDisambiguation",             true);
  lepDesc.add<bool>  ("doFiltering",                  true);
  desc.addOptional<edm::ParameterSetDescription>("leptonic", lepDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsDiagnosticAnalyzer);
