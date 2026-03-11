// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      HyddraSVsHadronicDiagnosticAnalyzer
//
// Description: Reads the 5 intermediate hadronic vertex collections produced by
//              HyddraSVsHadronicDiagnosticProducer and gen-matches each one to
//              trace where signal vertices are lost at each pipeline stage.
//              Writes three TTrees:
//                genFunnel    — one row per gen vertex per event
//                stageCounts  — one row per event (total + matched counts)
//                allStageVtx  — one row per vertex per stage (extended for
//                               hadronic: nTracks, nSignalTracks, hasSharedTrack)
//                hadronicConfig — single-row metadata tree
//
//              Compared to HyddraSVsDiagnosticAnalyzer (leptonic):
//                - No cleaningTracks tree (hadronic cleaning is vertex-level only)
//                - No seedTracks tree
//                - allStageVtx has nTracks/nSignalTracks/hasSharedTrack instead
//                  of per-track cosTheta/cosThetaCM vectors
//                - hadronicConfig tree stores the hadronic cut thresholds
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <limits>
#include <array>
#include <map>
#include <set>

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

// ROOT
#include "TTree.h"

// ─── Stage index constants ────────────────────────────────────────────────────
static constexpr size_t kHadSeed    = 0;
static constexpr size_t kHadMerged  = 1;
static constexpr size_t kHadCleaned = 2;
static constexpr size_t kHadDisambig= 3;
static constexpr size_t kHadFiltered= 4;
static constexpr size_t kHadNStages = 5;

// ─── Helper: result of gen-matching against one stage collection ──────────────
struct HadStageMatch {
  bool isGold   = false;
  bool isSilver = false;
  bool isBronze = false;
  const reco::Vertex* bestVertex = nullptr;
};

// ─────────────────────────────────────────────────────────────────────────────
class HyddraSVsHadronicDiagnosticAnalyzer
    : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit HyddraSVsHadronicDiagnosticAnalyzer(const edm::ParameterSet&);
  ~HyddraSVsHadronicDiagnosticAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  void clearBranches();

  // Gen matching helpers
  HadStageMatch findStageMatch(const GenVertex& genVtx,
                               const reco::VertexCollection& coll) const;
  float computeMatchRatio(const reco::Vertex& vtx) const;
  int   computeNSignalTracks(const reco::Vertex& vtx) const;
  float computeDxySignif(const reco::Vertex& vtx, const reco::Vertex& pv) const;
  void  fillRecoPropsAt(size_t stageIdx, const reco::Vertex& vtx,
                        const reco::Vertex& pv);

  // Per-stage helper: returns a bool vector indicating which vertices share
  // at least one track with another vertex in the same collection.
  std::vector<bool> computeSharedTrackFlags(const reco::VertexCollection& coll) const;

  // ── Configuration ──────────────────────────────────────────────────────────
  bool hasGenInfo_;
  bool isFullAOD_;
  double genMatchDeltaRCut_;

  // ── Tokens ─────────────────────────────────────────────────────────────────
  std::array<edm::EDGetTokenT<reco::VertexCollection>, kHadNStages> stageTokens_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<reco::TrackCollection>  tracksToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> packedGenToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // ── Handles ────────────────────────────────────────────────────────────────
  std::array<edm::Handle<reco::VertexCollection>, kHadNStages> stageHandles_;
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
  TTree* allStageVtxTree_;

  // ═══════════════════════════════════════════════════════════════════════════
  // genFunnel branches (vectors — one element per gen vertex)
  // ═══════════════════════════════════════════════════════════════════════════
  unsigned int genFunnel_n_;

  std::vector<float> genFunnel_dxy_;
  std::vector<float> genFunnel_pt_;
  std::vector<float> genFunnel_eta_;
  std::vector<float> genFunnel_mass_;
  std::vector<float> genFunnel_phi_;
  std::vector<bool>  genFunnel_hasTracks_;

  std::vector<bool> genFunnel_gold_seed_,    genFunnel_gold_merged_,
                    genFunnel_gold_cleaned_,  genFunnel_gold_disambig_,
                    genFunnel_gold_filtered_;
  std::vector<bool> genFunnel_silver_seed_,  genFunnel_silver_merged_,
                    genFunnel_silver_cleaned_,genFunnel_silver_disambig_,
                    genFunnel_silver_filtered_;
  std::vector<bool> genFunnel_bronze_seed_,  genFunnel_bronze_merged_,
                    genFunnel_bronze_cleaned_,genFunnel_bronze_disambig_,
                    genFunnel_bronze_filtered_;

  std::vector<int>  genFunnel_lostGoldAtStage_;
  std::vector<int>  genFunnel_lostSilverAtStage_;

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
  // allStageVtx branches — one TTree::Fill per vertex
  // Hadronic version: replaces per-track vectors with nTracks/nSignalTracks,
  // adds normChi2 (hadronic cleaning cut) and hasSharedTrack.
  // ═══════════════════════════════════════════════════════════════════════════
  int   stageVtx_stageIdx_;
  float stageVtx_mass_;
  float stageVtx_cosTheta_;
  float stageVtx_decayAngle_;
  float stageVtx_pOverE_;
  float stageVtx_dxySignif_;
  float stageVtx_normChi2_;
  int   stageVtx_nTracks_;
  int   stageVtx_nSignalTracks_;
  bool  stageVtx_hasSharedTrack_;
  bool  stageVtx_isGold_;
  bool  stageVtx_isSilver_;
  bool  stageVtx_isBronze_;

  // ═══════════════════════════════════════════════════════════════════════════
  // hadronicConfig branches (scalars — single Fill in beginJob)
  // ═══════════════════════════════════════════════════════════════════════════
  TTree* hadronicConfigTree_;
  int    cfg_minSize_;
  double cfg_minCosTheta_;
  double cfg_maxDecayAngle_;
  double cfg_minDxySignificance_;
  double cfg_minMass_;
  double cfg_minPOverE_;
  double cfg_maxNormChi2_;
  double cfg_seedCosThetaCut_;
};

// ─────────────────────────────────────────────────────────────────────────────
HyddraSVsHadronicDiagnosticAnalyzer::HyddraSVsHadronicDiagnosticAnalyzer(
    const edm::ParameterSet& iConfig) :
  hasGenInfo_(iConfig.getParameter<bool>("hasGenInfo")),
  isFullAOD_( iConfig.getParameter<bool>("isFullAOD")),
  genMatchDeltaRCut_(iConfig.getParameter<double>("genMatchDeltaRCut")),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder")))
{
  usesResource("TFileService");

  static const std::array<std::string, kHadNStages> kTagNames = {{
    "seedVertices", "mergedVertices", "cleanedVertices",
    "disambiguatedVertices", "filteredVertices"
  }};
  for (size_t i = 0; i < kHadNStages; ++i) {
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

  // Hadronic algorithm parameters — persisted as metadata in the output file.
  cfg_minSize_           = 5;
  cfg_minCosTheta_       = 0.0;
  cfg_maxDecayAngle_     = 0.9;
  cfg_minDxySignificance_= 40.0;
  cfg_minMass_           = 2.0;
  cfg_minPOverE_         = 0.6;
  cfg_maxNormChi2_       = 5.0;
  cfg_seedCosThetaCut_   = 0.0;
  if (iConfig.exists("hadronic")) {
    const edm::ParameterSet& had = iConfig.getParameterSet("hadronic");
    cfg_minSize_            = had.getParameter<int>   ("minSize");
    cfg_minCosTheta_        = had.getParameter<double>("minCosTheta");
    cfg_maxDecayAngle_      = had.getParameter<double>("maxDecayAngle");
    cfg_minDxySignificance_ = had.getParameter<double>("minDxySignificance");
    cfg_minMass_            = had.getParameter<double>("minMass");
    cfg_minPOverE_          = had.getParameter<double>("minPOverE");
    cfg_maxNormChi2_        = had.getParameter<double>("maxNormChi2");
    cfg_seedCosThetaCut_    = had.getParameter<double>("seedCosThetaCut");
  }
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsHadronicDiagnosticAnalyzer::beginJob() {
  edm::Service<TFileService> fs;

  // ── genFunnel ──────────────────────────────────────────────────────────────
  genFunnelTree_ = fs->make<TTree>("genFunnel", "Per-gen-vertex pipeline funnel (hadronic)");

  genFunnelTree_->Branch("GenFunnel_n",        &genFunnel_n_);
  genFunnelTree_->Branch("GenFunnel_dxy",       &genFunnel_dxy_);
  genFunnelTree_->Branch("GenFunnel_pt",        &genFunnel_pt_);
  genFunnelTree_->Branch("GenFunnel_eta",       &genFunnel_eta_);
  genFunnelTree_->Branch("GenFunnel_mass",      &genFunnel_mass_);
  genFunnelTree_->Branch("GenFunnel_phi",       &genFunnel_phi_);
  genFunnelTree_->Branch("GenFunnel_hasTracks", &genFunnel_hasTracks_);

  genFunnelTree_->Branch("GenFunnel_gold_seed",       &genFunnel_gold_seed_);
  genFunnelTree_->Branch("GenFunnel_gold_merged",     &genFunnel_gold_merged_);
  genFunnelTree_->Branch("GenFunnel_gold_cleaned",    &genFunnel_gold_cleaned_);
  genFunnelTree_->Branch("GenFunnel_gold_disambig",   &genFunnel_gold_disambig_);
  genFunnelTree_->Branch("GenFunnel_gold_filtered",   &genFunnel_gold_filtered_);

  genFunnelTree_->Branch("GenFunnel_silver_seed",     &genFunnel_silver_seed_);
  genFunnelTree_->Branch("GenFunnel_silver_merged",   &genFunnel_silver_merged_);
  genFunnelTree_->Branch("GenFunnel_silver_cleaned",  &genFunnel_silver_cleaned_);
  genFunnelTree_->Branch("GenFunnel_silver_disambig", &genFunnel_silver_disambig_);
  genFunnelTree_->Branch("GenFunnel_silver_filtered", &genFunnel_silver_filtered_);

  genFunnelTree_->Branch("GenFunnel_bronze_seed",     &genFunnel_bronze_seed_);
  genFunnelTree_->Branch("GenFunnel_bronze_merged",   &genFunnel_bronze_merged_);
  genFunnelTree_->Branch("GenFunnel_bronze_cleaned",  &genFunnel_bronze_cleaned_);
  genFunnelTree_->Branch("GenFunnel_bronze_disambig", &genFunnel_bronze_disambig_);
  genFunnelTree_->Branch("GenFunnel_bronze_filtered", &genFunnel_bronze_filtered_);

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
  stageCountsTree_ = fs->make<TTree>("stageCounts", "Per-event vertex counts per stage (hadronic)");

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

  // ── allStageVtx (hadronic) ─────────────────────────────────────────────────
  allStageVtxTree_ = fs->make<TTree>("allStageVtx",
      "All hadronic vertices at each stage with match quality and track purity");
  allStageVtxTree_->Branch("StageVtx_stageIdx",        &stageVtx_stageIdx_);
  allStageVtxTree_->Branch("StageVtx_mass",             &stageVtx_mass_);
  allStageVtxTree_->Branch("StageVtx_cosTheta",         &stageVtx_cosTheta_);
  allStageVtxTree_->Branch("StageVtx_decayAngle",       &stageVtx_decayAngle_);
  allStageVtxTree_->Branch("StageVtx_pOverE",           &stageVtx_pOverE_);
  allStageVtxTree_->Branch("StageVtx_dxySignif",        &stageVtx_dxySignif_);
  allStageVtxTree_->Branch("StageVtx_normChi2",         &stageVtx_normChi2_);
  allStageVtxTree_->Branch("StageVtx_nTracks",          &stageVtx_nTracks_);
  allStageVtxTree_->Branch("StageVtx_nSignalTracks",    &stageVtx_nSignalTracks_);
  allStageVtxTree_->Branch("StageVtx_hasSharedTrack",   &stageVtx_hasSharedTrack_);
  allStageVtxTree_->Branch("StageVtx_isGold",           &stageVtx_isGold_);
  allStageVtxTree_->Branch("StageVtx_isSilver",         &stageVtx_isSilver_);
  allStageVtxTree_->Branch("StageVtx_isBronze",         &stageVtx_isBronze_);

  // ── hadronicConfig ─────────────────────────────────────────────────────────
  hadronicConfigTree_ = fs->make<TTree>("hadronicConfig",
      "Hadronic HYDDRA configuration parameters (single row)");
  hadronicConfigTree_->Branch("minSize",           &cfg_minSize_);
  hadronicConfigTree_->Branch("minCosTheta",       &cfg_minCosTheta_);
  hadronicConfigTree_->Branch("maxDecayAngle",     &cfg_maxDecayAngle_);
  hadronicConfigTree_->Branch("minDxySignificance",&cfg_minDxySignificance_);
  hadronicConfigTree_->Branch("minMass",           &cfg_minMass_);
  hadronicConfigTree_->Branch("minPOverE",         &cfg_minPOverE_);
  hadronicConfigTree_->Branch("maxNormChi2",       &cfg_maxNormChi2_);
  hadronicConfigTree_->Branch("seedCosThetaCut",   &cfg_seedCosThetaCut_);
  hadronicConfigTree_->Fill();
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsHadronicDiagnosticAnalyzer::clearBranches() {
  genFunnel_n_ = 0;
  genFunnel_dxy_.clear();  genFunnel_pt_.clear();   genFunnel_eta_.clear();
  genFunnel_mass_.clear(); genFunnel_phi_.clear();  genFunnel_hasTracks_.clear();

  genFunnel_gold_seed_.clear();    genFunnel_gold_merged_.clear();
  genFunnel_gold_cleaned_.clear(); genFunnel_gold_disambig_.clear();
  genFunnel_gold_filtered_.clear();
  genFunnel_silver_seed_.clear();  genFunnel_silver_merged_.clear();
  genFunnel_silver_cleaned_.clear();genFunnel_silver_disambig_.clear();
  genFunnel_silver_filtered_.clear();
  genFunnel_bronze_seed_.clear();  genFunnel_bronze_merged_.clear();
  genFunnel_bronze_cleaned_.clear();genFunnel_bronze_disambig_.clear();
  genFunnel_bronze_filtered_.clear();

  genFunnel_lostGoldAtStage_.clear();
  genFunnel_lostSilverAtStage_.clear();

  genFunnel_mass_seed_.clear();     genFunnel_mass_merged_.clear();
  genFunnel_mass_cleaned_.clear();  genFunnel_mass_disambig_.clear();
  genFunnel_mass_filtered_.clear();
  genFunnel_dxySignif_seed_.clear(); genFunnel_dxySignif_merged_.clear();
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
  genFunnel_matchRatio_seed_.clear(); genFunnel_matchRatio_merged_.clear();
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
}

// ─────────────────────────────────────────────────────────────────────────────
HadStageMatch HyddraSVsHadronicDiagnosticAnalyzer::findStageMatch(
    const GenVertex& genVtx,
    const reco::VertexCollection& coll) const {

  HadStageMatch result;
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
float HyddraSVsHadronicDiagnosticAnalyzer::computeMatchRatio(
    const reco::Vertex& vtx) const {
  if (vtx.tracksSize() == 0) return 0.f;
  int nSignal = 0;
  for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
    reco::Track t(**it);
    if (TrackHelper::FindTrackIndex(t, signalTracks_) >= 0) nSignal++;
  }
  return float(nSignal) / float(vtx.tracksSize());
}

// ─────────────────────────────────────────────────────────────────────────────
int HyddraSVsHadronicDiagnosticAnalyzer::computeNSignalTracks(
    const reco::Vertex& vtx) const {
  int nSignal = 0;
  for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
    reco::Track t(**it);
    if (TrackHelper::FindTrackIndex(t, signalTracks_) >= 0) nSignal++;
  }
  return nSignal;
}

// ─────────────────────────────────────────────────────────────────────────────
float HyddraSVsHadronicDiagnosticAnalyzer::computeDxySignif(
    const reco::Vertex& vtx, const reco::Vertex& pv) const {
  float dxy    = float(VertexHelper::CalculateDxy(vtx, pv));
  float dxyErr = float(VertexHelper::CalculateDxyError(vtx, pv));
  return (dxyErr > 0) ? dxy / dxyErr : -1.f;
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsHadronicDiagnosticAnalyzer::fillRecoPropsAt(
    size_t s, const reco::Vertex& vtx, const reco::Vertex& pv) {

  std::array<std::vector<float>*, kHadNStages> massVecs = {{
      &genFunnel_mass_seed_,  &genFunnel_mass_merged_,
      &genFunnel_mass_cleaned_,&genFunnel_mass_disambig_,&genFunnel_mass_filtered_}};
  std::array<std::vector<float>*, kHadNStages> dxsVecs = {{
      &genFunnel_dxySignif_seed_,  &genFunnel_dxySignif_merged_,
      &genFunnel_dxySignif_cleaned_,&genFunnel_dxySignif_disambig_,&genFunnel_dxySignif_filtered_}};
  std::array<std::vector<float>*, kHadNStages> chi2Vecs = {{
      &genFunnel_normChi2_seed_,  &genFunnel_normChi2_merged_,
      &genFunnel_normChi2_cleaned_,&genFunnel_normChi2_disambig_,&genFunnel_normChi2_filtered_}};
  std::array<std::vector<int>*, kHadNStages> ntrkVecs = {{
      &genFunnel_nTracks_seed_,  &genFunnel_nTracks_merged_,
      &genFunnel_nTracks_cleaned_,&genFunnel_nTracks_disambig_,&genFunnel_nTracks_filtered_}};
  std::array<std::vector<float>*, kHadNStages> ctVecs = {{
      &genFunnel_cosTheta_seed_,  &genFunnel_cosTheta_merged_,
      &genFunnel_cosTheta_cleaned_,&genFunnel_cosTheta_disambig_,&genFunnel_cosTheta_filtered_}};
  std::array<std::vector<float>*, kHadNStages> mrVecs = {{
      &genFunnel_matchRatio_seed_,  &genFunnel_matchRatio_merged_,
      &genFunnel_matchRatio_cleaned_,&genFunnel_matchRatio_disambig_,&genFunnel_matchRatio_filtered_}};
  std::array<std::vector<float>*, kHadNStages> daVecs = {{
      &genFunnel_decayAngle_seed_,  &genFunnel_decayAngle_merged_,
      &genFunnel_decayAngle_cleaned_,&genFunnel_decayAngle_disambig_,&genFunnel_decayAngle_filtered_}};
  std::array<std::vector<float>*, kHadNStages> poeVecs = {{
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
std::vector<bool> HyddraSVsHadronicDiagnosticAnalyzer::computeSharedTrackFlags(
    const reco::VertexCollection& coll) const {

  std::vector<bool> flags(coll.size(), false);

  // Build track → vertex index map
  std::map<reco::TrackRef, std::vector<size_t>> trackToVtx;
  for (size_t i = 0; i < coll.size(); ++i) {
    for (auto it = coll[i].tracks_begin(); it != coll[i].tracks_end(); ++it) {
      trackToVtx[it->castTo<reco::TrackRef>()].push_back(i);
    }
  }

  // Any track appearing in more than one vertex → both vertices are shared
  for (const auto& [tref, indices] : trackToVtx) {
    if (indices.size() > 1) {
      for (size_t idx : indices) {
        flags[idx] = true;
      }
    }
  }

  return flags;
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsHadronicDiagnosticAnalyzer::analyze(
    const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  clearBranches();

  for (size_t i = 0; i < kHadNStages; ++i)
    iEvent.getByToken(stageTokens_[i], stageHandles_[i]);

  iEvent.getByToken(pvToken_, pvHandle_);
  iEvent.getByToken(tracksToken_, tracksHandle_);

  const reco::Vertex& pv = pvHandle_->at(0);
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);

  // ── stageCounts totals ─────────────────────────────────────────────────────
  auto safeSize = [&](size_t i) {
    return stageHandles_[i].isValid() ? unsigned(stageHandles_[i]->size()) : 0u;
  };
  stage_n_seed_     = safeSize(kHadSeed);
  stage_n_merged_   = safeSize(kHadMerged);
  stage_n_cleaned_  = safeSize(kHadCleaned);
  stage_n_disambig_ = safeSize(kHadDisambig);
  stage_n_filtered_ = safeSize(kHadFiltered);

  // ── Build gen vertices and signal tracks ───────────────────────────────────
  genVertices_.clear();
  signalTracks_.clear();

  if (hasGenInfo_) {
    iEvent.getByToken(genToken_, genHandle_);
    if (!isFullAOD_)
      iEvent.getByToken(packedGenToken_, packedGenHandle_);

    std::vector<reco::TransientTrack> ttracks;
    for (const auto& track : *tracksHandle_)
      ttracks.emplace_back(ttBuilder->build(track));

    GenVertices allSignalSVs(*genHandle_);
    DeltaRGenMatchHungarian<reco::TransientTrack> assigner(
        ttracks, allSignalSVs.getAllGenParticles());

    genVertices_ = GenVertices(assigner.GetPairedObjects().ConvertFromTTracks(), genMatchDeltaRCut_);
    allSignalSVs += genVertices_;
    genVertices_ = allSignalSVs;

    for (const auto& gv : genVertices_) {
      if (!gv.hasTracks()) continue;
      for (const auto& pair : gv.genMatches())
        signalTracks_.emplace_back(pair.GetObjectA());
    }
  }

  // ── stageCounts gold/silver/bronze per vertex ──────────────────────────────
  if (hasGenInfo_) {
    auto countMatches = [&](size_t si,
                            unsigned& nG, unsigned& nS, unsigned& nB) {
      if (!stageHandles_[si].isValid()) return;
      for (const auto& vtx : *stageHandles_[si]) {
        bool anyGold = false, anySilver = false, anyBronze = false;
        for (const auto& gv : genVertices_) {
          if (gv.isGold(vtx))        { anyGold = true; break; }
          if (gv.isSilver(vtx))      anySilver = true;
          else if (gv.isBronze(vtx)) anyBronze = true;
        }
        if (anyGold)        nG++;
        else if (anySilver) nS++;
        else if (anyBronze) nB++;
      }
    };

    countMatches(kHadSeed,    stage_nGold_seed_,    stage_nSilver_seed_,    stage_nBronze_seed_);
    countMatches(kHadMerged,  stage_nGold_merged_,  stage_nSilver_merged_,  stage_nBronze_merged_);
    countMatches(kHadCleaned, stage_nGold_cleaned_,  stage_nSilver_cleaned_, stage_nBronze_cleaned_);
    countMatches(kHadDisambig,stage_nGold_disambig_, stage_nSilver_disambig_,stage_nBronze_disambig_);
    countMatches(kHadFiltered,stage_nGold_filtered_, stage_nSilver_filtered_,stage_nBronze_filtered_);
  }

  // ── genFunnel: one entry per gen vertex ────────────────────────────────────
  if (hasGenInfo_) {
    for (const auto& gv : genVertices_) {
      genFunnel_dxy_       .push_back(float(gv.dxy()));
      genFunnel_pt_        .push_back(float(gv.pt()));
      genFunnel_eta_       .push_back(float(gv.eta()));
      genFunnel_mass_      .push_back(float(gv.mass()));
      genFunnel_phi_       .push_back(float(gv.phi()));
      genFunnel_hasTracks_ .push_back(bool(gv.hasTracks()));

      std::array<HadStageMatch, kHadNStages> m;
      for (size_t si = 0; si < kHadNStages; ++si) {
        if (stageHandles_[si].isValid())
          m[si] = findStageMatch(gv, *stageHandles_[si]);
      }

      genFunnel_gold_seed_    .push_back(m[kHadSeed].isGold);
      genFunnel_gold_merged_  .push_back(m[kHadMerged].isGold);
      genFunnel_gold_cleaned_ .push_back(m[kHadCleaned].isGold);
      genFunnel_gold_disambig_.push_back(m[kHadDisambig].isGold);
      genFunnel_gold_filtered_.push_back(m[kHadFiltered].isGold);

      genFunnel_silver_seed_    .push_back(m[kHadSeed].isSilver);
      genFunnel_silver_merged_  .push_back(m[kHadMerged].isSilver);
      genFunnel_silver_cleaned_ .push_back(m[kHadCleaned].isSilver);
      genFunnel_silver_disambig_.push_back(m[kHadDisambig].isSilver);
      genFunnel_silver_filtered_.push_back(m[kHadFiltered].isSilver);

      genFunnel_bronze_seed_    .push_back(m[kHadSeed].isBronze);
      genFunnel_bronze_merged_  .push_back(m[kHadMerged].isBronze);
      genFunnel_bronze_cleaned_ .push_back(m[kHadCleaned].isBronze);
      genFunnel_bronze_disambig_.push_back(m[kHadDisambig].isBronze);
      genFunnel_bronze_filtered_.push_back(m[kHadFiltered].isBronze);

      auto computeLostStage = [&](bool useGold) -> int {
        auto hasMatch = [&](size_t si) {
          return useGold ? m[si].isGold : (m[si].isGold || m[si].isSilver);
        };
        if (hasMatch(kHadFiltered)) return -1;
        for (int si = int(kHadNStages) - 2; si >= 0; --si) {
          if (hasMatch(size_t(si))) return si + 1;
        }
        return 0;
      };
      genFunnel_lostGoldAtStage_  .push_back(computeLostStage(true));
      genFunnel_lostSilverAtStage_.push_back(computeLostStage(false));

      // Reco props per stage (same sentinel pattern as leptonic analyzer)
      std::array<std::vector<float>*, kHadNStages> massVecs = {{
          &genFunnel_mass_seed_,  &genFunnel_mass_merged_,
          &genFunnel_mass_cleaned_,&genFunnel_mass_disambig_,&genFunnel_mass_filtered_}};
      std::array<std::vector<float>*, kHadNStages> dxsVecs = {{
          &genFunnel_dxySignif_seed_,  &genFunnel_dxySignif_merged_,
          &genFunnel_dxySignif_cleaned_,&genFunnel_dxySignif_disambig_,&genFunnel_dxySignif_filtered_}};
      std::array<std::vector<float>*, kHadNStages> chi2Vecs = {{
          &genFunnel_normChi2_seed_,  &genFunnel_normChi2_merged_,
          &genFunnel_normChi2_cleaned_,&genFunnel_normChi2_disambig_,&genFunnel_normChi2_filtered_}};
      std::array<std::vector<int>*, kHadNStages> ntrkVecs = {{
          &genFunnel_nTracks_seed_,  &genFunnel_nTracks_merged_,
          &genFunnel_nTracks_cleaned_,&genFunnel_nTracks_disambig_,&genFunnel_nTracks_filtered_}};
      std::array<std::vector<float>*, kHadNStages> ctVecs = {{
          &genFunnel_cosTheta_seed_,  &genFunnel_cosTheta_merged_,
          &genFunnel_cosTheta_cleaned_,&genFunnel_cosTheta_disambig_,&genFunnel_cosTheta_filtered_}};
      std::array<std::vector<float>*, kHadNStages> mrVecs = {{
          &genFunnel_matchRatio_seed_,  &genFunnel_matchRatio_merged_,
          &genFunnel_matchRatio_cleaned_,&genFunnel_matchRatio_disambig_,&genFunnel_matchRatio_filtered_}};
      std::array<std::vector<float>*, kHadNStages> daVecs = {{
          &genFunnel_decayAngle_seed_,  &genFunnel_decayAngle_merged_,
          &genFunnel_decayAngle_cleaned_,&genFunnel_decayAngle_disambig_,&genFunnel_decayAngle_filtered_}};
      std::array<std::vector<float>*, kHadNStages> poeVecs = {{
          &genFunnel_pOverE_seed_,  &genFunnel_pOverE_merged_,
          &genFunnel_pOverE_cleaned_,&genFunnel_pOverE_disambig_,&genFunnel_pOverE_filtered_}};

      for (size_t si = 0; si < kHadNStages; ++si) {
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

  // ── allStageVtx: one row per vertex per stage ─────────────────────────────
  for (size_t si = 0; si < kHadNStages; ++si) {
    if (!stageHandles_[si].isValid()) continue;

    const reco::VertexCollection& coll = *stageHandles_[si];
    std::vector<bool> sharedFlags = computeSharedTrackFlags(coll);

    stageVtx_stageIdx_ = int(si);
    for (size_t vi = 0; vi < coll.size(); ++vi) {
      const reco::Vertex& vtx = coll[vi];

      auto p4  = VertexHelper::GetVertex4Vector(vtx);
      double p = p4.P(), mv = p4.M();
      double E = std::sqrt(p*p + mv*mv);

      stageVtx_mass_          = float(mv);
      stageVtx_cosTheta_      = float(VertexHelper::CalculateCosTheta(pv, vtx));
      stageVtx_decayAngle_    = float(VertexHelper::CalculateDecayAngle(vtx));
      stageVtx_pOverE_        = (E > 1e-6) ? float(p / E) : -1.f;
      stageVtx_dxySignif_     = computeDxySignif(vtx, pv);
      stageVtx_normChi2_      = float(vtx.normalizedChi2());
      stageVtx_nTracks_       = int(vtx.tracksSize());
      stageVtx_nSignalTracks_ = hasGenInfo_ ? computeNSignalTracks(vtx) : -1;
      stageVtx_hasSharedTrack_= sharedFlags[vi];

      stageVtx_isGold_   = false;
      stageVtx_isSilver_ = false;
      stageVtx_isBronze_ = false;
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

  // ── Fill event-level trees ─────────────────────────────────────────────────
  genFunnelTree_->Fill();
  stageCountsTree_->Fill();
}

// ─────────────────────────────────────────────────────────────────────────────
void HyddraSVsHadronicDiagnosticAnalyzer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.add<bool>("hasGenInfo", true);
  desc.add<bool>("isFullAOD",  true);
  desc.add<double>("genMatchDeltaRCut", 0.02);
  desc.add<edm::InputTag>("seedVertices",
      edm::InputTag("hyddraSVsHadronicDiag", "hadronicSeeds"));
  desc.add<edm::InputTag>("mergedVertices",
      edm::InputTag("hyddraSVsHadronicDiag", "hadronicMerged"));
  desc.add<edm::InputTag>("cleanedVertices",
      edm::InputTag("hyddraSVsHadronicDiag", "hadronicCleaned"));
  desc.add<edm::InputTag>("disambiguatedVertices",
      edm::InputTag("hyddraSVsHadronicDiag", "hadronicDisambiguated"));
  desc.add<edm::InputTag>("filteredVertices",
      edm::InputTag("hyddraSVsHadronicDiag", "hadronicFiltered"));
  desc.add<edm::InputTag>("pvCollection",   edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("tracks",         edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("genParticles",   edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("packedGenParticles", edm::InputTag(""));

  edm::ParameterSetDescription hadDesc;
  hadDesc.add<double>("seedCosThetaCut",   0.0);
  hadDesc.add<double>("minMass",           2.0);
  hadDesc.add<double>("minPOverE",         0.6);
  hadDesc.add<double>("maxNormChi2",       5.0);
  hadDesc.add<double>("minDxySignificance",40.0);
  hadDesc.add<int>   ("minSize",           5);
  hadDesc.add<double>("minCosTheta",       0.0);
  hadDesc.add<double>("maxDecayAngle",     0.9);
  hadDesc.add<bool>  ("doMerging",         true);
  hadDesc.add<bool>  ("doCleaning",        true);
  hadDesc.add<bool>  ("doDisambiguation",  true);
  hadDesc.add<bool>  ("doFiltering",       true);
  hadDesc.add<bool>  ("useVertexSmoothing",false);
  desc.addOptional<edm::ParameterSetDescription>("hadronic", hadDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsHadronicDiagnosticAnalyzer);
