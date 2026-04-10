// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      HyddraSVsEXOAnalyzer
//
// Description: EDAnalyzer for the EXONanoAOD HYDDRA prototype format.
//
//   Reads inclusiveVertices (Tier 0) and isolatedVertices (Tier 1) produced
//   by HyddraSVsEXOProducer and writes a flat TTree with the agreed branch
//   structure. All kinematic variables are stored without signal-dependent cuts.
//
//   Gen matching (MC only) implements a hybrid approach:
//     - Gold:   both tracks match daughters of the same signal gen vertex (ΔR < genDRCut)
//     - Bronze: at least one track matches a signal gen daughter
//     - Per-leg ΔR and relPtDiff stored without cuts (LLPNanoAOD-style)
//     - Spatial: min3D and nearestGenVtxIdx (HYDDRA-style)
//   A separate HyddraGenSV collection stores truth-level vertex kinematics
//   including passSelection (both gen daughters have a reco match within tighter ΔR).
//
// Original Author:  Andres Abreu
//

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

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
#include "DataFormats/PatCandidates/interface/MET.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// ROOT
#include "TTree.h"

// ---------------------------------------------------------------------------
// Helpers (file-local)
// ---------------------------------------------------------------------------
namespace {

constexpr float INV = -999.f;

float deltaPhi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  while (dphi >  M_PI) dphi -= 2.f * M_PI;
  while (dphi < -M_PI) dphi += 2.f * M_PI;
  return dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2) {
  float de = eta1 - eta2, dp = deltaPhi(phi1, phi2);
  return std::sqrt(de * de + dp * dp);
}

// Pointing angle: cos(angle between PV→SV displacement and SV momentum)
float computeCosTheta(const reco::Vertex& sv, const reco::Vertex& pv,
                      float px, float py, float pz) {
  float dx = sv.x() - pv.x(), dy = sv.y() - pv.y(), dz = sv.z() - pv.z();
  float disp = std::sqrt(dx * dx + dy * dy + dz * dz);
  float p    = std::sqrt(px * px + py * py + pz * pz);
  if (disp < 1e-9f || p < 1e-9f) return INV;
  return (dx * px + dy * py + dz * pz) / (disp * p);
}

// CM-frame cosine of a single track wrt the boost direction.
// Track is treated as massless. Returns INV if boost is not well-defined.
float computeTrackCosThetaCM(float trk_px, float trk_py, float trk_pz,
                              float sys_px, float sys_py, float sys_pz, float sys_e) {
  if (sys_e < 1e-9f) return INV;
  float bx = sys_px / sys_e, by = sys_py / sys_e, bz = sys_pz / sys_e;
  float b2 = bx * bx + by * by + bz * bz;
  if (b2 >= 1.f || b2 < 1e-12f) return INV;
  float gamma   = 1.f / std::sqrt(1.f - b2);
  float trk_e   = std::sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);
  float bdotp   = bx * trk_px + by * trk_py + bz * trk_pz;
  float fac     = (gamma - 1.f) * bdotp / b2 - gamma * trk_e;
  float bp_x = trk_px + fac * bx;
  float bp_y = trk_py + fac * by;
  float bp_z = trk_pz + fac * bz;
  float bp_mag = std::sqrt(bp_x * bp_x + bp_y * bp_y + bp_z * bp_z);
  float b_mag  = std::sqrt(b2);
  if (bp_mag < 1e-9f) return INV;
  return (bp_x * bx + bp_y * by + bp_z * bz) / (bp_mag * b_mag);
}

// CM-frame decay angle of the system using the negatively-charged track.
// Matches VertexHelper::CalculateDecayAngle / processLLPNanoSV.py convention.
float computeDecayAngle(float px1, float py1, float pz1, int q1,
                         float px2, float py2, float pz2, int q2) {
  float sys_px = px1 + px2, sys_py = py1 + py2, sys_pz = pz1 + pz2;
  float sys_e  = std::sqrt(px1*px1+py1*py1+pz1*pz1)
               + std::sqrt(px2*px2+py2*py2+pz2*pz2);
  // Use negatively charged track; fall back to track 1 if charges equal
  bool useTrack2 = (q2 < q1);
  float tpx = useTrack2 ? px2 : px1;
  float tpy = useTrack2 ? py2 : py1;
  float tpz = useTrack2 ? pz2 : pz1;
  return computeTrackCosThetaCM(tpx, tpy, tpz, sys_px, sys_py, sys_pz, sys_e);
}

// Transverse displacement of SV from PV and its uncertainty (SV covariance only).
float computeDxy(const reco::Vertex& sv, const reco::Vertex& pv) {
  float dx = sv.x() - pv.x(), dy = sv.y() - pv.y();
  return std::sqrt(dx * dx + dy * dy);
}

float computeDxyErr(const reco::Vertex& sv, const reco::Vertex& pv) {
  float dx = sv.x() - pv.x(), dy = sv.y() - pv.y();
  float dxy = std::sqrt(dx * dx + dy * dy);
  if (dxy < 1e-9f) return 0.f;
  float ux = dx / dxy, uy = dy / dxy;
  // SV covariance contribution (indices 0=x, 1=y, 2=z)
  float sv_err2 = sv.covariance(0, 0) * ux * ux
                + sv.covariance(1, 1) * uy * uy
                + 2.f * sv.covariance(0, 1) * ux * uy;
  // PV covariance contribution
  float pv_err2 = pv.covariance(0, 0) * ux * ux
                + pv.covariance(1, 1) * uy * uy
                + 2.f * pv.covariance(0, 1) * ux * uy;
  return std::sqrt(std::max(0.f, sv_err2 + pv_err2));
}

} // namespace

// ---------------------------------------------------------------------------
// Gen vertex helper struct
// ---------------------------------------------------------------------------
struct HyddraGenVertex {
  float x, y, z, dxy;
  float pt, eta, phi, mass;
  float trk1Pt, trk1Eta, trk1Phi; int trk1PdgId;
  float trk2Pt, trk2Eta, trk2Phi; int trk2PdgId;
  int   motherPdgId;
  std::vector<size_t> dauIdx;  // indices of daughters in GenParticle collection
  bool passSelection = false;
  bool isReconstructed = false;
};

// ---------------------------------------------------------------------------
// Analyzer class
// ---------------------------------------------------------------------------
class HyddraSVsEXOAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit HyddraSVsEXOAnalyzer(const edm::ParameterSet&);
  ~HyddraSVsEXOAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}
  void clearBranches();

  // Gen matching helpers
  std::vector<HyddraGenVertex> buildGenVertices(
      const reco::GenParticleCollection& genParts) const;
  // Returns (genVtxIdx, deltaR, relPtDiff) for the best match within genDRCut_.
  // Returns (-1, INV, INV) if no match found.
  std::tuple<int, float, float> matchTrackToGen(
      const reco::Track& trk,
      const std::vector<HyddraGenVertex>& genVtxs,
      const reco::GenParticleCollection& genParts) const;
  // Check if any reco track in the event matches a gen daughter within passSelDR_.
  bool hasSoftMatch(size_t dauIdx,
                    const reco::GenParticleCollection& genParts,
                    const reco::TrackCollection& tracks) const;

  // Tokens
  edm::EDGetTokenT<reco::VertexCollection>   inclusiveToken_;
  edm::EDGetTokenT<std::vector<int>>         isoFlagsToken_;
  edm::EDGetTokenT<reco::VertexCollection>   pvToken_;
  edm::EDGetTokenT<reco::TrackCollection>    tracksToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<pat::METCollection>          metToken_;

  // Config
  bool   hasGenInfo_;
  int    motherPdgId_;
  float  genDRCut_;
  float  passSelDRCut_;

  // Output tree
  TTree* tree_ = nullptr;

  // ── Event scalars ──────────────────────────────────────────────────────────
  unsigned int run_, lumi_;
  unsigned long long event_;
  int   nRecoElectrons_;
  // ── MET ────────────────────────────────────────────────────────────────────
  float event_MET_;

  // ── HyddraSV branches (per SV) ────────────────────────────────────────────
  int nHyddraSV_;
  // Vertex position
  std::vector<float> sv_x_, sv_y_, sv_z_;
  std::vector<float> sv_xErr_, sv_yErr_, sv_zErr_;
  // Displacement (from PV — switch to BS in production)
  std::vector<float> sv_dxy_, sv_dxyErr_, sv_dxySig_;
  // Fit quality
  std::vector<float> sv_chi2_, sv_ndof_, sv_normChi2_;
  // 4-momentum (massless track sum)
  std::vector<float> sv_pt_, sv_eta_, sv_phi_, sv_mass_, sv_p_;
  // Topology
  std::vector<int>   sv_charge_;
  std::vector<float> sv_cosTheta_;    // pointing angle wrt PV
  std::vector<float> sv_decayAngle_;  // CM-frame decay angle (neg. track convention)
  std::vector<float> sv_dR_;          // ΔR between the two tracks
  std::vector<float> sv_beta_;        // p / sqrt(p² + m²) = p/E
  // Working point flag
  std::vector<bool>  sv_passesIsolation_;
  // Track 1 (leading pT)
  std::vector<float> sv_trk1Pt_, sv_trk1Eta_, sv_trk1Phi_;
  std::vector<int>   sv_trk1Charge_;
  std::vector<float> sv_trk1Dxy_, sv_trk1DxyErr_, sv_trk1DxySig_;
  std::vector<float> sv_trk1Dz_, sv_trk1DzErr_;
  std::vector<float> sv_trk1NormChi2_;
  std::vector<float> sv_trk1CosTheta_;    // alignment with vertex momentum
  std::vector<float> sv_trk1CosThetaCM_;  // CM-frame angle (cleaning variable)
  // Track 2 (subleading pT)
  std::vector<float> sv_trk2Pt_, sv_trk2Eta_, sv_trk2Phi_;
  std::vector<int>   sv_trk2Charge_;
  std::vector<float> sv_trk2Dxy_, sv_trk2DxyErr_, sv_trk2DxySig_;
  std::vector<float> sv_trk2Dz_, sv_trk2DzErr_;
  std::vector<float> sv_trk2NormChi2_;
  std::vector<float> sv_trk2CosTheta_;
  std::vector<float> sv_trk2CosThetaCM_;
  // ── Gen match branches (MC only, on HyddraSV) ─────────────────────────────
  std::vector<float> sv_trk1GenDR_, sv_trk2GenDR_;
  std::vector<float> sv_trk1GenRelPtDiff_, sv_trk2GenRelPtDiff_;
  std::vector<int>   sv_genVtxIdx_, sv_nearestGenVtxIdx_;
  std::vector<float> sv_min3D_;
  std::vector<bool>  sv_isGold_, sv_isBronze_;

  // ── HyddraGenSV branches (per gen vertex) ─────────────────────────────────
  int nHyddraGenSV_;
  std::vector<float> gv_x_, gv_y_, gv_z_, gv_dxy_;
  std::vector<float> gv_pt_, gv_eta_, gv_phi_, gv_mass_;
  std::vector<float> gv_trk1Pt_, gv_trk1Eta_, gv_trk1Phi_;
  std::vector<int>   gv_trk1PdgId_;
  std::vector<float> gv_trk2Pt_, gv_trk2Eta_, gv_trk2Phi_;
  std::vector<int>   gv_trk2PdgId_;
  std::vector<int>   gv_motherPdgId_;
  std::vector<bool>  gv_passSelection_;
  std::vector<bool>  gv_isReconstructed_;
};

// ---------------------------------------------------------------------------
HyddraSVsEXOAnalyzer::HyddraSVsEXOAnalyzer(const edm::ParameterSet& iConfig)
    : inclusiveToken_(consumes<reco::VertexCollection>(
          iConfig.getParameter<edm::InputTag>("inclusiveVertices")))
    , isoFlagsToken_(consumes<std::vector<int>>(
          iConfig.getParameter<edm::InputTag>("isolationFlags")))
    , pvToken_(consumes<reco::VertexCollection>(
          iConfig.getParameter<edm::InputTag>("pvCollection")))
    , tracksToken_(consumes<reco::TrackCollection>(
          iConfig.getParameter<edm::InputTag>("tracks")))
    , hasGenInfo_(iConfig.getParameter<bool>("hasGenInfo"))
    , motherPdgId_(iConfig.getParameter<int>("motherPdgId"))
    , genDRCut_(iConfig.getParameter<double>("genDRCut"))
    , passSelDRCut_(iConfig.getParameter<double>("passSelDRCut"))
{
  usesResource("TFileService");
  if (hasGenInfo_)
    genToken_ = consumes<reco::GenParticleCollection>(
        iConfig.getParameter<edm::InputTag>("genParticles"));
  metToken_ = consumes<pat::METCollection>(
      iConfig.getParameter<edm::InputTag>("MET"));
}

// ---------------------------------------------------------------------------
void HyddraSVsEXOAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Events", "HyddraSV EXONanoAOD prototype");

  // Event ID
  tree_->Branch("run",   &run_);
  tree_->Branch("lumi",  &lumi_);
  tree_->Branch("event", &event_);

  // Event-level quantities
  tree_->Branch("nRecoElectrons",      &nRecoElectrons_);
  tree_->Branch("Event_MET",     &event_MET_);

  // SV count
  tree_->Branch("nHyddraSV", &nHyddraSV_);

  // Vertex position
  tree_->Branch("HyddraSV_x",    &sv_x_);
  tree_->Branch("HyddraSV_y",    &sv_y_);
  tree_->Branch("HyddraSV_z",    &sv_z_);
  tree_->Branch("HyddraSV_xErr", &sv_xErr_);
  tree_->Branch("HyddraSV_yErr", &sv_yErr_);
  tree_->Branch("HyddraSV_zErr", &sv_zErr_);

  // Displacement
  tree_->Branch("HyddraSV_dxy",    &sv_dxy_);
  tree_->Branch("HyddraSV_dxyErr", &sv_dxyErr_);
  tree_->Branch("HyddraSV_dxySig", &sv_dxySig_);

  // Fit quality
  tree_->Branch("HyddraSV_chi2",     &sv_chi2_);
  tree_->Branch("HyddraSV_ndof",     &sv_ndof_);
  tree_->Branch("HyddraSV_normChi2", &sv_normChi2_);

  // 4-momentum
  tree_->Branch("HyddraSV_pt",   &sv_pt_);
  tree_->Branch("HyddraSV_eta",  &sv_eta_);
  tree_->Branch("HyddraSV_phi",  &sv_phi_);
  tree_->Branch("HyddraSV_mass", &sv_mass_);
  tree_->Branch("HyddraSV_p",    &sv_p_);

  // Topology
  tree_->Branch("HyddraSV_charge",     &sv_charge_);
  tree_->Branch("HyddraSV_cosTheta",   &sv_cosTheta_);
  tree_->Branch("HyddraSV_decayAngle", &sv_decayAngle_);
  tree_->Branch("HyddraSV_dR",   &sv_dR_);
  tree_->Branch("HyddraSV_beta", &sv_beta_);

  // Isolation flag
  tree_->Branch("HyddraSV_passesIsolation", &sv_passesIsolation_);

  // Track 1
  tree_->Branch("HyddraSV_trk1Pt",         &sv_trk1Pt_);
  tree_->Branch("HyddraSV_trk1Eta",        &sv_trk1Eta_);
  tree_->Branch("HyddraSV_trk1Phi",        &sv_trk1Phi_);
  tree_->Branch("HyddraSV_trk1Charge",     &sv_trk1Charge_);
  tree_->Branch("HyddraSV_trk1Dxy",        &sv_trk1Dxy_);
  tree_->Branch("HyddraSV_trk1DxyErr",     &sv_trk1DxyErr_);
  tree_->Branch("HyddraSV_trk1DxySig",     &sv_trk1DxySig_);
  tree_->Branch("HyddraSV_trk1Dz",         &sv_trk1Dz_);
  tree_->Branch("HyddraSV_trk1DzErr",      &sv_trk1DzErr_);
  tree_->Branch("HyddraSV_trk1NormChi2",   &sv_trk1NormChi2_);
  tree_->Branch("HyddraSV_trk1CosTheta",   &sv_trk1CosTheta_);
  tree_->Branch("HyddraSV_trk1CosThetaCM", &sv_trk1CosThetaCM_);

  // Track 2
  tree_->Branch("HyddraSV_trk2Pt",         &sv_trk2Pt_);
  tree_->Branch("HyddraSV_trk2Eta",        &sv_trk2Eta_);
  tree_->Branch("HyddraSV_trk2Phi",        &sv_trk2Phi_);
  tree_->Branch("HyddraSV_trk2Charge",     &sv_trk2Charge_);
  tree_->Branch("HyddraSV_trk2Dxy",        &sv_trk2Dxy_);
  tree_->Branch("HyddraSV_trk2DxyErr",     &sv_trk2DxyErr_);
  tree_->Branch("HyddraSV_trk2DxySig",     &sv_trk2DxySig_);
  tree_->Branch("HyddraSV_trk2Dz",         &sv_trk2Dz_);
  tree_->Branch("HyddraSV_trk2DzErr",      &sv_trk2DzErr_);
  tree_->Branch("HyddraSV_trk2NormChi2",   &sv_trk2NormChi2_);
  tree_->Branch("HyddraSV_trk2CosTheta",   &sv_trk2CosTheta_);
  tree_->Branch("HyddraSV_trk2CosThetaCM", &sv_trk2CosThetaCM_);

  // Gen match (MC only)
  if (hasGenInfo_) {
    tree_->Branch("HyddraSV_trk1GenDR",       &sv_trk1GenDR_);
    tree_->Branch("HyddraSV_trk2GenDR",       &sv_trk2GenDR_);
    tree_->Branch("HyddraSV_trk1GenRelPtDiff",&sv_trk1GenRelPtDiff_);
    tree_->Branch("HyddraSV_trk2GenRelPtDiff",&sv_trk2GenRelPtDiff_);
    tree_->Branch("HyddraSV_genVtxIdx",        &sv_genVtxIdx_);
    tree_->Branch("HyddraSV_nearestGenVtxIdx", &sv_nearestGenVtxIdx_);
    tree_->Branch("HyddraSV_min3D",            &sv_min3D_);
    tree_->Branch("HyddraSV_isGold",           &sv_isGold_);
    tree_->Branch("HyddraSV_isBronze",         &sv_isBronze_);

    // Gen vertex collection
    tree_->Branch("nHyddraGenSV",              &nHyddraGenSV_);
    tree_->Branch("HyddraGenSV_x",             &gv_x_);
    tree_->Branch("HyddraGenSV_y",             &gv_y_);
    tree_->Branch("HyddraGenSV_z",             &gv_z_);
    tree_->Branch("HyddraGenSV_dxy",           &gv_dxy_);
    tree_->Branch("HyddraGenSV_pt",            &gv_pt_);
    tree_->Branch("HyddraGenSV_eta",           &gv_eta_);
    tree_->Branch("HyddraGenSV_phi",           &gv_phi_);
    tree_->Branch("HyddraGenSV_mass",          &gv_mass_);
    tree_->Branch("HyddraGenSV_trk1Pt",        &gv_trk1Pt_);
    tree_->Branch("HyddraGenSV_trk1Eta",       &gv_trk1Eta_);
    tree_->Branch("HyddraGenSV_trk1Phi",       &gv_trk1Phi_);
    tree_->Branch("HyddraGenSV_trk1PdgId",     &gv_trk1PdgId_);
    tree_->Branch("HyddraGenSV_trk2Pt",        &gv_trk2Pt_);
    tree_->Branch("HyddraGenSV_trk2Eta",       &gv_trk2Eta_);
    tree_->Branch("HyddraGenSV_trk2Phi",       &gv_trk2Phi_);
    tree_->Branch("HyddraGenSV_trk2PdgId",     &gv_trk2PdgId_);
    tree_->Branch("HyddraGenSV_motherPdgId",   &gv_motherPdgId_);
    tree_->Branch("HyddraGenSV_passSelection", &gv_passSelection_);
    tree_->Branch("HyddraGenSV_isReconstructed",&gv_isReconstructed_);
  }
}

// ---------------------------------------------------------------------------
void HyddraSVsEXOAnalyzer::clearBranches() {
  event_MET_ = -1.f;
  sv_x_.clear(); sv_y_.clear(); sv_z_.clear();
  sv_xErr_.clear(); sv_yErr_.clear(); sv_zErr_.clear();
  sv_dxy_.clear(); sv_dxyErr_.clear(); sv_dxySig_.clear();
  sv_chi2_.clear(); sv_ndof_.clear(); sv_normChi2_.clear();
  sv_pt_.clear(); sv_eta_.clear(); sv_phi_.clear(); sv_mass_.clear(); sv_p_.clear();
  sv_charge_.clear();
  sv_cosTheta_.clear(); sv_decayAngle_.clear(); sv_dR_.clear(); sv_beta_.clear();
  sv_passesIsolation_.clear();
  sv_trk1Pt_.clear(); sv_trk1Eta_.clear(); sv_trk1Phi_.clear(); sv_trk1Charge_.clear();
  sv_trk1Dxy_.clear(); sv_trk1DxyErr_.clear(); sv_trk1DxySig_.clear();
  sv_trk1Dz_.clear(); sv_trk1DzErr_.clear(); sv_trk1NormChi2_.clear();
  sv_trk1CosTheta_.clear(); sv_trk1CosThetaCM_.clear();
  sv_trk2Pt_.clear(); sv_trk2Eta_.clear(); sv_trk2Phi_.clear(); sv_trk2Charge_.clear();
  sv_trk2Dxy_.clear(); sv_trk2DxyErr_.clear(); sv_trk2DxySig_.clear();
  sv_trk2Dz_.clear(); sv_trk2DzErr_.clear(); sv_trk2NormChi2_.clear();
  sv_trk2CosTheta_.clear(); sv_trk2CosThetaCM_.clear();
  sv_trk1GenDR_.clear(); sv_trk2GenDR_.clear();
  sv_trk1GenRelPtDiff_.clear(); sv_trk2GenRelPtDiff_.clear();
  sv_genVtxIdx_.clear(); sv_nearestGenVtxIdx_.clear(); sv_min3D_.clear();
  sv_isGold_.clear(); sv_isBronze_.clear();
  gv_x_.clear(); gv_y_.clear(); gv_z_.clear(); gv_dxy_.clear();
  gv_pt_.clear(); gv_eta_.clear(); gv_phi_.clear(); gv_mass_.clear();
  gv_trk1Pt_.clear(); gv_trk1Eta_.clear(); gv_trk1Phi_.clear(); gv_trk1PdgId_.clear();
  gv_trk2Pt_.clear(); gv_trk2Eta_.clear(); gv_trk2Phi_.clear(); gv_trk2PdgId_.clear();
  gv_motherPdgId_.clear(); gv_passSelection_.clear(); gv_isReconstructed_.clear();
}

// ---------------------------------------------------------------------------
// Build gen signal vertices: find status-1 charged daughters whose mother
// chain leads to a particle with |pdgId| == motherPdgId_.
// ---------------------------------------------------------------------------
std::vector<HyddraGenVertex>
HyddraSVsEXOAnalyzer::buildGenVertices(const reco::GenParticleCollection& genParts) const {
  // Map: signal-mother index -> list of daughter indices
  std::map<size_t, std::vector<size_t>> motherToDaughters;

  for (size_t i = 0; i < genParts.size(); ++i) {
    const auto& p = genParts[i];
    if (p.status() != 1 || p.charge() == 0) continue;

    // Walk the mother chain (up to 10 steps) to find the signal mother
    const reco::GenParticle* cur = &p;
    for (int step = 0; step < 10; ++step) {
      if (cur->numberOfMothers() == 0) break;
      const reco::GenParticle* mom =
          dynamic_cast<const reco::GenParticle*>(cur->mother(0));
      if (!mom) break;
      if (std::abs(mom->pdgId()) == motherPdgId_) {
        // Find the index of this mother in the collection
        size_t momIdx = mom - &genParts[0];
        motherToDaughters[momIdx].push_back(i);
        break;
      }
      // Skip intermediate same-particle propagations
      if (mom->pdgId() == p.pdgId()) { cur = mom; continue; }
      break;
    }
  }

  std::vector<HyddraGenVertex> result;
  for (auto& [momIdx, dauIdxs] : motherToDaughters) {
    if (dauIdxs.size() < 2) continue;

    // Sum 4-momenta of daughters (massless approximation for track-level)
    float px_tot = 0, py_tot = 0, pz_tot = 0, e_tot = 0;
    for (size_t di : dauIdxs) {
      const auto& d = genParts[di];
      px_tot += d.px(); py_tot += d.py(); pz_tot += d.pz();
      e_tot  += std::sqrt(d.px()*d.px() + d.py()*d.py() + d.pz()*d.pz()
                          + d.mass()*d.mass());
    }
    float p_tot = std::sqrt(px_tot*px_tot + py_tot*py_tot + pz_tot*pz_tot);

    HyddraGenVertex gv;
    // Decay position from first daughter's production vertex
    gv.x = genParts[dauIdxs[0]].vx();
    gv.y = genParts[dauIdxs[0]].vy();
    gv.z = genParts[dauIdxs[0]].vz();
    gv.dxy = std::sqrt(gv.x*gv.x + gv.y*gv.y);
    gv.pt  = std::sqrt(px_tot*px_tot + py_tot*py_tot);
    gv.eta = (p_tot > 1e-9f) ? std::atanh(pz_tot / p_tot) : 0.f;
    gv.phi = std::atan2(py_tot, px_tot);
    gv.mass = std::sqrt(std::max(0.f, e_tot*e_tot - p_tot*p_tot));
    gv.motherPdgId = genParts[momIdx].pdgId();
    gv.dauIdx = dauIdxs;

    // Sort daughters by pT descending for trk1/trk2 convention
    std::sort(gv.dauIdx.begin(), gv.dauIdx.end(), [&](size_t a, size_t b) {
      return genParts[a].pt() > genParts[b].pt();
    });

    const auto& d1 = genParts[gv.dauIdx[0]];
    const auto& d2 = genParts[gv.dauIdx[1]];
    gv.trk1Pt = d1.pt(); gv.trk1Eta = d1.eta(); gv.trk1Phi = d1.phi(); gv.trk1PdgId = d1.pdgId();
    gv.trk2Pt = d2.pt(); gv.trk2Eta = d2.eta(); gv.trk2Phi = d2.phi(); gv.trk2PdgId = d2.pdgId();

    result.push_back(gv);
  }
  return result;
}

// ---------------------------------------------------------------------------
// Match a reco track to the nearest signal gen daughter.
// Returns (genVtxIdx, deltaR, relPtDiff). (-1, INV, INV) if no match.
// ---------------------------------------------------------------------------
std::tuple<int, float, float>
HyddraSVsEXOAnalyzer::matchTrackToGen(const reco::Track& trk,
                                       const std::vector<HyddraGenVertex>& genVtxs,
                                       const reco::GenParticleCollection& genParts) const {
  int   bestGV  = -1;
  float bestDR  = std::numeric_limits<float>::max();
  float bestRPT = INV;

  for (int gi = 0; gi < (int)genVtxs.size(); ++gi) {
    for (size_t di : genVtxs[gi].dauIdx) {
      const auto& d = genParts[di];
      float dr = deltaR(trk.eta(), trk.phi(), d.eta(), d.phi());
      if (dr < bestDR) {
        bestDR  = dr;
        bestGV  = gi;
        bestRPT = (d.pt() > 1e-9f) ? (trk.pt() - d.pt()) / d.pt() : INV;
      }
    }
  }
  // Return uncapped DR for storage; caller checks against genDRCut_ for gold/bronze
  return { bestGV, bestDR, bestRPT };
}

// ---------------------------------------------------------------------------
// passSelection: does this gen daughter have any reco track within passSelDRCut_?
// ---------------------------------------------------------------------------
bool HyddraSVsEXOAnalyzer::hasSoftMatch(size_t dauIdx,
                                         const reco::GenParticleCollection& genParts,
                                         const reco::TrackCollection& tracks) const {
  const auto& d = genParts[dauIdx];
  for (const auto& trk : tracks) {
    if (deltaR(trk.eta(), trk.phi(), d.eta(), d.phi()) < passSelDRCut_)
      return true;
  }
  return false;
}

// ---------------------------------------------------------------------------
void HyddraSVsEXOAnalyzer::analyze(const edm::Event& iEvent,
                                    const edm::EventSetup& iSetup) {
  clearBranches();

  run_   = iEvent.id().run();
  lumi_  = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();

  // ── Inputs ──────────────────────────────────────────────────────────────
  auto inclusiveHandle = iEvent.getHandle(inclusiveToken_);
  auto isoFlagsHandle  = iEvent.getHandle(isoFlagsToken_);
  auto pvHandle        = iEvent.getHandle(pvToken_);
  auto tracksHandle    = iEvent.getHandle(tracksToken_);

  // ── Event-level quantities (always filled, even on early return) ──────────
  nRecoElectrons_ = tracksHandle.isValid() ? static_cast<int>(tracksHandle->size()) : 0;

  // ── MET ──────────────────────────────────────────────────────────────────
  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByToken(metToken_, metHandle);
  if (metHandle.isValid() && !metHandle->empty())
    event_MET_ = metHandle->at(0).pt();

  if (!inclusiveHandle.isValid() || pvHandle->empty()) {
    tree_->Fill();
    return;
  }
  const reco::Vertex& pv = pvHandle->at(0);

  // ── Gen vertices (MC only) ───────────────────────────────────────────────
  std::vector<HyddraGenVertex> genVtxs;
  const reco::GenParticleCollection* genParts = nullptr;
  edm::Handle<reco::GenParticleCollection> genHandle;
  if (hasGenInfo_) {
    iEvent.getByToken(genToken_, genHandle);
    if (genHandle.isValid()) {
      genParts = genHandle.product();
      genVtxs  = buildGenVertices(*genParts);
      // Compute passSelection for each gen vertex
      const reco::TrackCollection& tracks = *tracksHandle;
      for (auto& gv : genVtxs) {
        bool allMatched = true;
        for (size_t di : gv.dauIdx) {
          if (!hasSoftMatch(di, *genParts, tracks)) { allMatched = false; break; }
        }
        gv.passSelection = allMatched;
      }
    }
  }

  // ── Fill HyddraSV branches ───────────────────────────────────────────────
  int svIdx = 0;
  for (const auto& sv : *inclusiveHandle) {
    // Collect tracks from vertex, sort by pT descending
    std::vector<reco::TrackRef> trkRefs;
    for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it)
      trkRefs.push_back(it->castTo<reco::TrackRef>());
    if (trkRefs.size() != 2) { ++svIdx; continue; }
    std::sort(trkRefs.begin(), trkRefs.end(), [](const reco::TrackRef& a, const reco::TrackRef& b) {
      return a->pt() > b->pt();
    });
    const reco::Track& t1 = *trkRefs[0];
    const reco::Track& t2 = *trkRefs[1];

    // ── Isolation flag (from producer-computed flags, parallel to inclusiveVertices) ──
    bool passIso = (isoFlagsHandle.isValid() && svIdx < (int)isoFlagsHandle->size())
                   ? (*isoFlagsHandle)[svIdx] != 0 : false;

    // ── 4-momentum (massless track sum) ───────────────────────────────────
    float px = t1.px() + t2.px();
    float py = t1.py() + t2.py();
    float pz = t1.pz() + t2.pz();
    float e  = t1.p()  + t2.p();
    float p  = std::sqrt(px*px + py*py + pz*pz);
    float pt = std::sqrt(px*px + py*py);
    float eta = (p > 1e-9f) ? std::atanh(pz / p) : 0.f;
    float phi = std::atan2(py, px);
    float mass = std::sqrt(std::max(0.f, e*e - p*p));
    float beta_e = std::sqrt(p*p + mass*mass);
    float beta = (beta_e > 1e-9f) ? p / beta_e : INV;

    // ── Displacement ────────────────────────────────────────────────────
    float dxy    = computeDxy(sv, pv);
    float dxyErr = computeDxyErr(sv, pv);
    float dxySig = (dxyErr > 1e-9f) ? dxy / dxyErr : 0.f;

    // ── Topology ────────────────────────────────────────────────────────
    float cosTheta   = computeCosTheta(sv, pv, px, py, pz);
    float decayAngle = computeDecayAngle(t1.px(), t1.py(), t1.pz(), t1.charge(),
                                          t2.px(), t2.py(), t2.pz(), t2.charge());
    float dR = deltaR(t1.eta(), t1.phi(), t2.eta(), t2.phi());
    int   charge = t1.charge() + t2.charge();

    // ── Per-track variables ──────────────────────────────────────────────
    // Track cosTheta: angle between individual track and vertex momentum
    float t1ct = (pt > 1e-9f || std::abs(pz) > 1e-9f) ?
        (t1.px()*px + t1.py()*py + t1.pz()*pz) / (t1.p() * p) : INV;
    float t2ct = (p > 1e-9f) ?
        (t2.px()*px + t2.py()*py + t2.pz()*pz) / (t2.p() * p) : INV;
    // Track CM cosTheta
    float t1ctCM = computeTrackCosThetaCM(t1.px(), t1.py(), t1.pz(), px, py, pz, e);
    float t2ctCM = computeTrackCosThetaCM(t2.px(), t2.py(), t2.pz(), px, py, pz, e);

    // dxy/dz wrt PV
    float t1dxy = t1.dxy(pv.position()); float t1dxyErr = t1.dxyError(pv.position(), pv.covariance());
    float t2dxy = t2.dxy(pv.position()); float t2dxyErr = t2.dxyError(pv.position(), pv.covariance());
    float t1dz  = t1.dz(pv.position());  float t1dzErr  = t1.dzError();
    float t2dz  = t2.dz(pv.position());  float t2dzErr  = t2.dzError();

    // ── Fill SV branches ────────────────────────────────────────────────
    sv_x_.push_back(sv.x()); sv_y_.push_back(sv.y()); sv_z_.push_back(sv.z());
    sv_xErr_.push_back(sv.xError()); sv_yErr_.push_back(sv.yError()); sv_zErr_.push_back(sv.zError());
    sv_dxy_.push_back(dxy); sv_dxyErr_.push_back(dxyErr); sv_dxySig_.push_back(dxySig);
    sv_chi2_.push_back(sv.chi2()); sv_ndof_.push_back(sv.ndof());
    sv_normChi2_.push_back(sv.normalizedChi2());
    sv_pt_.push_back(pt); sv_eta_.push_back(eta); sv_phi_.push_back(phi);
    sv_mass_.push_back(mass); sv_p_.push_back(p);
    sv_charge_.push_back(charge);
    sv_cosTheta_.push_back(cosTheta); sv_decayAngle_.push_back(decayAngle);
    sv_dR_.push_back(dR); sv_beta_.push_back(beta);
    sv_passesIsolation_.push_back(passIso);
    sv_trk1Pt_.push_back(t1.pt()); sv_trk1Eta_.push_back(t1.eta()); sv_trk1Phi_.push_back(t1.phi());
    sv_trk1Charge_.push_back(t1.charge());
    sv_trk1Dxy_.push_back(t1dxy); sv_trk1DxyErr_.push_back(t1dxyErr);
    sv_trk1DxySig_.push_back((t1dxyErr > 1e-9f) ? t1dxy / t1dxyErr : 0.f);
    sv_trk1Dz_.push_back(t1dz); sv_trk1DzErr_.push_back(t1dzErr);
    sv_trk1NormChi2_.push_back(t1.normalizedChi2());
    sv_trk1CosTheta_.push_back(t1ct); sv_trk1CosThetaCM_.push_back(t1ctCM);
    sv_trk2Pt_.push_back(t2.pt()); sv_trk2Eta_.push_back(t2.eta()); sv_trk2Phi_.push_back(t2.phi());
    sv_trk2Charge_.push_back(t2.charge());
    sv_trk2Dxy_.push_back(t2dxy); sv_trk2DxyErr_.push_back(t2dxyErr);
    sv_trk2DxySig_.push_back((t2dxyErr > 1e-9f) ? t2dxy / t2dxyErr : 0.f);
    sv_trk2Dz_.push_back(t2dz); sv_trk2DzErr_.push_back(t2dzErr);
    sv_trk2NormChi2_.push_back(t2.normalizedChi2());
    sv_trk2CosTheta_.push_back(t2ct); sv_trk2CosThetaCM_.push_back(t2ctCM);

    ++svIdx;

    // ── Gen matching ────────────────────────────────────────────────────
    if (hasGenInfo_ && genParts) {
      auto [gv1, dr1, rpt1] = matchTrackToGen(t1, genVtxs, *genParts);
      auto [gv2, dr2, rpt2] = matchTrackToGen(t2, genVtxs, *genParts);

      sv_trk1GenDR_.push_back(dr1 < 900.f ? dr1 : INV);
      sv_trk2GenDR_.push_back(dr2 < 900.f ? dr2 : INV);
      sv_trk1GenRelPtDiff_.push_back(rpt1);
      sv_trk2GenRelPtDiff_.push_back(rpt2);

      bool t1matched = (gv1 >= 0 && dr1 < genDRCut_);
      bool t2matched = (gv2 >= 0 && dr2 < genDRCut_);
      bool isGold    = t1matched && t2matched && (gv1 == gv2);
      bool isBronze  = t1matched || t2matched;

      int goldGV = isGold ? gv1 : -1;
      if (isGold) genVtxs[goldGV].isReconstructed = true;

      sv_isGold_.push_back(isGold);
      sv_isBronze_.push_back(isBronze);
      sv_genVtxIdx_.push_back(goldGV);

      // Nearest gen vertex by 3D distance (independent of track matching)
      int   nearestIdx = -1;
      float minDist    = std::numeric_limits<float>::max();
      for (int gi = 0; gi < (int)genVtxs.size(); ++gi) {
        float dx = sv.x() - genVtxs[gi].x;
        float dy = sv.y() - genVtxs[gi].y;
        float dz = sv.z() - genVtxs[gi].z;
        float d  = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (d < minDist) { minDist = d; nearestIdx = gi; }
      }
      sv_nearestGenVtxIdx_.push_back(nearestIdx);
      sv_min3D_.push_back(nearestIdx >= 0 ? minDist : INV);
    }
  }

  nHyddraSV_ = static_cast<int>(sv_pt_.size());

  // ── Fill HyddraGenSV branches ─────────────────────────────────────────────
  if (hasGenInfo_) {
    for (const auto& gv : genVtxs) {
      gv_x_.push_back(gv.x); gv_y_.push_back(gv.y); gv_z_.push_back(gv.z);
      gv_dxy_.push_back(gv.dxy);
      gv_pt_.push_back(gv.pt); gv_eta_.push_back(gv.eta);
      gv_phi_.push_back(gv.phi); gv_mass_.push_back(gv.mass);
      gv_trk1Pt_.push_back(gv.trk1Pt); gv_trk1Eta_.push_back(gv.trk1Eta);
      gv_trk1Phi_.push_back(gv.trk1Phi); gv_trk1PdgId_.push_back(gv.trk1PdgId);
      gv_trk2Pt_.push_back(gv.trk2Pt); gv_trk2Eta_.push_back(gv.trk2Eta);
      gv_trk2Phi_.push_back(gv.trk2Phi); gv_trk2PdgId_.push_back(gv.trk2PdgId);
      gv_motherPdgId_.push_back(gv.motherPdgId);
      gv_passSelection_.push_back(gv.passSelection);
      gv_isReconstructed_.push_back(gv.isReconstructed);
    }
    nHyddraGenSV_ = static_cast<int>(gv_pt_.size());
  }

  tree_->Fill();
}

// ---------------------------------------------------------------------------
void HyddraSVsEXOAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inclusiveVertices",
      edm::InputTag("hyddraEXO", "inclusiveVertices"));
  desc.add<edm::InputTag>("isolationFlags",
      edm::InputTag("hyddraEXO", "isolationFlags"));
  desc.add<edm::InputTag>("pvCollection",
      edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("tracks",
      edm::InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"));
  desc.add<edm::InputTag>("genParticles",
      edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("MET",
      edm::InputTag("slimmedMETs"));
  desc.add<bool>("hasGenInfo",     true);
  desc.add<int> ("motherPdgId",   54);     // default: dark photon proxy
  desc.add<double>("genDRCut",    0.05);   // gold/bronze track matching threshold
  desc.add<double>("passSelDRCut",0.02);   // passSelection (acceptance) threshold
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsEXOAnalyzer);
