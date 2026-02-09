// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      LLPNanoSVAnalyzer
//
// Description: EDAnalyzer that reads LLPNanoAOD FlatTables (PatMuonVertex,
//              PatDSAMuonVertex, GenPart, Muon, DSAMuon) and writes a flat
//              TTree with the same branch structure as HyddraSVAnalyzer.
//              This allows direct comparison of SV reconstruction efficiency
//              between HYDDRA and the LLPNanoAOD dimuon vertex tools using
//              the existing efficiency plotting scripts.
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <limits>
#include <cmath>
#include <map>
#include <vector>
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

// NanoAOD FlatTable
#include "DataFormats/NanoAOD/interface/FlatTable.h"

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"

class LLPNanoSVAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit LLPNanoSVAnalyzer(const edm::ParameterSet&);
  ~LLPNanoSVAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  void clearBranches();

  // Helper to get a column from a FlatTable, returning an empty vector if not found
  template <typename T>
  std::vector<T> getColumn(const nanoaod::FlatTable& table, const std::string& name) const;

  // Specialization for bool columns stored as uint8
  std::vector<bool> getBoolColumn(const nanoaod::FlatTable& table, const std::string& name) const;

  // Gen vertex reconstruction from GenPart table
  struct GenMuonVertex {
    float x, y, z;
    float dxy;
    float pt, eta, phi, mass;
    int motherIndex;  // GenPart index of the mother
    std::vector<int> muonGenPartIndices;  // GenPart indices of the daughter muons
  };

  std::vector<GenMuonVertex> buildGenVertices(
    const std::vector<int>& genPdgId,
    const std::vector<int>& genMotherIdx,
    const std::vector<float>& genVx,
    const std::vector<float>& genVy,
    const std::vector<float>& genVz,
    const std::vector<float>& genPt,
    const std::vector<float>& genEta,
    const std::vector<float>& genPhi,
    const std::vector<float>& genMass,
    const std::vector<int>& genStatus) const;

  // Gold matching: trace reco vertex muon legs back to gen particles
  int findGoldMatch(
    int vtxIdx,
    const std::vector<int>& vtxOrigMuonIdx1,
    const std::vector<int>& vtxOrigMuonIdx2,
    const std::vector<bool>& vtxIsDSAMuon1,
    const std::vector<bool>& vtxIsDSAMuon2,
    const std::vector<int>& muonGenPartIdx,
    const std::vector<float>& muonPt,
    const std::vector<float>& muonEta,
    const std::vector<float>& muonPhi,
    const std::vector<float>& dsaMuonPt,
    const std::vector<float>& dsaMuonEta,
    const std::vector<float>& dsaMuonPhi,
    const std::vector<int>& genPdgId,
    const std::vector<int>& genMotherIdx,
    const std::vector<float>& genPt,
    const std::vector<float>& genEta,
    const std::vector<float>& genPhi,
    const std::vector<GenMuonVertex>& genVertices) const;

  // DSA muon to gen muon matching via deltaR
  int matchDSAToGen(
    float dsaEta, float dsaPhi, float dsaPt,
    const std::vector<float>& genPt,
    const std::vector<float>& genEta,
    const std::vector<float>& genPhi,
    const std::vector<int>& genPdgId) const;

  // Find which gen vertex a gen particle belongs to
  int findGenVertexForParticle(int genPartIdx, const std::vector<GenMuonVertex>& genVertices) const;

  // Spatial matching
  int findNearestGenVertex(float vx, float vy, float vz,
                           const std::vector<GenMuonVertex>& genVertices,
                           float& minDistance) const;

  // Configuration
  std::string collection_;
  int motherPdgId_;
  double dsaDeltaR_;
  double dsaRelPtDiff_;

  // Tokens
  edm::EDGetTokenT<nanoaod::FlatTable> genPartTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> muonTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> dsaMuonTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> vertexTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> refittedTracksTableToken_;

  // TTree
  TTree* tree_;

  // Branch variables - Reco vertex branches (matching HyddraSVAnalyzer names)
  unsigned int nVertices_;
  std::vector<bool> isLeptonic_;
  std::vector<unsigned int> nTracks_;
  std::vector<float> x_;
  std::vector<float> y_;
  std::vector<float> z_;
  std::vector<float> dxy_;
  std::vector<float> mass_;
  std::vector<float> pt_;
  std::vector<float> eta_;
  std::vector<float> phi_;
  std::vector<float> chi2_;
  std::vector<float> normalizedChi2_;
  std::vector<float> ndof_;

  // Gen matching branches (reco vertex)
  std::vector<int> genVertexIndex_;
  std::vector<int> nearestGenVertexIndex_;
  std::vector<float> min3D_;
  std::vector<bool> isBronze_;
  std::vector<bool> isSilver_;
  std::vector<bool> isGold_;

  // Gen vertex branches
  unsigned int genVertex_nTotal_;
  unsigned int genVertex_nMuon_;
  unsigned int genVertex_nElectron_;
  unsigned int genVertex_nHadronic_;
  std::vector<float> genVertex_dxy_;
  std::vector<float> genVertex_x_;
  std::vector<float> genVertex_y_;
  std::vector<float> genVertex_z_;
  std::vector<float> genVertex_pt_;
  std::vector<float> genVertex_eta_;
  std::vector<float> genVertex_phi_;
  std::vector<float> genVertex_mass_;
  std::vector<bool> genVertex_isMuon_;
  std::vector<bool> genVertex_isElectron_;
  std::vector<bool> genVertex_isHadronic_;
  std::vector<bool> genVertex_passSelection_;
};

// ============================================================================
// Constructor
// ============================================================================

LLPNanoSVAnalyzer::LLPNanoSVAnalyzer(const edm::ParameterSet& iConfig) :
  collection_(iConfig.getParameter<std::string>("collection")),
  motherPdgId_(iConfig.getParameter<int>("motherPdgId")),
  dsaDeltaR_(iConfig.getParameter<double>("dsaDeltaR")),
  dsaRelPtDiff_(iConfig.getParameter<double>("dsaRelPtDiff")),
  genPartTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("genPartTable"))),
  muonTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("muonTable"))),
  dsaMuonTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("dsaMuonTable"))),
  vertexTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("vertexTable"))),
  refittedTracksTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("refittedTracksTable")))
{
  usesResource("TFileService");
}

// ============================================================================
// beginJob - create TTree with branches matching HyddraSVAnalyzer
// ============================================================================

void LLPNanoSVAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "LLPNanoSV Tree");

  // Reco vertex branches
  tree_->Branch("HyddraSV_nVertices", &nVertices_);
  tree_->Branch("HyddraSV_isLeptonic", &isLeptonic_);
  tree_->Branch("HyddraSV_nTracks", &nTracks_);
  tree_->Branch("HyddraSV_x", &x_);
  tree_->Branch("HyddraSV_y", &y_);
  tree_->Branch("HyddraSV_z", &z_);
  tree_->Branch("HyddraSV_dxy", &dxy_);
  tree_->Branch("HyddraSV_mass", &mass_);
  tree_->Branch("HyddraSV_pt", &pt_);
  tree_->Branch("HyddraSV_eta", &eta_);
  tree_->Branch("HyddraSV_phi", &phi_);
  tree_->Branch("HyddraSV_chi2", &chi2_);
  tree_->Branch("HyddraSV_normalizedChi2", &normalizedChi2_);
  tree_->Branch("HyddraSV_ndof", &ndof_);

  // Gen matching branches
  tree_->Branch("HyddraSV_genVertexIndex", &genVertexIndex_);
  tree_->Branch("HyddraSV_nearestGenVertexIndex", &nearestGenVertexIndex_);
  tree_->Branch("HyddraSV_min3D", &min3D_);
  tree_->Branch("HyddraSV_isBronze", &isBronze_);
  tree_->Branch("HyddraSV_isSilver", &isSilver_);
  tree_->Branch("HyddraSV_isGold", &isGold_);

  // Gen vertex branches
  tree_->Branch("HyddraGenVertex_nTotal", &genVertex_nTotal_);
  tree_->Branch("HyddraGenVertex_nMuon", &genVertex_nMuon_);
  tree_->Branch("HyddraGenVertex_nElectron", &genVertex_nElectron_);
  tree_->Branch("HyddraGenVertex_nHadronic", &genVertex_nHadronic_);
  tree_->Branch("HyddraGenVertex_dxy", &genVertex_dxy_);
  tree_->Branch("HyddraGenVertex_x", &genVertex_x_);
  tree_->Branch("HyddraGenVertex_y", &genVertex_y_);
  tree_->Branch("HyddraGenVertex_z", &genVertex_z_);
  tree_->Branch("HyddraGenVertex_pt", &genVertex_pt_);
  tree_->Branch("HyddraGenVertex_eta", &genVertex_eta_);
  tree_->Branch("HyddraGenVertex_phi", &genVertex_phi_);
  tree_->Branch("HyddraGenVertex_mass", &genVertex_mass_);
  tree_->Branch("HyddraGenVertex_isMuon", &genVertex_isMuon_);
  tree_->Branch("HyddraGenVertex_isElectron", &genVertex_isElectron_);
  tree_->Branch("HyddraGenVertex_isHadronic", &genVertex_isHadronic_);
  tree_->Branch("HyddraGenVertex_passSelection", &genVertex_passSelection_);
}

// ============================================================================
// clearBranches
// ============================================================================

void LLPNanoSVAnalyzer::clearBranches() {
  nVertices_ = 0;
  isLeptonic_.clear();
  nTracks_.clear();
  x_.clear();
  y_.clear();
  z_.clear();
  dxy_.clear();
  mass_.clear();
  pt_.clear();
  eta_.clear();
  phi_.clear();
  chi2_.clear();
  normalizedChi2_.clear();
  ndof_.clear();

  genVertexIndex_.clear();
  nearestGenVertexIndex_.clear();
  min3D_.clear();
  isBronze_.clear();
  isSilver_.clear();
  isGold_.clear();

  genVertex_nTotal_ = 0;
  genVertex_nMuon_ = 0;
  genVertex_nElectron_ = 0;
  genVertex_nHadronic_ = 0;
  genVertex_dxy_.clear();
  genVertex_x_.clear();
  genVertex_y_.clear();
  genVertex_z_.clear();
  genVertex_pt_.clear();
  genVertex_eta_.clear();
  genVertex_phi_.clear();
  genVertex_mass_.clear();
  genVertex_isMuon_.clear();
  genVertex_isElectron_.clear();
  genVertex_isHadronic_.clear();
  genVertex_passSelection_.clear();
}

// ============================================================================
// FlatTable column accessors
// ============================================================================

template <typename T>
std::vector<T> LLPNanoSVAnalyzer::getColumn(const nanoaod::FlatTable& table, const std::string& name) const {
  int idx = table.columnIndex(name);
  if (idx < 0) return std::vector<T>();

  std::vector<T> result(table.size());
  for (unsigned int i = 0; i < table.size(); ++i) {
    result[i] = table.columValue<T>(idx, i);
  }
  return result;
}

std::vector<bool> LLPNanoSVAnalyzer::getBoolColumn(const nanoaod::FlatTable& table, const std::string& name) const {
  int idx = table.columnIndex(name);
  if (idx < 0) return std::vector<bool>();

  std::vector<bool> result(table.size());
  for (unsigned int i = 0; i < table.size(); ++i) {
    result[i] = static_cast<bool>(table.columValue<uint8_t>(idx, i));
  }
  return result;
}

// ============================================================================
// analyze - main event loop
// ============================================================================

void LLPNanoSVAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {

  clearBranches();

  // Get FlatTable handles
  edm::Handle<nanoaod::FlatTable> genPartHandle;
  edm::Handle<nanoaod::FlatTable> muonHandle;
  edm::Handle<nanoaod::FlatTable> dsaMuonHandle;
  edm::Handle<nanoaod::FlatTable> vertexHandle;
  edm::Handle<nanoaod::FlatTable> refTracksHandle;

  iEvent.getByToken(genPartTableToken_, genPartHandle);
  iEvent.getByToken(muonTableToken_, muonHandle);
  iEvent.getByToken(dsaMuonTableToken_, dsaMuonHandle);
  iEvent.getByToken(vertexTableToken_, vertexHandle);
  iEvent.getByToken(refittedTracksTableToken_, refTracksHandle);

  if (!vertexHandle.isValid()) {
    tree_->Fill();
    return;
  }

  const auto& vtxTable = *vertexHandle;
  const auto& genTable = *genPartHandle;
  const auto& muonTable = *muonHandle;

  // ---- Read GenPart columns ----
  auto genPdgId = getColumn<int>(genTable, "pdgId");
  auto genMotherIdx = getColumn<int>(genTable, "genPartIdxMother");
  auto genVx = getColumn<float>(genTable, "vx");
  auto genVy = getColumn<float>(genTable, "vy");
  auto genVz = getColumn<float>(genTable, "vz");
  auto genPt = getColumn<float>(genTable, "pt");
  auto genEta = getColumn<float>(genTable, "eta");
  auto genPhi = getColumn<float>(genTable, "phi");
  auto genMass = getColumn<float>(genTable, "mass");
  auto genStatus = getColumn<int>(genTable, "status");

  // ---- Read Muon columns ----
  auto muonGenPartIdx = getColumn<int>(muonTable, "genPartIdx");
  auto muonPt = getColumn<float>(muonTable, "pt");
  auto muonEta = getColumn<float>(muonTable, "eta");
  auto muonPhi = getColumn<float>(muonTable, "phi");

  // ---- Read DSAMuon columns ----
  std::vector<float> dsaMuonPt, dsaMuonEta, dsaMuonPhi;
  if (dsaMuonHandle.isValid()) {
    const auto& dsaTable = *dsaMuonHandle;
    dsaMuonPt = getColumn<float>(dsaTable, "pt");
    dsaMuonEta = getColumn<float>(dsaTable, "eta");
    dsaMuonPhi = getColumn<float>(dsaTable, "phi");
  }

  // ---- Read vertex columns ----
  auto vtxIsValid = getBoolColumn(vtxTable, "isValid");
  auto vtxVx = getColumn<float>(vtxTable, "vx");
  auto vtxVy = getColumn<float>(vtxTable, "vy");
  auto vtxVz = getColumn<float>(vtxTable, "vz");
  auto vtxVxy = getColumn<float>(vtxTable, "vxy");
  auto vtxChi2 = getColumn<float>(vtxTable, "chi2");
  auto vtxNormChi2 = getColumn<float>(vtxTable, "normChi2");
  auto vtxNdof = getColumn<float>(vtxTable, "ndof");
  auto vtxDR = getColumn<float>(vtxTable, "dR");
  auto vtxOrigMuonIdx1 = getColumn<int>(vtxTable, "originalMuonIdx1");
  auto vtxOrigMuonIdx2 = getColumn<int>(vtxTable, "originalMuonIdx2");
  auto vtxIsDSAMuon1 = getBoolColumn(vtxTable, "isDSAMuon1");
  auto vtxIsDSAMuon2 = getBoolColumn(vtxTable, "isDSAMuon2");

  // ---- Read refitted tracks columns ----
  std::vector<float> refPx, refPy, refPz, refPt, refEta, refPhi;
  std::vector<int> refIdx;
  if (refTracksHandle.isValid()) {
    const auto& refTable = *refTracksHandle;
    refPx = getColumn<float>(refTable, "px");
    refPy = getColumn<float>(refTable, "py");
    refPz = getColumn<float>(refTable, "pz");
    refPt = getColumn<float>(refTable, "pt");
    refEta = getColumn<float>(refTable, "eta");
    refPhi = getColumn<float>(refTable, "phi");
    refIdx = getColumn<int>(refTable, "idx");
  }

  // ---- Build gen muon vertices ----
  auto genVertices = buildGenVertices(genPdgId, genMotherIdx, genVx, genVy, genVz,
                                      genPt, genEta, genPhi, genMass, genStatus);

  // ---- Fill gen vertex branches ----
  genVertex_nTotal_ = static_cast<unsigned int>(genVertices.size());
  genVertex_nMuon_ = static_cast<unsigned int>(genVertices.size());  // All are muon vertices
  genVertex_nElectron_ = 0;
  genVertex_nHadronic_ = 0;

  for (const auto& gv : genVertices) {
    genVertex_dxy_.push_back(gv.dxy);
    genVertex_x_.push_back(gv.x);
    genVertex_y_.push_back(gv.y);
    genVertex_z_.push_back(gv.z);
    genVertex_pt_.push_back(gv.pt);
    genVertex_eta_.push_back(gv.eta);
    genVertex_phi_.push_back(gv.phi);
    genVertex_mass_.push_back(gv.mass);
    genVertex_isMuon_.push_back(true);
    genVertex_isElectron_.push_back(false);
    genVertex_isHadronic_.push_back(false);
    genVertex_passSelection_.push_back(true);
  }

  // ---- Process reco vertices ----
  unsigned int nValid = 0;

  for (unsigned int iv = 0; iv < vtxTable.size(); ++iv) {

    // Only keep valid vertices
    if (!vtxIsValid.empty() && !vtxIsValid[iv]) continue;

    nValid++;

    // All PatMuonVertex/PatDSAMuonVertex are dimuon (leptonic)
    isLeptonic_.push_back(true);
    nTracks_.push_back(2);

    x_.push_back(vtxVx[iv]);
    y_.push_back(vtxVy[iv]);
    z_.push_back(vtxVz[iv]);
    dxy_.push_back(vtxVxy[iv]);
    chi2_.push_back(vtxChi2[iv]);
    normalizedChi2_.push_back(vtxNormChi2[iv]);
    ndof_.push_back(vtxNdof[iv]);

    // Compute vertex mass/pt/eta/phi from refitted tracks
    float vtxMass = 0, vtxPt = 0, vtxEta = 0, vtxPhi = 0;

    // Find refitted tracks belonging to this vertex via the idx column
    // The refitted tracks table has an "idx" column that maps back to the vertex index
    TLorentzVector p4sum;
    bool foundTracks = false;

    if (!refIdx.empty()) {
      for (unsigned int it = 0; it < refIdx.size(); ++it) {
        if (refIdx[it] == static_cast<int>(iv)) {
          TLorentzVector p4;
          // Use muon mass hypothesis (0.10566 GeV)
          double energy = std::sqrt(refPx[it]*refPx[it] + refPy[it]*refPy[it] +
                                    refPz[it]*refPz[it] + 0.10566*0.10566);
          p4.SetPxPyPzE(refPx[it], refPy[it], refPz[it], energy);
          p4sum += p4;
          foundTracks = true;
        }
      }
    }

    if (foundTracks) {
      vtxMass = static_cast<float>(p4sum.M());
      vtxPt = static_cast<float>(p4sum.Pt());
      vtxEta = static_cast<float>(p4sum.Eta());
      vtxPhi = static_cast<float>(p4sum.Phi());
    }

    mass_.push_back(vtxMass);
    pt_.push_back(vtxPt);
    eta_.push_back(vtxEta);
    phi_.push_back(vtxPhi);

    // ---- Gold matching ----
    int goldGenVtxIdx = findGoldMatch(
      iv,
      vtxOrigMuonIdx1, vtxOrigMuonIdx2,
      vtxIsDSAMuon1, vtxIsDSAMuon2,
      muonGenPartIdx, muonPt, muonEta, muonPhi,
      dsaMuonPt, dsaMuonEta, dsaMuonPhi,
      genPdgId, genMotherIdx, genPt, genEta, genPhi,
      genVertices);

    bool isGoldMatch = (goldGenVtxIdx >= 0);
    genVertexIndex_.push_back(goldGenVtxIdx);
    isGold_.push_back(isGoldMatch);
    isSilver_.push_back(isGoldMatch);  // Same as gold (only one matching tier)
    isBronze_.push_back(isGoldMatch);  // Same as gold

    // Spatial matching to nearest gen vertex
    float minDist = -1.f;
    int nearestIdx = findNearestGenVertex(vtxVx[iv], vtxVy[iv], vtxVz[iv],
                                           genVertices, minDist);
    nearestGenVertexIndex_.push_back(nearestIdx);
    min3D_.push_back(minDist);
  }

  nVertices_ = nValid;

  tree_->Fill();
}

// ============================================================================
// buildGenVertices - reconstruct gen muon vertices from GenPart table
// ============================================================================

std::vector<LLPNanoSVAnalyzer::GenMuonVertex> LLPNanoSVAnalyzer::buildGenVertices(
  const std::vector<int>& genPdgId,
  const std::vector<int>& genMotherIdx,
  const std::vector<float>& genVx,
  const std::vector<float>& genVy,
  const std::vector<float>& genVz,
  const std::vector<float>& genPt,
  const std::vector<float>& genEta,
  const std::vector<float>& genPhi,
  const std::vector<float>& genMass,
  const std::vector<int>& genStatus) const
{
  // Find muons (|pdgId| == 13) whose mother has pdgId == motherPdgId_
  // Group by mother index -> each group = 1 gen vertex

  std::map<int, std::vector<int>> motherToMuons;

  for (unsigned int i = 0; i < genPdgId.size(); ++i) {
    if (std::abs(genPdgId[i]) != 13) continue;
    if (genStatus[i] != 1) continue;  // Only final-state muons

    int motherIdx = genMotherIdx[i];
    if (motherIdx < 0 || motherIdx >= static_cast<int>(genPdgId.size())) continue;

    // Walk up the chain to find the signal mother
    // Direct mother might be the muon itself (status copy) or an intermediate
    int currentMother = motherIdx;
    bool foundSignalMother = false;

    // Check up to 10 levels to avoid infinite loops
    for (int level = 0; level < 10; ++level) {
      if (currentMother < 0 || currentMother >= static_cast<int>(genPdgId.size())) break;

      if (genPdgId[currentMother] == motherPdgId_) {
        foundSignalMother = true;
        break;
      }

      // If mother is same pdgId muon (status copy), go up
      if (std::abs(genPdgId[currentMother]) == 13) {
        currentMother = genMotherIdx[currentMother];
        continue;
      }

      break;
    }

    if (foundSignalMother) {
      motherToMuons[currentMother].push_back(static_cast<int>(i));
    }
  }

  // Build vertices from grouped muons
  std::vector<GenMuonVertex> vertices;

  for (const auto& [motherIdx, muonIndices] : motherToMuons) {
    if (muonIndices.size() < 2) continue;  // Need at least 2 muons for a vertex

    GenMuonVertex gv;
    gv.motherIndex = motherIdx;
    gv.muonGenPartIndices = muonIndices;

    // Gen vertex position = decay point of the mother = production point of the muons
    int firstMuon = muonIndices[0];
    gv.x = genVx[firstMuon];
    gv.y = genVy[firstMuon];
    gv.z = genVz[firstMuon];
    gv.dxy = std::sqrt(gv.x * gv.x + gv.y * gv.y);

    // Compute 4-vector sum of muon daughters
    TLorentzVector p4sum;
    for (int idx : muonIndices) {
      TLorentzVector p4;
      p4.SetPtEtaPhiM(genPt[idx], genEta[idx], genPhi[idx], genMass[idx]);
      p4sum += p4;
    }

    gv.pt = static_cast<float>(p4sum.Pt());
    gv.eta = static_cast<float>(p4sum.Eta());
    gv.phi = static_cast<float>(p4sum.Phi());
    gv.mass = static_cast<float>(p4sum.M());

    vertices.push_back(gv);
  }

  return vertices;
}

// ============================================================================
// findGoldMatch - gold-equivalent matching for reco vertex
// ============================================================================

int LLPNanoSVAnalyzer::findGoldMatch(
  int vtxIdx,
  const std::vector<int>& vtxOrigMuonIdx1,
  const std::vector<int>& vtxOrigMuonIdx2,
  const std::vector<bool>& vtxIsDSAMuon1,
  const std::vector<bool>& vtxIsDSAMuon2,
  const std::vector<int>& muonGenPartIdx,
  const std::vector<float>& muonPt,
  const std::vector<float>& muonEta,
  const std::vector<float>& muonPhi,
  const std::vector<float>& dsaMuonPt,
  const std::vector<float>& dsaMuonEta,
  const std::vector<float>& dsaMuonPhi,
  const std::vector<int>& genPdgId,
  const std::vector<int>& genMotherIdx,
  const std::vector<float>& genPt,
  const std::vector<float>& genEta,
  const std::vector<float>& genPhi,
  const std::vector<GenMuonVertex>& genVertices) const
{
  if (genVertices.empty()) return -1;

  // Trace leg 1 to a gen vertex index
  int genVtxIdx1 = -1;
  {
    int origIdx1 = vtxOrigMuonIdx1[vtxIdx];
    bool isDSA1 = !vtxIsDSAMuon1.empty() && vtxIsDSAMuon1[vtxIdx];

    if (isDSA1) {
      // DSA leg: deltaR match to gen muons
      if (origIdx1 >= 0 && origIdx1 < static_cast<int>(dsaMuonPt.size())) {
        int genPartIdx = matchDSAToGen(dsaMuonEta[origIdx1], dsaMuonPhi[origIdx1], dsaMuonPt[origIdx1],
                                        genPt, genEta, genPhi, genPdgId);
        if (genPartIdx >= 0) {
          genVtxIdx1 = findGenVertexForParticle(genPartIdx, genVertices);
        }
      }
    } else {
      // PAT muon leg: use genPartIdx from Muon table
      if (origIdx1 >= 0 && origIdx1 < static_cast<int>(muonGenPartIdx.size())) {
        int gpIdx = muonGenPartIdx[origIdx1];
        if (gpIdx >= 0 && gpIdx < static_cast<int>(genPdgId.size())) {
          if (std::abs(genPdgId[gpIdx]) == 13) {
            genVtxIdx1 = findGenVertexForParticle(gpIdx, genVertices);
          }
        }
      }
    }
  }

  // Trace leg 2 to a gen vertex index
  int genVtxIdx2 = -1;
  {
    int origIdx2 = vtxOrigMuonIdx2[vtxIdx];
    bool isDSA2 = !vtxIsDSAMuon2.empty() && vtxIsDSAMuon2[vtxIdx];

    if (isDSA2) {
      if (origIdx2 >= 0 && origIdx2 < static_cast<int>(dsaMuonPt.size())) {
        int genPartIdx = matchDSAToGen(dsaMuonEta[origIdx2], dsaMuonPhi[origIdx2], dsaMuonPt[origIdx2],
                                        genPt, genEta, genPhi, genPdgId);
        if (genPartIdx >= 0) {
          genVtxIdx2 = findGenVertexForParticle(genPartIdx, genVertices);
        }
      }
    } else {
      if (origIdx2 >= 0 && origIdx2 < static_cast<int>(muonGenPartIdx.size())) {
        int gpIdx = muonGenPartIdx[origIdx2];
        if (gpIdx >= 0 && gpIdx < static_cast<int>(genPdgId.size())) {
          if (std::abs(genPdgId[gpIdx]) == 13) {
            genVtxIdx2 = findGenVertexForParticle(gpIdx, genVertices);
          }
        }
      }
    }
  }

  // Gold = both legs match to gen muons from the SAME gen vertex
  if (genVtxIdx1 >= 0 && genVtxIdx1 == genVtxIdx2)
    return genVtxIdx1;

  return -1;
}

// ============================================================================
// matchDSAToGen - deltaR match DSA muon to gen muons
// ============================================================================

int LLPNanoSVAnalyzer::matchDSAToGen(
  float dsaEta, float dsaPhi, float dsaPt,
  const std::vector<float>& genPt,
  const std::vector<float>& genEta,
  const std::vector<float>& genPhi,
  const std::vector<int>& genPdgId) const
{
  int bestIdx = -1;
  double bestDR = dsaDeltaR_;  // threshold

  for (unsigned int i = 0; i < genPdgId.size(); ++i) {
    if (std::abs(genPdgId[i]) != 13) continue;

    double deta = dsaEta - genEta[i];
    double dphi = dsaPhi - genPhi[i];
    // Normalize dphi to [-pi, pi]
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi < -M_PI) dphi += 2 * M_PI;
    double dr = std::sqrt(deta * deta + dphi * dphi);

    // Check relative pt difference
    double relPtDiff = std::abs(dsaPt - genPt[i]) / std::max(genPt[i], 0.001f);
    if (relPtDiff > dsaRelPtDiff_) continue;

    if (dr < bestDR) {
      bestDR = dr;
      bestIdx = static_cast<int>(i);
    }
  }

  return bestIdx;
}

// ============================================================================
// findGenVertexForParticle - find which gen vertex contains this gen particle
// ============================================================================

int LLPNanoSVAnalyzer::findGenVertexForParticle(int genPartIdx,
                                                  const std::vector<GenMuonVertex>& genVertices) const {
  for (unsigned int iv = 0; iv < genVertices.size(); ++iv) {
    for (int muIdx : genVertices[iv].muonGenPartIndices) {
      if (muIdx == genPartIdx) return static_cast<int>(iv);
    }
  }
  return -1;
}

// ============================================================================
// findNearestGenVertex - spatial matching
// ============================================================================

int LLPNanoSVAnalyzer::findNearestGenVertex(float vx, float vy, float vz,
                                              const std::vector<GenMuonVertex>& genVertices,
                                              float& minDistance) const {
  int minIdx = -1;
  minDistance = -1.f;
  float bestDist = std::numeric_limits<float>::max();

  for (unsigned int iv = 0; iv < genVertices.size(); ++iv) {
    float dx = vx - genVertices[iv].x;
    float dy = vy - genVertices[iv].y;
    float dz = vz - genVertices[iv].z;
    float dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (dist < bestDist) {
      bestDist = dist;
      minIdx = static_cast<int>(iv);
    }
  }

  if (minIdx >= 0)
    minDistance = bestDist;

  return minIdx;
}

// ============================================================================
// fillDescriptions
// ============================================================================

void LLPNanoSVAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("collection", "PatMuonVertex");
  desc.add<int>("motherPdgId", 54);
  desc.add<double>("dsaDeltaR", 0.3);
  desc.add<double>("dsaRelPtDiff", 0.5);
  desc.add<edm::InputTag>("genPartTable", edm::InputTag("genParticleTable"));
  desc.add<edm::InputTag>("muonTable", edm::InputTag("muonTable"));
  desc.add<edm::InputTag>("dsaMuonTable", edm::InputTag("dsaMuonTable"));
  desc.add<edm::InputTag>("vertexTable", edm::InputTag("patMuonVertexTable"));
  desc.add<edm::InputTag>("refittedTracksTable", edm::InputTag("patMuonVertexRefittedTracksTable"));

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LLPNanoSVAnalyzer);
