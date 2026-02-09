// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      LLPNanoSVAnalyzer
//
// Description: EDAnalyzer that reads LLPNanoAOD flat ROOT files directly
//              (Events TTree with branches like GenPart_pdgId,
//              PatMuonVertex_vx, etc.) and writes a flat TTree with the same
//              branch structure as HyddraSVAnalyzer.  This allows direct
//              comparison of SV reconstruction efficiency between HYDDRA and
//              the LLPNanoAOD dimuon vertex tools using the existing
//              efficiency plotting scripts.
//
//              Uses EmptySource + manual TTree reading because LLPNanoAOD
//              files may not be PoolSource-compatible across CMSSW versions.
//
// Original Author:  Andres Abreu
//

#include <cmath>
#include <limits>
#include <map>
#include <vector>
#include <string>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

// Maximum array sizes for NanoAOD variable-length branches
static constexpr unsigned int MAX_GENPART = 1000;
static constexpr unsigned int MAX_MUON = 100;
static constexpr unsigned int MAX_DSAMUON = 100;
static constexpr unsigned int MAX_VTX = 2000;
static constexpr unsigned int MAX_REFTRACKS = 4000;

class LLPNanoSVAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit LLPNanoSVAnalyzer(const edm::ParameterSet&);
  ~LLPNanoSVAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  void clearOutputBranches();
  void processFile(const std::string& filename);
  void processEvent();

  // Gen vertex reconstruction
  struct GenMuonVertex {
    float x, y, z, dxy;
    float pt, eta, phi, mass;
    int motherIndex;
    std::vector<int> muonGenPartIndices;
  };

  std::vector<GenMuonVertex> buildGenVertices() const;

  // Gold matching
  int findGoldMatch(unsigned int vtxIdx,
                    const std::vector<GenMuonVertex>& genVertices) const;

  // DSA-to-gen deltaR match
  int matchDSAToGen(float dsaEta, float dsaPhi, float dsaPt) const;

  // Find which gen vertex a gen particle belongs to
  int findGenVertexForParticle(int genPartIdx,
                               const std::vector<GenMuonVertex>& genVertices) const;

  // Spatial matching
  int findNearestGenVertex(float vx, float vy, float vz,
                           const std::vector<GenMuonVertex>& genVertices,
                           float& minDistance) const;

  // Configuration
  std::string collection_;
  int motherPdgId_;
  double dsaDeltaR_;
  double dsaRelPtDiff_;
  std::vector<std::string> inputFiles_;

  // ---- Input NanoAOD branch buffers (filled by SetBranchAddress) ----

  // GenPart
  Int_t in_nGenPart_;
  Int_t in_GenPart_pdgId_[MAX_GENPART];
  Short_t in_GenPart_genPartIdxMother_[MAX_GENPART];
  Float_t in_GenPart_vx_[MAX_GENPART];
  Float_t in_GenPart_vy_[MAX_GENPART];
  Float_t in_GenPart_vz_[MAX_GENPART];
  Float_t in_GenPart_pt_[MAX_GENPART];
  Float_t in_GenPart_eta_[MAX_GENPART];
  Float_t in_GenPart_phi_[MAX_GENPART];
  Float_t in_GenPart_mass_[MAX_GENPART];
  Int_t in_GenPart_status_[MAX_GENPART];

  // Muon
  Int_t in_nMuon_;
  Short_t in_Muon_genPartIdx_[MAX_MUON];
  Float_t in_Muon_pt_[MAX_MUON];
  Float_t in_Muon_eta_[MAX_MUON];
  Float_t in_Muon_phi_[MAX_MUON];

  // DSAMuon
  Int_t in_nDSAMuon_;
  Float_t in_DSAMuon_pt_[MAX_DSAMUON];
  Float_t in_DSAMuon_eta_[MAX_DSAMUON];
  Float_t in_DSAMuon_phi_[MAX_DSAMUON];

  // Vertex (PatMuonVertex or PatDSAMuonVertex)
  Int_t in_nVtx_;
  Float_t in_Vtx_isValid_[MAX_VTX];
  Float_t in_Vtx_vx_[MAX_VTX];
  Float_t in_Vtx_vy_[MAX_VTX];
  Float_t in_Vtx_vz_[MAX_VTX];
  Float_t in_Vtx_vxy_[MAX_VTX];
  Float_t in_Vtx_chi2_[MAX_VTX];
  Float_t in_Vtx_normChi2_[MAX_VTX];
  Float_t in_Vtx_ndof_[MAX_VTX];
  Float_t in_Vtx_originalMuonIdx1_[MAX_VTX];
  Float_t in_Vtx_originalMuonIdx2_[MAX_VTX];
  Float_t in_Vtx_isDSAMuon1_[MAX_VTX];
  Float_t in_Vtx_isDSAMuon2_[MAX_VTX];

  // RefittedTracks
  Int_t in_nRefTracks_;
  Float_t in_RefTracks_px_[MAX_REFTRACKS];
  Float_t in_RefTracks_py_[MAX_REFTRACKS];
  Float_t in_RefTracks_pz_[MAX_REFTRACKS];
  Float_t in_RefTracks_idx_[MAX_REFTRACKS];

  // ---- Output TTree and branches ----
  TTree* tree_;

  unsigned int nVertices_;
  std::vector<bool> isLeptonic_;
  std::vector<unsigned int> nTracks_;
  std::vector<float> x_, y_, z_, dxy_;
  std::vector<float> mass_, pt_, eta_, phi_;
  std::vector<float> chi2_, normalizedChi2_, ndof_;

  std::vector<int> genVertexIndex_;
  std::vector<int> nearestGenVertexIndex_;
  std::vector<float> min3D_;
  std::vector<bool> isBronze_, isSilver_, isGold_;

  unsigned int genVertex_nTotal_;
  unsigned int genVertex_nMuon_;
  unsigned int genVertex_nElectron_;
  unsigned int genVertex_nHadronic_;
  std::vector<float> genVertex_dxy_;
  std::vector<float> genVertex_x_, genVertex_y_, genVertex_z_;
  std::vector<float> genVertex_pt_, genVertex_eta_, genVertex_phi_, genVertex_mass_;
  std::vector<bool> genVertex_isMuon_, genVertex_isElectron_, genVertex_isHadronic_;
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
  inputFiles_(iConfig.getParameter<std::vector<std::string>>("inputFiles"))
{
  usesResource("TFileService");
}

// ============================================================================
// beginJob
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
// analyze - called once by EmptySource; loops over all input files internally
// ============================================================================

void LLPNanoSVAnalyzer::analyze(const edm::Event&, const edm::EventSetup&) {
  for (const auto& filename : inputFiles_) {
    processFile(filename);
  }
}

// ============================================================================
// processFile - open one NanoAOD file and loop over Events TTree entries
// ============================================================================

void LLPNanoSVAnalyzer::processFile(const std::string& filename) {

  TFile* inFile = TFile::Open(filename.c_str(), "READ");
  if (!inFile || inFile->IsZombie()) {
    edm::LogWarning("LLPNanoSVAnalyzer") << "Cannot open file: " << filename;
    return;
  }

  TTree* events = dynamic_cast<TTree*>(inFile->Get("Events"));
  if (!events) {
    edm::LogWarning("LLPNanoSVAnalyzer") << "No Events tree in: " << filename;
    inFile->Close();
    delete inFile;
    return;
  }

  // Build branch prefix from collection name
  std::string vtxPrefix = collection_ + "_";     // e.g. "PatMuonVertex_"
  std::string refPrefix = collection_ + "RefittedTracks_";
  std::string nVtxBranch = "n" + collection_;     // e.g. "nPatMuonVertex"
  std::string nRefBranch = "n" + collection_ + "RefittedTracks";

  // ---- Set branch addresses: GenPart ----
  events->SetBranchAddress("nGenPart", &in_nGenPart_);
  events->SetBranchAddress("GenPart_pdgId", in_GenPart_pdgId_);
  events->SetBranchAddress("GenPart_genPartIdxMother", in_GenPart_genPartIdxMother_);
  events->SetBranchAddress("GenPart_vx", in_GenPart_vx_);
  events->SetBranchAddress("GenPart_vy", in_GenPart_vy_);
  events->SetBranchAddress("GenPart_vz", in_GenPart_vz_);
  events->SetBranchAddress("GenPart_pt", in_GenPart_pt_);
  events->SetBranchAddress("GenPart_eta", in_GenPart_eta_);
  events->SetBranchAddress("GenPart_phi", in_GenPart_phi_);
  events->SetBranchAddress("GenPart_mass", in_GenPart_mass_);
  events->SetBranchAddress("GenPart_status", in_GenPart_status_);

  // ---- Muon ----
  events->SetBranchAddress("nMuon", &in_nMuon_);
  events->SetBranchAddress("Muon_genPartIdx", in_Muon_genPartIdx_);
  events->SetBranchAddress("Muon_pt", in_Muon_pt_);
  events->SetBranchAddress("Muon_eta", in_Muon_eta_);
  events->SetBranchAddress("Muon_phi", in_Muon_phi_);

  // ---- DSAMuon ----
  events->SetBranchAddress("nDSAMuon", &in_nDSAMuon_);
  events->SetBranchAddress("DSAMuon_pt", in_DSAMuon_pt_);
  events->SetBranchAddress("DSAMuon_eta", in_DSAMuon_eta_);
  events->SetBranchAddress("DSAMuon_phi", in_DSAMuon_phi_);

  // ---- Vertex table ----
  events->SetBranchAddress(nVtxBranch.c_str(), &in_nVtx_);
  events->SetBranchAddress((vtxPrefix + "isValid").c_str(), in_Vtx_isValid_);
  events->SetBranchAddress((vtxPrefix + "vx").c_str(), in_Vtx_vx_);
  events->SetBranchAddress((vtxPrefix + "vy").c_str(), in_Vtx_vy_);
  events->SetBranchAddress((vtxPrefix + "vz").c_str(), in_Vtx_vz_);
  events->SetBranchAddress((vtxPrefix + "vxy").c_str(), in_Vtx_vxy_);
  events->SetBranchAddress((vtxPrefix + "chi2").c_str(), in_Vtx_chi2_);
  events->SetBranchAddress((vtxPrefix + "normChi2").c_str(), in_Vtx_normChi2_);
  events->SetBranchAddress((vtxPrefix + "ndof").c_str(), in_Vtx_ndof_);
  events->SetBranchAddress((vtxPrefix + "originalMuonIdx1").c_str(), in_Vtx_originalMuonIdx1_);
  events->SetBranchAddress((vtxPrefix + "originalMuonIdx2").c_str(), in_Vtx_originalMuonIdx2_);

  // isDSAMuon1/2 exist for both PatMuonVertex and PatDSAMuonVertex
  events->SetBranchAddress((vtxPrefix + "isDSAMuon1").c_str(), in_Vtx_isDSAMuon1_);
  events->SetBranchAddress((vtxPrefix + "isDSAMuon2").c_str(), in_Vtx_isDSAMuon2_);

  // ---- RefittedTracks ----
  events->SetBranchAddress(nRefBranch.c_str(), &in_nRefTracks_);
  events->SetBranchAddress((refPrefix + "px").c_str(), in_RefTracks_px_);
  events->SetBranchAddress((refPrefix + "py").c_str(), in_RefTracks_py_);
  events->SetBranchAddress((refPrefix + "pz").c_str(), in_RefTracks_pz_);
  events->SetBranchAddress((refPrefix + "idx").c_str(), in_RefTracks_idx_);

  Long64_t nEntries = events->GetEntries();
  edm::LogInfo("LLPNanoSVAnalyzer") << "Processing " << filename
                                     << " (" << nEntries << " events)";

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    events->GetEntry(entry);
    processEvent();
  }

  inFile->Close();
  delete inFile;
}

// ============================================================================
// processEvent - process one NanoAOD event (buffers already filled by GetEntry)
// ============================================================================

void LLPNanoSVAnalyzer::processEvent() {

  clearOutputBranches();

  // ---- Build gen muon vertices ----
  auto genVertices = buildGenVertices();

  // ---- Fill gen vertex branches ----
  genVertex_nTotal_ = static_cast<unsigned int>(genVertices.size());
  genVertex_nMuon_ = genVertex_nTotal_;
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

  for (int iv = 0; iv < in_nVtx_; ++iv) {
    if (in_Vtx_isValid_[iv] < 0.5f) continue;

    nValid++;

    isLeptonic_.push_back(true);
    nTracks_.push_back(2);

    x_.push_back(in_Vtx_vx_[iv]);
    y_.push_back(in_Vtx_vy_[iv]);
    z_.push_back(in_Vtx_vz_[iv]);
    dxy_.push_back(in_Vtx_vxy_[iv]);
    chi2_.push_back(in_Vtx_chi2_[iv]);
    normalizedChi2_.push_back(in_Vtx_normChi2_[iv]);
    ndof_.push_back(in_Vtx_ndof_[iv]);

    // Compute vertex mass/pt/eta/phi from refitted tracks
    TLorentzVector p4sum;
    bool foundTracks = false;
    for (int it = 0; it < in_nRefTracks_; ++it) {
      if (static_cast<int>(in_RefTracks_idx_[it]) == iv) {
        double energy = std::sqrt(in_RefTracks_px_[it] * in_RefTracks_px_[it] +
                                  in_RefTracks_py_[it] * in_RefTracks_py_[it] +
                                  in_RefTracks_pz_[it] * in_RefTracks_pz_[it] +
                                  0.10566 * 0.10566);
        TLorentzVector p4;
        p4.SetPxPyPzE(in_RefTracks_px_[it], in_RefTracks_py_[it],
                       in_RefTracks_pz_[it], energy);
        p4sum += p4;
        foundTracks = true;
      }
    }

    mass_.push_back(foundTracks ? static_cast<float>(p4sum.M()) : 0.f);
    pt_.push_back(foundTracks ? static_cast<float>(p4sum.Pt()) : 0.f);
    eta_.push_back(foundTracks ? static_cast<float>(p4sum.Eta()) : 0.f);
    phi_.push_back(foundTracks ? static_cast<float>(p4sum.Phi()) : 0.f);

    // ---- Gold matching ----
    int goldGenVtxIdx = findGoldMatch(iv, genVertices);
    bool isGoldMatch = (goldGenVtxIdx >= 0);

    genVertexIndex_.push_back(goldGenVtxIdx);
    isGold_.push_back(isGoldMatch);
    isSilver_.push_back(isGoldMatch);
    isBronze_.push_back(isGoldMatch);

    // Spatial matching
    float minDist = -1.f;
    int nearestIdx = findNearestGenVertex(in_Vtx_vx_[iv], in_Vtx_vy_[iv],
                                           in_Vtx_vz_[iv], genVertices, minDist);
    nearestGenVertexIndex_.push_back(nearestIdx);
    min3D_.push_back(minDist);
  }

  nVertices_ = nValid;
  tree_->Fill();
}

// ============================================================================
// clearOutputBranches
// ============================================================================

void LLPNanoSVAnalyzer::clearOutputBranches() {
  nVertices_ = 0;
  isLeptonic_.clear(); nTracks_.clear();
  x_.clear(); y_.clear(); z_.clear(); dxy_.clear();
  mass_.clear(); pt_.clear(); eta_.clear(); phi_.clear();
  chi2_.clear(); normalizedChi2_.clear(); ndof_.clear();

  genVertexIndex_.clear(); nearestGenVertexIndex_.clear(); min3D_.clear();
  isBronze_.clear(); isSilver_.clear(); isGold_.clear();

  genVertex_nTotal_ = 0; genVertex_nMuon_ = 0;
  genVertex_nElectron_ = 0; genVertex_nHadronic_ = 0;
  genVertex_dxy_.clear();
  genVertex_x_.clear(); genVertex_y_.clear(); genVertex_z_.clear();
  genVertex_pt_.clear(); genVertex_eta_.clear();
  genVertex_phi_.clear(); genVertex_mass_.clear();
  genVertex_isMuon_.clear(); genVertex_isElectron_.clear();
  genVertex_isHadronic_.clear(); genVertex_passSelection_.clear();
}

// ============================================================================
// buildGenVertices - uses in_GenPart_* buffers
// ============================================================================

std::vector<LLPNanoSVAnalyzer::GenMuonVertex>
LLPNanoSVAnalyzer::buildGenVertices() const {

  std::map<int, std::vector<int>> motherToMuons;

  for (int i = 0; i < in_nGenPart_; ++i) {
    if (std::abs(in_GenPart_pdgId_[i]) != 13) continue;
    if (in_GenPart_status_[i] != 1) continue;

    int motherIdx = static_cast<int>(in_GenPart_genPartIdxMother_[i]);
    if (motherIdx < 0 || motherIdx >= in_nGenPart_) continue;

    // Walk up the chain to find the signal mother
    int currentMother = motherIdx;
    bool foundSignalMother = false;

    for (int level = 0; level < 10; ++level) {
      if (currentMother < 0 || currentMother >= in_nGenPart_)
        break;

      if (in_GenPart_pdgId_[currentMother] == motherPdgId_) {
        foundSignalMother = true;
        break;
      }

      // If mother is same-pdgId muon (status copy), keep going up
      if (std::abs(in_GenPart_pdgId_[currentMother]) == 13) {
        currentMother = static_cast<int>(in_GenPart_genPartIdxMother_[currentMother]);
        continue;
      }
      break;
    }

    if (foundSignalMother) {
      motherToMuons[currentMother].push_back(static_cast<int>(i));
    }
  }

  std::vector<GenMuonVertex> vertices;

  for (const auto& pair : motherToMuons) {
    const auto& muonIndices = pair.second;
    if (muonIndices.size() < 2) continue;

    GenMuonVertex gv;
    gv.motherIndex = pair.first;
    gv.muonGenPartIndices = muonIndices;

    int firstMuon = muonIndices[0];
    gv.x = in_GenPart_vx_[firstMuon];
    gv.y = in_GenPart_vy_[firstMuon];
    gv.z = in_GenPart_vz_[firstMuon];
    gv.dxy = std::sqrt(gv.x * gv.x + gv.y * gv.y);

    TLorentzVector p4sum;
    for (int idx : muonIndices) {
      TLorentzVector p4;
      p4.SetPtEtaPhiM(in_GenPart_pt_[idx], in_GenPart_eta_[idx],
                       in_GenPart_phi_[idx], in_GenPart_mass_[idx]);
      p4sum += p4;
    }

    gv.pt   = static_cast<float>(p4sum.Pt());
    gv.eta  = static_cast<float>(p4sum.Eta());
    gv.phi  = static_cast<float>(p4sum.Phi());
    gv.mass = static_cast<float>(p4sum.M());

    vertices.push_back(gv);
  }

  return vertices;
}

// ============================================================================
// findGoldMatch
// ============================================================================

int LLPNanoSVAnalyzer::findGoldMatch(
  unsigned int vtxIdx,
  const std::vector<GenMuonVertex>& genVertices) const
{
  if (genVertices.empty()) return -1;

  // -- Trace leg 1 --
  int genVtxIdx1 = -1;
  {
    int origIdx = static_cast<int>(in_Vtx_originalMuonIdx1_[vtxIdx]);
    bool isDSA = in_Vtx_isDSAMuon1_[vtxIdx] > 0.5f;

    if (isDSA) {
      if (origIdx >= 0 && origIdx < in_nDSAMuon_) {
        int gpIdx = matchDSAToGen(in_DSAMuon_eta_[origIdx],
                                   in_DSAMuon_phi_[origIdx],
                                   in_DSAMuon_pt_[origIdx]);
        if (gpIdx >= 0)
          genVtxIdx1 = findGenVertexForParticle(gpIdx, genVertices);
      }
    } else {
      if (origIdx >= 0 && origIdx < in_nMuon_) {
        int gpIdx = static_cast<int>(in_Muon_genPartIdx_[origIdx]);
        if (gpIdx >= 0 && gpIdx < in_nGenPart_ &&
            std::abs(in_GenPart_pdgId_[gpIdx]) == 13) {
          genVtxIdx1 = findGenVertexForParticle(gpIdx, genVertices);
        }
      }
    }
  }

  // -- Trace leg 2 --
  int genVtxIdx2 = -1;
  {
    int origIdx = static_cast<int>(in_Vtx_originalMuonIdx2_[vtxIdx]);
    bool isDSA = in_Vtx_isDSAMuon2_[vtxIdx] > 0.5f;

    if (isDSA) {
      if (origIdx >= 0 && origIdx < in_nDSAMuon_) {
        int gpIdx = matchDSAToGen(in_DSAMuon_eta_[origIdx],
                                   in_DSAMuon_phi_[origIdx],
                                   in_DSAMuon_pt_[origIdx]);
        if (gpIdx >= 0)
          genVtxIdx2 = findGenVertexForParticle(gpIdx, genVertices);
      }
    } else {
      if (origIdx >= 0 && origIdx < in_nMuon_) {
        int gpIdx = static_cast<int>(in_Muon_genPartIdx_[origIdx]);
        if (gpIdx >= 0 && gpIdx < in_nGenPart_ &&
            std::abs(in_GenPart_pdgId_[gpIdx]) == 13) {
          genVtxIdx2 = findGenVertexForParticle(gpIdx, genVertices);
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
// matchDSAToGen
// ============================================================================

int LLPNanoSVAnalyzer::matchDSAToGen(float dsaEta, float dsaPhi,
                                       float dsaPt) const {
  int bestIdx = -1;
  double bestDR = dsaDeltaR_;

  for (int i = 0; i < in_nGenPart_; ++i) {
    if (std::abs(in_GenPart_pdgId_[i]) != 13) continue;

    double deta = dsaEta - in_GenPart_eta_[i];
    double dphi = dsaPhi - in_GenPart_phi_[i];
    while (dphi > M_PI)  dphi -= 2 * M_PI;
    while (dphi < -M_PI) dphi += 2 * M_PI;
    double dr = std::sqrt(deta * deta + dphi * dphi);

    double relPtDiff = std::abs(dsaPt - in_GenPart_pt_[i]) /
                       std::max(double(in_GenPart_pt_[i]), 0.001);
    if (relPtDiff > dsaRelPtDiff_) continue;

    if (dr < bestDR) {
      bestDR = dr;
      bestIdx = static_cast<int>(i);
    }
  }

  return bestIdx;
}

// ============================================================================
// findGenVertexForParticle
// ============================================================================

int LLPNanoSVAnalyzer::findGenVertexForParticle(
  int genPartIdx,
  const std::vector<GenMuonVertex>& genVertices) const
{
  for (unsigned int iv = 0; iv < genVertices.size(); ++iv) {
    for (int muIdx : genVertices[iv].muonGenPartIndices) {
      if (muIdx == genPartIdx) return static_cast<int>(iv);
    }
  }
  return -1;
}

// ============================================================================
// findNearestGenVertex
// ============================================================================

int LLPNanoSVAnalyzer::findNearestGenVertex(
  float vx, float vy, float vz,
  const std::vector<GenMuonVertex>& genVertices,
  float& minDistance) const
{
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
  desc.add<std::vector<std::string>>("inputFiles", {});

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LLPNanoSVAnalyzer);
