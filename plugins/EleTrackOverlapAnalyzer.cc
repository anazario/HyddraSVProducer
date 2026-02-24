// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      EleTrackOverlapAnalyzer
//
// Description: Diagnostic analyzer to study track overlap between the three
//              MiniAOD track collections: packedPFCandidates, lostTracks, and
//              lostTracks:eleTracks.
//
//              - eleLost vs PF electrons (|pdgId|==11): characterizes the
//                GSF/KF dual-representation overlap; opposite-charge version
//                provided as cross-check.
//              - lostTracks vs all charged PF candidates: lostTracks is not
//                an electron-specific collection so the comparison is made
//                against all charged PF objects.
//              - lostTracks vs eleLost: checks whether the lostTracks-PF
//                overlap population is also present in eleLost.
//
//              Runs directly on MiniAOD — no geometry or B-field needed.
//
// Original Author:  Andres Abreu
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"

#include <cmath>
#include <limits>

namespace {
  float deltaPhi(float phi1, float phi2) {
    float dp = phi1 - phi2;
    while (dp >  M_PI) dp -= 2 * M_PI;
    while (dp < -M_PI) dp += 2 * M_PI;
    return dp;
  }
  float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = deltaPhi(phi1, phi2);
    return std::sqrt(deta * deta + dphi * dphi);
  }
}

class EleTrackOverlapAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit EleTrackOverlapAnalyzer(const edm::ParameterSet&);
  ~EleTrackOverlapAnalyzer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override {}
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> eleLostTracksToken_;

  // --- eleLost vs PF electrons (|pdgId|==11) ---
  TH1F* h_eleLost_pfEle_minDR_;
  TH1F* h_eleLost_pfEle_relPt_;
  TH2F* h2_eleLost_pfEle_DR_relPt_;
  // opposite-charge cross-check
  TH1F* h_eleLost_pfEle_minDR_opp_;
  TH2F* h2_eleLost_pfEle_DR_relPt_opp_;
  // pT of eleLost with no nearby same-charge PF electron
  TH1F* h_eleLost_unmatched_pt_;

  // --- lostTracks vs all charged PF candidates ---
  TH1F* h_lost_pfCharged_minDR_;
  TH2F* h2_lost_pfCharged_DR_relPt_;
  // opposite-charge cross-check
  TH1F* h_lost_pfCharged_minDR_opp_;
  TH2F* h2_lost_pfCharged_DR_relPt_opp_;

  // --- lostTracks vs eleLost ---
  TH1F* h_lost_eleLost_minDR_;
  TH2F* h2_lost_eleLost_DR_relPt_;
  // opposite-charge cross-check
  TH1F* h_lost_eleLost_minDR_opp_;
  TH2F* h2_lost_eleLost_DR_relPt_opp_;

  // --- multiplicity ---
  TH1F* h_nPFElectrons_;
  TH1F* h_nPFCharged_;
  TH1F* h_nEleLost_;
  TH1F* h_nLost_;

  float drNoMatchThreshold_;
};


EleTrackOverlapAnalyzer::EleTrackOverlapAnalyzer(const edm::ParameterSet& iConfig)
    : pfCandsToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      lostTracksToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("lostTracks"))),
      eleLostTracksToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("eleLostTracks"))),
      drNoMatchThreshold_(iConfig.getParameter<double>("drNoMatchThreshold")) {

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  // eleLost vs PF electrons
  h_eleLost_pfEle_minDR_ = fs->make<TH1F>(
      "h_eleLost_pfEle_minDR",
      ";#DeltaR(eleLost, nearest same-charge PF e^{#pm});Entries",
      200, 0., 1.0);
  h_eleLost_pfEle_relPt_ = fs->make<TH1F>(
      "h_eleLost_pfEle_relPt",
      ";|p_{T}^{eleLost} - p_{T}^{PF e}| / p_{T}^{PF e};Entries",
      200, 0., 2.0);
  h2_eleLost_pfEle_DR_relPt_ = fs->make<TH2F>(
      "h2_eleLost_pfEle_DR_relPt",
      ";#DeltaR(eleLost, nearest same-charge PF e^{#pm});|p_{T}^{eleLost} - p_{T}^{PF e}| / p_{T}^{PF e}",
      200, 0., 1.0, 200, 0., 2.0);
  h_eleLost_pfEle_minDR_opp_ = fs->make<TH1F>(
      "h_eleLost_pfEle_minDR_oppCharge",
      ";#DeltaR(eleLost, nearest opp-charge PF e^{#pm});Entries",
      200, 0., 1.0);
  h2_eleLost_pfEle_DR_relPt_opp_ = fs->make<TH2F>(
      "h2_eleLost_pfEle_DR_relPt_oppCharge",
      ";#DeltaR(eleLost, nearest opp-charge PF e^{#pm});|p_{T}^{eleLost} - p_{T}^{PF e}| / p_{T}^{PF e}",
      200, 0., 1.0, 200, 0., 2.0);
  h_eleLost_unmatched_pt_ = fs->make<TH1F>(
      "h_eleLost_unmatched_pt",
      ";p_{T}^{eleLost} [GeV] (no same-charge PF e^{#pm} within #DeltaR < threshold);Entries",
      100, 0., 50.);

  // lostTracks vs all charged PF candidates
  h_lost_pfCharged_minDR_ = fs->make<TH1F>(
      "h_lost_pfCharged_minDR",
      ";#DeltaR(lostTrack, nearest same-charge PF candidate);Entries",
      200, 0., 1.0);
  h2_lost_pfCharged_DR_relPt_ = fs->make<TH2F>(
      "h2_lost_pfCharged_DR_relPt",
      ";#DeltaR(lostTrack, nearest same-charge PF candidate);|p_{T}^{lost} - p_{T}^{PF}| / p_{T}^{PF}",
      200, 0., 1.0, 200, 0., 2.0);
  h_lost_pfCharged_minDR_opp_ = fs->make<TH1F>(
      "h_lost_pfCharged_minDR_oppCharge",
      ";#DeltaR(lostTrack, nearest opp-charge PF candidate);Entries",
      200, 0., 1.0);
  h2_lost_pfCharged_DR_relPt_opp_ = fs->make<TH2F>(
      "h2_lost_pfCharged_DR_relPt_oppCharge",
      ";#DeltaR(lostTrack, nearest opp-charge PF candidate);|p_{T}^{lost} - p_{T}^{PF}| / p_{T}^{PF}",
      200, 0., 1.0, 200, 0., 2.0);

  // lostTracks vs eleLost
  h_lost_eleLost_minDR_ = fs->make<TH1F>(
      "h_lost_eleLost_minDR",
      ";#DeltaR(lostTrack, nearest same-charge eleLost);Entries",
      200, 0., 1.0);
  h2_lost_eleLost_DR_relPt_ = fs->make<TH2F>(
      "h2_lost_eleLost_DR_relPt",
      ";#DeltaR(lostTrack, nearest same-charge eleLost);|p_{T}^{lost} - p_{T}^{eleLost}| / p_{T}^{eleLost}",
      200, 0., 1.0, 200, 0., 2.0);
  h_lost_eleLost_minDR_opp_ = fs->make<TH1F>(
      "h_lost_eleLost_minDR_oppCharge",
      ";#DeltaR(lostTrack, nearest opp-charge eleLost);Entries",
      200, 0., 1.0);
  h2_lost_eleLost_DR_relPt_opp_ = fs->make<TH2F>(
      "h2_lost_eleLost_DR_relPt_oppCharge",
      ";#DeltaR(lostTrack, nearest opp-charge eleLost);|p_{T}^{lost} - p_{T}^{eleLost}| / p_{T}^{eleLost}",
      200, 0., 1.0, 200, 0., 2.0);

  // multiplicity
  h_nPFElectrons_ = fs->make<TH1F>("h_nPFElectrons", ";N(PF e^{#pm});Events",       20,  0.,  20.);
  h_nPFCharged_   = fs->make<TH1F>("h_nPFCharged",   ";N(charged PF cands);Events", 200, 0., 200.);
  h_nEleLost_     = fs->make<TH1F>("h_nEleLost",     ";N(eleLost tracks);Events",    50, 0.,  50.);
  h_nLost_        = fs->make<TH1F>("h_nLost",        ";N(lostTracks);Events",        200, 0., 200.);
}


void EleTrackOverlapAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {

  edm::Handle<pat::PackedCandidateCollection> pfCandsHandle;
  iEvent.getByToken(pfCandsToken_, pfCandsHandle);

  edm::Handle<pat::PackedCandidateCollection> lostTracksHandle;
  iEvent.getByToken(lostTracksToken_, lostTracksHandle);

  edm::Handle<pat::PackedCandidateCollection> eleLostHandle;
  iEvent.getByToken(eleLostTracksToken_, eleLostHandle);

  // Separate PF into electrons and all charged
  std::vector<const pat::PackedCandidate*> pfElectrons, pfCharged;
  for (const auto& cand : *pfCandsHandle) {
    if (cand.charge() == 0) continue;
    pfCharged.push_back(&cand);
    if (std::abs(cand.pdgId()) == 11)
      pfElectrons.push_back(&cand);
  }

  h_nPFElectrons_->Fill(pfElectrons.size());
  h_nPFCharged_->Fill(pfCharged.size());
  h_nEleLost_->Fill(eleLostHandle->size());
  h_nLost_->Fill(lostTracksHandle->size());

  const float inf = std::numeric_limits<float>::max();

  // For each probe track, find nearest same- and opposite-charge match and fill histograms
  auto fillPair = [&](float eta, float phi, float pt, int charge,
                      const std::vector<const pat::PackedCandidate*>& coll,
                      TH1F* hDR_same, TH2F* h2_same,
                      TH1F* hDR_opp,  TH2F* h2_opp) {
    float minDR_same = inf, relPt_same = -1.f;
    float minDR_opp  = inf, relPt_opp  = -1.f;

    for (const auto* c : coll) {
      float dr  = deltaR(eta, phi, c->eta(), c->phi());
      float rpt = (c->pt() > 0.f) ? std::abs(pt - c->pt()) / c->pt() : -1.f;
      if (c->charge() == charge) {
        if (dr < minDR_same) { minDR_same = dr; relPt_same = rpt; }
      } else {
        if (dr < minDR_opp)  { minDR_opp  = dr; relPt_opp  = rpt; }
      }
    }

    if (minDR_same < inf) {
      hDR_same->Fill(minDR_same);
      if (relPt_same >= 0.f) h2_same->Fill(minDR_same, relPt_same);
    }
    if (minDR_opp < inf) {
      hDR_opp->Fill(minDR_opp);
      if (relPt_opp >= 0.f) h2_opp->Fill(minDR_opp, relPt_opp);
    }
    return minDR_same;
  };

  // --- eleLost vs PF electrons ---
  for (const auto& eleLost : *eleLostHandle) {
    float minDR_same = fillPair(eleLost.eta(), eleLost.phi(), eleLost.pt(), eleLost.charge(),
                                pfElectrons,
                                h_eleLost_pfEle_minDR_,    h2_eleLost_pfEle_DR_relPt_,
                                h_eleLost_pfEle_minDR_opp_, h2_eleLost_pfEle_DR_relPt_opp_);
    if (minDR_same >= drNoMatchThreshold_)
      h_eleLost_unmatched_pt_->Fill(eleLost.pt());
  }

  // --- lostTracks vs all charged PF candidates ---
  for (const auto& lost : *lostTracksHandle) {
    fillPair(lost.eta(), lost.phi(), lost.pt(), lost.charge(),
             pfCharged,
             h_lost_pfCharged_minDR_,     h2_lost_pfCharged_DR_relPt_,
             h_lost_pfCharged_minDR_opp_, h2_lost_pfCharged_DR_relPt_opp_);
  }

  // --- lostTracks vs eleLost ---
  std::vector<const pat::PackedCandidate*> eleLostVec;
  for (const auto& c : *eleLostHandle) eleLostVec.push_back(&c);

  for (const auto& lost : *lostTracksHandle) {
    fillPair(lost.eta(), lost.phi(), lost.pt(), lost.charge(),
             eleLostVec,
             h_lost_eleLost_minDR_,     h2_lost_eleLost_DR_relPt_,
             h_lost_eleLost_minDR_opp_, h2_lost_eleLost_DR_relPt_opp_);
  }
}


void EleTrackOverlapAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfCandidates",  edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("lostTracks",    edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("eleLostTracks", edm::InputTag("lostTracks", "eleTracks"));
  desc.add<double>("drNoMatchThreshold", 0.4);
  descriptions.add("eleTrackOverlapAnalyzer", desc);
}

DEFINE_FWK_MODULE(EleTrackOverlapAnalyzer);
