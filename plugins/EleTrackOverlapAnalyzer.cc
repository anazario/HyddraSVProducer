// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      EleTrackOverlapAnalyzer
//
// Description: Diagnostic analyzer to study the overlap between
//              packedPFCandidates and lostTracks:eleTracks in MiniAOD.
//              Compares each eleLost track against the nearest PF electron
//              candidate (same charge) in (eta, phi) space and fills
//              deltaR and relative pT difference plots.
//              Also checks lostTracks vs PF electrons for reference.
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

  // --- eleLost vs PF electrons ---
  TH1F* h_eleLost_pfEle_minDR_;          // min deltaR to nearest same-charge PF electron
  TH1F* h_eleLost_pfEle_relPt_;          // |pT_eleLost - pT_pfEle| / pT_pfEle, nearest match
  TH2F* h2_eleLost_pfEle_DR_relPt_;      // 2D: deltaR vs relative pT (nearest same-charge PF ele)

  // same plots but opposite charge (cross-check / background estimate)
  TH1F* h_eleLost_pfEle_minDR_opp_;
  TH2F* h2_eleLost_pfEle_DR_relPt_opp_;

  // pT of eleLost tracks with no PF electron match within DR < drNoMatchThreshold_
  TH1F* h_eleLost_unmatched_pt_;

  // --- lost vs PF electrons (reference) ---
  TH1F* h_lost_pfEle_minDR_;
  TH2F* h2_lost_pfEle_DR_relPt_;

  // --- multiplicity ---
  TH1F* h_nPFElectrons_;
  TH1F* h_nEleLost_;
  TH1F* h_nLost_;

  float drNoMatchThreshold_;  // eleLost tracks with min DR above this are "unmatched"
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
      200, 0., 1.0,
      200, 0., 2.0);

  // opposite charge cross-check
  h_eleLost_pfEle_minDR_opp_ = fs->make<TH1F>(
      "h_eleLost_pfEle_minDR_oppCharge",
      ";#DeltaR(eleLost, nearest opposite-charge PF e^{#pm});Entries",
      200, 0., 1.0);

  h2_eleLost_pfEle_DR_relPt_opp_ = fs->make<TH2F>(
      "h2_eleLost_pfEle_DR_relPt_oppCharge",
      ";#DeltaR(eleLost, nearest opp-charge PF e^{#pm});|p_{T}^{eleLost} - p_{T}^{PF e}| / p_{T}^{PF e}",
      200, 0., 1.0,
      200, 0., 2.0);

  // unmatched eleLost pT spectrum
  h_eleLost_unmatched_pt_ = fs->make<TH1F>(
      "h_eleLost_unmatched_pt",
      ";p_{T}^{eleLost} [GeV] (no PF e^{#pm} within #DeltaR < threshold);Entries",
      100, 0., 50.);

  // lost vs PF electrons
  h_lost_pfEle_minDR_ = fs->make<TH1F>(
      "h_lost_pfEle_minDR",
      ";#DeltaR(lostTrack, nearest same-charge PF e^{#pm});Entries",
      200, 0., 1.0);

  h2_lost_pfEle_DR_relPt_ = fs->make<TH2F>(
      "h2_lost_pfEle_DR_relPt",
      ";#DeltaR(lostTrack, nearest same-charge PF e^{#pm});|p_{T}^{lost} - p_{T}^{PF e}| / p_{T}^{PF e}",
      200, 0., 1.0,
      200, 0., 2.0);

  // multiplicity
  h_nPFElectrons_ = fs->make<TH1F>("h_nPFElectrons", ";N(PF e^{#pm});Events", 20, 0., 20.);
  h_nEleLost_     = fs->make<TH1F>("h_nEleLost",     ";N(eleLost tracks);Events",  50, 0., 50.);
  h_nLost_        = fs->make<TH1F>("h_nLost",        ";N(lostTracks);Events",     200, 0., 200.);
}


void EleTrackOverlapAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {

  edm::Handle<pat::PackedCandidateCollection> pfCandsHandle;
  iEvent.getByToken(pfCandsToken_, pfCandsHandle);

  edm::Handle<pat::PackedCandidateCollection> lostTracksHandle;
  iEvent.getByToken(lostTracksToken_, lostTracksHandle);

  edm::Handle<pat::PackedCandidateCollection> eleLostHandle;
  iEvent.getByToken(eleLostTracksToken_, eleLostHandle);

  // Collect PF electron candidates (|pdgId| == 11)
  std::vector<const pat::PackedCandidate*> pfElectrons;
  for (const auto& cand : *pfCandsHandle) {
    if (std::abs(cand.pdgId()) == 11)
      pfElectrons.push_back(&cand);
  }

  h_nPFElectrons_->Fill(pfElectrons.size());
  h_nEleLost_->Fill(eleLostHandle->size());
  h_nLost_->Fill(lostTracksHandle->size());

  const float inf = std::numeric_limits<float>::max();

  // --- eleLost vs PF electrons ---
  for (const auto& eleLost : *eleLostHandle) {

    float minDR_same = inf, relPt_same = -1.f;
    float minDR_opp  = inf, relPt_opp  = -1.f;

    for (const auto* pfEle : pfElectrons) {
      float dr = deltaR(eleLost.eta(), eleLost.phi(), pfEle->eta(), pfEle->phi());
      float rpt = (pfEle->pt() > 0.)
                  ? std::abs(eleLost.pt() - pfEle->pt()) / pfEle->pt()
                  : -1.f;

      if (eleLost.charge() == pfEle->charge()) {
        if (dr < minDR_same) { minDR_same = dr; relPt_same = rpt; }
      } else {
        if (dr < minDR_opp)  { minDR_opp  = dr; relPt_opp  = rpt; }
      }
    }

    if (minDR_same < inf) {
      h_eleLost_pfEle_minDR_->Fill(minDR_same);
      if (relPt_same >= 0.f) {
        h_eleLost_pfEle_relPt_->Fill(relPt_same);
        h2_eleLost_pfEle_DR_relPt_->Fill(minDR_same, relPt_same);
      }
    }

    if (minDR_opp < inf) {
      h_eleLost_pfEle_minDR_opp_->Fill(minDR_opp);
      if (relPt_opp >= 0.f)
        h2_eleLost_pfEle_DR_relPt_opp_->Fill(minDR_opp, relPt_opp);
    }

    // Unmatched: no same-charge PF electron in the event, or nearest one is far away
    if (minDR_same >= drNoMatchThreshold_)
      h_eleLost_unmatched_pt_->Fill(eleLost.pt());
  }

  // --- lostTracks vs PF electrons (reference) ---
  for (const auto& lost : *lostTracksHandle) {
    float minDR_same = inf, relPt_same = -1.f;

    for (const auto* pfEle : pfElectrons) {
      if (lost.charge() != pfEle->charge()) continue;
      float dr = deltaR(lost.eta(), lost.phi(), pfEle->eta(), pfEle->phi());
      if (dr < minDR_same) {
        minDR_same = dr;
        relPt_same = (pfEle->pt() > 0.)
                     ? std::abs(lost.pt() - pfEle->pt()) / pfEle->pt()
                     : -1.f;
      }
    }

    if (minDR_same < inf) {
      h_lost_pfEle_minDR_->Fill(minDR_same);
      if (relPt_same >= 0.f)
        h2_lost_pfEle_DR_relPt_->Fill(minDR_same, relPt_same);
    }
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
