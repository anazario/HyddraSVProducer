// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      MiniAODElectronTrackProducer
//
// Description: Extracts reco::Track collections from GED and low-pT electrons
//              in MiniAOD for use with HYDDRA displaced vertex reconstruction.
//
//              Implements the electron selection from the iDM analysis (AN-21-XX):
//                - GED electrons:   MVA90WP ID (configurable), pT > minPt, |eta| < maxAbsEta
//                - Low-pT electrons: IDmva > lowPtMVAThreshold, pT > minPt, |eta| < maxAbsEta
//                - Cross-cleaning:  low-pT electrons within deltaR < crossCleanDR of any
//                                   selected GED electron are removed
//                - HEM veto:        disabled by default (data artifact, not applicable to MC)
//
//              GSF tracks are sliced to reco::Track (base-class copy) for compatibility
//              with the HyddraSVsEXOProducer reco::TrackCollection interface.
//
//              Three output collections:
//                "gedElectronTracks"    - tracks from selected GED electrons
//                "lowPtElectronTracks"  - tracks from selected low-pT electrons (cross-cleaned)
//                "mergedElectronTracks" - union of the above (GED first, then low-pT)
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <vector>
#include <cmath>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

class MiniAODElectronTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit MiniAODElectronTrackProducer(const edm::ParameterSet&);
  ~MiniAODElectronTrackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  bool passGEDSelection(const pat::Electron& ele) const;
  bool passLowPtSelection(const pat::Electron& ele) const;
  bool inHEMRegion(const pat::Electron& ele) const;
  double deltaR(const pat::Electron& a, const pat::Electron& b) const;

  // Input tokens
  edm::EDGetTokenT<std::vector<pat::Electron>> gedElectronsToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> lowPtElectronsToken_;

  // Configuration
  double minPt_;
  double maxAbsEta_;
  bool   applyHEMVeto_;
  std::string gedMVALabel_;       // electronID label for GED MVA WP; "" disables cut
  std::string lowPtMVALabel_;     // userFloat key for low-pT MVA score; "" disables cut
  double lowPtMVAThreshold_;      // IDmva > this value required
  double crossCleanDR_;           // low-pT vs GED deltaR cross-cleaning threshold
};

// ---------------------------------------------------------------------------
MiniAODElectronTrackProducer::MiniAODElectronTrackProducer(const edm::ParameterSet& iConfig)
    : gedElectronsToken_(consumes<std::vector<pat::Electron>>(
          iConfig.getParameter<edm::InputTag>("gedElectrons"))),
      lowPtElectronsToken_(consumes<std::vector<pat::Electron>>(
          iConfig.getParameter<edm::InputTag>("lowPtElectrons"))),
      minPt_(iConfig.getParameter<double>("minPt")),
      maxAbsEta_(iConfig.getParameter<double>("maxAbsEta")),
      applyHEMVeto_(iConfig.getParameter<bool>("applyHEMVeto")),
      gedMVALabel_(iConfig.getParameter<std::string>("gedMVALabel")),
      lowPtMVALabel_(iConfig.getParameter<std::string>("lowPtMVALabel")),
      lowPtMVAThreshold_(iConfig.getParameter<double>("lowPtMVAThreshold")),
      crossCleanDR_(iConfig.getParameter<double>("crossCleanDR"))
{
  produces<reco::TrackCollection>("gedElectronTracks");
  produces<reco::TrackCollection>("lowPtElectronTracks");
  produces<reco::TrackCollection>("mergedElectronTracks");
}

// ---------------------------------------------------------------------------
bool MiniAODElectronTrackProducer::inHEMRegion(const pat::Electron& ele) const {
  // ECAL HEM 15/16 failure region (2018 data only)
  return (ele.eta() > -3.0 && ele.eta() < -1.3 &&
          ele.phi() > -1.57 && ele.phi() < -0.87);
}

bool MiniAODElectronTrackProducer::passGEDSelection(const pat::Electron& ele) const {
  if (ele.pt() < minPt_) return false;
  if (std::fabs(ele.eta()) > maxAbsEta_) return false;
  if (applyHEMVeto_ && inHEMRegion(ele)) return false;

  if (!gedMVALabel_.empty()) {
    // Reject if ID label not embedded or electron fails WP (electronID returns 1.0=pass, 0.0=fail)
    if (!ele.isElectronIDAvailable(gedMVALabel_) || ele.electronID(gedMVALabel_) < 0.5f)
      return false;
  }
  return true;
}

bool MiniAODElectronTrackProducer::passLowPtSelection(const pat::Electron& ele) const {
  if (ele.pt() < minPt_) return false;
  if (std::fabs(ele.eta()) > maxAbsEta_) return false;
  if (applyHEMVeto_ && inHEMRegion(ele)) return false;

  if (!lowPtMVALabel_.empty()) {
    // Scores are stored via electronID(), not userFloat(), in Run2 UL MiniAOD
    if (!ele.isElectronIDAvailable(lowPtMVALabel_) ||
        ele.electronID(lowPtMVALabel_) <= static_cast<float>(lowPtMVAThreshold_))
      return false;
  }
  return true;
}

double MiniAODElectronTrackProducer::deltaR(const pat::Electron& a, const pat::Electron& b) const {
  const double dEta = a.eta() - b.eta();
  double dPhi = std::fabs(a.phi() - b.phi());
  if (dPhi > M_PI) dPhi = 2.0 * M_PI - dPhi;
  return std::sqrt(dEta * dEta + dPhi * dPhi);
}

// ---------------------------------------------------------------------------
void MiniAODElectronTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<std::vector<pat::Electron>> gedHandle;
  iEvent.getByToken(gedElectronsToken_, gedHandle);

  edm::Handle<std::vector<pat::Electron>> lowPtHandle;
  iEvent.getByToken(lowPtElectronsToken_, lowPtHandle);

  auto gedTracks    = std::make_unique<reco::TrackCollection>();
  auto lowPtTracks  = std::make_unique<reco::TrackCollection>();
  auto mergedTracks = std::make_unique<reco::TrackCollection>();

  // ── GED electrons ──────────────────────────────────────────────────────────
  std::vector<const pat::Electron*> selectedGED;
  if (gedHandle.isValid()) {
    for (const auto& ele : *gedHandle) {
      if (!passGEDSelection(ele)) continue;
      if (ele.gsfTrack().isNull()) continue;
      selectedGED.push_back(&ele);
      // Slice GsfTrack → reco::Track (base-class copy; keeps kinematic parameters)
      gedTracks->push_back(static_cast<const reco::Track&>(*ele.gsfTrack()));
    }
  }

  // ── Low-pT electrons (cross-cleaned against selected GED) ──────────────────
  if (lowPtHandle.isValid()) {
    for (const auto& ele : *lowPtHandle) {
      if (!passLowPtSelection(ele)) continue;
      if (ele.gsfTrack().isNull()) continue;

      bool overlapsGED = false;
      for (const auto* gedEle : selectedGED) {
        if (deltaR(ele, *gedEle) < crossCleanDR_) { overlapsGED = true; break; }
      }
      if (overlapsGED) continue;

      lowPtTracks->push_back(static_cast<const reco::Track&>(*ele.gsfTrack()));
    }
  }

  // ── Merge ──────────────────────────────────────────────────────────────────
  mergedTracks->insert(mergedTracks->end(), gedTracks->begin(),   gedTracks->end());
  mergedTracks->insert(mergedTracks->end(), lowPtTracks->begin(), lowPtTracks->end());

  edm::LogInfo("MiniAODElectronTrackProducer")
      << "GED electron tracks: " << gedTracks->size()
      << "  low-pT electron tracks: " << lowPtTracks->size()
      << "  merged: " << mergedTracks->size();

  iEvent.put(std::move(gedTracks),    "gedElectronTracks");
  iEvent.put(std::move(lowPtTracks),  "lowPtElectronTracks");
  iEvent.put(std::move(mergedTracks), "mergedElectronTracks");
}

// ---------------------------------------------------------------------------
void MiniAODElectronTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("gedElectrons",   edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("lowPtElectrons", edm::InputTag("slimmedLowPtElectrons"));

  // Kinematic baseline (same for GED and low-pT per iDM analysis)
  desc.add<double>("minPt",      1.0);
  desc.add<double>("maxAbsEta",  2.4);

  // HEM veto: OFF by default (2018 data artifact; irrelevant for MC efficiency studies)
  desc.add<bool>("applyHEMVeto", false);

  // GED MVA ID: "mvaEleID-Fall17-noIso-V2-wp90" for Run2 UL after EGammaPostRecoTools.
  // Set to "" to disable the ID cut (e.g. while getting things running, or for Run3).
  desc.add<std::string>("gedMVALabel", "mvaEleID-Fall17-noIso-V2-wp90");

  // Low-pT MVA: standard userFloat key after re-running 2020Nov28 module is "ID".
  // Set to "" to disable the ID cut.
  desc.add<std::string>("lowPtMVALabel",     "ID");
  desc.add<double>     ("lowPtMVAThreshold", -0.25);

  // Cross-cleaning deltaR threshold (iDM analysis uses 0.05)
  desc.add<double>("crossCleanDR", 0.05);

  descriptions.add("miniAODElectronTrackProducer", desc);
}

DEFINE_FWK_MODULE(MiniAODElectronTrackProducer);
