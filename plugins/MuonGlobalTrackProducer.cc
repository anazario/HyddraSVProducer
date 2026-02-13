// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      MuonGlobalTrackProducer
//
// Description: Extracts tracks from the AOD muon collection into a
//              reco::TrackCollection for use with downstream vertex reconstruction.
//              Supports three modes:
//                "globalTrack" - extract only globalTrack (original behavior)
//                "bestTrack"   - use CMSSW's muonBestTrack()
//                "priority"    - try track types in configurable priority order
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class MuonGlobalTrackProducer : public edm::stream::EDProducer<> {

public:
  explicit MuonGlobalTrackProducer(const edm::ParameterSet&);
  ~MuonGlobalTrackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Get a track ref from a muon by type name
  static reco::TrackRef getTrackRef(const reco::Muon& muon, const std::string& type);

  // Select the best track for a muon based on configured mode
  reco::TrackRef selectTrack(const reco::Muon& muon) const;

  // Extract tracks from a muon collection
  void extractTracks(const reco::MuonCollection& muons,
                     reco::TrackCollection& outputTracks) const;

  // Input tokens
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::MuonCollection> displacedMuonsToken_;

  // Track selection configuration
  std::string mode_;
  std::vector<std::string> trackPriority_;

  // Valid track type names
  static const std::set<std::string> validTrackTypes_;
};

const std::set<std::string> MuonGlobalTrackProducer::validTrackTypes_ = {
    "globalTrack", "innerTrack", "outerTrack", "bestTrack"};

MuonGlobalTrackProducer::MuonGlobalTrackProducer(const edm::ParameterSet& iConfig)
    : muonsToken_(consumes<reco::MuonCollection>(
          iConfig.getParameter<edm::InputTag>("muons"))),
      displacedMuonsToken_(consumes<reco::MuonCollection>(
          iConfig.getParameter<edm::InputTag>("displacedMuons"))),
      mode_(iConfig.getParameter<std::string>("mode")),
      trackPriority_(iConfig.getParameter<std::vector<std::string>>("trackPriority")) {

  // Validate mode
  if (mode_ != "globalTrack" && mode_ != "bestTrack" && mode_ != "priority") {
    throw cms::Exception("InvalidConfiguration")
        << "Unknown mode: '" << mode_
        << "'. Must be one of: globalTrack, bestTrack, priority";
  }

  // Validate trackPriority entries
  if (mode_ == "priority") {
    if (trackPriority_.empty()) {
      throw cms::Exception("InvalidConfiguration")
          << "trackPriority must not be empty when mode is 'priority'";
    }
    for (const auto& type : trackPriority_) {
      if (validTrackTypes_.find(type) == validTrackTypes_.end()) {
        throw cms::Exception("InvalidConfiguration")
            << "Unknown track type in trackPriority: '" << type
            << "'. Valid types: globalTrack, innerTrack, outerTrack, bestTrack";
      }
    }
  }

  // Output collections (labels unchanged for backward compatibility)
  produces<reco::TrackCollection>("globalTracks");
  produces<reco::TrackCollection>("displacedGlobalTracks");
}

reco::TrackRef MuonGlobalTrackProducer::getTrackRef(
    const reco::Muon& muon, const std::string& type) {
  if (type == "globalTrack") return muon.globalTrack();
  if (type == "innerTrack")  return muon.innerTrack();
  if (type == "outerTrack")  return muon.outerTrack();
  if (type == "bestTrack")   return muon.muonBestTrack();
  return reco::TrackRef();
}

reco::TrackRef MuonGlobalTrackProducer::selectTrack(const reco::Muon& muon) const {
  if (mode_ == "globalTrack") {
    return muon.globalTrack();
  }

  if (mode_ == "bestTrack") {
    return muon.muonBestTrack();
  }

  // "priority" mode: try each type in order, take first non-null
  for (const auto& type : trackPriority_) {
    reco::TrackRef ref = getTrackRef(muon, type);
    if (ref.isNonnull()) return ref;
  }
  return reco::TrackRef();
}

void MuonGlobalTrackProducer::extractTracks(
    const reco::MuonCollection& muons,
    reco::TrackCollection& outputTracks) const {
  for (const auto& muon : muons) {
    reco::TrackRef trackRef = selectTrack(muon);
    if (trackRef.isNonnull()) {
      outputTracks.push_back(*trackRef);
    }
  }
}

void MuonGlobalTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get input collections
  edm::Handle<reco::MuonCollection> muonsHandle;
  iEvent.getByToken(muonsToken_, muonsHandle);

  edm::Handle<reco::MuonCollection> displacedMuonsHandle;
  iEvent.getByToken(displacedMuonsToken_, displacedMuonsHandle);

  // Create output collections
  auto globalTracks = std::make_unique<reco::TrackCollection>();
  auto displacedGlobalTracks = std::make_unique<reco::TrackCollection>();

  // Extract tracks from standard muons
  if (muonsHandle.isValid()) {
    extractTracks(*muonsHandle, *globalTracks);
  }

  // Extract tracks from displaced muons
  if (displacedMuonsHandle.isValid()) {
    extractTracks(*displacedMuonsHandle, *displacedGlobalTracks);
  }

  // Log some stats
  edm::LogInfo("MuonGlobalTrackProducer")
      << "Extracted " << globalTracks->size() << " tracks (mode=" << mode_ << ") from "
      << (muonsHandle.isValid() ? muonsHandle->size() : 0) << " muons";

  edm::LogInfo("MuonGlobalTrackProducer")
      << "Extracted " << displacedGlobalTracks->size() << " tracks (mode=" << mode_ << ") from "
      << (displacedMuonsHandle.isValid() ? displacedMuonsHandle->size() : 0) << " displaced muons";

  // Put collections into the event
  iEvent.put(std::move(globalTracks), "globalTracks");
  iEvent.put(std::move(displacedGlobalTracks), "displacedGlobalTracks");
}

void MuonGlobalTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Input collections (AOD defaults)
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("displacedMuons", edm::InputTag("displacedMuons"));

  // Track selection mode: "globalTrack", "bestTrack", or "priority"
  desc.add<std::string>("mode", "priority");

  // Priority order for "priority" mode (ignored in other modes)
  desc.add<std::vector<std::string>>("trackPriority",
      {"globalTrack", "innerTrack", "outerTrack"});

  descriptions.add("muonGlobalTrackProducer", desc);
}

DEFINE_FWK_MODULE(MuonGlobalTrackProducer);
