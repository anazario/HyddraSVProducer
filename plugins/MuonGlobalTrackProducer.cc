// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      MuonGlobalTrackProducer
//
// Description: Extracts global tracks from the AOD muon collection into a
//              reco::TrackCollection for use with downstream vertex reconstruction.
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <vector>

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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class MuonGlobalTrackProducer : public edm::stream::EDProducer<> {

public:
  explicit MuonGlobalTrackProducer(const edm::ParameterSet&);
  ~MuonGlobalTrackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::MuonCollection> displacedMuonsToken_;
};

MuonGlobalTrackProducer::MuonGlobalTrackProducer(const edm::ParameterSet& iConfig)
    : muonsToken_(consumes<reco::MuonCollection>(
          iConfig.getParameter<edm::InputTag>("muons"))),
      displacedMuonsToken_(consumes<reco::MuonCollection>(
          iConfig.getParameter<edm::InputTag>("displacedMuons"))) {

  // Output collections
  produces<reco::TrackCollection>("globalTracks");
  produces<reco::TrackCollection>("displacedGlobalTracks");
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

  // Extract global tracks from standard muons
  if (muonsHandle.isValid()) {
    for (const auto& muon : *muonsHandle) {
      // Check if muon has a valid global track
      if (muon.globalTrack().isNonnull()) {
        globalTracks->push_back(*muon.globalTrack());
      }
    }
  }

  // Extract global tracks from displaced muons
  if (displacedMuonsHandle.isValid()) {
    for (const auto& muon : *displacedMuonsHandle) {
      // Check if muon has a valid global track
      if (muon.globalTrack().isNonnull()) {
        displacedGlobalTracks->push_back(*muon.globalTrack());
      }
    }
  }

  // Log some stats
  edm::LogInfo("MuonGlobalTrackProducer")
      << "Extracted " << globalTracks->size() << " global tracks from "
      << (muonsHandle.isValid() ? muonsHandle->size() : 0) << " muons";

  edm::LogInfo("MuonGlobalTrackProducer")
      << "Extracted " << displacedGlobalTracks->size() << " global tracks from "
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

  descriptions.add("muonGlobalTrackProducer", desc);
}

DEFINE_FWK_MODULE(MuonGlobalTrackProducer);
