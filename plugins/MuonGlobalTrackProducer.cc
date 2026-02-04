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

  // Input token
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
};

MuonGlobalTrackProducer::MuonGlobalTrackProducer(const edm::ParameterSet& iConfig)
    : muonsToken_(consumes<reco::MuonCollection>(
          iConfig.getParameter<edm::InputTag>("muons"))) {

  // Output collection
  produces<reco::TrackCollection>("globalTracks");
}

void MuonGlobalTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get input collection
  edm::Handle<reco::MuonCollection> muonsHandle;
  iEvent.getByToken(muonsToken_, muonsHandle);

  // Create output collection
  auto globalTracks = std::make_unique<reco::TrackCollection>();

  // Extract global tracks from muons
  if (muonsHandle.isValid()) {
    for (const auto& muon : *muonsHandle) {
      // Check if muon has a valid global track
      if (muon.globalTrack().isNonnull()) {
        globalTracks->push_back(*muon.globalTrack());
      }
    }
  }

  // Log some stats
  edm::LogInfo("MuonGlobalTrackProducer")
      << "Extracted " << globalTracks->size() << " global tracks from "
      << (muonsHandle.isValid() ? muonsHandle->size() : 0) << " muons";

  // Put collection into the event
  iEvent.put(std::move(globalTracks), "globalTracks");
}

void MuonGlobalTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Input collection (AOD default)
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));

  descriptions.add("muonGlobalTrackProducer", desc);
}

DEFINE_FWK_MODULE(MuonGlobalTrackProducer);
