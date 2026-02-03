// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      MiniAODTrackProducer
//
// Description: Extracts reco::Track collections from MINIAOD packed candidates.
//              Produces separate and merged collections from packedPFCandidates
//              and lostTracks for use with downstream vertex reconstruction.
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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

class MiniAODTrackProducer : public edm::stream::EDProducer<> {

public:
  explicit MiniAODTrackProducer(const edm::ParameterSet&);
  ~MiniAODTrackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Helper to extract tracks from a PackedCandidateCollection
  void extractTracks(const pat::PackedCandidateCollection& candidates,
                     reco::TrackCollection& outputTracks) const;

  // Input tokens
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> eleLostTracksToken_;
};

MiniAODTrackProducer::MiniAODTrackProducer(const edm::ParameterSet& iConfig)
    : pfCandidatesToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      lostTracksToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("lostTracks"))),
      eleLostTracksToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("eleLostTracks"))) {

  // Individual collections
  produces<reco::TrackCollection>("pfCandidateTracks");
  produces<reco::TrackCollection>("lostTracks");
  produces<reco::TrackCollection>("eleLostTracks");

  // Merged collections
  produces<reco::TrackCollection>("merged");        // pfCandidates + lostTracks
  produces<reco::TrackCollection>("mergedWithEle"); // pfCandidates + eleLostTracks
}

void MiniAODTrackProducer::extractTracks(
    const pat::PackedCandidateCollection& candidates,
    reco::TrackCollection& outputTracks) const {

  for (const auto& cand : candidates) {
    // Must have a track
    if (!cand.hasTrackDetails()) {
      continue;
    }

    // Get the pseudo-track and add to output
    outputTracks.push_back(cand.pseudoTrack());
  }
}

void MiniAODTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get input collections
  edm::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken_, pfCandidatesHandle);

  edm::Handle<pat::PackedCandidateCollection> lostTracksHandle;
  iEvent.getByToken(lostTracksToken_, lostTracksHandle);

  edm::Handle<pat::PackedCandidateCollection> eleLostTracksHandle;
  iEvent.getByToken(eleLostTracksToken_, eleLostTracksHandle);

  // Create output collections
  auto pfCandidateTracks = std::make_unique<reco::TrackCollection>();
  auto lostTracks = std::make_unique<reco::TrackCollection>();
  auto eleLostTracks = std::make_unique<reco::TrackCollection>();
  auto mergedTracks = std::make_unique<reco::TrackCollection>();
  auto mergedTracksWithEle = std::make_unique<reco::TrackCollection>();

  // Extract tracks from each source
  if (pfCandidatesHandle.isValid()) {
    extractTracks(*pfCandidatesHandle, *pfCandidateTracks);
  }

  if (lostTracksHandle.isValid()) {
    extractTracks(*lostTracksHandle, *lostTracks);
  }

  if (eleLostTracksHandle.isValid()) {
    extractTracks(*eleLostTracksHandle, *eleLostTracks);
  }

  // Build merged collections
  // merged = pfCandidateTracks + lostTracks
  mergedTracks->reserve(pfCandidateTracks->size() + lostTracks->size());
  for (const auto& track : *pfCandidateTracks) {
    mergedTracks->push_back(track);
  }
  for (const auto& track : *lostTracks) {
    mergedTracks->push_back(track);
  }

  // mergedWithEle = pfCandidateTracks + eleLostTracks
  mergedTracksWithEle->reserve(pfCandidateTracks->size() + eleLostTracks->size());
  for (const auto& track : *pfCandidateTracks) {
    mergedTracksWithEle->push_back(track);
  }
  for (const auto& track : *eleLostTracks) {
    mergedTracksWithEle->push_back(track);
  }

  // Log some stats
  edm::LogInfo("MiniAODTrackProducer")
      << "Extracted tracks - PF: " << pfCandidateTracks->size()
      << ", Lost: " << lostTracks->size()
      << ", EleLost: " << eleLostTracks->size()
      << ", Merged: " << mergedTracks->size()
      << ", MergedWithEle: " << mergedTracksWithEle->size();

  // Put collections into the event
  iEvent.put(std::move(pfCandidateTracks), "pfCandidateTracks");
  iEvent.put(std::move(lostTracks), "lostTracks");
  iEvent.put(std::move(eleLostTracks), "eleLostTracks");
  iEvent.put(std::move(mergedTracks), "merged");
  iEvent.put(std::move(mergedTracksWithEle), "mergedWithEle");
}

void MiniAODTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Input collections (MINIAOD defaults)
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("eleLostTracks", edm::InputTag("lostTracks", "eleTracks"));

  descriptions.add("miniAODTrackProducer", desc);
}

DEFINE_FWK_MODULE(MiniAODTrackProducer);
