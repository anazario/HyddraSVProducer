// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      FilteredTrackProducer
//
// Description: Filters a track collection based on configurable quality cuts.
//              Produces a reco::TrackCollection for use with downstream analysis.
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

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class FilteredTrackProducer : public edm::stream::EDProducer<> {

public:
  explicit FilteredTrackProducer(const edm::ParameterSet&);
  ~FilteredTrackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // Cut parameters
  double minPt_;
  double maxPtResolution_;
  double maxNormalizedChi2_;
  double maxAbsSip2D_;
};

FilteredTrackProducer::FilteredTrackProducer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<reco::TrackCollection>(
          iConfig.getParameter<edm::InputTag>("tracks"))),
      pvToken_(consumes<reco::VertexCollection>(
          iConfig.getParameter<edm::InputTag>("pvCollection"))),
      ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      minPt_(iConfig.getParameter<double>("minPt")),
      maxPtResolution_(iConfig.getParameter<double>("maxPtResolution")),
      maxNormalizedChi2_(iConfig.getParameter<double>("maxNormalizedChi2")),
      maxAbsSip2D_(iConfig.getParameter<double>("maxAbsSip2D")) {

  produces<reco::TrackCollection>("filteredTracks");
}

void FilteredTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get input collections
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByToken(tracksToken_, tracksHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);

  // Create output collection
  auto filteredTracks = std::make_unique<reco::TrackCollection>();

  if (!tracksHandle.isValid() || !pvHandle.isValid() || pvHandle->empty()) {
    iEvent.put(std::move(filteredTracks), "filteredTracks");
    return;
  }

  const reco::Vertex& pv = pvHandle->at(0);

  size_t nInput = tracksHandle->size();
  size_t nFailPt = 0;
  size_t nFailPtRes = 0;
  size_t nFailChi2 = 0;
  size_t nFailSip2D = 0;

  for (const auto& track : *tracksHandle) {
    // pt cut
    if (track.pt() <= minPt_) {
      nFailPt++;
      continue;
    }

    // pt resolution cut
    double ptRes = track.ptError() / track.pt();
    if (ptRes >= maxPtResolution_) {
      nFailPtRes++;
      continue;
    }

    // normalized chi2 cut
    if (track.normalizedChi2() >= maxNormalizedChi2_) {
      nFailChi2++;
      continue;
    }

    // sip2D cut (requires transient track)
    reco::TransientTrack ttrack = ttBuilder->build(track);
    auto ip2dResult = IPTools::signedTransverseImpactParameter(ttrack, GlobalVector(track.px(), track.py(), track.pz()), pv);

    if (ip2dResult.first) {
      double sip2D = ip2dResult.second.significance();
      if (std::fabs(sip2D) >= maxAbsSip2D_) {
        nFailSip2D++;
        continue;
      }
    }

    // Track passed all cuts
    filteredTracks->push_back(track);
  }

  edm::LogInfo("FilteredTrackProducer")
      << "Filtered tracks: " << nInput << " -> " << filteredTracks->size()
      << " | Failed: pt=" << nFailPt
      << ", ptRes=" << nFailPtRes
      << ", chi2=" << nFailChi2
      << ", sip2D=" << nFailSip2D;

  iEvent.put(std::move(filteredTracks), "filteredTracks");
}

void FilteredTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));
  desc.add<double>("minPt", 5.0);
  desc.add<double>("maxPtResolution", 0.08);
  desc.add<double>("maxNormalizedChi2", 5.0);
  desc.add<double>("maxAbsSip2D", 5.0);

  descriptions.add("filteredTrackProducer", desc);
}

DEFINE_FWK_MODULE(FilteredTrackProducer);
