// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      HyddraSVsHadronicDiagnosticProducer
//
// Description: Diagnostic producer that runs the hadronic HYDDRA pipeline and
//              publishes intermediate vertex collections after each of the 5
//              stages (seeding, merging, cleaning, disambiguation, filtering).
//              Intended for signal efficiency studies.
//              Mirrors HyddraSVsDiagnosticProducer but uses HadronicHYDDRA.
//
// Original Author:  Andres Abreu
//

#include <memory>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// HYDDRA
#include "KUCMSNtupleizer/HyddraSVProducer/interface/HadronicHYDDRA.h"

class HyddraSVsHadronicDiagnosticProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraSVsHadronicDiagnosticProducer(const edm::ParameterSet&);
  ~HyddraSVsHadronicDiagnosticProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // HYDDRA hadronic reconstruction object
  HadronicHYDDRA hadronic_;
};

HyddraSVsHadronicDiagnosticProducer::HyddraSVsHadronicDiagnosticProducer(
    const edm::ParameterSet& iConfig) :
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  hadronic_(iConfig.getParameter<edm::ParameterSet>("hadronic"))
{
  produces<reco::VertexCollection>("hadronicSeeds");
  produces<reco::VertexCollection>("hadronicMerged");
  produces<reco::VertexCollection>("hadronicCleaned");
  produces<reco::VertexCollection>("hadronicDisambiguated");
  produces<reco::VertexCollection>("hadronicFiltered");
}

void HyddraSVsHadronicDiagnosticProducer::produce(edm::Event& iEvent,
                                                   const edm::EventSetup& iSetup) {
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByToken(tracksToken_, tracksHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);

  std::vector<reco::TrackRef> trackRefs;
  trackRefs.reserve(tracksHandle->size());
  for (size_t i = 0; i < tracksHandle->size(); ++i) {
    trackRefs.emplace_back(tracksHandle, i);
  }

  if (trackRefs.size() > 500) {
    edm::LogWarning("HyddraSVsHadronicDiagnosticProducer")
        << "Large input track collection (" << trackRefs.size()
        << " tracks). Vertex reconstruction may be very slow.";
  }

  const reco::Vertex& pv = pvHandle->at(0);

  auto snaps = hadronic_.run_with_snapshots(trackRefs, ttBuilder, pv);

  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterSeeding),        "hadronicSeeds");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterMerging),        "hadronicMerged");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterCleaning),       "hadronicCleaned");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterDisambiguation), "hadronicDisambiguated");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterFiltering),      "hadronicFiltered");
}

void HyddraSVsHadronicDiagnosticProducer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks",       edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));

  // Hadronic PSet — same defaults as HyddraSVsProducer
  edm::ParameterSetDescription hadronicDesc;
  hadronicDesc.add<double>("seedCosThetaCut",  0.0);
  hadronicDesc.add<double>("minMass",          2.0);
  hadronicDesc.add<double>("minPOverE",        0.6);
  hadronicDesc.add<double>("maxNormChi2",      5.0);
  hadronicDesc.add<double>("minDxySignificance", 40.0);
  hadronicDesc.add<int>   ("minSize",          5);
  hadronicDesc.add<double>("minCosTheta",      0.0);
  hadronicDesc.add<double>("maxDecayAngle",    0.9);
  hadronicDesc.add<bool>  ("doMerging",        true);
  hadronicDesc.add<bool>  ("doCleaning",       true);
  hadronicDesc.add<bool>  ("doDisambiguation", true);
  hadronicDesc.add<bool>  ("doFiltering",      true);
  hadronicDesc.add<bool>  ("useVertexSmoothing", false);
  desc.add<edm::ParameterSetDescription>("hadronic", hadronicDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsHadronicDiagnosticProducer);
