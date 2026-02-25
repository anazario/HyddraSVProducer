// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      HyddraSVsDiagnosticProducer
//
// Description: Diagnostic producer that runs the leptonic HYDDRA pipeline and
//              publishes intermediate vertex collections after each of the 5
//              stages (seeding, merging, cleaning, disambiguation, filtering).
//              Intended for signal efficiency studies; does NOT run hadronic.
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
#include "KUCMSNtupleizer/HyddraSVProducer/interface/LeptonicHYDDRA.h"

class HyddraSVsDiagnosticProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraSVsDiagnosticProducer(const edm::ParameterSet&);
  ~HyddraSVsDiagnosticProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // HYDDRA leptonic reconstruction object
  LeptonicHYDDRA leptonic_;
};

HyddraSVsDiagnosticProducer::HyddraSVsDiagnosticProducer(const edm::ParameterSet& iConfig) :
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  leptonic_(iConfig.getParameter<edm::ParameterSet>("leptonic"))
{
  produces<reco::VertexCollection>("leptonicSeeds");
  produces<reco::VertexCollection>("leptonicMerged");
  produces<reco::VertexCollection>("leptonicCleaned");
  produces<reco::VertexCollection>("leptonicDisambiguated");
  produces<reco::VertexCollection>("leptonicFiltered");
}

void HyddraSVsDiagnosticProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get inputs
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByToken(tracksToken_, tracksHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);

  // Build TrackRef vector from the input collection
  std::vector<reco::TrackRef> trackRefs;
  trackRefs.reserve(tracksHandle->size());
  for (size_t i = 0; i < tracksHandle->size(); ++i) {
    trackRefs.emplace_back(tracksHandle, i);
  }

  if (trackRefs.size() > 500) {
    edm::LogWarning("HyddraSVsDiagnosticProducer")
        << "Large input track collection (" << trackRefs.size()
        << " tracks). Vertex reconstruction may be very slow.";
  }

  const reco::Vertex& pv = pvHandle->at(0);

  // Run leptonic pipeline with intermediate snapshots
  auto snaps = leptonic_.run_with_snapshots(trackRefs, ttBuilder, pv);

  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterSeeding),       "leptonicSeeds");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterMerging),       "leptonicMerged");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterCleaning),      "leptonicCleaned");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterDisambiguation),"leptonicDisambiguated");
  iEvent.put(std::make_unique<reco::VertexCollection>(snaps.afterFiltering),     "leptonicFiltered");
}

void HyddraSVsDiagnosticProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks",       edm::InputTag("muonBestTrackProducer", "globalTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));

  // Leptonic PSet — same defaults as HyddraSVsProducer
  edm::ParameterSetDescription leptonicDesc;
  leptonicDesc.add<double>("seedCosThetaCut",              0.75);
  leptonicDesc.add<double>("minMass",                      2.0);
  leptonicDesc.add<double>("minPOverE",                    0.6);
  leptonicDesc.add<double>("maxNormChi2",                  5.0);
  leptonicDesc.add<double>("minDxySignificance",           25.0);
  leptonicDesc.add<double>("maxCompatibility",             1.5);
  leptonicDesc.add<double>("minCleanCosTheta",             0.5);
  leptonicDesc.add<bool>  ("useDiagonalCut",               false);
  leptonicDesc.add<double>("cleanCutSlope",                0.0);
  leptonicDesc.add<double>("minTrackCosTheta",             0.5);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Limit",     0.95);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Intercept", 1.8);
  leptonicDesc.add<double>("trackCosThetaCM_Slope",        -1.0);
  leptonicDesc.add<bool>  ("requireChargeNeutrality",      true);
  leptonicDesc.add<double>("minVtxCosTheta",               -1.0);
  leptonicDesc.add<bool>  ("useAbsVtxCosTheta",            false);
  leptonicDesc.add<double>("maxVtxDecayAngle",              1.0);
  leptonicDesc.add<bool>  ("useAbsVtxDecayAngle",          false);
  leptonicDesc.add<bool>  ("applyVtxDecayAngleCleaning",   false);
  leptonicDesc.add<bool>  ("applyVtxDecayAngleFiltering",  false);
  leptonicDesc.add<bool>  ("doMerging",                    true);
  leptonicDesc.add<bool>  ("doCleaning",                   true);
  leptonicDesc.add<bool>  ("doDisambiguation",             true);
  leptonicDesc.add<bool>  ("doFiltering",                  true);
  desc.add<edm::ParameterSetDescription>("leptonic", leptonicDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsDiagnosticProducer);
