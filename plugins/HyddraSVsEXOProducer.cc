// -*- C++ -*-
//
// Package:    HYDDRA
// Class:      HyddraSVsEXOProducer
//
// Description: Produces leptonic HYDDRA displaced vertex collections for
//              the EXONanoAOD prototype format. Runs a forked pipeline:
//
//   Tier 0 (inclusive): seeds -> disambiguation
//   Tier 1 (isolated):  seeds -> merging -> drop 3+ track -> disambiguation
//
//   Seeding is performed once and shared between both tiers. No signal-dependent
//   kinematic or angular cuts are applied; all variables are stored in the ntuple
//   for analyst-level selection.
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

class HyddraSVsEXOProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraSVsEXOProducer(const edm::ParameterSet&);
  ~HyddraSVsEXOProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // HYDDRA reconstruction object (leptonic only)
  LeptonicHYDDRA leptonic_;
};

HyddraSVsEXOProducer::HyddraSVsEXOProducer(const edm::ParameterSet& iConfig) :
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  leptonic_(iConfig.getParameter<edm::ParameterSet>("leptonic"))
{
  produces<reco::VertexCollection>("inclusiveVertices");
  produces<reco::VertexCollection>("isolatedVertices");
  produces<std::vector<int>>("isolationFlags");
}

void HyddraSVsEXOProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByToken(tracksToken_, tracksHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);

  std::vector<reco::TrackRef> trackRefs;
  trackRefs.reserve(tracksHandle->size());
  for (size_t i = 0; i < tracksHandle->size(); ++i)
    trackRefs.emplace_back(tracksHandle, i);

  if (trackRefs.size() > 500)
    edm::LogWarning("HyddraSVsEXOProducer") << "Large input track collection ("
      << trackRefs.size() << " tracks). Reconstruction may be slow.";

  const reco::Vertex& pv = pvHandle->at(0);

  leptonic_.run_forked(trackRefs, ttBuilder, pv);

  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.vertices()),         "inclusiveVertices");
  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.isolatedVertices()), "isolatedVertices");
  iEvent.put(std::make_unique<std::vector<int>>(leptonic_.computeIsolationFlags()),  "isolationFlags");
}

void HyddraSVsEXOProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks",       edm::InputTag("muonEnhancedTracks", "muonEnhancedTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));

  edm::ParameterSetDescription leptonicDesc;
  // Seeding: only chi2 gate; cosTheta cut disabled for universality
  leptonicDesc.add<double>("seedCosThetaCut",  -1.0);
  leptonicDesc.add<double>("maxNormChi2",        5.0);
  // Unused by the forked pipeline (doCleaning=false, doFiltering=false)
  // but required by LeptonicHYDDRA constructor
  leptonicDesc.add<double>("minMass",            0.0);
  leptonicDesc.add<double>("minPOverE",          0.0);
  leptonicDesc.add<double>("minDxySignificance", 0.0);
  leptonicDesc.add<double>("maxCompatibility",   999.0);
  leptonicDesc.add<double>("minCleanCosTheta",  -1.0);
  leptonicDesc.add<double>("maxCleanCosTheta",   1.0);
  leptonicDesc.add<bool>  ("invertCleanCosThetaCut",      false);
  leptonicDesc.add<bool>  ("useDiagonalCut",              false);
  leptonicDesc.add<double>("cleanCutSlope",               0.0);
  leptonicDesc.add<double>("minTrackCosTheta",           -1.0);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Limit",    1.0);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Intercept",999.0);
  leptonicDesc.add<double>("trackCosThetaCM_Slope",       0.0);
  leptonicDesc.add<bool>  ("requireChargeNeutrality",    false);
  leptonicDesc.add<double>("minVtxCosTheta",             -1.0);
  leptonicDesc.add<double>("maxVtxCosTheta",              1.0);
  leptonicDesc.add<bool>  ("useAbsVtxCosTheta",          false);
  leptonicDesc.add<double>("maxVtxDecayAngle",            1.0);
  leptonicDesc.add<bool>  ("useAbsVtxDecayAngle",        false);
  leptonicDesc.add<bool>  ("applyVtxDecayAngleCleaning", false);
  leptonicDesc.add<bool>  ("applyVtxDecayAngleFiltering",false);
  leptonicDesc.add<double>("minMassFilter",               0.0);
  leptonicDesc.add<double>("minBetaFilter",               0.0);
  // Pipeline stage flags: merging and disambiguation run; cleaning and filtering disabled
  leptonicDesc.add<bool>("doMerging",        true);
  leptonicDesc.add<bool>("doCleaning",       false);
  leptonicDesc.add<bool>("doDisambiguation", true);
  leptonicDesc.add<bool>("doFiltering",      false);
  leptonicDesc.add<bool>("useVertexSmoothing", false);
  desc.add<edm::ParameterSetDescription>("leptonic", leptonicDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsEXOProducer);
