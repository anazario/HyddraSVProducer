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
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Track/SC matching helpers (must come before HYDDRA headers)
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchTracksToSC.h"

// HYDDRA (transitively provides VertexHelper)
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
  edm::EDGetTokenT<reco::TrackCollection> muonTracksToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> mergedSCsToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  // HYDDRA leptonic reconstruction object
  LeptonicHYDDRA leptonic_;
};

HyddraSVsDiagnosticProducer::HyddraSVsDiagnosticProducer(const edm::ParameterSet& iConfig) :
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  muonTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("muonTracks"))),
  mergedSCsToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("mergedSCs"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  magneticFieldToken_(esConsumes()),
  caloGeometryToken_(esConsumes()),
  leptonic_(iConfig.getParameter<edm::ParameterSet>("leptonic"))
{
  trackAssociator_.useDefaultPropagator();
  edm::ParameterSet trackAssocPSet = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  edm::ConsumesCollector iC = consumesCollector();
  trackAssocParameters_.loadParameters(trackAssocPSet, iC);

  produces<reco::VertexCollection>("leptonicSeeds");
  produces<reco::VertexCollection>("leptonicMerged");
  produces<reco::VertexCollection>("leptonicCleaned");
  produces<reco::VertexCollection>("leptonicDisambiguated");
  produces<reco::VertexCollection>("leptonicFiltered");
  produces<reco::VertexCollection>("leptonicID");
}

void HyddraSVsDiagnosticProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get inputs
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByToken(tracksToken_, tracksHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::TrackCollection> muonTracksHandle;
  iEvent.getByToken(muonTracksToken_, muonTracksHandle);

  edm::Handle<reco::SuperClusterCollection> mergedSCsHandle;
  iEvent.getByToken(mergedSCsToken_, mergedSCsHandle);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

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

  // ID stage: build electron track collection via SC matching, then apply
  // muon ID (both tracks in muon collection) or electron ID (both SC-matched).
  reco::TrackCollection allFilteredTracks;
  for (const auto& vtx : snaps.afterFiltering)
    for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it)
      allFilteredTracks.emplace_back(**it);

  reco::TrackCollection electronTracks;
  if (mergedSCsHandle.isValid() && !mergedSCsHandle->empty()) {
    MatchTracksToSC<reco::Track> scAssigner(iEvent, iSetup, magfield, ecalGeometry,
                                            trackAssocParameters_, allFilteredTracks,
                                            *mergedSCsHandle);
    for (const auto& pair : scAssigner.GetMatchedTrackSCPairs())
      if (pair.GetDeltaR() < 0.04)
        electronTracks.emplace_back(pair.GetTrack());
  }

  auto leptonicID = std::make_unique<reco::VertexCollection>();
  for (const auto& vtx : snaps.afterFiltering) {
    if (vtx.tracksSize() != 2) continue;
    const bool passLooseMuonID     = VertexHelper::CountInstances(vtx, *muonTracksHandle) == 2;
    const bool passLooseElectronID = VertexHelper::CountInstances(vtx, electronTracks)    == 2;
    if (passLooseMuonID || passLooseElectronID)
      leptonicID->push_back(vtx);
  }
  iEvent.put(std::move(leptonicID), "leptonicID");
}

void HyddraSVsDiagnosticProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks",       edm::InputTag("muonBestTrackProducer", "globalTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("muonTracks",   edm::InputTag("muonEnhancedTracks", "muonGlobalTracks"));
  desc.add<edm::InputTag>("mergedSCs",    edm::InputTag("ecalTracks", "displacedElectronSCs"));

  edm::ParameterSetDescription trackAssocDesc;
  trackAssocDesc.setAllowAnything();
  desc.add<edm::ParameterSetDescription>("TrackAssociatorParameters", trackAssocDesc);

  // Leptonic PSet — same defaults as HyddraSVsProducer
  edm::ParameterSetDescription leptonicDesc;
  leptonicDesc.add<double>("seedCosThetaCut",              0.75);
  leptonicDesc.add<double>("minMass",                      2.0);
  leptonicDesc.add<double>("minPOverE",                    0.6);
  leptonicDesc.add<double>("maxNormChi2",                  5.0);
  leptonicDesc.add<double>("minDxySignificance",           25.0);
  leptonicDesc.add<double>("maxCompatibility",             1.5);
  leptonicDesc.add<double>("minCleanCosTheta",             0.5);
  leptonicDesc.add<double>("maxCleanCosTheta",             1.0);
  leptonicDesc.add<bool>  ("invertCleanCosThetaCut",       false);
  leptonicDesc.add<bool>  ("useDiagonalCut",               false);
  leptonicDesc.add<double>("cleanCutSlope",                0.0);
  leptonicDesc.add<double>("minTrackCosTheta",             0.5);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Limit",     0.95);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Intercept", 1.8);
  leptonicDesc.add<double>("trackCosThetaCM_Slope",        -1.0);
  leptonicDesc.add<bool>  ("requireChargeNeutrality",      true);
  leptonicDesc.add<double>("minVtxCosTheta",               -1.0);
  leptonicDesc.add<double>("maxVtxCosTheta",                1.0);
  leptonicDesc.add<bool>  ("useAbsVtxCosTheta",            false);
  leptonicDesc.add<double>("maxVtxDecayAngle",              1.0);
  leptonicDesc.add<bool>  ("useAbsVtxDecayAngle",          false);
  leptonicDesc.add<bool>  ("applyVtxDecayAngleCleaning",   false);
  leptonicDesc.add<bool>  ("applyVtxDecayAngleFiltering",  false);
  leptonicDesc.add<double>("minMassFilter",                 0.0);
  leptonicDesc.add<double>("minBetaFilter",                 0.0);
  leptonicDesc.add<bool>  ("doMerging",                    true);
  leptonicDesc.add<bool>  ("doCleaning",                   true);
  leptonicDesc.add<bool>  ("doDisambiguation",             true);
  leptonicDesc.add<bool>  ("doFiltering",                  true);
  leptonicDesc.add<bool>  ("useVertexSmoothing",           false);
  desc.add<edm::ParameterSetDescription>("leptonic", leptonicDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsDiagnosticProducer);
