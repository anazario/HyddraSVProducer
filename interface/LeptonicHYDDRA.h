#pragma once

#include "HYDDRABase.h"

// Leptonic displaced vertex reconstruction.
// Targets 2-track vertices (e.g. displaced dileptons).
// Disambiguation keeps the vertex with the smallest CM cosTheta.
class LeptonicHYDDRA : public HYDDRABase<LeptonicHYDDRA> {

 public:

  LeptonicHYDDRA(const edm::ParameterSet& pset) : HYDDRABase(pset) {
    // Cleaning thresholds
    maxCompatibility_ = pset.getParameter<double>("maxCompatibility");
    minCleanCosTheta_ = pset.getParameter<double>("minCleanCosTheta");

    // Post-cleaning kinematic cut
    minPostCosTheta_ = pset.getParameter<double>("minPostCosTheta");

    // Final filtering cuts (post-disambiguation, 2-track only)
    minTrackCosTheta_         = pset.getParameter<double>("minTrackCosTheta");
    maxTrackCosThetaCM_Limit_ = pset.getParameter<double>("maxTrackCosThetaCM_Limit");
    maxTrackCosThetaCM_Slope_ = pset.getParameter<double>("maxTrackCosThetaCM_Slope");
  }

  // Stage 1: Seeding
  void seedingImpl(const std::vector<reco::TrackRef>& tracks) {
    this->generateSeeds(tracks);
  }

  // Stage 2: Merging
  void mergingImpl() {
    this->mergeVertices();
  }

  // Stage 3: Cleaning
  // For 2-track vertices, only apply kinematic cuts.
  // For larger vertices, remove tracks with poor compatibility or alignment,
  // then apply kinematic cuts to whatever remains.
  void cleaningImpl() {
    if (this->empty() || !primaryVertex_ || !primaryVertex_->isValid()) return;

    HYDDRA_DBG("[Leptonic] Cleaning " << this->size() << " vertices...\n");

    TrackVertexSetCollection cleanedVertices;

    for (const auto& vertex : *this) {
      // Two-track vertices skip track removal
      if (vertex.size() <= 2) {
	if (!failsKinematics(vertex)) cleanedVertices.add(vertex);
	continue;
      }

      // Scan for bad tracks before modifying
      std::vector<reco::TrackRef> tracksToRemove;
      for (const auto& trackRef : vertex.tracks()) {
	double compatibility = vertex.compatibility(trackRef);
	double cosTheta      = vertex.trackCosTheta(*primaryVertex_, trackRef);

	if (compatibility > maxCompatibility_ || cosTheta < minCleanCosTheta_) {
	  tracksToRemove.push_back(trackRef);
	}
      }

      // No bad tracks — keep original if it passes kinematics
      if (tracksToRemove.empty()) {
	if (!failsKinematics(vertex)) cleanedVertices.add(vertex);
	continue;
      }

      // Remove bad tracks from a working copy
      TrackVertexSet workingVertex(vertex);
      for (const auto& trackRef : tracksToRemove) {
	if (workingVertex.size() > 2) {
	  workingVertex.removeTrack(trackRef);
	}
      }

      if (workingVertex.isValid() && !failsKinematics(workingVertex)) {
	cleanedVertices.add(workingVertex);
      }
    }

    this->clear();
    for (const auto& v : cleanedVertices) {
      this->add(v);
    }
  }

  // Stage 4: Disambiguation
  // For each track shared between vertices, keep the vertex with the
  // smallest CM cosTheta (most isotropic decay) and reject the rest.
  void disambiguationImpl() {
    if (this->empty()) return;

    // Map each track to the vertices that contain it
    std::map<reco::TrackRef, std::vector<const TrackVertexSet*>> trackToVertexMap;
    for (const auto& vertex : *this) {
      for (const auto& track : vertex) {
	trackToVertexMap[track].push_back(&vertex);
      }
    }

    // Identify vertices to reject
    std::set<const TrackVertexSet*> rejectedPtrs;

    for (const auto& [track, vtxPtrList] : trackToVertexMap) {
      std::vector<const TrackVertexSet*> candidates;
      for (const auto* vtxPtr : vtxPtrList) {
	if (rejectedPtrs.find(vtxPtr) == rejectedPtrs.end()) {
	  candidates.push_back(vtxPtr);
	}
      }

      if (candidates.size() > 1) {
	double minCosTheta = std::numeric_limits<double>::max();
	const TrackVertexSet* bestVertex = nullptr;

	for (const auto* vtxPtr : candidates) {
	  double cosTheta = VertexHelper::CalculateCMCosTheta(*vtxPtr, *track);
	  if (cosTheta < minCosTheta) {
	    minCosTheta = cosTheta;
	    bestVertex = vtxPtr;
	  }
	}

	for (const auto* vtxPtr : candidates) {
	  if (vtxPtr != bestVertex) {
	    rejectedPtrs.insert(vtxPtr);
	  }
	}
      }
    }

    // Rebuild collection without rejected vertices
    if (!rejectedPtrs.empty()) {
      HYDDRA_DBG("[Leptonic] Disambiguation dropping " << rejectedPtrs.size() << " vertices.\n");

      TrackVertexSetCollection survivors;
      for (const auto& vtx : *this) {
	if (rejectedPtrs.find(&vtx) == rejectedPtrs.end()) {
	  survivors.add(vtx);
	}
      }
      this->clear();
      for (const auto& v : survivors) {
	this->add(v);
      }
    }
  }

  // Stage 5: Filtering
  // Keep only 2-track, charge-neutral vertices that pass angular and
  // displacement significance cuts.
  void filteringImpl() {
    if (this->empty() || !primaryVertex_ || !primaryVertex_->isValid()) return;

    HYDDRA_DBG("[Leptonic] Final Filtering on " << this->size() << " vertices...\n");

    TrackVertexSetCollection finalVertices;

    for (const auto& vertex : *this) {
      if (vertex.size() != 2) continue;

      // Track angular cuts
      bool passAngular = true;
      for (const auto& trackRef : vertex) {
	const double trackCosTheta   = vertex.trackCosTheta(*primaryVertex_, trackRef);
	const double trackCosThetaCM = vertex.trackDecayAngleCM(trackRef);

	if (trackCosTheta < minTrackCosTheta_ ||
	    std::fabs(trackCosThetaCM) > (maxTrackCosThetaCM_Slope_ - trackCosTheta) ||
	    std::fabs(trackCosThetaCM) > maxTrackCosThetaCM_Limit_) {
	  passAngular = false;
	  break;
	}
      }
      if (!passAngular) continue;

      // Charge neutrality
      if (VertexHelper::CalculateTotalCharge(vertex) != 0) continue;

      // Displacement significance
      const double dxy      = VertexHelper::CalculateDxy(vertex, *primaryVertex_);
      const double dxyError = VertexHelper::CalculateDxyError(vertex, *primaryVertex_);
      if (dxyError <= 0) continue;

      if (dxy / dxyError > minDxySignificance_) {
	finalVertices.add(vertex);
      }
    }

    this->clear();
    for (const auto& v : finalVertices) {
      this->add(v);
    }

    HYDDRA_DBG("[Leptonic] Final count: " << this->size() << "\n");
  }

 private:

  // Cleaning thresholds
  double maxCompatibility_;
  double minCleanCosTheta_;

  // Post-cleaning kinematic cut (leptonic-specific)
  double minPostCosTheta_;

  // Final filtering thresholds
  double minTrackCosTheta_;
  double maxTrackCosThetaCM_Limit_;
  double maxTrackCosThetaCM_Slope_;

  // Returns true if the vertex fails kinematic selection.
  // Checked after cleaning and before disambiguation.
  bool failsKinematics(const TrackVertexSet& v) const {
    if (v.vertex().normalizedChi2() > maxNormChi2_) return true;
    if (v.cosTheta(*primaryVertex_) < minPostCosTheta_) return true;

    auto p4 = VertexHelper::GetVertex4Vector(v);
    double mass = p4.M();
    if (mass < minMass_) return true;

    double p = p4.P();
    double E = std::sqrt(p * p + mass * mass);
    double pOverE = (E > 1e-6) ? p / E : 0.0;
    if (pOverE < minPOverE_) return true;

    return false;
  }
};
