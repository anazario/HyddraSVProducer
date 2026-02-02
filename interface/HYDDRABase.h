#pragma once

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackVertexSetCollection.h"
#include "TrackHelper.h"

// Debug output toggle. Comment out to disable.
#define HYDDRA_DEBUG_ENABLED

#ifdef HYDDRA_DEBUG_ENABLED
    #define HYDDRA_DBG(x) std::cout << x
#else
    #define HYDDRA_DBG(x)
#endif

// CRTP base class for HYDDRA vertex reconstruction.
// Provides the 5-stage pipeline (seed, merge, clean, disambiguate, filter)
// and common implementations for seeding and merging.
// Derived classes implement the stage-specific logic via *Impl() methods.
template <class Derived>
class HYDDRABase : public TrackVertexSetCollection {

 protected:

  // Shared cut thresholds (parsed from ParameterSet)
  double seedCosThetaCut_;
  double minMass_;
  double minPOverE_;
  double maxNormChi2_;
  double minDxySignificance_;

  // Event context (set per-event in run_reconstruction)
  const TransientTrackBuilder* ttBuilder_ = nullptr;
  const reco::Vertex* primaryVertex_ = nullptr;

 public:

  HYDDRABase(const edm::ParameterSet& pset) {
    seedCosThetaCut_    = pset.getParameter<double>("seedCosThetaCut");
    minMass_            = pset.getParameter<double>("minMass");
    minPOverE_          = pset.getParameter<double>("minPOverE");
    maxNormChi2_        = pset.getParameter<double>("maxNormChi2");
    minDxySignificance_ = pset.getParameter<double>("minDxySignificance");
  }

  // Main entry point. Runs the full 5-stage pipeline on the given tracks.
  void run_reconstruction(const std::vector<reco::TrackRef>& tracks,
			  const TransientTrackBuilder* builder,
			  const reco::Vertex& pv) {
    ttBuilder_ = builder;
    primaryVertex_ = &pv;
    this->clear();

    Derived& self = static_cast<Derived&>(*this);
    self.seedingImpl(tracks);
    self.mergingImpl();
    self.cleaningImpl();
    self.disambiguationImpl();
    self.filteringImpl();
  }

 protected:

  // Form all valid 2-track seed vertices from the input tracks.
  // Seeds are filtered by chi2 validity and cosTheta alignment with the PV.
  void generateSeeds(const std::vector<reco::TrackRef>& tracks) {
    if (tracks.size() < 2) return;

    size_t nPairs = 0;
    size_t nOverlapping = 0;
    size_t nInvalid = 0;
    size_t nFailedCosTheta = 0;

    auto end = tracks.end();
    auto endm1 = end - 1;

    for (auto x = tracks.begin(); x != endm1; ++x) {
      for (auto y = x + 1; y != end; ++y) {
	nPairs++;

	if (TrackHelper::OverlappingTrack(**x, **y, ttBuilder_)) {
	  nOverlapping++;
	  continue;
	}

	TrackVertexSet seed({*x, *y}, ttBuilder_);

	if (!isValidVertex(seed)) {
	  nInvalid++;
	  continue;
	}

	// CosTheta alignment cut (check before insert for speed)
	if (primaryVertex_ && primaryVertex_->isValid()) {
	  if (seed.cosTheta(*primaryVertex_) < seedCosThetaCut_) {
	    nFailedCosTheta++;
	    continue;
	  }
	}

	this->insert(seed);
      }
    }

    size_t nPassedChi2 = this->size() + nFailedCosTheta;
    HYDDRA_DBG("[HYDDRA] Seeds: " << tracks.size() << " tracks, "
               << nPairs << " pairs, "
               << nOverlapping << " overlapping, "
               << nInvalid << " failed chi2, "
               << nPassedChi2 << " passed chi2, "
               << this->size() << " passed cosTheta\n");
  }

  // Iteratively merge vertices that share tracks and are spatially close
  // (distance significance < 4), until no more merges are possible.
  void mergeVertices() {
    if (this->empty()) return;

    int iteration = 0;
    bool madeChange = true;

    while (madeChange) {
      madeChange = false;

      TrackVertexSetCollection deadList;
      TrackVertexSetCollection mergedVertices;
      int nMerged = 0, nMergeFailed = 0;

      // Snapshot to vector for safe pairwise iteration
      std::vector<TrackVertexSet> vtxVec(this->begin(), this->end());

      for (size_t i = 0; i < vtxVec.size(); ++i) {
	if (deadList.contains(vtxVec[i])) continue;

	for (size_t j = i + 1; j < vtxVec.size(); ++j) {
	  if (deadList.contains(vtxVec[j]) || deadList.contains(vtxVec[i])) continue;
	  if ((vtxVec[i] & vtxVec[j]) == 0) continue;

	  if (vtxVec[i].distanceSignificance(vtxVec[j]) < 4) {
	    TrackVertexSet mergedVertex(vtxVec[i] + vtxVec[j]);

	    if (isValidVertex(mergedVertex)) {
	      deadList.insert(vtxVec[i]);
	      deadList.insert(vtxVec[j]);
	      mergedVertices.add(mergedVertex);
	      madeChange = true;
	      nMerged++;
	    } else {
	      nMergeFailed++;
	    }
	  }
	}
      }

      // Build updated collection: survivors + newly merged
      TrackVertexSetCollection updated;
      for (const auto& v : vtxVec) {
	if (!deadList.contains(v)) {
	  updated.insert(v);
	}
      }
      updated += mergedVertices;

      iteration++;
      HYDDRA_DBG("[HYDDRA] Merge iter " << iteration << ": "
		 << this->size() << " -> " << updated.size()
		 << " | merged=" << nMerged
		 << " failed=" << nMergeFailed << "\n");

      this->clear();
      this->insert(updated.begin(), updated.end());
    }

    HYDDRA_DBG("[HYDDRA] Merge result: " << this->size() << " vertices, "
	       << iteration << " iterations\n");
  }

  // A vertex is valid if the fit converged and normChi2 < 5.
  bool isValidVertex(const TrackVertexSet& set) const {
    return set.isValid() && set.normChi2() < 5;
  }
};
