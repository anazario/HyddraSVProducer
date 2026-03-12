#pragma once

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackVertexSetCollection.h"
#include "TrackHelper.h"

// Debug output toggle. Comment out to disable.
//#define HYDDRA_DEBUG_ENABLED

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

  // Stage enable flags
  bool doMerging_;
  bool doCleaning_;
  bool doDisambiguation_;
  bool doFiltering_;

  // Output option: refit final vertices with Kalman smoother to obtain
  // vertex-constrained track parameters stored in the output reco::Vertex.
  bool useVertexSmoothing_;

  // Event context (set per-event in run_reconstruction)
  const TransientTrackBuilder* ttBuilder_ = nullptr;
  const reco::Vertex* primaryVertex_ = nullptr;

  // Master list to track all vertices ever created (prevents duplicates across iterations)
  TrackVertexSetCollection masterList_;

  // Diagnostic counters
  struct DiagnosticCounters {
    size_t nSeeds = 0;
    size_t nAfterMerge = 0;
    size_t nAfterClean = 0;
    size_t nAfterDisambiguate = 0;
    size_t nAfterFilter = 0;
    size_t nMergeIterations = 0;
    size_t nTotalMerges = 0;
    size_t nMergeFailed = 0;
    size_t nMergeAlreadyExists = 0;
    // Forked pipeline counters
    size_t nTier0 = 0;
    size_t nTier1 = 0;
    size_t nDroppedMultiTrack = 0;
  };
  DiagnosticCounters diagnostics_;

  // Forked pipeline intermediate collections
  TrackVertexSetCollection seeds_;
  TrackVertexSetCollection tier0_;
  TrackVertexSetCollection tier1_;

 public:

  HYDDRABase(const edm::ParameterSet& pset) {
    seedCosThetaCut_    = pset.getParameter<double>("seedCosThetaCut");
    minMass_            = pset.getParameter<double>("minMass");
    minPOverE_          = pset.getParameter<double>("minPOverE");
    maxNormChi2_        = pset.getParameter<double>("maxNormChi2");
    minDxySignificance_ = pset.getParameter<double>("minDxySignificance");
    doMerging_            = pset.getParameter<bool>("doMerging");
    doCleaning_           = pset.getParameter<bool>("doCleaning");
    doDisambiguation_     = pset.getParameter<bool>("doDisambiguation");
    doFiltering_          = pset.getParameter<bool>("doFiltering");
    useVertexSmoothing_   = pset.getParameter<bool>("useVertexSmoothing");
  }

  // Access diagnostics
  const DiagnosticCounters& getDiagnostics() const { return diagnostics_; }

  // Snapshots of each intermediate stage for diagnostic purposes.
  struct PipelineSnapshots {
    reco::VertexCollection afterSeeding;
    reco::VertexCollection afterMerging;
    reco::VertexCollection afterCleaning;
    reco::VertexCollection afterDisambiguation;
    reco::VertexCollection afterFiltering;
  };

  // Diagnostic entry point. Runs the full pipeline and returns a snapshot of
  // the vertex collection after each stage. The final state is identical to
  // what run_reconstruction() would produce.
  PipelineSnapshots run_with_snapshots(const std::vector<reco::TrackRef>& tracks,
                                       const TransientTrackBuilder* builder,
                                       const reco::Vertex& pv) {
    ttBuilder_ = builder;
    primaryVertex_ = &pv;
    this->clear();
    masterList_.clear();
    diagnostics_ = DiagnosticCounters();

    Derived& self = static_cast<Derived&>(*this);
    PipelineSnapshots snaps;

    self.seedingImpl(tracks);
    diagnostics_.nSeeds = this->size();
    snaps.afterSeeding = this->vertices();

    if (doMerging_) self.mergingImpl();
    diagnostics_.nAfterMerge = this->size();
    snaps.afterMerging = this->vertices();

    if (doCleaning_) self.cleaningImpl();
    diagnostics_.nAfterClean = this->size();
    snaps.afterCleaning = this->vertices();

    if (doDisambiguation_) self.disambiguationImpl();
    diagnostics_.nAfterDisambiguate = this->size();
    snaps.afterDisambiguation = this->vertices();

    if (doFiltering_) self.filteringImpl();
    diagnostics_.nAfterFilter = this->size();
    snaps.afterFiltering = this->vertices();

    HYDDRA_DBG("[HYDDRA] Snapshot pipeline summary: seeds=" << diagnostics_.nSeeds
               << " -> merge=" << diagnostics_.nAfterMerge
               << " -> clean=" << diagnostics_.nAfterClean
               << " -> disambig=" << diagnostics_.nAfterDisambiguate
               << " -> filter=" << diagnostics_.nAfterFilter << "\n");

    return snaps;
  }

  // Accessor for Tier 1 (isolated) vertices from the forked pipeline.
  reco::VertexCollection isolatedVertices() const { return tier1_.vertices(); }

  // Forked pipeline entry point for EXONanoAOD.
  // Runs seeding once, then branches into two disambiguation paths:
  //   Tier 0 (inclusive): seeds -> disambiguation
  //   Tier 1 (isolated):  seeds -> merging -> drop 3+ track -> disambiguation
  // After the call, *this contains Tier 0 (inclusive) vertices.
  // Tier 1 vertices are accessible via isolatedVertices().
  void run_forked(const std::vector<reco::TrackRef>& tracks,
                  const TransientTrackBuilder* builder,
                  const reco::Vertex& pv) {
    ttBuilder_ = builder;
    primaryVertex_ = &pv;
    this->clear();
    masterList_.clear();
    diagnostics_ = DiagnosticCounters();

    Derived& self = static_cast<Derived&>(*this);

    // Stage 1: Seeding — computed once, shared by both tiers
    self.seedingImpl(tracks);
    diagnostics_.nSeeds = this->size();
    seeds_ = static_cast<const TrackVertexSetCollection&>(*this);

    // Tier 0: seeds -> disambiguation (most inclusive)
    if (doDisambiguation_) self.disambiguationImpl();
    diagnostics_.nTier0 = this->size();
    tier0_ = static_cast<const TrackVertexSetCollection&>(*this);

    // Tier 1: seeds -> merging -> drop 3+ track -> disambiguation (isolated)
    static_cast<std::set<TrackVertexSet>&>(*this) = seeds_;
    masterList_ = seeds_;
    if (doMerging_) self.mergingImpl();
    diagnostics_.nAfterMerge = this->size();
    dropMultiTrack();
    diagnostics_.nDroppedMultiTrack = diagnostics_.nAfterMerge - this->size();
    if (doDisambiguation_) self.disambiguationImpl();
    diagnostics_.nTier1 = this->size();
    tier1_ = static_cast<const TrackVertexSetCollection&>(*this);

    // Restore *this to Tier 0 as the primary output collection
    static_cast<std::set<TrackVertexSet>&>(*this) = tier0_;

    HYDDRA_DBG("[HYDDRA] Forked pipeline: seeds=" << diagnostics_.nSeeds
               << " | tier0=" << diagnostics_.nTier0
               << " | merged=" << diagnostics_.nAfterMerge
               << " dropped=" << diagnostics_.nDroppedMultiTrack
               << " tier1=" << diagnostics_.nTier1 << "\n");
  }

  // Main entry point. Runs the full 5-stage pipeline on the given tracks.
  void run_reconstruction(const std::vector<reco::TrackRef>& tracks,
			  const TransientTrackBuilder* builder,
			  const reco::Vertex& pv) {
    ttBuilder_ = builder;
    primaryVertex_ = &pv;
    this->clear();
    masterList_.clear();
    diagnostics_ = DiagnosticCounters();

    Derived& self = static_cast<Derived&>(*this);

    self.seedingImpl(tracks);
    diagnostics_.nSeeds = this->size();

    if (doMerging_) self.mergingImpl();
    diagnostics_.nAfterMerge = this->size();

    if (doCleaning_) self.cleaningImpl();
    diagnostics_.nAfterClean = this->size();

    if (doDisambiguation_) self.disambiguationImpl();
    diagnostics_.nAfterDisambiguate = this->size();

    if (doFiltering_) self.filteringImpl();
    diagnostics_.nAfterFilter = this->size();

    HYDDRA_DBG("[HYDDRA] Pipeline summary: seeds=" << diagnostics_.nSeeds
               << " -> merge=" << diagnostics_.nAfterMerge
               << " -> clean=" << diagnostics_.nAfterClean
               << " -> disambig=" << diagnostics_.nAfterDisambiguate
               << " -> filter=" << diagnostics_.nAfterFilter << "\n");
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
	
	TrackVertexSet seed({*x, *y}, ttBuilder_, useVertexSmoothing_);
	
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
	
	// Use add() to properly handle subsets
	this->add(seed);
	masterList_.add(seed);
      }
    }
    
    HYDDRA_DBG("[HYDDRA] Seeds: " << tracks.size() << " tracks, "
               << nPairs << " pairs, "
               << nOverlapping << " overlapping, "
               << nInvalid << " failed chi2, "
               << (this->size() + nFailedCosTheta) << " passed chi2, "
               << this->size() << " passed cosTheta\n");
  }
  
  // Iteratively merge vertices that share tracks and are spatially close
  // (distance significance < 4), until no more merges are possible.
  void mergeVertices() {
    if (this->empty()) return;
    
    int iteration = 0;
    bool madeChange = true;
    
    while (madeChange && iteration < 30) {
      madeChange = false;
      
      TrackVertexSetCollection ignoreList;
      TrackVertexSetCollection mergedList;
      int nMerged = 0, nMergeFailed = 0, nAlreadyExists = 0;
      
      // Snapshot to vector for safe pairwise iteration
      std::vector<TrackVertexSet> vtxVec(this->begin(), this->end());
      
      for (size_t i = 0; i < vtxVec.size(); ++i) {
	if (ignoreList.contains(vtxVec[i])) continue;
	
	for (size_t j = i + 1; j < vtxVec.size(); ++j) {
	  if (ignoreList.contains(vtxVec[j]) || ignoreList.contains(vtxVec[i])) continue;
	  if ((vtxVec[i] & vtxVec[j]) == 0) continue;
	  
	  if (vtxVec[i].distanceSignificance(vtxVec[j]) < 4) {
	    TrackVertexSet mergedVertex(vtxVec[i] + vtxVec[j]);
	    
	    if (isValidVertex(mergedVertex) && masterList_.doesNotContain(mergedVertex)) {
	      ignoreList.add(vtxVec[i]);
	      ignoreList.add(vtxVec[j]);
	      masterList_.add(mergedVertex);
	      mergedList.add(mergedVertex);
	      madeChange = true;
	      nMerged++;
	    } else if (masterList_.contains(mergedVertex)) {
	      // Merged vertex already exists - still mark parents for removal
	      ignoreList.add(vtxVec[i]);
	      ignoreList.add(vtxVec[j]);
	      nAlreadyExists++;
	    } else {
	      nMergeFailed++;
	    }
	  }
	}
      }
      
      // Build updated collection: survivors + newly merged
      TrackVertexSetCollection updated;
      for (const auto& v : vtxVec) {
	if (!ignoreList.contains(v)) {
	  updated.add(v);
	}
      }
      updated += mergedList;
      
      iteration++;
      diagnostics_.nTotalMerges += nMerged;
      diagnostics_.nMergeFailed += nMergeFailed;
      diagnostics_.nMergeAlreadyExists += nAlreadyExists;
      
      HYDDRA_DBG("[HYDDRA] Merge iter " << iteration << ": "
		 << this->size() << " -> " << updated.size()
		 << " | merged=" << nMerged
		 << " failed=" << nMergeFailed
		 << " exists=" << nAlreadyExists << "\n");
      
      this->clear();
      for (const auto& v : updated) {
	this->add(v);
      }
    }
    
    diagnostics_.nMergeIterations = iteration;
    HYDDRA_DBG("[HYDDRA] Merge result: " << this->size() << " vertices, "
	       << iteration << " iterations\n");
  }
  
  // A vertex is valid if the fit converged and normChi2 < 5.
  bool isValidVertex(const TrackVertexSet& set) const {
    return set.isValid() && set.normChi2() < 5;
  }

  // Remove any vertex with more than 2 tracks from the working collection.
  // Used in the Tier 1 path after merging to enforce the dilepton-only constraint.
  void dropMultiTrack() {
    TrackVertexSetCollection filtered;
    for (const auto& v : *this)
      if (v.size() == 2) filtered.add(v);
    static_cast<std::set<TrackVertexSet>&>(*this) = filtered;
  }
};
