#pragma once

#include <unordered_map>
#include <limits>
#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "KUCMSNtupleizer/HyddraSVProducer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/HyddraSVProducer/interface/TrackHelper.h"

class TrackVertexSet : public std::set<reco::TrackRef> {
 public:

  // Constructor with initializer list
  TrackVertexSet() = default;
  TrackVertexSet(const std::vector<reco::TrackRef> &init, const TransientTrackBuilder* ttBuilder);
  TrackVertexSet(std::initializer_list<reco::TrackRef> init, const TransientTrackBuilder* ttBuilder);
  
  // Copy constructor - explicitly inherit from base class
  TrackVertexSet(const TrackVertexSet& other);
  
  // Destructor
  ~TrackVertexSet() = default;
  
  TransientVertex tvertex() const { return vertex_; }
  reco::Vertex vertex() const {return reco::Vertex(*this);}
  
  bool isValid() const {return vertex_.isValid();}
  double pt() const {return VertexHelper::GetVertex4Vector(*this).pt();}
  double eta() const {return VertexHelper::GetVertex4Vector(*this).eta();}
  double phi() const {return VertexHelper::GetVertex4Vector(*this).phi();}
  double mass() const {return VertexHelper::GetVertex4Vector(*this).M();}
  double decayAngle() const {return VertexHelper::CalculateDecayAngle(*this);}
  double normChi2() const {return vertex_.normalisedChiSquared();}
  double pValue() const {return calculateChiSquaredPValue(vertex_.totalChiSquared(), vertex_.degreesOfFreedom());}
  double dxySeparation(const TrackVertexSet& other) const;
  double dxySeparationError(const TrackVertexSet& other) const;
  double dxySeparationSignificance(const TrackVertexSet& other) const;
  double distance(const TrackVertexSet& other) const;
  double distanceError(const TrackVertexSet& other) const;
  double distanceSignificance(const TrackVertexSet& other) const;
  double compatibility(const reco::TrackRef &track) const;
  double dxy(const reco::Vertex &primaryVertex)	const {return VertexHelper::CalculateDxy(*this, primaryVertex);}
  double dxyError(const reco::Vertex &primaryVertex) const {return VertexHelper::CalculateDxyError(*this, primaryVertex);}
  double cosTheta(const reco::Vertex &primaryVertex) const { return VertexHelper::CalculateCosTheta(primaryVertex, *this); }
  double trackCosTheta(const reco::Vertex &primaryVertex, const reco::TrackRef &track) const{ return TrackHelper::CalculateCosTheta(primaryVertex, *this, *track); }
  double trackDecayAngleCM(const reco::TrackRef &track) const { return VertexHelper::CalculateCMCosTheta(*this, *track); }
  double shiftDzAfterTrackRemoval(const reco::TrackRef &track) const;
  double shift3DAfterTrackRemoval(const reco::TrackRef &track) const;
  GlobalPoint position() const {return vertex_.position();}
  std::vector<reco::Track> trackList() const;
  std::vector<reco::TrackRef> tracks() const;
  std::vector<reco::TrackRef> commonTracks(const TrackVertexSet& other) const;
  
  void printTrackInfo() const;

  // element manipulation
  void clearTracks();
  void addTrack(const reco::TrackRef& track);
  void removeTrack(const reco::TrackRef& track);
  bool contains(const reco::TrackRef& track) const;
  
  // interaction with other TrackVertexSets
  bool isSmallerThan(const TrackVertexSet &other) const {return this->size() < other.size();}
  bool isLargerThan(const TrackVertexSet &other) const {return this->size() > other.size();}
  TrackVertexSet merge(const TrackVertexSet &other) const;

  TrackVertexSet& operator=(const TrackVertexSet& other); 
  TrackVertexSet operator+(const TrackVertexSet& other) const;
  size_t operator&(const TrackVertexSet& other) const;
  bool operator<(const TrackVertexSet& other) const;
  bool operator==(const TrackVertexSet& other) const;
  bool operator!=(const TrackVertexSet& other) const;
  bool operator|=(const TrackVertexSet& other) const;

  operator reco::Vertex() const;
  
 private:
  const TransientTrackBuilder* ttBuilder_;
  std::unique_ptr<KalmanVertexFitter> fitter_; 
  TransientVertex vertex_;
  
  void fit();
  std::vector<reco::TransientTrack> convertTracks() const;
  double calculateChiSquaredPValue(double chiSquaredValue, int degreesOfFreedom) const;
};


