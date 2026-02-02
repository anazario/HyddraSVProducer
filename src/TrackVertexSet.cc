#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackVertexSet.h"

TrackVertexSet::TrackVertexSet(const std::vector<reco::TrackRef> &init, const TransientTrackBuilder* ttBuilder) :
  std::set<reco::TrackRef>(init.begin(), init.end()),
  ttBuilder_(ttBuilder),
  fitter_(std::make_unique<KalmanVertexFitter>()) {
  fit();
}

TrackVertexSet::TrackVertexSet(std::initializer_list<reco::TrackRef> init, const TransientTrackBuilder* ttBuilder) :
  std::set<reco::TrackRef>(init),
  ttBuilder_(ttBuilder),
  fitter_(std::make_unique<KalmanVertexFitter>()) {
  fit();
}

TrackVertexSet::TrackVertexSet(const TrackVertexSet& other) :
  std::set<reco::TrackRef>(other),
  ttBuilder_(other.ttBuilder_),
  fitter_(std::make_unique<KalmanVertexFitter>())
{fit();}

// Transverse distance between two vertices
double TrackVertexSet::dxySeparation(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparation(this->vertex_, other.vertex_);
}

// Error in transverse distance between two vertices
double TrackVertexSet::dxySeparationError(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparationError(this->vertex_, other.vertex_);
}

// Significance in transverse distance between two vertices
double TrackVertexSet::dxySeparationSignificance(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparationSignificance(this->vertex_, other.vertex_);
}

// Distance between two vertices
double TrackVertexSet::distance(const TrackVertexSet& other) const {
  return VertexHelper::GetDistance(this->vertex_, other.vertex_);
}

// Error in the distance between two vertices
double TrackVertexSet::distanceError(const TrackVertexSet& other) const {
  return VertexHelper::GetDistanceError(this->vertex_, other.vertex_);
}

// Significance of the distance between two vertices: distance/distanceError
double TrackVertexSet::distanceSignificance(const TrackVertexSet& other) const {
  return VertexHelper::GetDistanceSignificance(this->vertex_, other.vertex_);
}

// Uses KalmanVertexTrackCompatibilityEstimator to check compatibility between this vertex and one of its tracks
double TrackVertexSet::compatibility(const reco::TrackRef &track) const {

  if(!contains(track)) {
    std::cout << "Warning in TrackVertexSet::compatibility: Given track is not in vertex!" << std::endl;
    return 999.;
  }
  const KalmanVertexTrackCompatibilityEstimator<5> estimator;
  auto result(estimator.estimate(*this, ttBuilder_->build(*track)));

  if(!result.first)
    throw std::runtime_error("Track compatibility estimation failed");

  return std::sqrt(result.second);
}

double TrackVertexSet::shiftDzAfterTrackRemoval(const reco::TrackRef &track) const {

  TrackVertexSet thisSetMinusTrack = *this;
  thisSetMinusTrack.removeTrack(track);

  return thisSetMinusTrack.isValid()? fabs(this->position().z() - thisSetMinusTrack.position().z()) : -1.;
}

double TrackVertexSet::shift3DAfterTrackRemoval(const reco::TrackRef &track) const {

  TrackVertexSet thisSetMinusTrack = *this;
  thisSetMinusTrack.removeTrack(track);

  return thisSetMinusTrack.isValid()? distance(thisSetMinusTrack) : -1.;
}

std::vector<reco::Track> TrackVertexSet::trackList() const {
  std::vector<reco::TrackRef> returnTracks(this->begin(), this->end());
  std::vector<reco::Track> tracks;
  for(const auto& trackRef : returnTracks)
    tracks.emplace_back(*trackRef);
  return tracks;
}

std::vector<reco::TrackRef> TrackVertexSet::tracks() const {

  std::vector<reco::TrackRef> returnTracks(this->begin(), this->end());
  std::sort(returnTracks.begin(), returnTracks.end(),
    [](const reco::TrackRef& a, const reco::TrackRef& b) {
      return a->pt() < b->pt();
    });
  return returnTracks;
}

// Function to return common tracks between two sets
std::vector<reco::TrackRef> TrackVertexSet::commonTracks(const TrackVertexSet& other) const {
  std::vector<reco::TrackRef> common;
  for (const auto& track : other) {
    if (this->contains(track)) {
      common.emplace_back(track);
    }
  }
  return common;
}

void TrackVertexSet::printTrackInfo() const {

  int count(0);
  std::cout << "\nThis set (normChi2 = " << normChi2() << ") has " << this->size() << " tracks: " << std::endl;
  for (const auto& track : *this) {
    std::cout << "\ttrack " << count << std::endl;
    std::cout << "\t\tpt: " << track->pt() << std::endl;
    std::cout << "\t\teta: " << track->eta() << std::endl;
    std::cout << "\t\tphi: " << track->phi() << std::endl;
    std::cout << "\t\tcompatibility: " << this->compatibility(track) << std::endl;
    count++;
  }
}

// Merge two TrackVertexSets together
TrackVertexSet TrackVertexSet::merge(const TrackVertexSet &other) const {
  return *this + other;
}

// Add a track to the set
void TrackVertexSet::addTrack(const reco::TrackRef& track) {
  this->insert(track);
  fit();
}

// Check if a track exists in the set
bool TrackVertexSet::contains(const reco::TrackRef& track) const {
  return this->find(track) != this->end();
}

// Remove a track from the set
void TrackVertexSet::removeTrack(const reco::TrackRef& track) {
  this->erase(track);
  if(this->size() > 1)
    fit();
  else {
    this->clear();
    vertex_ = TransientVertex();
  }
}

// Clear all tracks from the set
void TrackVertexSet::clearTracks() {
  this->clear();
}

TrackVertexSet& TrackVertexSet::operator=(const TrackVertexSet& other) {
  if (this != &other) {
    std::set<reco::TrackRef>::operator=(other);
    ttBuilder_ = other.ttBuilder_;
    fitter_ = std::make_unique<KalmanVertexFitter>(*other.fitter_);
    fit();
  }
  return *this;
}

// Overloaded operator to combine two sets into one
TrackVertexSet TrackVertexSet::operator+(const TrackVertexSet& other) const {
  TrackVertexSet result = *this; // Start with the current set
  for(const auto& track : other) {
    result.addTrack(track);
  }

  return result;
}

// Overloaded operator to count overlapping tracks between two sets
size_t TrackVertexSet::operator&(const TrackVertexSet& other) const {
  size_t count = 0;
  for (const auto& track : other) {
    if (this->contains(track)) {
      ++count;
    }
  }
  return count;
}

// Overloaded comparison operator for maintaining uniqueness in a std::set<TrackVertexSet>
bool TrackVertexSet::operator<(const TrackVertexSet& other) const {
    // Compare sizes first
    if (this->size() != other.size())
        return this->size() < other.size();

    // Use std::lexicographical_compare on the set's iterators directly
    return std::lexicographical_compare(this->begin(), this->end(),
                                        other.begin(), other.end());
}

// Overloaded equality operator to compare elements
bool TrackVertexSet::operator==(const TrackVertexSet& other) const {
  if(this->size() != other.size()) {
    return false;
  }
  return std::equal(this->begin(), this->end(), other.begin());
}

bool TrackVertexSet::operator!=(const TrackVertexSet& other) const {
  return(!(*this == other));
}

// Overloaded subset operator to check if all tracks in *this are in other
bool TrackVertexSet::operator|=(const TrackVertexSet& other) const {
  if (this->size() > other.size()) {
    return false;
  }
  for (const auto& track : *this) {
    if (!other.contains(track)) {
      return false;
    }
  }
  return true;
}

// private methods
void TrackVertexSet::fit() {
  vertex_ = this->size() < 2? TransientVertex() : fitter_->vertex(convertTracks());
}

std::vector<reco::TransientTrack> TrackVertexSet::convertTracks() const {

  std::vector<reco::TransientTrack> ttracks;
  for(const auto& track : *this) {
    ttracks.emplace_back(ttBuilder_->build(*track));
  }
  return ttracks;
}

double TrackVertexSet::calculateChiSquaredPValue(double chiSquaredValue, int degreesOfFreedom) const {

  // Define the chi-squared distribution with the specified degrees of freedom
  boost::math::chi_squared chiSqDist(degreesOfFreedom);

  // Compute the p-value as 1 - CDF(chiSquaredValue)
  double pValue = 1.0 - boost::math::cdf(chiSqDist, chiSquaredValue);

  return pValue;
}

// Direct casting from TrackVertexSet to reco::Vertex
TrackVertexSet::operator reco::Vertex() const {

  if(!vertex_.isValid())
    return reco::Vertex();

  reco::Vertex recoVertex(VertexHelper::ConvertFitVertex(vertex_));

  // Check for rare edge case where reco::Vertex can fail matrix inversion
  //leading to an exception when checking track compatibility
  int ifail(0);
  if(this->isValid()) {
    GlobalError err(recoVertex.covariance());
    err.matrix().Inverse(ifail);
  }

  if(!recoVertex.isValid() || ifail != 0)
    return reco::Vertex();

  for(const auto &trackRef : *this) {
    recoVertex.add(reco::TrackBaseRef(trackRef), recoVertex.trackWeight(reco::TrackBaseRef(trackRef)));
  }

  return recoVertex;
}
