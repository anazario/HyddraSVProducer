#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackVertexSetCollection.h"

bool TrackVertexSetCollection::contains(const TrackVertexSet &set) const {
  return this->find(set) != this->end();
}

bool TrackVertexSetCollection::isVertexUnique(const TrackVertexSet &set) const {

  for (auto it = this->begin(); it != this->end(); ++it) {
    if(*it != set && (*it & set) > 0)
      return false;
  }

  return true;
}

TrackVertexSetCollection TrackVertexSetCollection::testUniqueness(const TrackVertexSet &set) const {

  TrackVertexSetCollection retCollection;
  if(!isVertexUnique(set)) {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if(*it != set && (*it & set) > 0)
	retCollection.add(*it);
    }
  }

  return retCollection;
}

// Add a set to the collection, avoiding subsets and replacing subsets with larger sets
void TrackVertexSetCollection::add(const TrackVertexSet &set) {

  if (!set.isValid())
    return;

  // Iterate over existing sets to check for subsets
  for (auto it = this->begin(); it != this->end();) {
    if (*it |= set) { // If the new set is a superset of an existing set
      it = this->erase(it); // Remove the smaller subset
    } else if (set |= *it) { // If an existing set is a superset of the new set
      return; // Do not insert the new set
    } else {
      ++it;
    }
  }

  // Insert the new set if it is not a subset or superset
  this->insert(set);
}

// Subtraction operator: returns a new collection with unique elements from lhs not in rhs
TrackVertexSetCollection TrackVertexSetCollection::operator-(const TrackVertexSetCollection& rhs) const {
  TrackVertexSetCollection result = *this; // Copy current object
  result -= rhs;                           // Use -= to handle subtraction logic
  return result;                           // Return the result
}

// Subtraction assignment operator as a member function
TrackVertexSetCollection& TrackVertexSetCollection::operator-=(const TrackVertexSetCollection& rhs) {
  TrackVertexSetCollection result;

  for (const auto& set : *this) {
    if (rhs.doesNotContain(set)) { // Keep only elements not in rhs
      result.add(set);           // Use add to ensure subset handling
    }
  }

  *this = std::move(result); // Replace the current object with the result
  return *this;
}

// Addition operator: combines both collections, ensuring no subsets
TrackVertexSetCollection TrackVertexSetCollection::operator+(const TrackVertexSetCollection& rhs) const {
  TrackVertexSetCollection result = *this; // Copy current object
  result += rhs;                           // Use += to handle addition logic
  return result;                           // Return the result
}

TrackVertexSetCollection& TrackVertexSetCollection::operator+=(const TrackVertexSetCollection& rhs) {
  for (const auto& set : rhs) {
    this->add(set); // Use add to ensure no subsets are inserted
  }
  return *this;
}

reco::VertexCollection TrackVertexSetCollection::vertices() const {

  reco::VertexCollection uniqueVertices;
  for(const auto &set : *this) {
    if(set.isValid())
      uniqueVertices.emplace_back(set);
  }

  return uniqueVertices;
}

reco::TrackCollection TrackVertexSetCollection::tracks() const {

  reco::TrackCollection tracks;
  for(const auto &trackRef : this->completeTrackSet())
    tracks.emplace_back(*trackRef);

  return tracks;
}

std::set<reco::TrackRef> TrackVertexSetCollection::completeTrackSet() const {

  std::set<reco::TrackRef> allTracks;
  for(const auto &vertex : *this) {
    for(const auto &track : vertex)
      allTracks.insert(track);
  }

  return allTracks;
}

std::set<reco::TrackRef> TrackVertexSetCollection::overlappingTracks() const {
  // Step 1: Track occurrence count
  std::map<reco::TrackRef, int> trackCount;

  // Step 2: Iterate over all TrackVertexSets to count each track
  for (const auto& vertexSet : *this) {
    for (const auto& track : vertexSet.tracks()) {
      trackCount[track]++;
    }
  }

  // Step 3: Collect tracks that appear more than once
  std::set<reco::TrackRef> overlappingTracks;
  for (const auto& trackEntry : trackCount) {
    if (trackEntry.second > 1) {
      overlappingTracks.insert(trackEntry.first);
    }
  }

  return overlappingTracks;
}

bool TrackVertexSetCollection::hasExclusiveVertices() const {

  std::set<reco::TrackRef> encounteredTracks;

  // Iterate through each TrackVertexSet in the collection
  for (const auto& trackVertexSet : *this) {
    // Iterate through each track in the TrackVertexSet
    for (const auto& track : trackVertexSet) {
      // Check if the track is already in the encountered set
      if (encounteredTracks.find(track) != encounteredTracks.end()) {
	return false; // Duplicate found
      }
      // Add the track to the encountered set
      encounteredTracks.insert(track);
    }
  }

  return true;
}

size_t TrackVertexSetCollection::countOverlaps() const {

  size_t overlaps = 0;

  // Convert to vector for random access iteration
  std::vector<std::reference_wrapper<const TrackVertexSet>> sets(this->begin(), this->end());

  // Use indices for clearer bounds and potentially better optimization
  for (size_t i = 0; i < sets.size() - 1; ++i) {
    for (size_t j = i + 1; j < sets.size(); ++j) {
      if ((sets[i].get() & sets[j].get()) > 0) {
	overlaps++;
      }
    }
  }
  return overlaps;
}
