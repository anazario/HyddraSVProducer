#pragma once

#include "TrackVertexSet.h"

class TrackVertexSetCollection : public std::set<TrackVertexSet> {

public:
  TrackVertexSetCollection() = default;

  TrackVertexSetCollection(const TrackVertexSetCollection &other) : std::set<TrackVertexSet>(other) {}

  virtual ~TrackVertexSetCollection() = default;

  bool contains(const TrackVertexSet &set) const;
  bool doesNotContain(const TrackVertexSet &set) const {return !contains(set); }
  bool isVertexUnique(const TrackVertexSet &set) const;
  TrackVertexSetCollection testUniqueness(const TrackVertexSet &set) const;

  TrackVertexSetCollection operator+(const TrackVertexSetCollection& rhs) const;
  TrackVertexSetCollection operator-(const TrackVertexSetCollection& rhs) const;
  TrackVertexSetCollection& operator+=(const TrackVertexSetCollection& rhs);
  TrackVertexSetCollection& operator-=(const TrackVertexSetCollection& rhs);

  void add(const TrackVertexSet &set);

  reco::VertexCollection vertices() const;
  reco::TrackCollection tracks() const;
  std::set<reco::TrackRef> completeTrackSet() const;
  std::set<reco::TrackRef> overlappingTracks() const;
  bool hasExclusiveVertices() const;
  size_t countOverlaps() const;
};
