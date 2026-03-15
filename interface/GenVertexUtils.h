#pragma once

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"

namespace HyddraUtils {

/// Return all stable charged hadrons from packed gen particles that descend
/// from the same Z boson as the given gen vertex.  Used in miniAOD mode where
/// prunedGenParticles does not contain final-state hadrons.
inline reco::GenParticleCollection
getStableChargedDaughtersFromPacked(
    const GenVertex& genVertex,
    const std::vector<pat::PackedGenParticle>& packed)
{
  const reco::Candidate* genZ = genVertex.genPair().first.mother();
  reco::GenParticleCollection result;
  if (!genZ) return result;

  for (const auto& p : packed) {
    if (p.status() != 1 || p.charge() == 0) continue;
    const reco::Candidate* mom = p.mother(0);
    if (!mom) continue;

    bool isFromSameZ = false;
    const reco::Candidate* prev = nullptr;
    while (mom && mom != prev) {
      if (mom->pdgId() == 23) {
        if (mom == genZ) isFromSameZ = true;
        break;
      }
      prev = mom;
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }
    if (isFromSameZ) {
      result.emplace_back(reco::GenParticle(
          p.charge(), p.p4(), p.vertex(), p.pdgId(), p.status(), true));
    }
  }
  return result;
}

} // namespace HyddraUtils
