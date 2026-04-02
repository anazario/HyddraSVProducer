#!/usr/bin/env python3
"""
inspectGenParticles.py - Inspect gen particle content of a MiniAOD file.

Finds all particles that decay (directly or indirectly) to electrons and
prints their PDG IDs, status codes, and decay chains. Run this once on a
signal MiniAOD file to identify the correct motherPdgId for gen matching.

Requires CMSSW FWLite (run inside a cmsenv).

Usage:
    python scripts/inspectGenParticles.py file.root
    python scripts/inspectGenParticles.py file.root --max-events 5
    python scripts/inspectGenParticles.py file.root --all-particles
"""

import sys
import argparse

# FWLite setup
import ROOT
ROOT.gSystem.Load('libFWCoreFWLite.so')
ROOT.gSystem.Load('libDataFormatsFWLite.so')
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Events, Handle

# ---------------------------------------------------------------------------

def ancestors(p):
    """Return list of (pdgId, status) for all mothers going up the chain."""
    chain = []
    visited = set()
    queue = [p]
    while queue:
        current = queue.pop(0)
        for i in range(current.numberOfMothers()):
            m = current.mother(i)
            key = (m.pdgId(), id(m))
            if key not in visited:
                visited.add(key)
                chain.append((m.pdgId(), m.status(), m.pt()))
                queue.append(m)
    return chain

def has_electron_daughter(p, depth=0, max_depth=10):
    """Return True if p has an electron (pdgId 11 or -11) anywhere in its decay tree."""
    if depth > max_depth:
        return False
    for i in range(p.numberOfDaughters()):
        d = p.daughter(i)
        if abs(d.pdgId()) == 11:
            return True
        if has_electron_daughter(d, depth + 1, max_depth):
            return True
    return False

def print_decay_tree(p, indent=0, max_depth=6):
    """Print the decay tree rooted at p."""
    prefix = "  " * indent
    nd = p.numberOfDaughters()
    print(f"{prefix}pdgId={p.pdgId():8d}  status={p.status()}  "
          f"pt={p.pt():.3f}  eta={p.eta():.3f}  nDau={nd}")
    if indent < max_depth:
        for i in range(nd):
            print_decay_tree(p.daughter(i), indent + 1, max_depth)

# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='Input MiniAOD ROOT file (use file: prefix or plain path)')
    parser.add_argument('--max-events', type=int, default=3,
                        help='Number of events to inspect (default: 3)')
    parser.add_argument('--all-particles', action='store_true',
                        help='Print all gen particles, not just those with electron daughters')
    parser.add_argument('--gen-collection', default='prunedGenParticles',
                        help='Gen particle collection label (default: prunedGenParticles)')
    args = parser.parse_args()

    if args.input.startswith(('file:', 'root://')):
        fn = args.input
    else:
        fn = f'file:{args.input}'

    events = Events(fn)
    handle = Handle('std::vector<reco::GenParticle>')
    label  = args.gen_collection

    print(f"\nInspecting {args.max_events} event(s) from: {fn}")
    print(f"Gen collection: {label}\n")

    # Collect unique PDG IDs of particles with electron daughters across all events
    candidates = {}  # pdgId -> list of (status, pt) examples

    for i_ev, event in enumerate(events):
        if i_ev >= args.max_events:
            break

        event.getByLabel(label, handle)
        gens = handle.product()

        print(f"{'='*60}")
        print(f"Event {i_ev}  —  {gens.size()} gen particles")
        print(f"{'='*60}")

        if args.all_particles:
            print(f"\n{'PDG ID':>10}  {'Status':>8}  {'pT':>10}  {'eta':>8}  {'nMothers':>10}  {'nDaughters':>12}")
            print("-" * 65)
            for j in range(gens.size()):
                p = gens[j]
                print(f"{p.pdgId():>10}  {p.status():>8}  {p.pt():>10.3f}  "
                      f"{p.eta():>8.3f}  {p.numberOfMothers():>10}  {p.numberOfDaughters():>12}")
        else:
            # Find particles with electron daughters (direct BSM candidates)
            print("\nParticles with electron daughters:\n")
            found_any = False
            for j in range(gens.size()):
                p = gens[j]
                if has_electron_daughter(p):
                    found_any = True
                    pdg = p.pdgId()
                    if pdg not in candidates:
                        candidates[pdg] = []
                    candidates[pdg].append((p.status(), p.pt()))
                    print_decay_tree(p)
                    print()
            if not found_any:
                print("  (none found — try --all-particles to see everything)\n")

        print()

    if not args.all_particles and candidates:
        print(f"{'='*60}")
        print("SUMMARY — unique PDG IDs of particles with electron daughters:")
        print(f"{'='*60}")
        for pdg, examples in sorted(candidates.items(), key=lambda x: abs(x[0])):
            statuses = sorted(set(s for s, _ in examples))
            print(f"  pdgId = {pdg:8d}   status codes seen: {statuses}   "
                  f"occurrences: {len(examples)}")
        print()
        print("→ Set motherPdgId to the BSM particle's PDG ID in your config.")
        print("  Typical iDM candidates: dark photon (A'), mediator, or excited DM state.\n")

if __name__ == '__main__':
    main()
