#!/usr/bin/env python3
"""
HyddraSV Diagnostic Printout

Reads an HyddraSVAnalyzer ntuple and prints a comprehensive text-based
diagnostic summary of the signal content — no plots, just terminal output.

Sections:
  1. File & Tree Info
  2. Gen Vertex Summary
  3. Gen Vertex Kinematics
  4. Reco SV Summary
  5. Reco SV Kinematics
  6. Matching Quality
  7. Signal Efficiency
  8. Track-Level Signal
  9. ID Pass Rates

Usage:
    python hyddraSVDiagnostic.py -i hyddraSV_ntuple.root
    python hyddraSVDiagnostic.py -i hyddraSV_ntuple.root --tree-path hyddraSVAnalyzer/tree
"""

import argparse
import numpy as np
import uproot


# ============================================================================
# Helpers
# ============================================================================

def safe_pct(num, den):
    """Return percentage, or 0 if denominator is zero."""
    return 100.0 * num / den if den > 0 else 0.0


def stat_line(arr, label, unit=''):
    """Return a formatted min/max/mean/median line for an array."""
    if len(arr) == 0:
        return f"  {label:20s}  (no entries)"
    u = f' {unit}' if unit else ''
    return (f"  {label:20s}  min={np.min(arr):10.4f}{u}  max={np.max(arr):10.4f}{u}  "
            f"mean={np.mean(arr):10.4f}{u}  median={np.median(arr):10.4f}{u}")


def flatten(jagged):
    """Flatten a numpy object-array of variable-length arrays."""
    return np.concatenate(jagged) if len(jagged) > 0 else np.array([])


def section(title):
    """Print a section header."""
    width = 72
    print()
    print('=' * width)
    print(f"  {title}")
    print('=' * width)


def has_branch(data, name):
    """Check whether a branch was loaded."""
    return name in data


# ============================================================================
# Loading
# ============================================================================

def load_data(filename, tree_path=None):
    """Load all relevant branches from a HyddraSVAnalyzer ROOT file."""
    file = uproot.open(filename)

    if tree_path:
        tree = file[tree_path]
        found_path = tree_path
    else:
        possible_paths = ['hyddraSVAnalyzer/tree', 'llpNanoSVAnalyzer/tree', 'tree']
        tree = None
        found_path = None
        for path in possible_paths:
            if path in file:
                tree = file[path]
                found_path = path
                break
        if tree is None:
            print(f"  Available keys: {file.keys()}")
            raise KeyError("Could not find tree. Try specifying --tree-path")

    available = tree.keys()

    # All branches we might want
    wanted = [
        # Gen vertex branches
        'HyddraGenVertex_nTotal', 'HyddraGenVertex_nElectron',
        'HyddraGenVertex_nMuon', 'HyddraGenVertex_nHadronic',
        'HyddraGenVertex_dxy', 'HyddraGenVertex_pt',
        'HyddraGenVertex_eta', 'HyddraGenVertex_phi',
        'HyddraGenVertex_mass',
        'HyddraGenVertex_isElectron', 'HyddraGenVertex_isMuon',
        'HyddraGenVertex_isHadronic',
        'HyddraGenVertex_passSelection', 'HyddraGenVertex_passSelectionAndCuts',
        'HyddraGenVertex_nTracks',
        # Reco SV branches
        'HyddraSV_nVertices', 'HyddraSV_isLeptonic', 'HyddraSV_nTracks',
        'HyddraSV_dxy', 'HyddraSV_mass', 'HyddraSV_pt',
        'HyddraSV_eta', 'HyddraSV_phi',
        'HyddraSV_chi2', 'HyddraSV_normalizedChi2',
        'HyddraSV_passLooseMuonID', 'HyddraSV_passLooseElectronID',
        'HyddraSV_passLooseID',
        'HyddraSV_genVertexIndex', 'HyddraSV_nearestGenVertexIndex',
        'HyddraSV_isBronze', 'HyddraSV_isSilver', 'HyddraSV_isGold',
        'HyddraSV_matchRatio', 'HyddraSV_min3D',
        # Track-level branches (may not be present)
        'HyddraSV_trackIsSignal', 'HyddraSV_trackIsSignalElectron',
        'HyddraSV_trackIsSignalMuon',
    ]

    branches_to_load = [b for b in wanted if b in available]
    missing = [b for b in wanted if b not in available]

    data = tree.arrays(branches_to_load, library='np')
    return data, found_path, available, missing


# ============================================================================
# Diagnostic sections
# ============================================================================

def print_file_info(filename, tree_path, available, missing):
    """Section 1: File & Tree Info."""
    section("1. File & Tree Info")
    print(f"  File:       {filename}")
    print(f"  Tree path:  {tree_path}")
    print(f"  Branches available: {len(available)}")
    if missing:
        print(f"  Branches missing (skipped): {len(missing)}")
        for b in missing:
            print(f"    - {b}")


def print_gen_vertex_summary(data):
    """Section 2: Gen Vertex Summary."""
    section("2. Gen Vertex Summary")

    if not has_branch(data, 'HyddraGenVertex_dxy'):
        print("  (no gen vertex branches found)")
        return

    n_events = len(data['HyddraGenVertex_dxy'])
    print(f"  Events: {n_events}")

    gen_dxy = data['HyddraGenVertex_dxy']
    n_total = sum(len(evt) for evt in gen_dxy)
    print(f"  Total gen vertices: {n_total}")

    # Per-event count stats
    counts = np.array([len(evt) for evt in gen_dxy])
    print(f"  Per-event count:  min={int(counts.min())}  max={int(counts.max())}  "
          f"mean={counts.mean():.2f}")

    # Breakdown by type
    if has_branch(data, 'HyddraGenVertex_isElectron'):
        n_ele = int(sum(np.sum(evt) for evt in data['HyddraGenVertex_isElectron']))
        n_mu = int(sum(np.sum(evt) for evt in data['HyddraGenVertex_isMuon']))
        n_had = int(sum(np.sum(evt) for evt in data['HyddraGenVertex_isHadronic']))
        print(f"  Electron:  {n_ele:8d}  ({safe_pct(n_ele, n_total):5.1f}%)")
        print(f"  Muon:      {n_mu:8d}  ({safe_pct(n_mu, n_total):5.1f}%)")
        print(f"  Hadronic:  {n_had:8d}  ({safe_pct(n_had, n_total):5.1f}%)")

    # passSelection / passSelectionAndCuts
    if has_branch(data, 'HyddraGenVertex_passSelection'):
        n_sel = int(sum(np.sum(evt) for evt in data['HyddraGenVertex_passSelection']))
        print(f"  passSelection:         {n_sel:8d}  ({safe_pct(n_sel, n_total):5.1f}%)")
    if has_branch(data, 'HyddraGenVertex_passSelectionAndCuts'):
        n_selc = int(sum(np.sum(evt) for evt in data['HyddraGenVertex_passSelectionAndCuts']))
        print(f"  passSelectionAndCuts:  {n_selc:8d}  ({safe_pct(n_selc, n_total):5.1f}%)")

    # nTotal/nElectron/nMuon/nHadronic (event-level counters)
    for branch, label in [
        ('HyddraGenVertex_nTotal', 'nTotal'),
        ('HyddraGenVertex_nElectron', 'nElectron'),
        ('HyddraGenVertex_nMuon', 'nMuon'),
        ('HyddraGenVertex_nHadronic', 'nHadronic'),
    ]:
        if has_branch(data, branch):
            arr = data[branch]
            print(f"  Event-level {label:12s}  "
                  f"min={int(np.min(arr))}  max={int(np.max(arr))}  mean={np.mean(arr):.2f}")


def print_gen_vertex_kinematics(data):
    """Section 3: Gen Vertex Kinematics."""
    section("3. Gen Vertex Kinematics")

    if not has_branch(data, 'HyddraGenVertex_dxy'):
        print("  (no gen vertex branches found)")
        return

    for branch, label, unit in [
        ('HyddraGenVertex_dxy', 'dxy', 'cm'),
        ('HyddraGenVertex_pt', 'pt', 'GeV'),
        ('HyddraGenVertex_eta', 'eta', ''),
        ('HyddraGenVertex_mass', 'mass', 'GeV'),
    ]:
        if has_branch(data, branch):
            arr = flatten(data[branch])
            print(stat_line(arr, label, unit))

    if has_branch(data, 'HyddraGenVertex_nTracks'):
        arr = flatten(data['HyddraGenVertex_nTracks'])
        if len(arr) > 0:
            print(f"  {'nTracks':20s}  min={int(np.min(arr)):10d}       max={int(np.max(arr)):10d}       "
                  f"mean={np.mean(arr):10.2f}       median={np.median(arr):10.1f}")


def print_reco_sv_summary(data):
    """Section 4: Reco SV Summary."""
    section("4. Reco SV Summary")

    if not has_branch(data, 'HyddraSV_dxy'):
        print("  (no reco SV branches found)")
        return

    n_events = len(data['HyddraSV_dxy'])
    sv_counts = np.array([len(evt) for evt in data['HyddraSV_dxy']])
    n_total = int(sv_counts.sum())
    print(f"  Events:     {n_events}")
    print(f"  Total SVs:  {n_total}")
    print(f"  SVs/event:  min={int(sv_counts.min())}  max={int(sv_counts.max())}  "
          f"mean={sv_counts.mean():.2f}")

    if has_branch(data, 'HyddraSV_isLeptonic'):
        is_lep = flatten(data['HyddraSV_isLeptonic']).astype(bool)
        n_lep = int(np.sum(is_lep))
        n_had = n_total - n_lep
        print(f"  Leptonic:   {n_lep:8d}  ({safe_pct(n_lep, n_total):5.1f}%)")
        print(f"  Hadronic:   {n_had:8d}  ({safe_pct(n_had, n_total):5.1f}%)")

    if has_branch(data, 'HyddraSV_nTracks'):
        nt = flatten(data['HyddraSV_nTracks'])
        print(f"  nTracks/SV: min={int(np.min(nt))}  max={int(np.max(nt))}  "
              f"mean={np.mean(nt):.2f}  median={np.median(nt):.0f}")

    if has_branch(data, 'HyddraSV_nVertices'):
        nv = data['HyddraSV_nVertices']
        print(f"  nVertices (event-level): min={int(np.min(nv))}  max={int(np.max(nv))}  "
              f"mean={np.mean(nv):.2f}")


def print_reco_sv_kinematics(data):
    """Section 5: Reco SV Kinematics."""
    section("5. Reco SV Kinematics")

    if not has_branch(data, 'HyddraSV_dxy'):
        print("  (no reco SV branches found)")
        return

    for branch, label, unit in [
        ('HyddraSV_dxy', 'dxy', 'cm'),
        ('HyddraSV_pt', 'pt', 'GeV'),
        ('HyddraSV_eta', 'eta', ''),
        ('HyddraSV_mass', 'mass', 'GeV'),
        ('HyddraSV_chi2', 'chi2', ''),
        ('HyddraSV_normalizedChi2', 'normalizedChi2', ''),
        ('HyddraSV_min3D', 'min3D', 'cm'),
    ]:
        if has_branch(data, branch):
            arr = flatten(data[branch])
            print(stat_line(arr, label, unit))


def print_matching_quality(data):
    """Section 6: Matching Quality."""
    section("6. Matching Quality")

    if not has_branch(data, 'HyddraSV_isGold'):
        print("  (no matching branches found)")
        return

    is_gold = flatten(data['HyddraSV_isGold']).astype(bool)
    is_silver = flatten(data['HyddraSV_isSilver']).astype(bool)
    is_bronze = flatten(data['HyddraSV_isBronze']).astype(bool)
    n_total = len(is_gold)

    # Leptonic breakdown
    if has_branch(data, 'HyddraSV_isLeptonic'):
        is_lep = flatten(data['HyddraSV_isLeptonic']).astype(bool)
        n_lep = int(np.sum(is_lep))
        n_gold_lep = int(np.sum(is_gold & is_lep))
        n_silver_lep = int(np.sum(is_silver & is_lep))
        n_bronze_lep = int(np.sum(is_bronze & is_lep))
        n_none_lep = n_lep - n_gold_lep - n_silver_lep - n_bronze_lep

        print(f"  Leptonic SVs ({n_lep} total):")
        print(f"    Gold:    {n_gold_lep:8d}  ({safe_pct(n_gold_lep, n_lep):5.1f}%)")
        print(f"    Silver:  {n_silver_lep:8d}  ({safe_pct(n_silver_lep, n_lep):5.1f}%)")
        print(f"    Bronze:  {n_bronze_lep:8d}  ({safe_pct(n_bronze_lep, n_lep):5.1f}%)")
        print(f"    None:    {n_none_lep:8d}  ({safe_pct(n_none_lep, n_lep):5.1f}%)")

        # Hadronic matchRatio
        is_had = ~is_lep
        n_had = int(np.sum(is_had))
        if n_had > 0 and has_branch(data, 'HyddraSV_matchRatio'):
            mr = flatten(data['HyddraSV_matchRatio'])
            mr_had = mr[is_had]
            valid = mr_had >= 0
            mr_valid = mr_had[valid]
            print(f"\n  Hadronic SVs ({n_had} total):")
            if len(mr_valid) > 0:
                print(stat_line(mr_valid, 'matchRatio'))
                n_r50 = int(np.sum(mr_valid >= 0.5))
                n_rany = int(np.sum(mr_valid > 0))
                print(f"    matchRatio >= 50%: {n_r50:8d}  ({safe_pct(n_r50, len(mr_valid)):5.1f}%)")
                print(f"    matchRatio >  0%:  {n_rany:8d}  ({safe_pct(n_rany, len(mr_valid)):5.1f}%)")
                print(f"    matchRatio == 0:   {int(len(mr_valid) - n_rany):8d}  "
                      f"({safe_pct(len(mr_valid) - n_rany, len(mr_valid)):5.1f}%)")
            else:
                print("    (no valid matchRatio entries)")
    else:
        # No leptonic flag — just print overall
        n_gold = int(np.sum(is_gold))
        n_silver = int(np.sum(is_silver))
        n_bronze = int(np.sum(is_bronze))
        print(f"  All SVs ({n_total} total):")
        print(f"    Gold:    {n_gold:8d}  ({safe_pct(n_gold, n_total):5.1f}%)")
        print(f"    Silver:  {n_silver:8d}  ({safe_pct(n_silver, n_total):5.1f}%)")
        print(f"    Bronze:  {n_bronze:8d}  ({safe_pct(n_bronze, n_total):5.1f}%)")


def print_signal_efficiency(data):
    """Section 7: Signal Efficiency (gen-to-reco matching)."""
    section("7. Signal Efficiency")

    if not has_branch(data, 'HyddraGenVertex_dxy') or not has_branch(data, 'HyddraSV_genVertexIndex'):
        print("  (requires both gen vertex and reco SV branches)")
        return

    n_events = len(data['HyddraGenVertex_dxy'])

    # Accumulate per-gen-vertex matching results
    gen_dxy_all, gen_type_all = [], []  # type: 0=electron, 1=muon, 2=hadronic
    gen_gold_match, gen_mr50_match, gen_mrany_match = [], [], []

    has_electron = has_branch(data, 'HyddraGenVertex_isElectron')
    has_muon = has_branch(data, 'HyddraGenVertex_isMuon')
    has_hadronic = has_branch(data, 'HyddraGenVertex_isHadronic')
    has_mr = has_branch(data, 'HyddraSV_matchRatio')

    for evt in range(n_events):
        gen_dxy = data['HyddraGenVertex_dxy'][evt]
        n_gen = len(gen_dxy)

        reco_gvi = data['HyddraSV_genVertexIndex'][evt]
        reco_gold = data['HyddraSV_isGold'][evt]

        # Build gold-matched gen indices
        gold_set = set()
        for ri in range(len(reco_gvi)):
            if int(reco_gvi[ri]) >= 0 and reco_gold[ri]:
                gold_set.add(int(reco_gvi[ri]))

        # Build matchRatio-matched gen indices
        mr50_set, mrany_set = set(), set()
        if has_mr and has_branch(data, 'HyddraSV_nearestGenVertexIndex'):
            reco_ngvi = data['HyddraSV_nearestGenVertexIndex'][evt]
            reco_mr = data['HyddraSV_matchRatio'][evt]
            for ri in range(len(reco_ngvi)):
                ngi = int(reco_ngvi[ri])
                if ngi >= 0:
                    if reco_mr[ri] >= 0.5:
                        mr50_set.add(ngi)
                    if reco_mr[ri] > 0:
                        mrany_set.add(ngi)

        for i in range(n_gen):
            gen_dxy_all.append(gen_dxy[i])

            gtype = -1
            if has_electron and data['HyddraGenVertex_isElectron'][evt][i]:
                gtype = 0
            elif has_muon and data['HyddraGenVertex_isMuon'][evt][i]:
                gtype = 1
            elif has_hadronic and data['HyddraGenVertex_isHadronic'][evt][i]:
                gtype = 2
            gen_type_all.append(gtype)

            gen_gold_match.append(i in gold_set)
            gen_mr50_match.append(i in mr50_set)
            gen_mrany_match.append(i in mrany_set)

    gen_dxy_all = np.array(gen_dxy_all)
    gen_type_all = np.array(gen_type_all)
    gen_gold_match = np.array(gen_gold_match)
    gen_mr50_match = np.array(gen_mr50_match)
    gen_mrany_match = np.array(gen_mrany_match)

    n_total = len(gen_dxy_all)

    def _eff_row(label, mask, match_arr):
        n = int(np.sum(mask))
        if n == 0:
            print(f"  {label:25s}  (no entries)")
            return
        n_pass = int(np.sum(match_arr[mask]))
        print(f"  {label:25s}  {n_pass:6d} / {n:6d}  = {safe_pct(n_pass, n):5.1f}%")

    print(f"  Total gen vertices: {n_total}")

    print("\n  --- Leptonic (isGold matching) ---")
    is_lep = (gen_type_all == 0) | (gen_type_all == 1)
    _eff_row("Leptonic overall", is_lep, gen_gold_match)
    _eff_row("  Electron", gen_type_all == 0, gen_gold_match)
    _eff_row("  Muon", gen_type_all == 1, gen_gold_match)

    print("\n  --- Hadronic (matchRatio) ---")
    is_had = gen_type_all == 2
    _eff_row("matchRatio >= 50%", is_had, gen_mr50_match)
    _eff_row("matchRatio >  0%", is_had, gen_mrany_match)


def print_track_signal(data):
    """Section 8: Track-Level Signal fractions."""
    section("8. Track-Level Signal")

    track_branches = [
        ('HyddraSV_trackIsSignal', 'isSignal'),
        ('HyddraSV_trackIsSignalElectron', 'isSignalElectron'),
        ('HyddraSV_trackIsSignalMuon', 'isSignalMuon'),
    ]

    found_any = False
    for branch, label in track_branches:
        if not has_branch(data, branch):
            continue
        found_any = True
        all_vals = flatten(data[branch])
        # Each entry is itself a per-track array
        if len(all_vals) > 0 and hasattr(all_vals[0], '__len__'):
            flat = np.concatenate(all_vals)
        else:
            flat = all_vals

        n = len(flat)
        if n == 0:
            print(f"  {label:25s}  (no entries)")
            continue
        n_sig = int(np.sum(flat.astype(bool)))
        print(f"  {label:25s}  {n_sig:8d} / {n:8d}  = {safe_pct(n_sig, n):5.1f}%")

    if not found_any:
        print("  (no track-level signal branches found)")


def print_id_pass_rates(data):
    """Section 9: ID Pass Rates (signal vs background SVs)."""
    section("9. ID Pass Rates")

    id_branches = [
        ('HyddraSV_passLooseMuonID', 'passLooseMuonID'),
        ('HyddraSV_passLooseElectronID', 'passLooseElectronID'),
        ('HyddraSV_passLooseID', 'passLooseID'),
    ]

    has_any_id = any(has_branch(data, b) for b, _ in id_branches)
    if not has_any_id:
        print("  (no ID branches found)")
        return

    # Determine signal vs background using isGold
    has_gold = has_branch(data, 'HyddraSV_isGold')
    has_lep = has_branch(data, 'HyddraSV_isLeptonic')
    has_mr = has_branch(data, 'HyddraSV_matchRatio')

    if has_gold:
        is_gold = flatten(data['HyddraSV_isGold']).astype(bool)
    if has_lep:
        is_lep = flatten(data['HyddraSV_isLeptonic']).astype(bool)
    if has_mr:
        mr = flatten(data['HyddraSV_matchRatio'])

    # Define "signal" as gold (leptonic) or matchRatio>0 (hadronic)
    n_total = len(flatten(data['HyddraSV_dxy']))
    if has_gold and has_lep and has_mr:
        is_signal = (is_gold & is_lep) | (~is_lep & (mr > 0))
    elif has_gold:
        is_signal = is_gold
    else:
        is_signal = None

    print(f"  {'ID':30s}  {'All':>18s}  {'Signal':>18s}  {'Background':>18s}")
    print(f"  {'-'*30}  {'-'*18}  {'-'*18}  {'-'*18}")

    for branch, label in id_branches:
        if not has_branch(data, branch):
            continue

        vals = flatten(data[branch]).astype(bool)
        n = len(vals)
        n_pass = int(np.sum(vals))
        all_str = f"{n_pass}/{n} ({safe_pct(n_pass, n):.1f}%)"

        if is_signal is not None:
            sig = vals[is_signal]
            bkg = vals[~is_signal]
            n_sig = len(sig)
            n_bkg = len(bkg)
            sig_pass = int(np.sum(sig))
            bkg_pass = int(np.sum(bkg))
            sig_str = f"{sig_pass}/{n_sig} ({safe_pct(sig_pass, n_sig):.1f}%)"
            bkg_str = f"{bkg_pass}/{n_bkg} ({safe_pct(bkg_pass, n_bkg):.1f}%)"
        else:
            sig_str = "N/A"
            bkg_str = "N/A"

        print(f"  {label:30s}  {all_str:>18s}  {sig_str:>18s}  {bkg_str:>18s}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='HyddraSV diagnostic printout — text-based summary of ntuple content')
    parser.add_argument('-i', '--input', required=True,
                        help='Input ROOT file (from HyddraSVAnalyzer)')
    parser.add_argument('--tree-path', default=None,
                        help='Path to tree in ROOT file (default: auto-detect)')
    args = parser.parse_args()

    print(f"\n{'*' * 72}")
    print(f"  HyddraSV Diagnostic Report")
    print(f"{'*' * 72}")

    data, tree_path, available, missing = load_data(args.input, args.tree_path)

    print_file_info(args.input, tree_path, available, missing)
    print_gen_vertex_summary(data)
    print_gen_vertex_kinematics(data)
    print_reco_sv_summary(data)
    print_reco_sv_kinematics(data)
    print_matching_quality(data)
    print_signal_efficiency(data)
    print_track_signal(data)
    print_id_pass_rates(data)

    print(f"\n{'*' * 72}")
    print(f"  Done.")
    print(f"{'*' * 72}\n")


if __name__ == "__main__":
    main()
