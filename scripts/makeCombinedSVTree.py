#!/usr/bin/env python3
"""
makeCombinedSVTree.py — Extract post-ID SVs from all SMS-GlGl diagnostic
files and merge them into a single flat ROOT TTree for interactive
exploration with svExplorer.py.

One row per SV surviving to the ID stage (stageIdx == 5).  Both hadronic
and leptonic SV types are included and tagged per row.  Raw signal
discriminants (matchRatio, min3D for hadronic; isGold for leptonic) are
stored so that svExplorer.py can apply configurable thresholds at run time.

If a StageVtx_eventIdx (or StageVtx_event) branch is found, per-event
hadronic/leptonic SV counts are computed and stored as evt_nHadSV /
evt_nLepSV.  Otherwise those branches are set to -1.

Usage:
    python3 scripts/makeCombinedSVTree.py \\
        --files "path/to/diag_*.root" --output combined_id.root [--workers 4]
"""

from __future__ import annotations

import argparse
import glob as glob_module
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

# ── Constants ──────────────────────────────────────────────────────────────────

_ID_STAGE = 5   # index of "id" in [seed, merged, cleaned, disambig, filtered, id]
_CAT_INT  = {"compressed": 0, "threshold": 1, "uncompressed": 2}

_SMS_RE = re.compile(
    r'mGl-(?P<mGl>\d+).*?mN2-(?P<mN2>\d+).*?mN1-(?P<mN1>\d+)'
    r'.*?N2ctau-(?P<ctau>[\dp]+)',
    re.IGNORECASE,
)


def _ctau_float(s: str) -> float:
    return float(s.replace('p', '.'))


def _parse_params(path: str) -> dict | None:
    m = _SMS_RE.search(os.path.basename(path))
    if m is None:
        return None
    mGl  = int(m.group('mGl'))
    mN2  = int(m.group('mN2'))
    mN1  = int(m.group('mN1'))
    ctau = _ctau_float(m.group('ctau'))
    dm   = mN2 - mN1
    if   dm <= 50:  cat = 'compressed'
    elif dm <= 100: cat = 'threshold'
    else:           cat = 'uncompressed'
    return dict(mGl=mGl, mN2=mN2, mN1=mN1, ctau=ctau, dm=dm, category=cat)


# ── Worker ─────────────────────────────────────────────────────────────────────

def _extract_svs(task: tuple) -> dict | None:
    """
    Open one file, extract hadronic + leptonic SVs at the ID stage.
    Returns a dict of {branch_name: np.ndarray} or None on failure.
    """
    path, params, run_lep, run_had = task
    import uproot
    import awkward as ak

    cat_int = _CAT_INT[params['category']]

    def _get(sv, key, default, dtype):
        if key in sv.fields:
            return ak.to_numpy(sv[key]).astype(dtype)
        return np.full(len(ak.to_numpy(sv[next(iter(sv.fields))])), default, dtype=dtype)

    def _find_event_idx(sv):
        for k in ("StageVtx_eventIdx", "StageVtx_event", "StageVtx_evtIdx"):
            if k in sv.fields:
                return ak.to_numpy(sv[k]).astype(np.int64)
        return None

    chunks = []

    try:
        with uproot.open(path) as f:

            # ── Hadronic ──────────────────────────────────────────────────────
            if run_had:
                had_base = "hyddraSVsHadronicDiagAnalyzer"
                try:
                    f[f"{had_base}/hadronicConfig"]
                    had_ok = True
                except KeyError:
                    had_ok = False

                if had_ok:
                    sv   = f[f"{had_base}/allStageVtx"].arrays(library="ak")
                    mask = ak.to_numpy(sv["StageVtx_stageIdx"]).astype(int) == _ID_STAGE
                    n    = int(mask.sum())
                    if n > 0:
                        mass = _get(sv, "StageVtx_mass",    np.nan, np.float32)[mask]
                        nt   = _get(sv, "StageVtx_nTracks", -1,     np.int16  )[mask]
                        mon  = np.where(nt > 0,
                                        mass / np.where(nt > 0, nt, 1).astype(np.float32),
                                        np.float32('nan'))
                        evt_idx = _find_event_idx(sv)
                        x = _get(sv, "StageVtx_x", np.nan, np.float32)[mask]
                        y = _get(sv, "StageVtx_y", np.nan, np.float32)[mask]
                        z = _get(sv, "StageVtx_z", np.nan, np.float32)[mask]
                        chunks.append({
                            "svType":          np.zeros(n, np.int8),
                            "category":        np.full(n, cat_int,            np.int8),
                            "mGl":             np.full(n, params['mGl'],      np.int16),
                            "mN2":             np.full(n, params['mN2'],      np.int16),
                            "mN1":             np.full(n, params['mN1'],      np.int16),
                            "deltaM":          np.full(n, params['dm'],       np.int16),
                            "ctau":            np.full(n, params['ctau'],     np.float32),
                            "matchRatio":      _get(sv, "StageVtx_matchRatio", np.nan, np.float32)[mask],
                            "min3D":           _get(sv, "StageVtx_min3D",     np.nan, np.float32)[mask],
                            "isGold":          np.zeros(n, bool),
                            "isMuon":          np.zeros(n, bool),
                            "isElectron":      np.zeros(n, bool),
                            "nTracks":         nt,
                            "massOverNTracks": mon,
                            "cosTheta":        _get(sv, "StageVtx_cosTheta",  np.nan, np.float32)[mask],
                            "decayAngle":      _get(sv, "StageVtx_decayAngle",np.nan, np.float32)[mask],
                            "pOverE":          _get(sv, "StageVtx_pOverE",    np.nan, np.float32)[mask],
                            "dxySignif":       _get(sv, "StageVtx_dxySignif", np.nan, np.float32)[mask],
                            "mass":            mass,
                            "dxy":             np.sqrt(x**2 + y**2),
                            "z":               z,
                            "evtIdx":          evt_idx[mask].astype(np.int64)
                                               if evt_idx is not None
                                               else np.full(n, -1, np.int64),
                        })

            # ── Leptonic ──────────────────────────────────────────────────────
            if run_lep:
                lep_base = "hyddraSVsDiagAnalyzer"
                try:
                    f[f"{lep_base}/leptonicConfig"]
                    lep_ok = True
                except KeyError:
                    lep_ok = False

                if lep_ok:
                    sv   = f[f"{lep_base}/allStageVtx"].arrays(library="ak")
                    mask = ak.to_numpy(sv["StageVtx_stageIdx"]).astype(int) == _ID_STAGE
                    n    = int(mask.sum())
                    if n > 0:
                        evt_idx = _find_event_idx(sv)
                        x = _get(sv, "StageVtx_x", np.nan, np.float32)[mask]
                        y = _get(sv, "StageVtx_y", np.nan, np.float32)[mask]
                        z = _get(sv, "StageVtx_z", np.nan, np.float32)[mask]
                        chunks.append({
                            "svType":          np.ones(n, np.int8),
                            "category":        np.full(n, cat_int,            np.int8),
                            "mGl":             np.full(n, params['mGl'],      np.int16),
                            "mN2":             np.full(n, params['mN2'],      np.int16),
                            "mN1":             np.full(n, params['mN1'],      np.int16),
                            "deltaM":          np.full(n, params['dm'],       np.int16),
                            "ctau":            np.full(n, params['ctau'],     np.float32),
                            "matchRatio":      np.full(n, np.nan,             np.float32),
                            "min3D":           np.full(n, np.nan,             np.float32),
                            "isGold":          _get(sv, "StageVtx_isGold",               False, bool)[mask],
                            "isMuon":          _get(sv, "StageVtx_passLooseMuonID",      False, bool)[mask],
                            "isElectron":      _get(sv, "StageVtx_passLooseElectronID",  False, bool)[mask],
                            "nTracks":         np.full(n, -1,                 np.int16),
                            "massOverNTracks": np.full(n, np.nan,             np.float32),
                            "cosTheta":        _get(sv, "StageVtx_cosTheta",  np.nan, np.float32)[mask],
                            "decayAngle":      _get(sv, "StageVtx_decayAngle",np.nan, np.float32)[mask],
                            "pOverE":          _get(sv, "StageVtx_pOverE",    np.nan, np.float32)[mask],
                            "dxySignif":       _get(sv, "StageVtx_dxySignif", np.nan, np.float32)[mask],
                            "mass":            _get(sv, "StageVtx_mass",      np.nan, np.float32)[mask],
                            "dxy":             np.sqrt(x**2 + y**2),
                            "z":               z,
                            "evtIdx":          evt_idx[mask].astype(np.int64)
                                               if evt_idx is not None
                                               else np.full(n, -1, np.int64),
                        })

    except Exception as e:
        print(f"  [warn] {os.path.basename(path)}: {e}", file=sys.stderr)
        return None

    if not chunks:
        return None
    combined = {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
    return combined


# ── Event-level counts ─────────────────────────────────────────────────────────

def _add_event_counts(data: dict) -> dict:
    """
    Compute per-event hadronic/leptonic SV counts and denormalise back to
    each row.  Requires evtIdx != -1; otherwise sets counts to -1.
    """
    import pandas as pd
    evt_idx = data["evtIdx"]
    n = len(evt_idx)
    if np.all(evt_idx == -1):
        data["evt_nHadSV"] = np.full(n, -1, np.int16)
        data["evt_nLepSV"] = np.full(n, -1, np.int16)
        return data

    sv_type = data["svType"]
    tmp = pd.DataFrame({
        "evtIdx": evt_idx,
        "isHad":  sv_type == 0,
        "isLep":  sv_type == 1,
    })
    counts = tmp.groupby("evtIdx").agg(
        evt_nHadSV=("isHad", "sum"),
        evt_nLepSV=("isLep", "sum"),
    ).astype(np.int16)
    merged = tmp[["evtIdx"]].join(counts, on="evtIdx")
    data["evt_nHadSV"] = merged["evt_nHadSV"].to_numpy()
    data["evt_nLepSV"] = merged["evt_nLepSV"].to_numpy()
    return data


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Merge post-ID SVs from all SMS-GlGl files into one ROOT TTree."
    )
    parser.add_argument("--files",   "-f", required=True,
                        help="Glob pattern for input ROOT files (quote it!)")
    parser.add_argument("--output",  "-o", default="combined_id.root")
    parser.add_argument("--workers", "-j", type=int, default=4)
    parser.add_argument("--mode", choices=["both", "hadronic", "leptonic"],
                        default="both")
    args = parser.parse_args()

    paths = sorted(glob_module.glob(args.files))
    if not paths:
        print(f"No files matched: {args.files}", file=sys.stderr)
        sys.exit(1)

    run_lep = args.mode in ("both", "leptonic")
    run_had = args.mode in ("both", "hadronic")

    tasks = []
    for p in paths:
        params = _parse_params(p)
        if params is None:
            print(f"  [skip] cannot parse params from: {os.path.basename(p)}")
            continue
        tasks.append((p, params, run_lep, run_had))

    print(f"Processing {len(tasks)} files with {args.workers} workers...")
    chunks = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_extract_svs, t): t[0] for t in tasks}
        for i, fut in enumerate(as_completed(futures), 1):
            p      = futures[fut]
            result = fut.result()
            if result is not None:
                chunks.append(result)
            print(f"  [{i:2d}/{len(tasks)}] {os.path.basename(p)}", flush=True)

    if not chunks:
        print("No SVs extracted.", file=sys.stderr)
        sys.exit(1)

    print("Merging...")
    data = {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
    data = _add_event_counts(data)

    n_had = int(np.sum(data["svType"] == 0))
    n_lep = int(np.sum(data["svType"] == 1))
    has_evt = not np.all(data["evtIdx"] == -1)
    print(f"Total SVs : {len(data['svType']):,}")
    print(f"  Hadronic: {n_had:,}")
    print(f"  Leptonic: {n_lep:,}")
    print(f"  Event grouping: {'yes' if has_evt else 'no (no evtIdx branch found)'}")

    import uproot
    print(f"Writing → {args.output}")
    with uproot.recreate(args.output) as out_f:
        out_f["svtree"] = data
    print("Done.")


if __name__ == "__main__":
    main()
