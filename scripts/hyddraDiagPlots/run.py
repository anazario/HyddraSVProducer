#!/usr/bin/env python3
"""
run.py — entry point for hyddraDiagPlots.

Usage (single file):
    python3 scripts/hyddraDiagPlots/run.py --signal signal.root

Usage (multiple files / glob):
    python3 scripts/hyddraDiagPlots/run.py --signal "diag_*.root" -j 4

Usage (with background file):
    python3 scripts/hyddraDiagPlots/run.py --signal signal.root --background bkg.root

Output:
  Single-file mode  : plots written at root level of --output file,
                      organised under seeding/merging/cleaning/disambiguation/filtering/summary
  Multi-file mode   : per-file stem subdirectories, each with the 6 stage subdirs

Cut values for plot annotations are read automatically from the leptonicConfig
tree stored inside each signal ROOT file.  No manual flags are required.
"""

import argparse
import contextlib
import glob as glob_module
import multiprocessing as mp
import os
import sys
import tempfile

import awkward as ak
import uproot
import ROOT
from tqdm import tqdm

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

# Import package modules
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

from hyddraDiagPlots.src    import config, loader
from hyddraDiagPlots.stages import (
    summary, seeding, merging, cleaning, disambiguation, filtering
)

_STAGE_KEYS   = ['seed',    'merged',  'cleaned',  'disambig',       'filtered']
_STAGE_LABELS = ['Seeding', 'Merging', 'Cleaning', 'Disambiguation', 'Filtering']


# ── Efficiency summary helpers ────────────────────────────────────────────────

def _compute_summary(stem, gf_sig, sc_sig):
    """Compute per-file signal efficiency summary from loaded data."""
    n_gen = int(ak.sum(gf_sig['GenFunnel_n']))
    rows = []
    for key, label in zip(_STAGE_KEYS, _STAGE_LABELS):
        n_gold_gen   = int(ak.sum(ak.flatten(gf_sig[f'GenFunnel_gold_{key}'])))
        n_silver_gen = int(ak.sum(ak.flatten(gf_sig[f'GenFunnel_silver_{key}'])))
        rows.append({
            'label':       label,
            'n_reco':      sc_sig.get(f'Stage_n_{key}',     0),
            'n_gold_reco': sc_sig.get(f'Stage_nGold_{key}', 0),
            'n_gold_gen':  n_gold_gen,
            'n_gs_gen':    n_gold_gen + n_silver_gen,
        })
    return {'stem': stem, 'n_gen': n_gen, 'rows': rows}


def _print_summaries(summaries):
    """Print a clean per-file efficiency table for all processed files."""
    # Sort by file stem for deterministic output order
    summaries = sorted(summaries, key=lambda s: s['stem'])

    col_w = [16, 10, 11, 12, 16]  # Stage, Reco vtx, Gold reco, Eff(gold), Eff(gold+sil)
    divider = '  ' + '-' * (sum(col_w) + len(col_w) * 3 - 1)

    hdr = (f"  {'Stage':<{col_w[0]}} {'Reco vtx':>{col_w[1]}} "
           f"{'Gold reco':>{col_w[2]}} {'Eff (gold)':>{col_w[3]}} "
           f"{'Eff (gld+sil)':>{col_w[4]}}")

    print(f"\n{'=' * 72}")
    print(f"  HYDDRA diagnostic summary")
    print(f"{'=' * 72}")

    for s in summaries:
        n_gen = s['n_gen']
        print(f"\n  File : {s['stem']}")
        print(f"  Gen signal vertices : {n_gen}")
        print(divider)
        print(hdr)
        print(divider)
        for row in s['rows']:
            if n_gen > 0:
                gold_eff = f"{row['n_gold_gen'] / n_gen * 100:.1f}%"
                gs_eff   = f"{row['n_gs_gen']   / n_gen * 100:.1f}%"
            else:
                gold_eff = gs_eff = 'N/A'
            print(f"  {row['label']:<{col_w[0]}} "
                  f"{row['n_reco']:>{col_w[1]}} "
                  f"{row['n_gold_reco']:>{col_w[2]}} "
                  f"{gold_eff:>{col_w[3]}} "
                  f"{gs_eff:>{col_w[4]}}")
        print(divider)

    print()


# ── Core plot runner ──────────────────────────────────────────────────────────

def _run_all_plots(out_file, sig_path, bkg_path):
    """Load data, write all stage plots into out_file, return efficiency summary."""
    stem = os.path.splitext(os.path.basename(sig_path))[0]

    with uproot.open(sig_path) as sig_f:
        gf_sig = loader.load_gen_funnel(sig_f)
        sc_sig = loader.load_stage_counts(sig_f)
        ct_sig = loader.load_cleaning_tracks(sig_f)
        sv_sig = loader.load_all_stage_vtx(sig_f)
        st_sig = loader.load_seed_tracks(sig_f)
        fv_sig = loader.load_filtering_vtx(sig_f)
        cfg    = loader.load_leptonic_config(sig_f)

    if cfg is None:
        print(f"  [{stem}] Warning: leptonicConfig not found; cut annotations will use defaults",
              file=sys.stderr)

    sc_bkg, sv_bkg = None, None
    if bkg_path:
        with uproot.open(bkg_path) as bkg_f:
            sc_bkg = loader.load_stage_counts(bkg_f)
            sv_bkg = loader.load_all_stage_vtx(bkg_f)

    stage_dirs = {d: out_file.mkdir(d) for d in config.STAGE_DIRS}

    stages = [
        ("seeding",        lambda: seeding.make_plots(      stage_dirs["seeding"],        gf_sig, sv_sig, sv_bkg, st_sig)),
        ("merging",        lambda: merging.make_plots(       stage_dirs["merging"],        gf_sig, sv_sig, sv_bkg)),
        ("cleaning",       lambda: cleaning.make_plots(      stage_dirs["cleaning"],       gf_sig, sv_sig, sv_bkg, ct_sig, cfg)),
        ("disambiguation", lambda: disambiguation.make_plots(stage_dirs["disambiguation"], gf_sig, sv_sig, sv_bkg)),
        ("filtering",      lambda: filtering.make_plots(     stage_dirs["filtering"],      gf_sig, sv_sig, sv_bkg, fv_sig, cfg)),
        ("summary",        lambda: summary.make_plots(       stage_dirs["summary"],        gf_sig, sc_sig, sc_bkg)),
    ]

    for stage_name, fn in tqdm(stages, desc=f"  {stem}", unit="stage", leave=False):
        fn()

    return _compute_summary(stem, gf_sig, sc_sig)


# ── Multi-file helpers ────────────────────────────────────────────────────────

def _worker(args_tuple):
    """Multiprocessing worker: process one signal file into a temp ROOT file."""
    sig_path, bkg_path, tmp_path = args_tuple
    with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
        tmp_file = ROOT.TFile(tmp_path, "RECREATE")
        summ = _run_all_plots(tmp_file, sig_path, bkg_path)
        tmp_file.Close()
    return tmp_path, summ


def _copy_to_dir(src_path, tdir):
    """Recursively copy all objects and subdirectories from src ROOT file into tdir."""
    src = ROOT.TFile.Open(src_path, "READ")
    if not src or src.IsZombie():
        print(f"  Warning: could not open temp file {src_path}", file=sys.stderr)
        return

    def _copy_dir(src_tdir, dst_tdir):
        dst_tdir.cd()
        for key in src_tdir.GetListOfKeys():
            obj = key.ReadObj()
            cls_name = key.GetClassName()
            if "TDirectory" in cls_name:
                sub = dst_tdir.mkdir(key.GetName())
                _copy_dir(obj, sub)
            else:
                dst_tdir.cd()
                obj.Write(key.GetName())

    _copy_dir(src, tdir)
    src.Close()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="HYDDRA diagnostic plots — full per-stage analysis"
    )
    parser.add_argument("--signal", required=True, nargs="+",
                        help="Signal ROOT file(s) or glob pattern(s)")
    parser.add_argument("--background", default=None,
                        help="Background ROOT file (optional, applied to all signal files)")
    parser.add_argument("--output", default="hyddra_diag_plots.root",
                        help="Output ROOT file (default: hyddra_diag_plots.root)")
    parser.add_argument("--jobs", "-j", type=int, default=0,
                        help="Max parallel workers for multi-file mode (0 = cpu_count)")
    args = parser.parse_args()

    # Expand glob patterns and deduplicate
    sig_files = []
    for pat in args.signal:
        expanded = sorted(glob_module.glob(pat))
        sig_files.extend(expanded if expanded else [pat])
    sig_files = list(dict.fromkeys(sig_files))

    if not sig_files:
        print("Error: no signal files found.", file=sys.stderr)
        sys.exit(1)

    print(f"[hyddraDiagPlots] {len(sig_files)} signal file(s)  |  "
          f"background: {args.background or 'none'}  |  output: {args.output}")

    summaries = []

    if len(sig_files) == 1:
        # ── Single-file mode ──────────────────────────────────────────────────
        out_file = ROOT.TFile(args.output, "RECREATE")
        summ = _run_all_plots(out_file, sig_files[0], args.background)
        out_file.Close()
        summaries.append(summ)

    else:
        # ── Multi-file mode: parallel workers → temp files → merge ────────────
        n_workers = args.jobs if args.jobs > 0 else min(len(sig_files), mp.cpu_count())

        tmpdir     = tempfile.mkdtemp(prefix="hyddra_diag_")
        work_items = []
        for sig_path in sig_files:
            stem    = os.path.splitext(os.path.basename(sig_path))[0]
            tmp_out = os.path.join(tmpdir, f"{stem}.root")
            work_items.append((sig_path, args.background, tmp_out))

        with mp.Pool(n_workers) as pool:
            with tqdm(total=len(work_items), desc="Processing", unit="file") as pbar:
                for _, summ in pool.imap_unordered(_worker, work_items):
                    summaries.append(summ)
                    pbar.update()

        # Merge: each file gets its own top-level subdirectory
        out_file = ROOT.TFile(args.output, "RECREATE")
        for sig_path, _, tmp_path in work_items:
            stem = os.path.splitext(os.path.basename(sig_path))[0]
            tdir = out_file.mkdir(stem)
            _copy_to_dir(tmp_path, tdir)
            os.unlink(tmp_path)
        out_file.Close()
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass

    _print_summaries(summaries)
    print(f"[hyddraDiagPlots] Plots saved to {args.output}")


if __name__ == "__main__":
    main()
