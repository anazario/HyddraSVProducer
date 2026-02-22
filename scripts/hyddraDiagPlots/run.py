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
"""

import argparse
import glob as glob_module
import multiprocessing as mp
import os
import sys
import tempfile

import uproot
import ROOT

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


# ── Core plot runner ──────────────────────────────────────────────────────────

def _run_all_plots(out_file, sig_path, bkg_path, max_compat, min_cos_theta,
                   use_diagonal_cut=False, clean_cut_slope=0.0):
    """Load data and write all stage plots into out_file."""
    stem = os.path.splitext(os.path.basename(sig_path))[0]
    print(f"  [{stem}] Loading signal data...")

    with uproot.open(sig_path) as sig_f:
        gf_sig = loader.load_gen_funnel(sig_f)
        sc_sig = loader.load_stage_counts(sig_f)
        ct_sig = loader.load_cleaning_tracks(sig_f)
        sv_sig = loader.load_all_stage_vtx(sig_f)
        st_sig = loader.load_seed_tracks(sig_f)

    sc_bkg, sv_bkg = None, None
    if bkg_path:
        print(f"  [{stem}] Loading background data...")
        with uproot.open(bkg_path) as bkg_f:
            sc_bkg = loader.load_stage_counts(bkg_f)
            sv_bkg = loader.load_all_stage_vtx(bkg_f)

    # Create one TDirectory per stage
    stage_dirs = {d: out_file.mkdir(d) for d in config.STAGE_DIRS}

    print(f"  [{stem}] Seeding plots...")
    seeding.make_plots(stage_dirs["seeding"], gf_sig, sv_sig, sv_bkg, st_sig)

    print(f"  [{stem}] Merging plots...")
    merging.make_plots(stage_dirs["merging"], gf_sig, sv_sig, sv_bkg)

    print(f"  [{stem}] Cleaning plots...")
    cleaning.make_plots(stage_dirs["cleaning"], gf_sig, sv_sig, sv_bkg,
                        ct_sig, max_compat, min_cos_theta, use_diagonal_cut, clean_cut_slope)

    print(f"  [{stem}] Disambiguation plots...")
    disambiguation.make_plots(stage_dirs["disambiguation"], gf_sig, sv_sig, sv_bkg)

    print(f"  [{stem}] Filtering plots...")
    filtering.make_plots(stage_dirs["filtering"], gf_sig, sv_sig, sv_bkg)

    print(f"  [{stem}] Summary plots...")
    summary.make_plots(stage_dirs["summary"], gf_sig, sc_sig, sc_bkg)

    print(f"  [{stem}] Done.")


# ── Multi-file helpers ────────────────────────────────────────────────────────

def _worker(args_tuple):
    """Multiprocessing worker: process one signal file into a temp ROOT file."""
    sig_path, bkg_path, tmp_path, max_compat, min_cos_theta, use_diagonal_cut, clean_cut_slope = args_tuple
    stem = os.path.splitext(os.path.basename(sig_path))[0]
    print(f"[{stem}] Worker started (PID {os.getpid()})")
    tmp_file = ROOT.TFile(tmp_path, "RECREATE")
    _run_all_plots(tmp_file, sig_path, bkg_path, max_compat, min_cos_theta,
                   use_diagonal_cut, clean_cut_slope)
    tmp_file.Close()
    return tmp_path


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
    parser.add_argument("--max-compat",       type=float, default=1.5,
                        help="maxCompatibility cut value shown on 2D cleaning plot (default: 1.5)")
    parser.add_argument("--min-cos-theta",    type=float, default=0.5,
                        help="minCleanCosTheta / diagonal intercept shown on 2D cleaning plot (default: 0.5)")
    parser.add_argument("--use-diagonal-cut", action="store_true",
                        help="Draw diagonal cut line instead of square cuts on 2D cleaning plot")
    parser.add_argument("--clean-cut-slope",  type=float, default=0.0,
                        help="Slope of the diagonal cut (default: 0.0)")
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

    print(f"[hyddraDiagPlots] Signal files ({len(sig_files)}):")
    for f in sig_files:
        print(f"  {f}")
    print(f"[hyddraDiagPlots] Background: {args.background or 'none'}")
    print(f"[hyddraDiagPlots] Output:     {args.output}")

    if len(sig_files) == 1:
        # ── Single-file mode: plots directly in stage subdirs at root level ──
        out_file = ROOT.TFile(args.output, "RECREATE")
        _run_all_plots(out_file, sig_files[0], args.background,
                       args.max_compat, args.min_cos_theta,
                       args.use_diagonal_cut, args.clean_cut_slope)
        out_file.Close()

    else:
        # ── Multi-file mode: parallel workers → temp files → merge ────────────
        n_workers = args.jobs if args.jobs > 0 else min(len(sig_files), mp.cpu_count())
        print(f"[hyddraDiagPlots] Workers: {n_workers}")

        tmpdir     = tempfile.mkdtemp(prefix="hyddra_diag_")
        work_items = []
        for sig_path in sig_files:
            stem    = os.path.splitext(os.path.basename(sig_path))[0]
            tmp_out = os.path.join(tmpdir, f"{stem}.root")
            work_items.append((sig_path, args.background, tmp_out,
                               args.max_compat, args.min_cos_theta,
                               args.use_diagonal_cut, args.clean_cut_slope))

        print(f"\n[hyddraDiagPlots] Launching {len(work_items)} workers...")
        with mp.Pool(n_workers) as pool:
            pool.map(_worker, work_items)

        # Merge: each file gets its own top-level subdirectory
        print(f"\n[hyddraDiagPlots] Merging into {args.output}...")
        out_file = ROOT.TFile(args.output, "RECREATE")
        for sig_path, _, tmp_path, _, _ in work_items:
            stem = os.path.splitext(os.path.basename(sig_path))[0]
            print(f"  Copying {stem}...")
            tdir = out_file.mkdir(stem)
            _copy_to_dir(tmp_path, tdir)
            os.unlink(tmp_path)
        out_file.Close()
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass

    print(f"\n[hyddraDiagPlots] All plots saved to {args.output}")


if __name__ == "__main__":
    main()
