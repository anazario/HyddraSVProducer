#!/usr/bin/env python3
"""
run.py — entry point for hyddraDiagPlots.

Usage (single file, auto-detect mode):
    python3 scripts/hyddraDiagPlots/run.py --signal signal.root

Usage (explicit mode):
    python3 scripts/hyddraDiagPlots/run.py --signal signal.root --mode leptonic
    python3 scripts/hyddraDiagPlots/run.py --signal signal.root --mode hadronic
    python3 scripts/hyddraDiagPlots/run.py --signal signal.root --mode both

Usage (multiple files / glob):
    python3 scripts/hyddraDiagPlots/run.py --signal "diag_*.root" -j 4

Usage (with background file):
    python3 scripts/hyddraDiagPlots/run.py --signal signal.root --background bkg.root

Mode detection:
  auto (default) — inspects each ROOT file for leptonicConfig / hadronicConfig
                   trees and runs whichever mode(s) are present.
  leptonic       — only leptonic plots.  Warns and prompts if data is absent.
  hadronic       — only hadronic plots.  Warns and prompts if data is absent.
  both           — both modes.  Warns for any absent mode.

Output structure:
  Plots are organised under <mode>/seeding/, <mode>/merging/, etc., where
  <mode> is 'leptonic' or 'hadronic'.  In multi-file mode each file stem
  gets its own top-level subdirectory first.

Cut values for plot annotations are read automatically from the config
tree stored inside each signal ROOT file.  No manual flags are required.
"""

import argparse
import contextlib
import glob as glob_module
import multiprocessing as mp
import os
import sys
import tempfile

import numpy as np
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
from hyddraDiagPlots.src.plotter import had_signal_mask
from hyddraDiagPlots.stages import (
    summary, seeding, merging, cleaning, disambiguation, filtering,
    id_stage,
    hadronic_seeding, hadronic_merging, hadronic_cleaning,
    hadronic_disambiguation, hadronic_filtering, hadronic_summary, hadronic_id,
)

_STAGE_KEYS   = ['seed',    'merged',  'cleaned',  'disambig',       'filtered',  'id']
_STAGE_LABELS = ['Seeding', 'Merging', 'Cleaning', 'Disambiguation', 'Filtering', 'ID']


# ── Mode detection ─────────────────────────────────────────────────────────────

def _detect_modes(sig_path):
    """Return (has_leptonic, has_hadronic) for the given ROOT file."""
    with uproot.open(sig_path) as f:
        return loader.has_leptonic_data(f), loader.has_hadronic_data(f)


def _resolve_modes(requested_mode, sig_path):
    """
    Return the set of modes to actually run, after checking availability.
    Warns and prompts the user when a requested mode is absent.
    Returns a set of strings from {'leptonic', 'hadronic'}, or None to abort.
    """
    has_lep, has_had = _detect_modes(sig_path)
    available = set()
    if has_lep: available.add('leptonic')
    if has_had: available.add('hadronic')

    if requested_mode == 'auto':
        if not available:
            print(f"  Warning: neither leptonicConfig nor hadronicConfig found in {sig_path}.")
            print("  No plots will be produced for this file.")
            return set()
        return available

    wanted = {'leptonic', 'hadronic'} if requested_mode == 'both' else {requested_mode}
    missing = wanted - available

    if not missing:
        return wanted

    for mode in sorted(missing):
        print(f"\n  WARNING: requested mode '{mode}' is not available in:")
        print(f"    {sig_path}")
        print(f"  (no {mode}Config tree found)")

    can_continue = available & wanted
    if not can_continue:
        print("  No requested modes are available.  Skipping this file.")
        return set()

    # y/n prompt — default to 'y' when running non-interactively
    if sys.stdin.isatty():
        ans = input(f"\n  Continue with available mode(s) {sorted(can_continue)}? [Y/n] ").strip().lower()
        if ans and ans not in ('y', 'yes'):
            print("  Skipping this file.")
            return set()
    else:
        print(f"  (non-interactive) Continuing with available mode(s): {sorted(can_continue)}")

    return can_continue


# ── Efficiency summary helpers ─────────────────────────────────────────────────

def _compute_summary(stem, gf_sig, sc_sig, mode_prefix="", active_stage_keys=None):
    """Compute per-file signal efficiency summary from loaded data.

    Only leptonic gen vertices (electrons + muons) are counted; hadronic gen
    vertices are never reconstructed by the leptonic pipeline and would
    artificially inflate the denominator if included.

    active_stage_keys restricts which stages appear in the summary table; pass
    a subset of _STAGE_KEYS to skip stages absent from older files.
    """
    import numpy as np
    if active_stage_keys is None:
        active_stage_keys = _STAGE_KEYS
    active_stage_labels = [_STAGE_LABELS[_STAGE_KEYS.index(k)] for k in active_stage_keys]

    is_elec = ak.to_numpy(ak.flatten(gf_sig['GenFunnel_isElectron'])).astype(bool)
    is_muon = ak.to_numpy(ak.flatten(gf_sig['GenFunnel_isMuon'])).astype(bool)
    is_lep  = is_elec | is_muon

    n_gen_elec = int(np.sum(is_elec))
    n_gen_muon = int(np.sum(is_muon))
    n_gen      = n_gen_elec + n_gen_muon

    rows = []
    for key, label in zip(active_stage_keys, active_stage_labels):
        gold_flat   = ak.to_numpy(ak.flatten(gf_sig[f'GenFunnel_gold_{key}'])).astype(bool)
        silver_flat = ak.to_numpy(ak.flatten(gf_sig[f'GenFunnel_silver_{key}'])).astype(bool)

        n_gold_gen_elec = int(np.sum(gold_flat & is_elec))
        n_gold_gen_muon = int(np.sum(gold_flat & is_muon))
        n_gold_gen      = n_gold_gen_elec + n_gold_gen_muon
        n_silver_gen    = int(np.sum(silver_flat & is_lep))

        n_reco        = sc_sig.get(f'Stage_n_{key}',       0)
        n_gold_reco   = sc_sig.get(f'Stage_nGold_{key}',   0)
        n_silver_reco = sc_sig.get(f'Stage_nSilver_{key}', 0)
        n_bronze_reco = sc_sig.get(f'Stage_nBronze_{key}', 0)
        rows.append({
            'label':            label,
            'n_reco':           n_reco,
            'n_gold_reco':      n_gold_reco,
            'n_nonsig_reco':    max(0, n_reco - n_gold_reco - n_silver_reco - n_bronze_reco),
            'n_gold_gen':       n_gold_gen,
            'n_gold_gen_elec':  n_gold_gen_elec,
            'n_gold_gen_muon':  n_gold_gen_muon,
            'n_gs_gen':         n_gold_gen + n_silver_gen,
            'n_dup_gold':       max(0, n_gold_reco - n_gold_gen),
        })
    return {
        'stem':       stem,
        'mode':       mode_prefix,
        'n_gen':      n_gen,
        'n_gen_elec': n_gen_elec,
        'n_gen_muon': n_gen_muon,
        'rows':       rows,
    }


def _print_summaries(summaries):
    """Print a clean per-file efficiency table for all processed files."""
    summaries = sorted(summaries, key=lambda s: (s['stem'], s['mode']))

    any_dups = any(row['n_dup_gold'] > 0
                   for s in summaries for row in s['rows'])

    if any_dups:
        col_w = [16, 10, 11, 12, 9, 11, 9, 9, 13]
        divider = '  ' + '-' * (sum(col_w) + len(col_w) * 3 - 1)
        hdr = (f"  {'Stage':<{col_w[0]}} {'Reco vtx':>{col_w[1]}} "
               f"{'Gold reco':>{col_w[2]}} {'Non-signal':>{col_w[3]}} "
               f"{'Dup gold':>{col_w[4]}} {'Eff (gold)':>{col_w[5]}} "
               f"{'Eff (ele)':>{col_w[6]}} {'Eff (mu)':>{col_w[7]}} "
               f"{'% NS removed':>{col_w[8]}}")
    else:
        col_w = [16, 10, 11, 12, 11, 9, 9, 13]
        divider = '  ' + '-' * (sum(col_w) + len(col_w) * 3 - 1)
        hdr = (f"  {'Stage':<{col_w[0]}} {'Reco vtx':>{col_w[1]}} "
               f"{'Gold reco':>{col_w[2]}} {'Non-signal':>{col_w[3]}} "
               f"{'Eff (gold)':>{col_w[4]}} "
               f"{'Eff (ele)':>{col_w[5]}} {'Eff (mu)':>{col_w[6]}} "
               f"{'% NS removed':>{col_w[7]}}")

    print(f"\n{'=' * 79}")
    print(f"  HYDDRA diagnostic summary")
    print(f"{'=' * 79}")

    for s in summaries:
        n_gen      = s['n_gen']
        n_gen_elec = s.get('n_gen_elec', None)
        n_gen_muon = s.get('n_gen_muon', None)
        mode_tag = f" [{s['mode']}]" if s['mode'] else ""
        print(f"\n  File : {s['stem']}{mode_tag}")
        if n_gen_elec is not None and n_gen_muon is not None:
            print(f"  Gen signal vertices : {n_gen}  "
                  f"(electrons: {n_gen_elec}, muons: {n_gen_muon})")
        else:
            print(f"  Gen signal vertices : {n_gen}")
        print(divider)
        print(hdr)
        print(divider)
        prev_nonsig = None
        for row in s['rows']:
            gold_eff = f"{row['n_gold_gen'] / n_gen * 100:.1f}%" if n_gen > 0 else 'N/A'
            ele_eff  = (f"{row['n_gold_gen_elec'] / n_gen_elec * 100:.1f}%"
                        if (n_gen_elec is not None and n_gen_elec > 0) else 'N/A')
            mu_eff   = (f"{row['n_gold_gen_muon'] / n_gen_muon * 100:.1f}%"
                        if (n_gen_muon is not None and n_gen_muon > 0) else 'N/A')
            if prev_nonsig is None:
                ns_removed = '---'
            elif prev_nonsig > 0:
                ns_removed = f"{(prev_nonsig - row['n_nonsig_reco']) / prev_nonsig * 100:.1f}%"
            else:
                ns_removed = 'N/A'
            prev_nonsig = row['n_nonsig_reco']
            if any_dups:
                dup_str = str(row['n_dup_gold']) if row['n_dup_gold'] > 0 else '-'
                print(f"  {row['label']:<{col_w[0]}} "
                      f"{row['n_reco']:>{col_w[1]}} "
                      f"{row['n_gold_reco']:>{col_w[2]}} "
                      f"{row['n_nonsig_reco']:>{col_w[3]}} "
                      f"{dup_str:>{col_w[4]}} "
                      f"{gold_eff:>{col_w[5]}} "
                      f"{ele_eff:>{col_w[6]}} "
                      f"{mu_eff:>{col_w[7]}} "
                      f"{ns_removed:>{col_w[8]}}")
            else:
                print(f"  {row['label']:<{col_w[0]}} "
                      f"{row['n_reco']:>{col_w[1]}} "
                      f"{row['n_gold_reco']:>{col_w[2]}} "
                      f"{row['n_nonsig_reco']:>{col_w[3]}} "
                      f"{gold_eff:>{col_w[4]}} "
                      f"{ele_eff:>{col_w[5]}} "
                      f"{mu_eff:>{col_w[6]}} "
                      f"{ns_removed:>{col_w[7]}}")
        print(divider)
        if s.get('cutflow'):
            print()
            print(s['cutflow'])

    if any_dups:
        print(f"\n  WARNING: duplicate gold matches detected.")
        print(f"  'Eff (gold)' counts unique gen vertices, not reco SVs.")

    print()


# ── ID stage summary helpers ───────────────────────────────────────────────────

def _compute_id_summary_leptonic(gf_sig, sc_sig, sv_sig):
    """Compute leptonic ID-stage stats.  Returns None if id data is absent."""
    if sc_sig.get('Stage_n_id', 0) == 0 and 'GenFunnel_gold_id' not in gf_sig.fields:
        return None

    is_muon = ak.to_numpy(ak.flatten(gf_sig['GenFunnel_isMuon'    ])).astype(bool)
    is_elec = ak.to_numpy(ak.flatten(gf_sig['GenFunnel_isElectron'])).astype(bool)
    n_gen_muon = int(np.sum(is_muon))
    n_gen_elec = int(np.sum(is_elec))
    n_gen      = n_gen_muon + n_gen_elec

    gold_id = (ak.to_numpy(ak.flatten(gf_sig['GenFunnel_gold_id'])).astype(bool)
               if 'GenFunnel_gold_id' in gf_sig.fields
               else np.zeros(len(is_muon), dtype=bool))
    n_gold_gen_muon = int(np.sum(gold_id & is_muon))
    n_gold_gen_elec = int(np.sum(gold_id & is_elec))
    n_gold_gen      = n_gold_gen_muon + n_gold_gen_elec

    n_id_reco     = sc_sig.get('Stage_n_id',       0)
    n_gold_reco   = sc_sig.get('Stage_nGold_id',   0)
    n_silver_reco = sc_sig.get('Stage_nSilver_id', 0)
    n_bronze_reco = sc_sig.get('Stage_nBronze_id', 0)
    n_nonsig_reco = max(0, n_id_reco - n_gold_reco - n_silver_reco - n_bronze_reco)
    n_filtered    = sc_sig.get('Stage_n_filtered', 0)

    # Flavor breakdown from allStageVtx (may be absent in older files)
    flavor = None
    sidx_id = config.STAGE_IDX.get('id', 5)
    if (sv_sig is not None and len(sv_sig) > 0
            and 'StageVtx_passLooseMuonID' in sv_sig.fields):
        mask_id = ak.to_numpy(sv_sig['StageVtx_stageIdx']).astype(int) == sidx_id
        if np.any(mask_id):
            is_gold  = ak.to_numpy(sv_sig['StageVtx_isGold'             ])[mask_id].astype(bool)
            pass_mu  = ak.to_numpy(sv_sig['StageVtx_passLooseMuonID'    ])[mask_id].astype(bool)
            pass_el  = ak.to_numpy(sv_sig['StageVtx_passLooseElectronID'])[mask_id].astype(bool)
            flavor = {
                'sig_muon': int(np.sum( is_gold & pass_mu)),
                'sig_elec': int(np.sum( is_gold & pass_el)),
                'ns_muon':  int(np.sum(~is_gold & pass_mu)),
                'ns_elec':  int(np.sum(~is_gold & pass_el)),
            }

    return {
        'n_gen': n_gen, 'n_gen_muon': n_gen_muon, 'n_gen_elec': n_gen_elec,
        'n_gold_gen': n_gold_gen, 'n_gold_gen_muon': n_gold_gen_muon,
        'n_gold_gen_elec': n_gold_gen_elec,
        'n_id_reco': n_id_reco, 'n_gold_reco': n_gold_reco,
        'n_silver_reco': n_silver_reco, 'n_bronze_reco': n_bronze_reco,
        'n_nonsig_reco': n_nonsig_reco, 'n_filtered': n_filtered,
        'flavor': flavor,
    }


def _compute_id_summary_hadronic(gf_sig, sc_sig, sv_sig):
    """Compute hadronic ID-stage stats.  Returns None if id data is absent."""
    if sc_sig.get('Stage_n_id', 0) == 0 and 'GenFunnel_matchRatio_id' not in gf_sig.fields:
        return None

    is_hadronic = ak.to_numpy(ak.flatten(gf_sig['GenFunnel_isHadronic'])).astype(bool)
    n_gen = int(np.sum(is_hadronic))

    mr_id   = (ak.to_numpy(ak.flatten(gf_sig['GenFunnel_matchRatio_id']))[is_hadronic]
               if 'GenFunnel_matchRatio_id' in gf_sig.fields
               else np.full(n_gen, -1.0))
    n_loose = int(np.sum(mr_id > 0))
    n_tight = int(np.sum(mr_id >= 0.5))

    n_id_reco  = sc_sig.get('Stage_n_id',       0)
    n_filtered = sc_sig.get('Stage_n_filtered', 0)

    sidx_id = config.STAGE_IDX.get('id', 5)
    n_nonsig_reco = None
    if sv_sig is not None and len(sv_sig) > 0 and 'StageVtx_matchRatio' in sv_sig.fields:
        mask_id = ak.to_numpy(sv_sig['StageVtx_stageIdx']).astype(int) == sidx_id
        if np.any(mask_id):
            mr = ak.to_numpy(sv_sig['StageVtx_matchRatio'])[mask_id]
            n_nonsig_reco = int(np.sum(mr <= 0))

    return {
        'n_gen': n_gen, 'n_loose': n_loose, 'n_tight': n_tight,
        'n_id_reco': n_id_reco, 'n_filtered': n_filtered,
        'n_nonsig_reco': n_nonsig_reco,
    }


def _print_id_summaries(summaries):
    """Print per-file ID stage summary for all processed files."""
    id_summs = [(s['stem'], s['mode'], s.get('id_summary'))
                for s in sorted(summaries, key=lambda s: (s['stem'], s['mode']))
                if s.get('id_summary') is not None]
    if not id_summs:
        return

    print(f"\n{'=' * 79}")
    print(f"  HYDDRA ID stage summary")
    print(f"{'=' * 79}")

    for stem, mode, d in id_summs:
        divider = '  ' + '-' * 75
        print(f"\n  File : {stem} [{mode}]")

        if mode == 'leptonic':
            n_gen = d['n_gen']
            print(f"  Gen signal      : {n_gen}"
                  f"  (muon: {d['n_gen_muon']}, electron: {d['n_gen_elec']})")
            print(divider)

            ns_pct = (f"{d['n_nonsig_reco'] / d['n_filtered'] * 100:.1f}%"
                      if d['n_filtered'] > 0 else 'N/A')
            print(f"  {'Reco vertices':<22}  Pre-ID: {d['n_filtered']:<8} "
                  f"Post-ID: {d['n_id_reco']:<8}"
                  f"(non-signal removed: {ns_pct})")
            print(f"  {'Gold reco at ID':<22}: {d['n_gold_reco']}"
                  f"  |  Silver: {d['n_silver_reco']}"
                  f"  |  Bronze: {d['n_bronze_reco']}"
                  f"  |  Non-signal: {d['n_nonsig_reco']}")

            eff     = (f"{d['n_gold_gen']      / n_gen                * 100:.1f}%"
                       if n_gen > 0 else 'N/A')
            eff_mu  = (f"{d['n_gold_gen_muon'] / d['n_gen_muon']      * 100:.1f}%"
                       if d['n_gen_muon'] > 0 else 'N/A')
            eff_el  = (f"{d['n_gold_gen_elec'] / d['n_gen_elec']      * 100:.1f}%"
                       if d['n_gen_elec'] > 0 else 'N/A')
            print(f"  {'Gen ID efficiency':<22}: {eff}"
                  f"  (muon: {eff_mu},  electron: {eff_el})")

            if d['flavor']:
                f = d['flavor']
                print(f"  {'Flavor breakdown':<22}  "
                      f"Muon ID:     signal (gold) {f['sig_muon']:<6}  "
                      f"non-signal {f['ns_muon']}")
                print(f"  {'':<22}  "
                      f"Electron ID: signal (gold) {f['sig_elec']:<6}  "
                      f"non-signal {f['ns_elec']}")

        else:  # hadronic
            n_gen = d['n_gen']
            print(f"  Gen signal      : {n_gen}")
            print(divider)

            ns_str = (str(d['n_nonsig_reco']) if d['n_nonsig_reco'] is not None else 'N/A')
            ns_pct = ('N/A' if d['n_nonsig_reco'] is None or d['n_filtered'] == 0
                      else f"{d['n_nonsig_reco'] / d['n_filtered'] * 100:.1f}%")
            print(f"  {'Reco vertices':<22}  Pre-ID: {d['n_filtered']:<8} "
                  f"Post-ID: {d['n_id_reco']:<8}"
                  f"(non-signal: {ns_str},  removed: {ns_pct})")

            eff_loose = (f"{d['n_loose'] / n_gen * 100:.1f}%" if n_gen > 0 else 'N/A')
            eff_tight = (f"{d['n_tight'] / n_gen * 100:.1f}%" if n_gen > 0 else 'N/A')
            print(f"  {'Gen efficiency':<22}  "
                  f"Loose (matchRatio > 0):    {d['n_loose']:<6} / {n_gen}  ({eff_loose})")
            print(f"  {'':<22}  "
                  f"Tight (matchRatio >= 0.5): {d['n_tight']:<6} / {n_gen}  ({eff_tight})")

        print(divider)

    print()


# ── Mode-specific plot runners ─────────────────────────────────────────────────

def _run_leptonic_plots(mode_dir, sig_path, bkg_path):
    """Run the full leptonic diagnostic plot suite into mode_dir."""
    stem = os.path.splitext(os.path.basename(sig_path))[0]

    with uproot.open(sig_path) as sig_f:
        has_id = loader.has_id_stage(sig_f, base=loader._BASE)
        gf_sig = loader.load_gen_funnel(sig_f)
        sc_sig = loader.load_stage_counts(sig_f)
        ct_sig = loader.load_cleaning_tracks(sig_f)
        sv_sig = loader.load_all_stage_vtx(sig_f)
        st_sig = loader.load_seed_tracks(sig_f)
        cfg    = loader.load_leptonic_config(sig_f)

    if not has_id:
        print(f"  [{stem}/leptonic] Note: ID stage not found in file "
              "(produced before ID stage was added); skipping ID plots.",
              file=sys.stderr)

    if cfg is None:
        print(f"  [{stem}/leptonic] Warning: leptonicConfig not found; "
              "cut annotations will use defaults", file=sys.stderr)

    sc_bkg, sv_bkg = None, None
    if bkg_path:
        with uproot.open(bkg_path) as bkg_f:
            sc_bkg = loader.load_stage_counts(bkg_f)
            sv_bkg = loader.load_all_stage_vtx(bkg_f)

    stage_dirs = {d: mode_dir.mkdir(d) for d in config.STAGE_DIRS}

    stages = [
        ("seeding",        lambda: seeding.make_plots(       stage_dirs["seeding"],        gf_sig, sv_sig, sv_bkg, st_sig)),
        ("merging",        lambda: merging.make_plots(        stage_dirs["merging"],        gf_sig, sv_sig, sv_bkg)),
        ("cleaning",       lambda: cleaning.make_plots(       stage_dirs["cleaning"],       gf_sig, sv_sig, sv_bkg, ct_sig, cfg)),
        ("disambiguation", lambda: disambiguation.make_plots( stage_dirs["disambiguation"], gf_sig, sv_sig, sv_bkg)),
        ("filtering",      lambda: filtering.make_plots(      stage_dirs["filtering"],      gf_sig, sv_sig, sv_bkg, cfg)),
        ("summary",        lambda: summary.make_plots(        stage_dirs["summary"],        gf_sig, sc_sig, sc_bkg)),
    ]
    if has_id:
        stages.insert(-1, ("id", lambda: id_stage.make_plots(stage_dirs["id"], gf_sig, sv_sig, sv_bkg)))

    active_stage_keys = _STAGE_KEYS if has_id else [k for k in _STAGE_KEYS if k != 'id']

    cutflow_str = None
    with tqdm(stages, desc=f"  {stem} [leptonic]", unit="stage", leave=True) as pbar:
        for stage_name, fn in pbar:
            pbar.set_postfix_str(stage_name)
            with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
                result = fn()
            if result is not None:
                cutflow_str = result
        pbar.set_postfix_str("done")

    summ = _compute_summary(stem, gf_sig, sc_sig, mode_prefix="leptonic",
                            active_stage_keys=active_stage_keys)
    summ['cutflow']    = cutflow_str
    summ['id_summary'] = _compute_id_summary_leptonic(gf_sig, sc_sig, sv_sig) if has_id else None
    return summ


def _compute_summary_hadronic(stem, gf_sig, sc_sig, sv_sig=None, active_stage_keys=None):
    """Compute per-file hadronic signal efficiency summary using matchRatio."""
    if active_stage_keys is None:
        active_stage_keys = _STAGE_KEYS
    active_stage_labels = [_STAGE_LABELS[_STAGE_KEYS.index(k)] for k in active_stage_keys]

    is_hadronic = ak.to_numpy(ak.flatten(gf_sig["GenFunnel_isHadronic"])).astype(bool)
    n_gen       = int(np.sum(is_hadronic))

    # Per-stage non-signal reco count from allStageVtx
    nonsig_by_stage = {}
    if sv_sig is not None and len(sv_sig) > 0:
        sv_stage  = ak.to_numpy(sv_sig["StageVtx_stageIdx"])
        sv_is_sig = had_signal_mask(sv_sig)
        for si, key in enumerate(active_stage_keys):
            mask = sv_stage == si
            nonsig_by_stage[key] = int(np.sum(mask & ~sv_is_sig))

    rows = []
    for key, label in zip(active_stage_keys, active_stage_labels):
        mr       = ak.to_numpy(ak.flatten(gf_sig[f"GenFunnel_matchRatio_{key}"]))[is_hadronic]
        n_loose  = int(np.sum(mr > 0))
        n_tight  = int(np.sum(mr >= 0.5))
        n_reco   = sc_sig.get(f"Stage_n_{key}", 0)
        rows.append({
            'label':        label,
            'n_reco':       n_reco,
            'n_loose':      n_loose,
            'n_tight':      n_tight,
            'n_nonsig_reco': nonsig_by_stage.get(key, None),
        })
    return {'stem': stem, 'mode': 'hadronic', 'n_gen': n_gen, 'rows': rows}


def _print_summaries_hadronic(hadronic_summaries):
    """Print hadronic efficiency table."""
    if not hadronic_summaries:
        return
    col_w = [16, 10, 11, 11, 12, 12, 13]
    divider = '  ' + '-' * (sum(col_w) + len(col_w) * 3 - 1)
    hdr = (f"  {'Stage':<{col_w[0]}} {'Reco vtx':>{col_w[1]}} "
           f"{'Loose gen':>{col_w[2]}} {'Tight gen':>{col_w[3]}} "
           f"{'Eff (loose)':>{col_w[4]}} {'Eff (tight)':>{col_w[5]}} "
           f"{'% NS removed':>{col_w[6]}}")
    for s in hadronic_summaries:
        n_gen = s['n_gen']
        print(f"\n  File : {s['stem']} [hadronic]")
        print(f"  Gen signal vertices : {n_gen}")
        print(divider)
        print(hdr)
        print(divider)
        prev_nonsig = None
        for row in s['rows']:
            loose_eff = f"{row['n_loose'] / n_gen * 100:.1f}%" if n_gen > 0 else 'N/A'
            tight_eff = f"{row['n_tight'] / n_gen * 100:.1f}%" if n_gen > 0 else 'N/A'
            n_ns = row.get('n_nonsig_reco', None)
            if n_ns is None or prev_nonsig is None:
                ns_removed = '---'
            elif prev_nonsig > 0:
                ns_removed = f"{(prev_nonsig - n_ns) / prev_nonsig * 100:.1f}%"
            else:
                ns_removed = 'N/A'
            prev_nonsig = n_ns
            print(f"  {row['label']:<{col_w[0]}} "
                  f"{row['n_reco']:>{col_w[1]}} "
                  f"{row['n_loose']:>{col_w[2]}} "
                  f"{row['n_tight']:>{col_w[3]}} "
                  f"{loose_eff:>{col_w[4]}} "
                  f"{tight_eff:>{col_w[5]}} "
                  f"{ns_removed:>{col_w[6]}}")
        print(divider)
        if s.get('cutflow'):
            print()
            print(s['cutflow'])


def _run_hadronic_plots(mode_dir, sig_path, bkg_path):
    """Run the full hadronic diagnostic plot suite into mode_dir."""
    stem = os.path.splitext(os.path.basename(sig_path))[0]

    with uproot.open(sig_path) as sig_f:
        has_id = loader.has_id_stage(sig_f, base=loader._HADRONIC_BASE)
        gf_sig = loader.load_gen_funnel(sig_f,               base=loader._HADRONIC_BASE)
        sc_sig = loader.load_stage_counts_hadronic(sig_f)
        sv_sig = loader.load_all_stage_vtx(sig_f,            base=loader._HADRONIC_BASE)
        cfg    = loader.load_hadronic_config(sig_f)

    if not has_id:
        print(f"  [{stem}/hadronic] Note: ID stage not found in file "
              "(produced before ID stage was added); skipping ID plots.",
              file=sys.stderr)

    if cfg is None:
        print(f"  [{stem}/hadronic] Warning: hadronicConfig not found; "
              "cut annotations will use defaults", file=sys.stderr)

    sc_bkg, sv_bkg = None, None
    if bkg_path:
        with uproot.open(bkg_path) as bkg_f:
            sc_bkg = loader.load_stage_counts_hadronic(bkg_f)
            sv_bkg = loader.load_all_stage_vtx(bkg_f, base=loader._HADRONIC_BASE)

    stage_dirs = {d: mode_dir.mkdir(d) for d in config.STAGE_DIRS}

    stages = [
        ("seeding",        lambda: hadronic_seeding.make_plots(       stage_dirs["seeding"],        gf_sig, sv_sig, sv_bkg)),
        ("merging",        lambda: hadronic_merging.make_plots(        stage_dirs["merging"],        gf_sig, sv_sig, sv_bkg)),
        ("cleaning",       lambda: hadronic_cleaning.make_plots(       stage_dirs["cleaning"],       gf_sig, sv_sig, sv_bkg, cfg)),
        ("disambiguation", lambda: hadronic_disambiguation.make_plots( stage_dirs["disambiguation"], gf_sig, sv_sig, sv_bkg)),
        ("filtering",      lambda: hadronic_filtering.make_plots(      stage_dirs["filtering"],      gf_sig, sv_sig, sv_bkg, cfg)),
        ("summary",        lambda: hadronic_summary.make_plots(        stage_dirs["summary"],        gf_sig, sc_sig, sc_bkg)),
    ]
    if has_id:
        stages.insert(-1, ("id", lambda: hadronic_id.make_plots(stage_dirs["id"], gf_sig, sv_sig, sv_bkg)))

    active_stage_keys = _STAGE_KEYS if has_id else [k for k in _STAGE_KEYS if k != 'id']

    cutflow_str = None
    with tqdm(stages, desc=f"  {stem} [hadronic]", unit="stage", leave=True) as pbar:
        for stage_name, fn in pbar:
            pbar.set_postfix_str(stage_name)
            with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
                result = fn()
            if result is not None:
                cutflow_str = result
        pbar.set_postfix_str("done")

    summ = _compute_summary_hadronic(stem, gf_sig, sc_sig, sv_sig,
                                     active_stage_keys=active_stage_keys)
    summ['cutflow']    = cutflow_str
    summ['id_summary'] = _compute_id_summary_hadronic(gf_sig, sc_sig, sv_sig) if has_id else None
    return summ


# ── Core plot runner ───────────────────────────────────────────────────────────

def _run_all_plots(out_file, sig_path, bkg_path, requested_mode):
    """
    Determine which modes to run, create mode subdirectories inside out_file,
    and dispatch to the appropriate plot runners.  Returns a list of summaries.
    """
    active_modes = _resolve_modes(requested_mode, sig_path)
    if not active_modes:
        return []

    summaries = []
    for mode in sorted(active_modes):
        mode_dir = out_file.mkdir(mode)
        if mode == 'leptonic':
            summ = _run_leptonic_plots(mode_dir, sig_path, bkg_path)
        else:
            summ = _run_hadronic_plots(mode_dir, sig_path, bkg_path)
        summaries.append(summ)

    return summaries


# ── Multi-file helpers ─────────────────────────────────────────────────────────

def _worker(args_tuple):
    """Multiprocessing worker: process one signal file into a temp ROOT file."""
    sig_path, bkg_path, tmp_path, requested_mode = args_tuple
    with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
        tmp_file = ROOT.TFile(tmp_path, "RECREATE")
        summs = _run_all_plots(tmp_file, sig_path, bkg_path, requested_mode)
        tmp_file.Close()
    return tmp_path, summs


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


# ── Main ──────────────────────────────────────────────────────────────────────

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
    parser.add_argument("--mode", default="auto",
                        choices=["auto", "leptonic", "hadronic", "both"],
                        help="Which HYDDRA variant to plot.  'auto' detects from file contents "
                             "(default).  'leptonic'/'hadronic' run only that variant — a "
                             "warning and y/n prompt is issued if the data is absent.  "
                             "'both' runs both variants.")
    parser.add_argument("--jobs", "-j", type=int, default=0,
                        help="Max parallel workers for multi-file mode (0 = cpu_count)")
    args = parser.parse_args()

    # Expand glob patterns and deduplicate
    sig_files = []
    for pat in args.signal:
        expanded = sorted(glob_module.glob(pat))
        if not expanded:
            print(f"Error: no files matched pattern: {pat}", file=sys.stderr)
            sys.exit(1)
        sig_files.extend(expanded)
    sig_files = list(dict.fromkeys(sig_files))

    if not sig_files:
        print("Error: no signal files found.", file=sys.stderr)
        sys.exit(1)

    print(f"[hyddraDiagPlots] {len(sig_files)} signal file(s)  |  "
          f"mode: {args.mode}  |  "
          f"background: {args.background or 'none'}  |  output: {args.output}")

    summaries = []

    if len(sig_files) == 1:
        # ── Single-file mode ──────────────────────────────────────────────────
        out_file = ROOT.TFile(args.output, "RECREATE")
        summs = _run_all_plots(out_file, sig_files[0], args.background, args.mode)
        out_file.Close()
        summaries.extend(summs)

    else:
        # ── Multi-file mode: parallel workers → temp files → merge ────────────
        n_workers = args.jobs if args.jobs > 0 else min(len(sig_files), mp.cpu_count())

        tmpdir     = tempfile.mkdtemp(prefix="hyddra_diag_")
        work_items = []
        for sig_path in sig_files:
            stem    = os.path.splitext(os.path.basename(sig_path))[0]
            tmp_out = os.path.join(tmpdir, f"{stem}.root")
            work_items.append((sig_path, args.background, tmp_out, args.mode))

        with mp.Pool(n_workers) as pool:
            with tqdm(total=len(work_items), desc="Processing", unit="file") as pbar:
                for _, summs in pool.imap_unordered(_worker, work_items):
                    summaries.extend(summs)
                    pbar.update()

        # Merge: each file gets its own top-level subdirectory
        out_file = ROOT.TFile(args.output, "RECREATE")
        for sig_path, _, tmp_path, _ in work_items:
            stem = os.path.splitext(os.path.basename(sig_path))[0]
            tdir = out_file.mkdir(stem)
            _copy_to_dir(tmp_path, tdir)
            os.unlink(tmp_path)
        out_file.Close()
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass

    lep_summs = [s for s in summaries if s['mode'] == 'leptonic']
    had_summs = [s for s in summaries if s['mode'] == 'hadronic']
    if lep_summs:
        _print_summaries(lep_summs)
    if had_summs:
        _print_summaries_hadronic(had_summs)
    _print_id_summaries(summaries)
    print(f"[hyddraDiagPlots] Plots saved to {args.output}")


if __name__ == "__main__":
    main()
