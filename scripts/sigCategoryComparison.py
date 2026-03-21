#!/usr/bin/env python3
"""
sigCategoryComparison.py — kinematic comparison across SMS-GlGl signal categories.

Classifies diagnostic ROOT files automatically by parsing mN2 and mN1 masses from
their filenames:
  ΔM = mN2 - mN1 ≤  50 GeV  → compressed
  ΔM             = 100 GeV   → threshold
  ΔM             > 100 GeV   → uncompressed

Usage:
    python3 scripts/sigCategoryComparison.py --files "path/diag_*.root" \\
            --output comparison.root [--workers 4] [--mode auto]

Memory design: each file is processed by a worker that opens it, extracts small
per-file statistics (efficiency arrays, normalized histograms), discards the full
awkward arrays, and returns only numpy arrays.  At most --workers files are held
in memory simultaneously; the main process only ever holds the collected small
stats dicts.

Output ROOT structure:
  leptonic/
    combined/  muon/  electron/
      compressed/    eff_vs_dxy_{filtered,id}  eff_vs_pt_{filtered,id}
      threshold/     (same)
      uncompressed/  (same)
      categories/    eff_vs_dxy_categories  eff_vs_pt_categories
                     stage_summary
                     reco_distributions/reco_{obs}_categories
  hadronic/
    tight/  loose/   (same layout, no flavor split)
"""

from __future__ import annotations

import argparse
import glob as glob_module
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

# ROOT and uproot are imported lazily below so that spawned worker processes
# do not trigger ROOT's global initialisation before they need it.

# ── Path setup (needed in both main and workers) ───────────────────────────────

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# ── Constants (no ROOT dependency) ────────────────────────────────────────────

CATEGORIES = ["compressed", "threshold", "uncompressed"]

# These mirror hyddraDiagPlots.src.config — duplicated here so the module can
# be imported without triggering ROOT in worker processes.
_STAGE_KEYS  = ["seed", "merged", "cleaned", "disambig", "filtered", "id"]
_STAGE_NAMES = ["Seeding", "Merging", "Cleaning", "Disambiguation", "Filtering", "ID"]

_GEN_DXY_BINS = list(np.concatenate([
    np.linspace(0,  5, 11),
    np.linspace(6, 20,  8),
    np.linspace(25, 50, 6),
]))
_GEN_PT_BINS = list(np.linspace(0, 100, 21))

_HAD_MIN3D_CUT = 0.05  # cm

_RECO_OBS = {
    "cosTheta":   {"label": "cos#theta (wrt PV)",       "bins": list(np.linspace(-1,   1, 51)), "log_y": True},
    "decayAngle": {"label": "cos#theta* (decay angle)", "bins": list(np.linspace(-1,   1, 51)), "log_y": True},
    "pOverE":     {"label": "p/E",                      "bins": list(np.linspace( 0,   1, 51)), "log_y": True},
    "dxySignif":  {"label": "dxy Significance",         "bins": list(np.linspace( 0, 150, 76)), "log_y": True},
    "mass":       {"label": "Invariant mass (GeV)",     "bins": list(np.linspace( 0, 100, 51)), "log_y": True},
}
_HAD_RECO_OBS = {
    **_RECO_OBS,
    "nTracks": {"label": "Number of tracks", "bins": list(range(2, 33)), "log_y": True},
}

# Leptonic flavor configs: (dir_key, genFunnel_field, legend_label)
_FLAVOR_CONFIGS = [
    ("combined", None,         "Combined"),
    ("muon",     "isMuon",     "Muon"),
    ("electron", "isElectron", "Electron"),
]

# Hadronic tiers: (dir_key, loose_bool, legend_label)
_HAD_TIERS = [
    ("tight", False, "Tight (matchRatio #geq 0.5)"),
    ("loose", True,  "Loose (matchRatio > 0)"),
]

# ── Filename parsing ───────────────────────────────────────────────────────────

_SMS_RE = re.compile(
    r'mGl-(?P<mGl>\d+).*?mN2-(?P<mN2>\d+).*?mN1-(?P<mN1>\d+)'
    r'.*?N2ctau-(?P<ctau>[\dp]+)',
    re.IGNORECASE,
)


def _ctau_float(s: str) -> float:
    return float(s.replace('p', '.'))


def parse_sms_params(path: str) -> dict | None:
    """Return parameter dict for an SMS-GlGl file, or None if no match."""
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
    return dict(mGl=mGl, mN2=mN2, mN1=mN1, ctau=ctau, dm=dm, category=cat,
                label=f"mGl={mGl}, #DeltaM={dm}, c#tau={ctau}")


def classify_files(paths: list[str]) -> dict[str, list]:
    """Return {category: [(path, params), ...]} sorted by (dm, mGl, ctau)."""
    groups = {cat: [] for cat in CATEGORIES}
    unclassified = []
    for p in paths:
        params = parse_sms_params(p)
        if params is None:
            unclassified.append(p)
        else:
            groups[params['category']].append((p, params))
    for cat in CATEGORIES:
        groups[cat].sort(key=lambda x: (x[1]['dm'], x[1]['mGl'], x[1]['ctau']))
    print("\n[classify]")
    for cat, files in groups.items():
        print(f"  {cat:15s}: {len(files):2d} file(s)")
    if unclassified:
        print(f"  unclassified   : {len(unclassified)} file(s)")
        for p in unclassified:
            print(f"    {os.path.basename(p)}")
    return groups


# ── Per-file stat extraction (worker-safe, no ROOT) ────────────────────────────

def _lep_eff_bins(gf, var_key: str, bins_a: np.ndarray, stage_key: str,
                  flavor_field: str | None) -> np.ndarray | None:
    """Leptonic gold efficiency per bin (NaN where denom=0), or None if stage absent."""
    import awkward as ak
    gold_key = f"GenFunnel_gold_{stage_key}"
    if gold_key not in gf.fields:
        return None
    has_t = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    var   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{var_key}"]))[has_t]
    gold  = ak.to_numpy(ak.flatten(gf[gold_key]))[has_t].astype(bool)
    if flavor_field is not None:
        fm   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{flavor_field}"]))[has_t].astype(bool)
        var, gold = var[fm], gold[fm]
    denom = np.histogram(var,       bins=bins_a)[0].astype(float)
    numer = np.histogram(var[gold], bins=bins_a)[0].astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(denom > 0, numer / denom, np.nan)


def _lep_stage_eff(gf, stage_key: str, flavor_field: str | None) -> float:
    """Single leptonic gold efficiency scalar at one stage."""
    import awkward as ak
    gold_key = f"GenFunnel_gold_{stage_key}"
    if gold_key not in gf.fields:
        return float('nan')
    has_t = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    gold  = ak.to_numpy(ak.flatten(gf[gold_key]))[has_t].astype(bool)
    if flavor_field is not None:
        fm      = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{flavor_field}"]))[has_t].astype(bool)
        n_denom = int(np.sum(fm))
        return float(np.sum(gold[fm])) / n_denom if n_denom > 0 else float('nan')
    n_denom = int(np.sum(has_t))
    return float(np.sum(gold)) / n_denom if n_denom > 0 else float('nan')


def _had_eff_bins(gf, var_key: str, bins_a: np.ndarray, stage_key: str,
                  loose: bool) -> np.ndarray | None:
    """Hadronic signal efficiency per bin (NaN where denom=0), or None if stage absent."""
    import awkward as ak
    mr_key = f"GenFunnel_matchRatio_{stage_key}"
    if mr_key not in gf.fields:
        return None
    is_had = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    var    = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{var_key}"]))[is_had]
    mr     = ak.to_numpy(ak.flatten(gf[mr_key]))[is_had].astype(float)
    sig    = (mr > 0) if loose else (mr >= 0.5)
    min3d_key = f"GenFunnel_min3D_{stage_key}"
    if min3d_key in gf.fields:
        m3d = ak.to_numpy(ak.flatten(gf[min3d_key]))[is_had].astype(float)
        sig = sig & (m3d >= 0) & (m3d < _HAD_MIN3D_CUT)
    denom = np.histogram(var,      bins=bins_a)[0].astype(float)
    numer = np.histogram(var[sig], bins=bins_a)[0].astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(denom > 0, numer / denom, np.nan)


def _had_stage_eff(gf, stage_key: str, loose: bool) -> float:
    """Single hadronic signal efficiency scalar at one stage."""
    import awkward as ak
    mr_key = f"GenFunnel_matchRatio_{stage_key}"
    if mr_key not in gf.fields:
        return float('nan')
    is_had  = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    n_denom = int(np.sum(is_had))
    if n_denom == 0:
        return float('nan')
    mr  = ak.to_numpy(ak.flatten(gf[mr_key]))[is_had].astype(float)
    sig = (mr > 0) if loose else (mr >= 0.5)
    min3d_key = f"GenFunnel_min3D_{stage_key}"
    if min3d_key in gf.fields:
        m3d = ak.to_numpy(ak.flatten(gf[min3d_key]))[is_had].astype(float)
        sig = sig & (m3d >= 0) & (m3d < _HAD_MIN3D_CUT)
    return float(np.sum(sig)) / n_denom


def _norm_hist(vals: np.ndarray, bins_a: np.ndarray) -> np.ndarray | None:
    cnts, _ = np.histogram(vals, bins=bins_a)
    tot = float(cnts.sum())
    return cnts.astype(float) / tot if tot > 0 else None


# ── Worker function ────────────────────────────────────────────────────────────

def _compute_file_stats(task: tuple) -> dict:
    """
    Worker: open one diagnostic ROOT file, compute all summary statistics,
    discard the full arrays, return a small dict of numpy arrays.

    This function is designed to be safe for multiprocessing (no ROOT drawing).
    """
    path, params, run_lep, run_had = task
    import uproot
    import awkward as ak

    stats: dict = {'params': params, 'has_lep': False, 'has_had': False}

    bins_dxy = np.array(_GEN_DXY_BINS, dtype=float)
    bins_pt  = np.array(_GEN_PT_BINS,  dtype=float)

    try:
        with uproot.open(path) as f:

            # ── leptonic ──────────────────────────────────────────────────────
            lep_base = "hyddraSVsDiagAnalyzer"
            if run_lep:
                try:
                    f[f"{lep_base}/leptonicConfig"]
                    has_lep = True
                except KeyError:
                    has_lep = False

                if has_lep:
                    gf = f[f"{lep_base}/genFunnel"].arrays(library="ak")
                    sv = f[f"{lep_base}/allStageVtx"].arrays(library="ak")

                    lep_eff_dxy:  dict = {}
                    lep_eff_pt:   dict = {}
                    lep_stage_eff: dict = {}
                    lep_reco:     dict = {}

                    for ff in [None, 'isMuon', 'isElectron']:
                        for sk in _STAGE_KEYS:
                            lep_eff_dxy[(sk, ff)]  = _lep_eff_bins(gf, 'dxy', bins_dxy, sk, ff)
                            lep_eff_pt[(sk, ff)]   = _lep_eff_bins(gf, 'pt',  bins_pt,  sk, ff)
                            lep_stage_eff[(sk, ff)] = _lep_stage_eff(gf, sk, ff)

                    sidx_map = {sk: i for i, sk in enumerate(_STAGE_KEYS)}
                    stg_idx  = ak.to_numpy(sv["StageVtx_stageIdx"]).astype(int)
                    sig_mask = ak.to_numpy(sv["StageVtx_isGold"]).astype(bool)

                    for sk, sidx in sidx_map.items():
                        sel = (stg_idx == sidx) & sig_mask
                        for obs_key, obs_cfg in _RECO_OBS.items():
                            branch = f"StageVtx_{obs_key}"
                            if branch not in sv.fields or not np.any(sel):
                                lep_reco[(sk, obs_key)] = None
                            else:
                                vals = ak.to_numpy(sv[branch]).astype(float)[sel]
                                lep_reco[(sk, obs_key)] = _norm_hist(
                                    vals, np.array(obs_cfg['bins'], dtype=float))

                    stats.update(has_lep=True,
                                 lep_eff_dxy=lep_eff_dxy, lep_eff_pt=lep_eff_pt,
                                 lep_stage_eff=lep_stage_eff, lep_reco=lep_reco)
                    del gf, sv  # free large arrays immediately

            # ── hadronic ──────────────────────────────────────────────────────
            had_base = "hyddraSVsHadronicDiagAnalyzer"
            if run_had:
                try:
                    f[f"{had_base}/hadronicConfig"]
                    has_had = True
                except KeyError:
                    has_had = False

                if has_had:
                    gf = f[f"{had_base}/genFunnel"].arrays(library="ak")
                    sv = f[f"{had_base}/allStageVtx"].arrays(library="ak")

                    had_eff_dxy:  dict = {}
                    had_eff_pt:   dict = {}
                    had_stage_eff: dict = {}
                    had_reco:     dict = {}

                    for loose in [False, True]:
                        for sk in _STAGE_KEYS:
                            had_eff_dxy[(sk, loose)]   = _had_eff_bins(gf, 'dxy', bins_dxy, sk, loose)
                            had_eff_pt[(sk, loose)]    = _had_eff_bins(gf, 'pt',  bins_pt,  sk, loose)
                            had_stage_eff[(sk, loose)] = _had_stage_eff(gf, sk, loose)

                    sidx_map  = {sk: i for i, sk in enumerate(_STAGE_KEYS)}
                    stg_idx   = ak.to_numpy(sv["StageVtx_stageIdx"]).astype(int)
                    if "StageVtx_matchRatio" in sv.fields:
                        mr      = ak.to_numpy(sv["StageVtx_matchRatio"]).astype(float)
                        sig_had = mr > 0
                        if "StageVtx_min3D" in sv.fields:
                            m3d     = ak.to_numpy(sv["StageVtx_min3D"]).astype(float)
                            sig_had = sig_had & (m3d >= 0) & (m3d < _HAD_MIN3D_CUT)
                    else:
                        sig_had = np.zeros(len(stg_idx), dtype=bool)

                    for sk, sidx in sidx_map.items():
                        sel = (stg_idx == sidx) & sig_had
                        for obs_key, obs_cfg in _HAD_RECO_OBS.items():
                            branch = f"StageVtx_{obs_key}"
                            if branch not in sv.fields or not np.any(sel):
                                had_reco[(sk, obs_key)] = None
                            else:
                                vals = ak.to_numpy(sv[branch]).astype(float)[sel]
                                had_reco[(sk, obs_key)] = _norm_hist(
                                    vals, np.array(obs_cfg['bins'], dtype=float))

                    stats.update(has_had=True,
                                 had_eff_dxy=had_eff_dxy, had_eff_pt=had_eff_pt,
                                 had_stage_eff=had_stage_eff, had_reco=had_reco)
                    del gf, sv

    except Exception as e:
        print(f"  [worker] {os.path.basename(path)}: {e}", file=sys.stderr)

    return stats


# ── Parallel stats computation ─────────────────────────────────────────────────

def compute_all_stats(groups: dict, run_lep: bool, run_had: bool,
                      n_workers: int) -> dict[str, list]:
    """
    Dispatch one worker per file using ProcessPoolExecutor.
    Returns {category: [(params, stats_dict), ...]} in the original sort order.
    """
    # Build flat task list preserving (category, index) for reassembly
    tasks, task_meta = [], []
    for cat in CATEGORIES:
        for idx, (path, params) in enumerate(groups[cat]):
            tasks.append((path, params, run_lep, run_had))
            task_meta.append((cat, idx))

    # Pre-allocate result slots so we can restore original order
    result_slots: dict[str, list] = {cat: [None] * len(groups[cat]) for cat in CATEGORIES}

    n_total = len(tasks)
    n_done  = 0

    print(f"\n  Processing {n_total} files with {n_workers} worker(s)...")
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        future_to_meta = {
            pool.submit(_compute_file_stats, task): meta
            for task, meta in zip(tasks, task_meta)
        }
        for future in as_completed(future_to_meta):
            cat, idx = future_to_meta[future]
            try:
                stats = future.result()
            except Exception as e:
                print(f"  [worker] future failed: {e}", file=sys.stderr)
                stats = {'params': groups[cat][idx][1],
                         'has_lep': False, 'has_had': False}
            result_slots[cat][idx] = (stats['params'], stats)
            n_done += 1
            print(f"  [{n_done}/{n_total}] {os.path.basename(groups[cat][idx][0])} done",
                  flush=True)

    # Filter out any None slots (shouldn't happen, but defensive)
    return {cat: [x for x in slots if x is not None]
            for cat, slots in result_slots.items()}


# ── Band computation ───────────────────────────────────────────────────────────

def build_eff_band(cat_stats: list, accessor, bins: list):
    """
    Aggregate per-file efficiency into median + IQR band.

    cat_stats : [(params, stats_dict), ...]
    accessor  : callable(stats_dict) → 1D numpy array or None
    bins      : list of bin edges

    Returns (x_centers, per_file_list, params_list, median, q25, q75) or None.
    """
    bins_a  = np.array(bins, dtype=float)
    x_cents = 0.5 * (bins_a[:-1] + bins_a[1:])
    per_file, params = [], []
    for p, s in cat_stats:
        arr = accessor(s)
        if arr is not None:
            per_file.append(arr)
            params.append(p)
    if not per_file:
        return None
    mat    = np.array(per_file)
    median = np.nanmedian(mat,         axis=0)
    q25    = np.nanpercentile(mat, 25, axis=0)
    q75    = np.nanpercentile(mat, 75, axis=0)
    return x_cents, per_file, params, median, q25, q75


# ── ROOT drawing helpers (main process only) ───────────────────────────────────

def _init_root():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    return ROOT


def _get_cat_colors(ROOT):
    return {
        "compressed":   ROOT.kRed   + 1,
        "threshold":    ROOT.kGreen + 2,
        "uncompressed": ROOT.kBlue  + 2,
    }


def _get_cat_palettes(ROOT):
    return {
        "compressed":   [ROOT.kRed-7,   ROOT.kRed-4,   ROOT.kRed-1,
                         ROOT.kRed+1,   ROOT.kRed+2,   ROOT.kRed+3],
        "threshold":    [ROOT.kGreen-9, ROOT.kGreen-6, ROOT.kGreen,
                         ROOT.kGreen+1, ROOT.kGreen+2, ROOT.kGreen+3],
        "uncompressed": [ROOT.kAzure-9, ROOT.kAzure-7, ROOT.kAzure+1,
                         ROOT.kAzure+3, ROOT.kBlue+1,  ROOT.kBlue+3],
    }


_CAT_LSTYLE = {"compressed": 1, "threshold": 2, "uncompressed": 3}


def _make_canvas(ROOT, name, w=800, h=600):
    c = ROOT.TCanvas(name, name, w, h)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)
    return c


def _axis_hist(ROOT, name, x_label, y_label, bins_a, y_min=0.0, y_max=1.2):
    n = len(bins_a) - 1
    h = ROOT.TH1F(f"h_ax_{name}", f";{x_label};{y_label}", n, bins_a)
    h.SetMinimum(y_min); h.SetMaximum(y_max)
    h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
    h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
    h.GetXaxis().SetLabelSize(0.040); h.GetYaxis().SetLabelSize(0.040)
    h.GetXaxis().SetTitleOffset(1.2); h.GetYaxis().SetTitleOffset(1.3)
    h.SetStats(0)
    return h


def _iqr_band(ROOT, x, q25, q75, color, alpha=0.30):
    """TGraphAsymmErrors IQR band (draw with 'E3')."""
    med  = 0.5 * (q25 + q75)
    el   = np.where(np.isnan(q25), 0.0, np.abs(med - q25))
    eh   = np.where(np.isnan(q75), 0.0, np.abs(q75 - med))
    med  = np.where(np.isnan(med), 0.0, med)
    n    = len(x)
    g    = ROOT.TGraphAsymmErrors(n,
                                   x.astype(float), med.astype(float),
                                   np.zeros(n), np.zeros(n),
                                   el.astype(float), eh.astype(float))
    g.SetFillColor(ROOT.TColor.GetColorTransparent(color, alpha))
    g.SetLineWidth(0)
    return g


def _draw_axis_grid(ROOT, h, logy=False):
    """Vertical grid lines at bin boundaries (mirrors hyddraDiagPlots style)."""
    lines = []
    n = h.GetNbinsX()
    y_lo = h.GetMinimum()
    y_hi = h.GetMaximum()
    for i in range(1, n + 1):
        x = h.GetBinLowEdge(i)
        ln = ROOT.TLine(x, y_lo, x, y_hi)
        ln.SetLineColor(ROOT.kGray + 1)
        ln.SetLineStyle(3)
        ln.SetLineWidth(1)
        ln.Draw()
        lines.append(ln)
    return lines


def _cms_label(ROOT):
    from hyddraDiagPlots.src.style import draw_cms_label
    draw_cms_label()


# ── Plot: per-category efficiency band ────────────────────────────────────────

def plot_eff_band(ROOT, cat_colors, cat_palettes,
                  tdir, cat, band, bins, x_label, y_label, canvas_name):
    if band is None:
        return
    x_cents, per_file, params_list, median, q25, q75 = band
    bins_a    = np.array(bins, dtype=float)
    cat_color = cat_colors[cat]
    palette   = cat_palettes[cat]

    c = _make_canvas(ROOT, canvas_name)
    h_ax = _axis_hist(ROOT, canvas_name, x_label, y_label, bins_a)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    g_band = _iqr_band(ROOT, x_cents, q25, q75, cat_color)
    g_band.Draw("E3 SAME"); c._g_band = g_band

    indiv = []
    for i, (eff_arr, p) in enumerate(zip(per_file, params_list)):
        col   = palette[i % len(palette)]
        lsty  = 1 if p['ctau'] < 0.3 else 2
        eff_c = np.where(np.isnan(eff_arr), -0.05, eff_arr)
        g = ROOT.TGraph(len(x_cents), x_cents.astype(float), eff_c.astype(float))
        g.SetTitle(p['label'])
        g.SetLineColor(col); g.SetLineWidth(1); g.SetLineStyle(lsty)
        g.SetMarkerSize(0)
        g.Draw("L SAME")
        indiv.append(g)
    c._indiv = indiv

    med_c = np.where(np.isnan(median), -0.05, median)
    g_med = ROOT.TGraph(len(x_cents), x_cents.astype(float), med_c.astype(float))
    g_med.SetLineColor(cat_color); g_med.SetLineWidth(3); g_med.SetMarkerSize(0)
    g_med.Draw("L SAME"); c._g_med = g_med

    c._grid = _draw_axis_grid(ROOT, h_ax)

    leg = ROOT.TLegend(0.50, 0.76, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.030)
    leg.AddEntry(g_med,  f"{cat.capitalize()} median ({len(per_file)} samples)", "l")
    leg.AddEntry(g_band, "IQR (25th–75th percentile)", "f")
    leg.Draw(); c._leg = leg

    note = ROOT.TLatex()
    note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.026)
    note.DrawLatex(0.16, 0.84, "Solid: c#tau = 0.1 mm   Dashed: c#tau = 0.5 mm")
    c._note = note

    _cms_label(ROOT)
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Plot: cross-category overlay ───────────────────────────────────────────────

def plot_category_overlay(ROOT, cat_colors,
                           tdir, bands_by_cat, bins,
                           x_label, y_label, canvas_name):
    bins_a = np.array(bins, dtype=float)
    c      = _make_canvas(ROOT, canvas_name)
    h_ax   = _axis_hist(ROOT, canvas_name, x_label, y_label, bins_a)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    graphs, gbands = [], []
    for cat in CATEGORIES:
        band = bands_by_cat.get(cat)
        if band is None:
            continue
        x_c, _, _, median, q25, q75 = band
        col    = cat_colors[cat]
        g_band = _iqr_band(ROOT, x_c, q25, q75, col, alpha=0.20)
        g_band.Draw("E3 SAME"); gbands.append(g_band)

        med_c = np.where(np.isnan(median), -0.05, median)
        g = ROOT.TGraph(len(x_c), x_c.astype(float), med_c.astype(float))
        g.SetTitle(cat.capitalize())
        g.SetLineColor(col); g.SetLineWidth(3)
        g.SetLineStyle(_CAT_LSTYLE[cat]); g.SetMarkerSize(0)
        g.Draw("L SAME")
        graphs.append(g)

    c._graphs = graphs; c._gbands = gbands
    c._grid   = _draw_axis_grid(ROOT, h_ax)

    leg = ROOT.TLegend(0.52, 0.72, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.034)
    for g in graphs:
        leg.AddEntry(g, g.GetTitle(), "l")
    leg.Draw(); c._leg = leg

    _cms_label(ROOT)
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Plot: stage-by-stage efficiency summary ────────────────────────────────────

def plot_stage_summary(ROOT, cat_colors,
                        tdir, all_stats_by_cat, stage_eff_accessor,
                        canvas_name, y_label):
    """
    One line per category: median efficiency at each algorithm stage.
    stage_eff_accessor(stats_dict, stage_key) → float
    """
    n_stages = len(_STAGE_KEYS)
    x        = np.arange(1, n_stages + 1, dtype=float)

    c = ROOT.TCanvas(canvas_name, canvas_name, 900, 600)
    c.SetLeftMargin(0.12); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)

    h_ax = ROOT.TH1F(f"h_ax_{canvas_name}", f";Algorithm Stage;{y_label}",
                     n_stages, 0.5, n_stages + 0.5)
    for i, name in enumerate(_STAGE_NAMES):
        h_ax.GetXaxis().SetBinLabel(i + 1, name)
    h_ax.SetMinimum(0.0); h_ax.SetMaximum(1.2)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.044); h_ax.GetYaxis().SetLabelSize(0.040)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    graphs, gbands = [], []
    for i_cat, cat in enumerate(CATEGORIES):
        cat_stats = all_stats_by_cat.get(cat, [])
        if not cat_stats:
            continue
        rows = [[stage_eff_accessor(s, sk) for sk in _STAGE_KEYS]
                for _, s in cat_stats]
        if not rows:
            continue
        mat    = np.array(rows, dtype=float)
        median = np.nanmedian(mat,         axis=0)
        q25    = np.nanpercentile(mat, 25, axis=0)
        q75    = np.nanpercentile(mat, 75, axis=0)
        col    = cat_colors[cat]

        g_band = _iqr_band(ROOT, x, q25, q75, col, alpha=0.25)
        g_band.Draw("E3 SAME"); gbands.append(g_band)

        med_c = np.where(np.isnan(median), 0.0, median)
        g = ROOT.TGraph(n_stages, x.astype(float), med_c.astype(float))
        g.SetTitle(cat.capitalize())
        g.SetLineColor(col); g.SetLineWidth(3)
        g.SetMarkerStyle(20 + i_cat); g.SetMarkerColor(col); g.SetMarkerSize(1.0)
        g.Draw("LP SAME")
        graphs.append(g)

    c._graphs = graphs; c._gbands = gbands

    leg = ROOT.TLegend(0.60, 0.72, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.034)
    for g in graphs:
        leg.AddEntry(g, g.GetTitle(), "lp")
    leg.Draw(); c._leg = leg

    _cms_label(ROOT)
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Plot: normalized reco distributions per category ──────────────────────────

def plot_reco_dist_categories(ROOT, cat_colors,
                               tdir, all_stats_by_cat, obs_key, obs_cfg,
                               reco_accessor, canvas_name):
    """
    Normalized signal-vertex distributions overlaid per category (median + IQR band).
    reco_accessor(stats_dict) → normalized 1D array or None
    """
    bins_a = np.array(obs_cfg['bins'], dtype=float)
    n_bins = len(bins_a) - 1
    x_c    = 0.5 * (bins_a[:-1] + bins_a[1:])
    logy   = obs_cfg.get('log_y', False)

    c = _make_canvas(ROOT, canvas_name)
    c.SetLogy(logy)
    h_ax = _axis_hist(ROOT, canvas_name, obs_cfg['label'],
                      "Normalised (signal, unit area)",
                      bins_a, y_min=1e-4 if logy else 0.0, y_max=2.0)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    keep  = []
    drawn = False

    for cat in CATEGORIES:
        cat_stats = all_stats_by_cat.get(cat, [])
        per_file  = [reco_accessor(s) for _, s in cat_stats]
        per_file  = [a for a in per_file if a is not None]
        if not per_file:
            continue
        mat    = np.array(per_file, dtype=float)
        median = np.nanmedian(mat,         axis=0)
        q25    = np.nanpercentile(mat, 25, axis=0)
        q75    = np.nanpercentile(mat, 75, axis=0)
        col    = cat_colors[cat]

        g_band = _iqr_band(ROOT, x_c, q25, q75, col, alpha=0.20)
        g_band.Draw("E3 SAME"); keep.append(g_band)

        h = ROOT.TH1F(f"h_{canvas_name}_{cat}", "", n_bins, bins_a)
        h.SetDirectory(0)
        for i, v in enumerate(median, 1):
            h.SetBinContent(i, float(v) if not np.isnan(v) else 0.0)
        h.SetLineColor(col); h.SetLineWidth(2)
        h.SetLineStyle(_CAT_LSTYLE[cat])
        h.SetFillStyle(0); h.SetStats(0)
        h.SetTitle(cat.capitalize())
        h.Draw("HIST SAME"); keep.append(h)
        drawn = True

    if not drawn:
        print(f"    [{canvas_name}] no data — skipping")
        return

    c._keep = keep
    c._grid = _draw_axis_grid(ROOT, h_ax, logy=logy)

    leg = ROOT.TLegend(0.52, 0.72, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.034)
    for obj in keep:
        if isinstance(obj, ROOT.TH1):
            leg.AddEntry(obj, obj.GetTitle(), "l")
    leg.Draw(); c._leg = leg

    _cms_label(ROOT)
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Final stage detection ──────────────────────────────────────────────────────

def _best_final_stage(all_stats_by_cat: dict, leptonic: bool) -> str:
    key = 'lep_eff_dxy' if leptonic else 'had_eff_dxy'
    for cat_stats in all_stats_by_cat.values():
        for _, s in cat_stats:
            for (sk, _other), arr in s.get(key, {}).items():
                if sk == 'id' and arr is not None and not np.all(np.isnan(arr)):
                    return 'id'
    return 'filtered'


# ── Leptonic orchestration ─────────────────────────────────────────────────────

def run_leptonic(ROOT, cat_colors, cat_palettes, all_stats, out_file):
    print("\n[leptonic]")
    final_stage = _best_final_stage(all_stats, leptonic=True)
    print(f"  Final stage: {final_stage}")

    lep_tdir = out_file.mkdir("leptonic")

    for flav_key, flav_field, flav_label in _FLAVOR_CONFIGS:
        print(f"\n  [leptonic/{flav_key}]")
        tflav = lep_tdir.mkdir(flav_key)

        bands_dxy_final: dict = {}
        bands_pt_final:  dict = {}

        for cat in CATEGORIES:
            cat_stats = all_stats.get(cat, [])
            if not cat_stats:
                continue
            tcat = tflav.mkdir(cat)

            stages_to_plot = (["filtered", "id"] if final_stage == "id"
                              else ["filtered"])

            for stage_key in stages_to_plot:
                ff = flav_field
                sk = stage_key

                band_dxy = build_eff_band(
                    cat_stats,
                    lambda s, _sk=sk, _ff=ff: s.get('lep_eff_dxy', {}).get((_sk, _ff)),
                    _GEN_DXY_BINS,
                )
                band_pt = build_eff_band(
                    cat_stats,
                    lambda s, _sk=sk, _ff=ff: s.get('lep_eff_pt', {}).get((_sk, _ff)),
                    _GEN_PT_BINS,
                )
                plot_eff_band(ROOT, cat_colors, cat_palettes,
                              tcat, cat, band_dxy, _GEN_DXY_BINS,
                              "Gen dxy (cm)", f"Gold Efficiency [{sk}]",
                              f"eff_vs_dxy_{sk}")
                plot_eff_band(ROOT, cat_colors, cat_palettes,
                              tcat, cat, band_pt, _GEN_PT_BINS,
                              "Gen p_{T} (GeV)", f"Gold Efficiency [{sk}]",
                              f"eff_vs_pt_{sk}")

                if stage_key == final_stage:
                    bands_dxy_final[cat] = band_dxy
                    bands_pt_final[cat]  = band_pt

        tover = tflav.mkdir("categories")

        plot_category_overlay(ROOT, cat_colors, tover,
                              bands_dxy_final, _GEN_DXY_BINS,
                              "Gen dxy (cm)",
                              f"Gold Efficiency [{final_stage}] ({flav_label})",
                              "eff_vs_dxy_categories")
        plot_category_overlay(ROOT, cat_colors, tover,
                              bands_pt_final, _GEN_PT_BINS,
                              "Gen p_{T} (GeV)",
                              f"Gold Efficiency [{final_stage}] ({flav_label})",
                              "eff_vs_pt_categories")

        ff = flav_field
        plot_stage_summary(
            ROOT, cat_colors, tover, all_stats,
            stage_eff_accessor=lambda s, sk, _ff=ff:
                s.get('lep_stage_eff', {}).get((sk, _ff), float('nan')),
            canvas_name=f"stage_summary_{flav_key}",
            y_label=f"Gold Efficiency ({flav_label})",
        )

        treco = tover.mkdir("reco_distributions")
        for obs_key in _RECO_OBS:
            plot_reco_dist_categories(
                ROOT, cat_colors, treco, all_stats, obs_key, _RECO_OBS[obs_key],
                reco_accessor=lambda s, _sk=final_stage, _ok=obs_key:
                    s.get('lep_reco', {}).get((_sk, _ok)),
                canvas_name=f"reco_{obs_key}_categories",
            )


# ── Hadronic orchestration ─────────────────────────────────────────────────────

def run_hadronic(ROOT, cat_colors, cat_palettes, all_stats, out_file):
    print("\n[hadronic]")
    final_stage = _best_final_stage(all_stats, leptonic=False)
    print(f"  Final stage: {final_stage}")

    had_tdir = out_file.mkdir("hadronic")

    for tier_key, loose, tier_label in _HAD_TIERS:
        print(f"\n  [hadronic/{tier_key}]")
        ttier = had_tdir.mkdir(tier_key)

        bands_dxy_final: dict = {}
        bands_pt_final:  dict = {}

        for cat in CATEGORIES:
            cat_stats = all_stats.get(cat, [])
            if not cat_stats:
                continue
            tcat = ttier.mkdir(cat)

            start_stage   = "merged"
            stages_to_plot = ([start_stage, final_stage]
                              if final_stage != start_stage else [start_stage])

            for stage_key in stages_to_plot:
                _l = loose
                sk = stage_key

                band_dxy = build_eff_band(
                    cat_stats,
                    lambda s, _sk=sk, __l=_l: s.get('had_eff_dxy', {}).get((_sk, __l)),
                    _GEN_DXY_BINS,
                )
                band_pt = build_eff_band(
                    cat_stats,
                    lambda s, _sk=sk, __l=_l: s.get('had_eff_pt', {}).get((_sk, __l)),
                    _GEN_PT_BINS,
                )
                plot_eff_band(ROOT, cat_colors, cat_palettes,
                              tcat, cat, band_dxy, _GEN_DXY_BINS,
                              "Gen dxy (cm)", f"Efficiency [{sk}, {tier_key}]",
                              f"eff_vs_dxy_{sk}")
                plot_eff_band(ROOT, cat_colors, cat_palettes,
                              tcat, cat, band_pt, _GEN_PT_BINS,
                              "Gen p_{T} (GeV)", f"Efficiency [{sk}, {tier_key}]",
                              f"eff_vs_pt_{sk}")

                if stage_key == final_stage:
                    bands_dxy_final[cat] = band_dxy
                    bands_pt_final[cat]  = band_pt

        tover = ttier.mkdir("categories")

        plot_category_overlay(ROOT, cat_colors, tover,
                              bands_dxy_final, _GEN_DXY_BINS,
                              "Gen dxy (cm)",
                              f"Efficiency [{final_stage}, {tier_key}]",
                              "eff_vs_dxy_categories")
        plot_category_overlay(ROOT, cat_colors, tover,
                              bands_pt_final, _GEN_PT_BINS,
                              "Gen p_{T} (GeV)",
                              f"Efficiency [{final_stage}, {tier_key}]",
                              "eff_vs_pt_categories")

        _l = loose
        plot_stage_summary(
            ROOT, cat_colors, tover, all_stats,
            stage_eff_accessor=lambda s, sk, __l=_l:
                s.get('had_stage_eff', {}).get((sk, __l), float('nan')),
            canvas_name=f"stage_summary_{tier_key}",
            y_label=f"Efficiency ({tier_label})",
        )

        treco = ttier.mkdir("reco_distributions")
        for obs_key in _HAD_RECO_OBS:
            plot_reco_dist_categories(
                ROOT, cat_colors, treco, all_stats, obs_key, _HAD_RECO_OBS[obs_key],
                reco_accessor=lambda s, _sk=final_stage, _ok=obs_key:
                    s.get('had_reco', {}).get((_sk, _ok)),
                canvas_name=f"reco_{obs_key}_categories",
            )


# ── CLI ────────────────────────────────────────────────────────────────────────

def _parse_args():
    ap = argparse.ArgumentParser(
        description="Cross-category comparison for SMS-GlGl HYDDRA diagnostic files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--files", nargs="+", required=True,
                    help="Diagnostic ROOT files (globs accepted).")
    ap.add_argument("--output", default="category_comparison.root",
                    help="Output ROOT file (default: category_comparison.root).")
    ap.add_argument("--mode", choices=["auto", "leptonic", "hadronic", "both"],
                    default="auto",
                    help="HYDDRA mode to analyse (default: auto-detect).")
    ap.add_argument("--workers", type=int, default=4,
                    help="Number of parallel file-processing workers (default: 4). "
                         "Use 1 for sequential processing.")
    return ap.parse_args()


def main():
    args = _parse_args()

    # Expand globs
    paths = []
    for pattern in args.files:
        expanded = sorted(glob_module.glob(pattern))
        if expanded:
            paths.extend(expanded)
        elif os.path.exists(pattern):
            paths.append(pattern)
        else:
            print(f"  Warning: no files matched '{pattern}'")
    if not paths:
        print("Error: no input files found.")
        sys.exit(1)
    print(f"Found {len(paths)} input file(s).")

    groups = classify_files(paths)

    # Detect modes from first classified file
    sample_path = next(
        (p for g in groups.values() for p, _ in g if os.path.exists(p)), None
    )
    if sample_path is None:
        print("Error: no valid classified files.")
        sys.exit(1)

    import uproot
    run_lep = run_had = False
    if args.mode == "auto":
        with uproot.open(sample_path) as f:
            try:
                f["hyddraSVsDiagAnalyzer/leptonicConfig"]
                run_lep = True
            except KeyError:
                pass
            try:
                f["hyddraSVsHadronicDiagAnalyzer/hadronicConfig"]
                run_had = True
            except KeyError:
                pass
        print(f"Auto-detected: leptonic={run_lep}, hadronic={run_had}")
    elif args.mode == "leptonic":
        run_lep = True
    elif args.mode == "hadronic":
        run_had = True
    else:
        run_lep = run_had = True

    # ── Parallel stats extraction (no ROOT in workers) ──────────────────────
    n_workers = min(args.workers, len(paths))
    all_stats = compute_all_stats(groups, run_lep, run_had, n_workers)

    # ── Drawing (ROOT in main process only) ─────────────────────────────────
    ROOT       = _init_root()
    cat_colors  = _get_cat_colors(ROOT)
    cat_palettes = _get_cat_palettes(ROOT)

    out = ROOT.TFile(args.output, "RECREATE")

    if run_lep:
        lep_has_data = any(
            s.get('has_lep', False)
            for cat_list in all_stats.values()
            for _, s in cat_list
        )
        if lep_has_data:
            run_leptonic(ROOT, cat_colors, cat_palettes, all_stats, out)
        else:
            print("\n[leptonic] No leptonic data found in any file — skipping.")

    if run_had:
        had_has_data = any(
            s.get('has_had', False)
            for cat_list in all_stats.values()
            for _, s in cat_list
        )
        if had_has_data:
            run_hadronic(ROOT, cat_colors, cat_palettes, all_stats, out)
        else:
            print("\n[hadronic] No hadronic data found in any file — skipping.")

    out.Close()
    print(f"\nDone. Output: {args.output}")


if __name__ == "__main__":
    main()
