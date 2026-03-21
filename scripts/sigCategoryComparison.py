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
            --output comparison.root

    python3 scripts/sigCategoryComparison.py \\
            --files diag1.root diag2.root ... --output out.root \\
            --mode leptonic   # leptonic / hadronic / both (default: auto)

Output ROOT structure:
  leptonic/
    combined/  muon/  electron/
      compressed/    eff_vs_dxy_{filtered,id}   eff_vs_pt_{filtered,id}
      threshold/     (same)
      uncompressed/  (same)
      categories/    eff_vs_dxy_categories  eff_vs_pt_categories
                     stage_summary
                     reco_{obs}_categories   (at final stage)
  hadronic/
    tight/  loose/
      (same layout, no flavor split)
"""

import argparse
import glob as glob_module
import os
import re
import sys

import numpy as np
import awkward as ak
import uproot
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from hyddraDiagPlots.src import loader
from hyddraDiagPlots.src.config import (
    STAGE_KEYS, STAGE_NAMES,
    GEN_DXY_BINS, GEN_PT_BINS,
    RECO_OBSERVABLES, HADRONIC_RECO_OBSERVABLES,
    HAD_MIN3D_CUT,
)
from hyddraDiagPlots.src.style  import draw_cms_label, draw_axis_grid
from hyddraDiagPlots.src.plotter import had_signal_mask

# ── Constants ──────────────────────────────────────────────────────────────────

CATEGORIES = ["compressed", "threshold", "uncompressed"]

_CAT_COLOR = {
    "compressed":   ROOT.kRed   + 1,
    "threshold":    ROOT.kGreen + 2,
    "uncompressed": ROOT.kBlue  + 2,
}

# Sequential palette (up to 6 shades) for individual curves within each category.
# Ordered light → dark so low-ΔM samples are lighter.
_CAT_PALETTE = {
    "compressed":   [ROOT.kRed-7,   ROOT.kRed-4,   ROOT.kRed-1,
                     ROOT.kRed+1,   ROOT.kRed+2,   ROOT.kRed+3],
    "threshold":    [ROOT.kGreen-9, ROOT.kGreen-6, ROOT.kGreen,
                     ROOT.kGreen+1, ROOT.kGreen+2, ROOT.kGreen+3],
    "uncompressed": [ROOT.kAzure-9, ROOT.kAzure-7, ROOT.kAzure+1,
                     ROOT.kAzure+3, ROOT.kBlue+1,  ROOT.kBlue+3],
}

# Line style for cross-category overlay
_CAT_LSTYLE = {"compressed": 1, "threshold": 2, "uncompressed": 3}

# Leptonic flavor configurations: (dir_key, genFunnel_field, legend_label)
FLAVOR_CONFIGS = [
    ("combined", None,          "Combined"),
    ("muon",     "isMuon",      "Muon"),
    ("electron", "isElectron",  "Electron"),
]

# Hadronic matchRatio tiers: (dir_key, loose_bool, legend_label)
HADRONIC_TIERS = [
    ("tight", False, "Tight (matchRatio #geq 0.5)"),
    ("loose", True,  "Loose (matchRatio > 0)"),
]

# ── Filename parsing ───────────────────────────────────────────────────────────

_SMS_RE = re.compile(
    r'mGl-(?P<mGl>\d+).*?mN2-(?P<mN2>\d+).*?mN1-(?P<mN1>\d+)'
    r'.*?N2ctau-(?P<ctau>[\dp]+)',
    re.IGNORECASE,
)


def _ctau_float(s):
    return float(s.replace('p', '.'))


def parse_sms_params(path):
    """Return parameter dict for an SMS-GlGl file, or None if pattern not matched."""
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


def classify_files(paths):
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


# ── Data loading ───────────────────────────────────────────────────────────────

def _load_leptonic(path):
    """Return (gf, sv, cfg) for leptonic data, or None."""
    try:
        with uproot.open(path) as f:
            if not loader.has_leptonic_data(f):
                return None
            return (loader.load_gen_funnel(f),
                    loader.load_all_stage_vtx(f),
                    loader.load_leptonic_config(f))
    except Exception as e:
        print(f"  [load] Warning: {os.path.basename(path)}: {e}")
        return None


def _load_hadronic(path):
    """Return (gf, sv, cfg) for hadronic data, or None."""
    base = loader._HADRONIC_BASE
    try:
        with uproot.open(path) as f:
            if not loader.has_hadronic_data(f):
                return None
            return (loader.load_gen_funnel(f, base=base),
                    loader.load_all_stage_vtx(f, base=base),
                    loader.load_hadronic_config(f))
    except Exception as e:
        print(f"  [load] Warning: {os.path.basename(path)}: {e}")
        return None


def load_all_data(groups, mode_fn):
    """Load all files for all categories.  Returns {cat: [(params, gf, sv, cfg)]}."""
    all_data = {cat: [] for cat in CATEGORIES}
    for cat, files in groups.items():
        print(f"  {cat} ({len(files)} files)...", end="", flush=True)
        n_ok = 0
        for path, params in files:
            data = mode_fn(path)
            if data is not None:
                gf, sv, cfg = data
                all_data[cat].append((params, gf, sv, cfg))
                n_ok += 1
        print(f" {n_ok} loaded")
    return all_data


# ── Per-file efficiency (numpy) ────────────────────────────────────────────────

def _lep_eff_per_bin(gf, var_key, bins_a, stage_key, flavor_field):
    """
    Leptonic gold efficiency per bin.  Returns 1D array (NaN where denom=0) or None.
    """
    gold_key = f"GenFunnel_gold_{stage_key}"
    if gold_key not in gf.fields:
        return None
    has_t = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    var   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{var_key}"]))[has_t]
    gold  = ak.to_numpy(ak.flatten(gf[gold_key]))[has_t].astype(bool)
    if flavor_field is not None:
        fmask = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{flavor_field}"]))[has_t].astype(bool)
        var, gold = var[fmask], gold[fmask]
    denom = np.histogram(var,       bins=bins_a)[0].astype(float)
    numer = np.histogram(var[gold], bins=bins_a)[0].astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(denom > 0, numer / denom, np.nan)


def _had_eff_per_bin(gf, var_key, bins_a, stage_key, loose):
    """
    Hadronic signal efficiency per bin.  Returns 1D array (NaN where denom=0) or None.
    """
    mr_key = f"GenFunnel_matchRatio_{stage_key}"
    if mr_key not in gf.fields:
        return None
    is_had = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    var    = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{var_key}"]))[is_had]
    sig    = had_signal_mask(gf, stage_key, loose=loose)[is_had]
    denom  = np.histogram(var,      bins=bins_a)[0].astype(float)
    numer  = np.histogram(var[sig], bins=bins_a)[0].astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(denom > 0, numer / denom, np.nan)


def _lep_stage_eff(gf, stage_key, flavor_field):
    """Single leptonic gold efficiency value at one stage (scalar or NaN)."""
    gold_key = f"GenFunnel_gold_{stage_key}"
    if gold_key not in gf.fields:
        return np.nan
    has_t = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    gold  = ak.to_numpy(ak.flatten(gf[gold_key]))[has_t].astype(bool)
    if flavor_field is not None:
        fmask  = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{flavor_field}"]))[has_t].astype(bool)
        n_denom = int(np.sum(fmask))
        if n_denom == 0:
            return np.nan
        return float(np.sum(gold[fmask])) / n_denom
    n_denom = int(np.sum(has_t))
    return float(np.sum(gold)) / n_denom if n_denom > 0 else np.nan


def _had_stage_eff(gf, stage_key, loose):
    """Single hadronic signal efficiency value at one stage (scalar or NaN)."""
    mr_key = f"GenFunnel_matchRatio_{stage_key}"
    if mr_key not in gf.fields:
        return np.nan
    is_had  = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    n_denom = int(np.sum(is_had))
    if n_denom == 0:
        return np.nan
    sig = had_signal_mask(gf, stage_key, loose=loose)[is_had]
    return float(np.sum(sig)) / n_denom


# ── Band computation ───────────────────────────────────────────────────────────

def build_eff_band(cat_data, per_file_fn, bins):
    """
    Compute per-file efficiency arrays and derive median + IQR band.

    per_file_fn(gf, bins_array) → 1D numpy array or None

    Returns (x_centers, per_file_list, params_list, median, q25, q75) or None.
    """
    bins_a   = np.array(bins, dtype=float)
    x_cents  = 0.5 * (bins_a[:-1] + bins_a[1:])
    per_file, params = [], []
    for p, gf, sv, cfg in cat_data:
        arr = per_file_fn(gf, bins_a)
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


# ── ROOT helpers ───────────────────────────────────────────────────────────────

def _make_canvas(name, w=800, h=600):
    c = ROOT.TCanvas(name, name, w, h)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)
    return c


def _axis_hist(name, x_label, y_label, bins_a, y_min=0.0, y_max=1.2):
    n = len(bins_a) - 1
    h = ROOT.TH1F(f"h_ax_{name}", f";{x_label};{y_label}", n, bins_a)
    h.SetMinimum(y_min); h.SetMaximum(y_max)
    h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
    h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
    h.GetXaxis().SetLabelSize(0.040); h.GetYaxis().SetLabelSize(0.040)
    h.GetXaxis().SetTitleOffset(1.2); h.GetYaxis().SetTitleOffset(1.3)
    h.SetStats(0)
    return h


def _iqr_band_graph(x, q25, q75, cat_color, alpha=0.30):
    """TGraphAsymmErrors filled band spanning q25 – q75 (draw with 'E3')."""
    med  = 0.5 * (q25 + q75)
    el   = np.where(np.isnan(q25), 0.0, np.abs(med - q25))
    eh   = np.where(np.isnan(q75), 0.0, np.abs(q75 - med))
    med  = np.where(np.isnan(med), 0.0, med)
    n    = len(x)
    g    = ROOT.TGraphAsymmErrors(n,
                                   x.astype(float), med.astype(float),
                                   np.zeros(n), np.zeros(n),
                                   el.astype(float), eh.astype(float))
    g.SetFillColor(ROOT.TColor.GetColorTransparent(cat_color, alpha))
    g.SetLineWidth(0)
    return g


# ── Plot: per-category efficiency band ────────────────────────────────────────

def plot_eff_band(tdir, cat, band, bins, x_label, y_label, canvas_name):
    """
    Individual thin curves + IQR band + median for a single category.
    Individual curves: colored by position in sorted sample list,
                       solid = ctau 0.1, dashed = ctau 0.5.
    """
    if band is None:
        return
    x_cents, per_file, params_list, median, q25, q75 = band
    bins_a    = np.array(bins, dtype=float)
    cat_color = _CAT_COLOR[cat]
    palette   = _CAT_PALETTE[cat]

    c = _make_canvas(canvas_name)
    h_ax = _axis_hist(canvas_name, x_label, y_label, bins_a)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    # IQR band (draw first so curves go on top)
    g_band = _iqr_band_graph(x_cents, q25, q75, cat_color)
    g_band.Draw("E3 SAME"); c._g_band = g_band

    # Individual thin curves
    indiv = []
    for i, (eff_arr, p) in enumerate(zip(per_file, params_list)):
        col  = palette[i % len(palette)]
        lsty = 1 if p['ctau'] < 0.3 else 2
        # Push NaN below axis so ROOT draws a gap automatically
        eff_c = np.where(np.isnan(eff_arr), -0.05, eff_arr)
        g = ROOT.TGraph(len(x_cents), x_cents.astype(float), eff_c.astype(float))
        g.SetTitle(p['label'])
        g.SetLineColor(col); g.SetLineWidth(1); g.SetLineStyle(lsty)
        g.SetMarkerSize(0)
        g.Draw("L SAME")
        indiv.append(g)
    c._indiv = indiv

    # Median (drawn last, on top)
    med_c = np.where(np.isnan(median), -0.05, median)
    g_med = ROOT.TGraph(len(x_cents), x_cents.astype(float), med_c.astype(float))
    g_med.SetLineColor(cat_color); g_med.SetLineWidth(3)
    g_med.SetMarkerSize(0)
    g_med.Draw("L SAME"); c._g_med = g_med

    c._grid = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.50, 0.76, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.030)
    leg.AddEntry(g_med,  f"{cat.capitalize()} median ({len(per_file)} samples)", "l")
    leg.AddEntry(g_band, "IQR (25th–75th percentile)", "f")
    leg.Draw(); c._leg = leg

    note = ROOT.TLatex()
    note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.026)
    note.DrawLatex(0.16, 0.84, "Solid: c#tau = 0.1 mm   Dashed: c#tau = 0.5 mm")
    c._note = note

    draw_cms_label()
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Plot: cross-category efficiency overlay ────────────────────────────────────

def plot_category_overlay(tdir, bands_by_cat, bins, x_label, y_label, canvas_name):
    """Three median curves (one per category) with IQR bands on one canvas."""
    bins_a = np.array(bins, dtype=float)
    c      = _make_canvas(canvas_name)
    h_ax   = _axis_hist(canvas_name, x_label, y_label, bins_a)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    graphs, gbands = [], []
    for cat in CATEGORIES:
        band = bands_by_cat.get(cat)
        if band is None:
            continue
        x_c, _, _, median, q25, q75 = band
        col    = _CAT_COLOR[cat]
        g_band = _iqr_band_graph(x_c, q25, q75, col, alpha=0.20)
        g_band.Draw("E3 SAME"); gbands.append(g_band)

        med_c = np.where(np.isnan(median), -0.05, median)
        g = ROOT.TGraph(len(x_c), x_c.astype(float), med_c.astype(float))
        g.SetTitle(cat.capitalize())
        g.SetLineColor(col); g.SetLineWidth(3)
        g.SetLineStyle(_CAT_LSTYLE[cat]); g.SetMarkerSize(0)
        g.Draw("L SAME")
        graphs.append(g)

    c._graphs = graphs; c._gbands = gbands
    c._grid   = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.52, 0.72, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.034)
    for g in graphs:
        leg.AddEntry(g, g.GetTitle(), "l")
    leg.Draw(); c._leg = leg

    draw_cms_label()
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Plot: stage-by-stage efficiency summary ────────────────────────────────────

def plot_stage_summary(tdir, all_data, stage_eff_fn, canvas_name, y_label):
    """
    One line per category showing median gold efficiency at each algorithm stage.
    stage_eff_fn(gf, stage_key) → float  (different for leptonic vs hadronic)
    """
    n_stages = len(STAGE_KEYS)
    x        = np.arange(1, n_stages + 1, dtype=float)

    c = ROOT.TCanvas(canvas_name, canvas_name, 900, 600)
    c.SetLeftMargin(0.12); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)

    h_ax = ROOT.TH1F(f"h_ax_{canvas_name}", f";Algorithm Stage;{y_label}",
                     n_stages, 0.5, n_stages + 0.5)
    for i, name in enumerate(STAGE_NAMES):
        h_ax.GetXaxis().SetBinLabel(i + 1, name)
    h_ax.SetMinimum(0.0); h_ax.SetMaximum(1.2)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.044); h_ax.GetYaxis().SetLabelSize(0.040)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    graphs, gbands = [], []
    for i_cat, cat in enumerate(CATEGORIES):
        cat_data = all_data.get(cat, [])
        if not cat_data:
            continue
        rows = [[stage_eff_fn(gf, sk) for sk in STAGE_KEYS]
                for _, gf, sv, cfg in cat_data]
        if not rows:
            continue
        mat    = np.array(rows)
        median = np.nanmedian(mat,         axis=0)
        q25    = np.nanpercentile(mat, 25, axis=0)
        q75    = np.nanpercentile(mat, 75, axis=0)
        col    = _CAT_COLOR[cat]

        g_band = _iqr_band_graph(x, q25, q75, col, alpha=0.25)
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

    draw_cms_label()
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Plot: normalized reco distributions per category ──────────────────────────

def plot_reco_dist_categories(tdir, all_data, stage_key, obs_key, obs_cfg,
                              sig_mask_fn, canvas_name):
    """
    Normalized signal-vertex distributions at stage_key, overlaid per category.
    Median histogram + IQR band per category, drawn on one canvas.

    sig_mask_fn(sv) → boolean 1D array selecting signal vertices.
    """
    sidx   = STAGE_KEYS.index(stage_key) if stage_key in STAGE_KEYS else -1
    bins_a = np.array(obs_cfg["bins"], dtype=float)
    n_bins = len(bins_a) - 1
    x_c    = 0.5 * (bins_a[:-1] + bins_a[1:])
    logy   = obs_cfg.get("log_y", False)

    c = _make_canvas(canvas_name)
    c.SetLogy(logy)
    h_ax = _axis_hist(canvas_name, obs_cfg["label"],
                      "Normalised (signal, unit area)",
                      bins_a,
                      y_min=1e-4 if logy else 0.0,
                      y_max=2.0)
    h_ax.Draw("AXIS"); c._h_ax = h_ax

    keep_alive = []   # prevent ROOT GC
    any_drawn  = False

    for cat in CATEGORIES:
        cat_data  = all_data.get(cat, [])
        per_file  = []
        branch    = f"StageVtx_{obs_key}"
        for _, gf, sv, cfg in cat_data:
            stg_mask = ak.to_numpy(sv["StageVtx_stageIdx"]).astype(int) == sidx
            sig_mask = sig_mask_fn(sv)
            sel      = stg_mask & sig_mask
            if branch not in sv.fields or not np.any(sel):
                continue
            vals  = ak.to_numpy(sv[branch]).astype(float)[sel]
            cnts, _ = np.histogram(vals, bins=bins_a)
            tot = float(cnts.sum())
            if tot > 0:
                per_file.append(cnts / tot)

        if not per_file:
            continue
        mat    = np.array(per_file)
        median = np.nanmedian(mat,         axis=0)
        q25    = np.nanpercentile(mat, 25, axis=0)
        q75    = np.nanpercentile(mat, 75, axis=0)
        col    = _CAT_COLOR[cat]

        # IQR band
        g_band = _iqr_band_graph(x_c, q25, q75, col, alpha=0.20)
        g_band.Draw("E3 SAME"); keep_alive.append(g_band)

        # Median histogram
        h = ROOT.TH1F(f"h_{canvas_name}_{cat}", "", n_bins, bins_a)
        h.SetDirectory(0)
        for i, v in enumerate(median, 1):
            h.SetBinContent(i, float(v) if not np.isnan(v) else 0.0)
        h.SetLineColor(col); h.SetLineWidth(2)
        h.SetLineStyle(_CAT_LSTYLE[cat])
        h.SetFillStyle(0); h.SetStats(0)
        h.SetTitle(cat.capitalize())
        h.Draw("HIST SAME"); keep_alive.append(h)
        any_drawn = True

    if not any_drawn:
        print(f"    [{canvas_name}] no data — skipping")
        return

    c._keep = keep_alive
    c._grid = draw_axis_grid(h_ax, logy=logy)

    leg = ROOT.TLegend(0.52, 0.72, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.034)
    for obj in keep_alive:
        if isinstance(obj, ROOT.TH1):
            leg.AddEntry(obj, obj.GetTitle(), "l")
    leg.Draw(); c._leg = leg

    draw_cms_label()
    c.Update()
    tdir.cd(); c.Write()
    print(f"    [{canvas_name}] done")


# ── Detect best final stage ────────────────────────────────────────────────────

def _best_final_stage(all_data, leptonic=True):
    """Return 'id' if any file has the id stage, else 'filtered'."""
    key = "GenFunnel_gold_id" if leptonic else "GenFunnel_matchRatio_id"
    for cat_data in all_data.values():
        for _, gf, _, _ in cat_data:
            if key in gf.fields:
                return "id"
    return "filtered"


# ── Leptonic orchestration ─────────────────────────────────────────────────────

def run_leptonic(all_data, out_file):
    print("\n[leptonic]")
    final_stage = _best_final_stage(all_data, leptonic=True)
    print(f"  Final stage: {final_stage}")

    lep_tdir = out_file.mkdir("leptonic")

    lep_sig_fn = lambda sv: ak.to_numpy(sv["StageVtx_isGold"]).astype(bool)

    for flav_key, flav_field, flav_label in FLAVOR_CONFIGS:
        print(f"\n  [leptonic/{flav_key}]")
        tflav = lep_tdir.mkdir(flav_key)

        # -- Per-category efficiency bands ------------------------------------
        bands_dxy_final = {}
        bands_pt_final  = {}

        for cat in CATEGORIES:
            cat_data = all_data.get(cat, [])
            if not cat_data:
                continue
            tcat = tflav.mkdir(cat)

            for stage_key in ["filtered", final_stage]:
                if stage_key == "filtered" and final_stage == "filtered":
                    pass  # only one stage to show
                elif stage_key == final_stage and stage_key == "filtered":
                    continue  # already done above

                ff = flav_field  # capture for closure
                sk = stage_key

                band_dxy = build_eff_band(
                    cat_data,
                    lambda gf, ba, _sk=sk, _ff=ff:
                        _lep_eff_per_bin(gf, "dxy", ba, _sk, _ff),
                    GEN_DXY_BINS,
                )
                band_pt = build_eff_band(
                    cat_data,
                    lambda gf, ba, _sk=sk, _ff=ff:
                        _lep_eff_per_bin(gf, "pt", ba, _sk, _ff),
                    GEN_PT_BINS,
                )
                plot_eff_band(tcat, cat, band_dxy, GEN_DXY_BINS,
                              "Gen dxy (cm)", f"Gold Efficiency [{sk}]",
                              f"eff_vs_dxy_{sk}")
                plot_eff_band(tcat, cat, band_pt, GEN_PT_BINS,
                              "Gen p_{T} (GeV)", f"Gold Efficiency [{sk}]",
                              f"eff_vs_pt_{sk}")

                if stage_key == final_stage:
                    bands_dxy_final[cat] = band_dxy
                    bands_pt_final[cat]  = band_pt

        # -- Cross-category overlays ------------------------------------------
        tover = tflav.mkdir("categories")

        plot_category_overlay(
            tover, bands_dxy_final, GEN_DXY_BINS,
            "Gen dxy (cm)",
            f"Gold Efficiency [{final_stage}] ({flav_label})",
            "eff_vs_dxy_categories",
        )
        plot_category_overlay(
            tover, bands_pt_final, GEN_PT_BINS,
            "Gen p_{T} (GeV)",
            f"Gold Efficiency [{final_stage}] ({flav_label})",
            "eff_vs_pt_categories",
        )

        # -- Stage summary ----------------------------------------------------
        ff = flav_field
        plot_stage_summary(
            tover, all_data,
            stage_eff_fn=lambda gf, sk, _ff=ff: _lep_stage_eff(gf, sk, _ff),
            canvas_name=f"stage_summary_{flav_key}",
            y_label=f"Gold Efficiency ({flav_label})",
        )

        # -- Reco distributions at final stage --------------------------------
        treco = tover.mkdir("reco_distributions")
        sidx  = STAGE_KEYS.index(final_stage) if final_stage in STAGE_KEYS else -1
        for obs_key, obs_cfg in RECO_OBSERVABLES.items():
            plot_reco_dist_categories(
                treco, all_data, final_stage, obs_key, obs_cfg,
                lep_sig_fn, f"reco_{obs_key}_categories",
            )


# ── Hadronic orchestration ─────────────────────────────────────────────────────

def run_hadronic(all_data, out_file):
    print("\n[hadronic]")
    final_stage = _best_final_stage(all_data, leptonic=False)
    print(f"  Final stage: {final_stage}")

    had_tdir = out_file.mkdir("hadronic")

    def had_sig_fn(sv):
        mr    = ak.to_numpy(sv["StageVtx_matchRatio"]).astype(float)
        mask  = mr > 0
        if "StageVtx_min3D" in sv.fields:
            min3d = ak.to_numpy(sv["StageVtx_min3D"]).astype(float)
            mask  = mask & (min3d >= 0) & (min3d < HAD_MIN3D_CUT)
        return mask

    for tier_key, loose, tier_label in HADRONIC_TIERS:
        print(f"\n  [hadronic/{tier_key}]")
        ttier = had_tdir.mkdir(tier_key)

        bands_dxy_final = {}
        bands_pt_final  = {}

        for cat in CATEGORIES:
            cat_data = all_data.get(cat, [])
            if not cat_data:
                continue
            tcat = ttier.mkdir(cat)

            start_stage = "merged"  # hadronic SVs first form at merging
            for stage_key in [start_stage, final_stage]:
                if stage_key == start_stage and final_stage == start_stage:
                    pass
                elif stage_key == final_stage and stage_key == start_stage:
                    continue

                _l  = loose
                sk  = stage_key
                band_dxy = build_eff_band(
                    cat_data,
                    lambda gf, ba, _sk=sk, __l=_l:
                        _had_eff_per_bin(gf, "dxy", ba, _sk, __l),
                    GEN_DXY_BINS,
                )
                band_pt = build_eff_band(
                    cat_data,
                    lambda gf, ba, _sk=sk, __l=_l:
                        _had_eff_per_bin(gf, "pt", ba, _sk, __l),
                    GEN_PT_BINS,
                )
                plot_eff_band(tcat, cat, band_dxy, GEN_DXY_BINS,
                              "Gen dxy (cm)", f"Efficiency [{sk}, {tier_key}]",
                              f"eff_vs_dxy_{sk}")
                plot_eff_band(tcat, cat, band_pt, GEN_PT_BINS,
                              "Gen p_{T} (GeV)", f"Efficiency [{sk}, {tier_key}]",
                              f"eff_vs_pt_{sk}")

                if stage_key == final_stage:
                    bands_dxy_final[cat] = band_dxy
                    bands_pt_final[cat]  = band_pt

        tover = ttier.mkdir("categories")

        plot_category_overlay(
            tover, bands_dxy_final, GEN_DXY_BINS,
            "Gen dxy (cm)",
            f"Efficiency [{final_stage}, {tier_key}]",
            "eff_vs_dxy_categories",
        )
        plot_category_overlay(
            tover, bands_pt_final, GEN_PT_BINS,
            "Gen p_{T} (GeV)",
            f"Efficiency [{final_stage}, {tier_key}]",
            "eff_vs_pt_categories",
        )

        _l = loose
        plot_stage_summary(
            tover, all_data,
            stage_eff_fn=lambda gf, sk, __l=_l: _had_stage_eff(gf, sk, __l),
            canvas_name=f"stage_summary_{tier_key}",
            y_label=f"Efficiency ({tier_label})",
        )

        treco = tover.mkdir("reco_distributions")
        for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
            plot_reco_dist_categories(
                treco, all_data, final_stage, obs_key, obs_cfg,
                had_sig_fn, f"reco_{obs_key}_categories",
            )


# ── CLI ────────────────────────────────────────────────────────────────────────

def _parse_args():
    ap = argparse.ArgumentParser(
        description="Cross-category kinematic comparison for SMS-GlGl diagnostic files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--files", nargs="+", required=True,
                    help="Diagnostic ROOT files (globs OK).")
    ap.add_argument("--output", default="category_comparison.root",
                    help="Output ROOT file (default: category_comparison.root).")
    ap.add_argument("--mode", choices=["auto", "leptonic", "hadronic", "both"],
                    default="auto",
                    help="HYDDRA mode to analyse (default: auto-detect from first file).")
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

    # Detect modes from first available file
    sample_path = next(
        (p for g in groups.values() for p, _ in g if os.path.exists(p)), None
    )
    if sample_path is None:
        print("Error: no valid classified files found.")
        sys.exit(1)

    run_lep = run_had = False
    if args.mode == "auto":
        with uproot.open(sample_path) as f:
            run_lep = loader.has_leptonic_data(f)
            run_had = loader.has_hadronic_data(f)
        print(f"Auto-detected: leptonic={run_lep}, hadronic={run_had}")
    elif args.mode == "leptonic":
        run_lep = True
    elif args.mode == "hadronic":
        run_had = True
    else:
        run_lep = run_had = True

    out = ROOT.TFile(args.output, "RECREATE")

    if run_lep:
        print("\n--- Loading leptonic data ---")
        lep_data = load_all_data(groups, _load_leptonic)
        run_leptonic(lep_data, out)

    if run_had:
        print("\n--- Loading hadronic data ---")
        had_data = load_all_data(groups, _load_hadronic)
        run_hadronic(had_data, out)

    out.Close()
    print(f"\nDone. Output: {args.output}")


if __name__ == "__main__":
    main()
