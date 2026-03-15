"""
stages/hadronic_summary.py — cross-stage summary plots for Hadronic HYDDRA.

Signal categorisation uses matchRatio thresholds applied at plot time:
  loose : matchRatio > 0   (at least one signal track)
  tight : matchRatio >= 0.5 (majority signal tracks)

Plots:
  - yield_flow_unnorm / yield_flow_norm
  - efficiency_funnel
  - eff_vs_dxy / eff_vs_pt
  - loss_stage_distribution
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config import (
    STAGE_NAMES, STAGE_KEYS,
    COLOR_BKG, COLOR_GOLD, COLOR_SILVER, COLORS_STAGE, MARKERS,
    GEN_DXY_BINS, GEN_PT_BINS,
)
from ..src.style  import draw_cms_label, make_canvas, draw_grid_lines, draw_axis_grid
from ..src.plotter import setup_axis_hist, make_tgraph, add_legend

_COLOR_TIGHT = COLOR_GOLD
_COLOR_LOOSE = COLOR_SILVER


# ── Yield flow ────────────────────────────────────────────────────────────────

def plot_yield_flow(tdir, gf, sc_sig, sc_bkg=None):
    """
    Total vertex yield per stage (Stage_n_*), plus gen-level loose/tight signal
    counts derived from GenFunnel_matchRatio_*.
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool) if gf is not None and len(gf) > 0 else None

    for normalized in [False, True]:
        suffix = "norm" if normalized else "unnorm"
        canvas = make_canvas(
            f"yield_flow_{suffix}",
            f"Hadronic HYDDRA Yield Flow ({'Normalized' if normalized else 'Raw'})",
            logy=(not normalized),
        )

        total_reco = [sc_sig.get(f"Stage_n_{s}", 0) for s in STAGE_KEYS]

        if has_tracks is not None:
            loose_count = []
            tight_count = []
            for s in STAGE_KEYS:
                mr = ak.to_numpy(ak.flatten(gf[f"GenFunnel_matchRatio_{s}"]))[has_tracks]
                loose_count.append(int(np.sum(mr > 0)))
                tight_count.append(int(np.sum(mr >= 0.5)))
        else:
            loose_count = [0] * len(STAGE_KEYS)
            tight_count = [0] * len(STAGE_KEYS)

        if sc_bkg:
            bkg_total = [sc_bkg.get(f"Stage_n_{s}", 0) for s in STAGE_KEYS]
        else:
            bkg_total = [max(0, total_reco[i] - loose_count[i]) for i in range(len(STAGE_KEYS))]

        norm_tight = tight_count[0] if (normalized and tight_count[0] > 0) else 1
        norm_loose = loose_count[0] if (normalized and loose_count[0] > 0) else 1
        norm_bkg   = bkg_total[0]   if (normalized and bkg_total[0]   > 0) else 1

        if normalized:
            tight_y = [v / norm_tight for v in tight_count]
            loose_y = [v / norm_loose for v in loose_count]
            bkg_y   = [v / norm_bkg   for v in bkg_total]
            y_title = "Survival Fraction (relative to Seeding)"
            min_y, max_y = 1e-6, 1.2
        else:
            tight_y = [float(v) for v in tight_count]
            loose_y = [float(v) for v in loose_count]
            bkg_y   = [float(v) for v in bkg_total]
            y_title = "Vertex Yield"
            all_vals = tight_y + loose_y + bkg_y
            nonzero  = [v for v in all_vals if v > 0]
            min_y    = max(0.5, min(nonzero) / 2) if nonzero else 0.5
            max_y    = max(nonzero) * 5 if nonzero else 10

        h = setup_axis_hist(canvas, STAGE_NAMES, y_title, min_y, max_y)
        h.Draw("AXIS")

        x        = np.arange(1, len(STAGE_NAMES) + 1, dtype=float)
        bkg_label = "Background (total)" if sc_bkg else "Non-signal SVs (signal file)"
        g_tight = make_tgraph(x, tight_y, _COLOR_TIGHT, 20, "Signal tight (matchRatio #geq 0.5)", 1)
        g_loose = make_tgraph(x, loose_y, _COLOR_LOOSE, 24, "Signal loose (matchRatio > 0)",      1)
        g_bkg   = make_tgraph(x, bkg_y,   COLOR_BKG,   21, bkg_label,                             2)
        graphs  = [g_tight, g_loose, g_bkg]

        mg = ROOT.TMultiGraph()
        for g in graphs:
            mg.Add(g)
        mg.Draw("lp")
        canvas._mg = mg
        canvas._graphs = graphs

        add_legend(canvas, graphs, x1=0.50, y1=0.72, x2=0.88, y2=0.88)
        draw_grid_lines(canvas, len(STAGE_NAMES), min_y, max_y, logy=not normalized)
        draw_cms_label()
        canvas.Update()
        tdir.cd()
        canvas.Write()
        print(f"    [yield_flow_{suffix}] done")


# ── Efficiency funnel ─────────────────────────────────────────────────────────

def plot_efficiency_funnel(tdir, gf):
    """
    Efficiency per stage using matchRatio thresholds:
      tight : matchRatio >= 0.5
      loose : matchRatio > 0 (excl. tight)
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    n_denom    = int(np.sum(has_tracks))
    if n_denom == 0:
        print("    [efficiency_funnel] No gen vertices with tracks — skipping")
        return

    canvas = make_canvas("efficiency_funnel",
                         "Hadronic HYDDRA — Signal Efficiency per Stage", logy=False)
    canvas.SetGridx(0)

    n_stages   = len(STAGE_KEYS)
    eff_tight  = np.zeros(n_stages)
    eff_loose  = np.zeros(n_stages)

    for i, s in enumerate(STAGE_KEYS):
        mr = ak.to_numpy(ak.flatten(gf[f"GenFunnel_matchRatio_{s}"]))[has_tracks]
        eff_tight[i] = np.sum(mr >= 0.5)        / n_denom
        eff_loose[i] = np.sum((mr > 0) & (mr < 0.5)) / n_denom

    h_tight = ROOT.TH1F("h_fgt", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)
    h_loose = ROOT.TH1F("h_fgl", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)

    for i, name in enumerate(STAGE_NAMES):
        for h in (h_tight, h_loose):
            h.GetXaxis().SetBinLabel(i + 1, name)
        h_tight.SetBinContent(i + 1, eff_tight[i])
        h_loose.SetBinContent(i + 1, eff_loose[i])

    h_tight.SetFillColor(_COLOR_TIGHT); h_tight.SetLineColor(_COLOR_TIGHT)
    h_loose.SetFillColor(_COLOR_LOOSE); h_loose.SetLineColor(_COLOR_LOOSE)

    max_y = max((eff_tight + eff_loose).max(), 0.01) * 1.3
    h_tight.SetMinimum(1e-4)
    h_tight.SetMaximum(max_y)
    h_tight.GetYaxis().SetTitle("Efficiency")
    h_tight.GetXaxis().CenterTitle(True)
    h_tight.GetYaxis().CenterTitle(True)
    h_tight.SetStats(0)
    h_tight.Draw("HIST")

    hs = ROOT.THStack("hs_funnel", "")
    hs.Add(h_loose)
    hs.Add(h_tight)
    hs.Draw("HIST SAME")
    draw_grid_lines(canvas, n_stages, 1e-4, max_y, logy=False)

    leg = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0)
    leg.AddEntry(h_tight, "Tight (matchRatio #geq 0.5)",            "f")
    leg.AddEntry(h_loose, "Loose excl. (0 < matchRatio < 0.5)",     "f")
    leg.Draw()
    canvas._leg = leg; canvas._hs = hs
    canvas._htight = h_tight; canvas._hloose = h_loose

    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    print("    [efficiency_funnel] done")


# ── Efficiency vs gen variable ─────────────────────────────────────────────────

def plot_efficiency_vs_var(tdir, gf, var_key, var_label, bins, canvas_name):
    """
    Tight and loose signal efficiency vs a gen variable (dxy, pT).
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    var_vals   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{var_key}"]))[has_tracks]
    n_bins_ax  = len(bins) - 1

    canvas = make_canvas(canvas_name,
                         f"Hadronic HYDDRA — Signal efficiency vs {var_label}", logy=False)
    canvas.SetGridx(0)

    h_denom = ROOT.TH1F(f"h_denom_{var_key}", "", n_bins_ax, np.array(bins, dtype=float))
    for v in var_vals:
        h_denom.Fill(v)

    graphs, hists = [], []
    for i, s in enumerate(STAGE_KEYS):
        mr = ak.to_numpy(ak.flatten(gf[f"GenFunnel_matchRatio_{s}"]))[has_tracks]
        for threshold, label_suffix, line_style in [
            (0.5, "tight", 1),
            (0.0, "loose", 2),
        ]:
            sig_mask = mr >= 0.5 if threshold == 0.5 else mr > 0
            h_num = ROOT.TH1F(f"h_num_{s}_{var_key}_{label_suffix}", "",
                              n_bins_ax, np.array(bins, dtype=float))
            for v in var_vals[sig_mask]:
                h_num.Fill(v)
            eff = ROOT.TEfficiency(h_num, h_denom)
            eff.SetStatisticOption(ROOT.TEfficiency.kFCP)
            g = eff.CreateGraph()
            g.SetTitle(f"{STAGE_NAMES[i]} ({label_suffix})")
            g.SetMarkerColor(COLORS_STAGE[i]); g.SetLineColor(COLORS_STAGE[i])
            g.SetMarkerStyle(MARKERS[i]);      g.SetLineWidth(2)
            g.SetLineStyle(line_style)
            graphs.append(g)
            hists.append((h_num, eff))

    h_ax = ROOT.TH1F(f"h_ax_{var_key}", f";{var_label};Signal Efficiency",
                     n_bins_ax, np.array(bins, dtype=float))
    h_ax.SetMinimum(1e-4); h_ax.SetMaximum(1.1)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")
    canvas._h_ax = h_ax; canvas._hists = hists; canvas._graphs = graphs

    for g in graphs:
        g.Draw("P SAME")

    canvas._grid_lines = draw_axis_grid(h_ax, logy=False)
    add_legend(canvas, graphs, x1=0.55, y1=0.45, x2=0.88, y2=0.88)
    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    print(f"    [eff_vs_{var_key}] done")


# ── Loss stage distribution ───────────────────────────────────────────────────

def plot_loss_stage_distribution(tdir, gf):
    """
    Bar chart: how many gen vertices (with tracks) first lose their signal match
    at each pipeline stage.  A vertex is 'lost' at stage S if matchRatio_{S} <= 0
    (or -1) but matchRatio_{S-1} > 0.  Vertices never matched at any stage are
    counted as 'Never found'.  Vertices with a loose match at the final stage
    are counted as 'Survived'.
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)

    # Build a [n_gen_vtx, n_stages] array of booleans: True = loose match
    matched_at = np.zeros((int(np.sum(has_tracks)), len(STAGE_KEYS)), dtype=bool)
    for i, s in enumerate(STAGE_KEYS):
        mr = ak.to_numpy(ak.flatten(gf[f"GenFunnel_matchRatio_{s}"]))[has_tracks]
        matched_at[:, i] = mr > 0

    # Classify each gen vertex
    # -1 = survived, 0 = never found, 1..N = lost at stage index (1-based)
    n_gen    = matched_at.shape[0]
    lost_idx = np.full(n_gen, -1, dtype=int)  # default: survived

    never_found = ~matched_at[:, 0]  # not matched even at seeding
    lost_idx[never_found] = 0

    for i in range(1, len(STAGE_KEYS)):
        lost_here = matched_at[:, i - 1] & ~matched_at[:, i] & (lost_idx == -1)
        lost_idx[lost_here] = i

    labels = ["Never found"] + [f"Lost at {STAGE_NAMES[i]}" for i in range(1, len(STAGE_KEYS))]
    n_cats  = len(labels)
    counts  = np.array([np.sum(lost_idx == v) for v in range(0, n_cats)], dtype=float)
    n_surv  = int(np.sum(lost_idx == -1))

    canvas = make_canvas("loss_stage_distribution",
                         "Hadronic HYDDRA — Signal vertex loss per stage", logy=False)
    canvas.SetGridx(0)

    h = ROOT.TH1F("h_had_loss", ";Loss Category;Count", n_cats, 0.5, n_cats + 0.5)
    for i, (c, lbl) in enumerate(zip(counts, labels)):
        h.SetBinContent(i + 1, c)
        h.GetXaxis().SetBinLabel(i + 1, lbl)
    h.SetFillColor(ROOT.kOrange + 1)
    h.SetLineColor(ROOT.kOrange + 1)
    h.GetYaxis().SetTitle("Gen vertex count")
    h.GetXaxis().CenterTitle(True)
    h.GetYaxis().CenterTitle(True)
    h.SetStats(0)
    max_y_val = max(counts.max() * 1.4, 1.0)
    h.SetMinimum(1e-4)
    h.SetMaximum(max_y_val)
    h.Draw("HIST")
    draw_grid_lines(canvas, n_cats, 1e-4, max_y_val, logy=False)

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.035)
    latex.DrawLatex(0.18, 0.84,
                    f"Survived to filter: {n_surv}/{n_gen}")
    canvas._h = h
    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    print("    [loss_stage_distribution] done")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sc_sig, sc_bkg=None):
    print("  [had/summary] Yield flow...")
    plot_yield_flow(tdir, gf, sc_sig, sc_bkg)
    if gf is not None and len(gf) > 0:
        print("  [had/summary] Efficiency funnel...")
        plot_efficiency_funnel(tdir, gf)
        print("  [had/summary] Efficiency vs gen dxy...")
        plot_efficiency_vs_var(tdir, gf, "dxy", "Gen dxy (cm)", GEN_DXY_BINS, "eff_vs_dxy")
        print("  [had/summary] Efficiency vs gen pt...")
        plot_efficiency_vs_var(tdir, gf, "pt",  "Gen p_{T} (GeV)", GEN_PT_BINS, "eff_vs_pt")
        print("  [had/summary] Loss stage distribution...")
        plot_loss_stage_distribution(tdir, gf)
