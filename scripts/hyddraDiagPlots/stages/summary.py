"""
stages/summary.py — cross-stage summary plots for hyddraDiagPlots.

Plots:
  - yield_flow_unnorm / yield_flow_norm
  - efficiency_funnel
  - eff_vs_dxy / eff_vs_pt
  - silver_to_gold_recovery
  - loss_stage_distribution
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config import (
    STAGE_NAMES, STAGE_KEYS,
    COLOR_BKG, COLORS_STAGE, MARKERS,
    GEN_DXY_BINS, GEN_PT_BINS,
)
from ..src.style  import draw_cms_label, make_canvas, draw_grid_lines
from ..src.plotter import setup_axis_hist, make_tgraph, add_legend

_COLOR_SIG = ROOT.kBlue + 2


# ── Yield flow ────────────────────────────────────────────────────────────────

def plot_yield_flow(tdir, sc_sig, sc_bkg=None):
    for normalized in [False, True]:
        suffix = "norm" if normalized else "unnorm"
        canvas = make_canvas(
            f"yield_flow_{suffix}",
            f"HYDDRA Leptonic Yield Flow ({'Normalized' if normalized else 'Raw'})",
            logy=(not normalized),
        )

        sig_gold = [sc_sig[f"Stage_nGold_{s}"] for s in STAGE_KEYS]
        norm_sig = sig_gold[0] if (normalized and sig_gold[0] > 0) else 1

        if sc_bkg:
            bkg_total = [sc_bkg[f"Stage_n_{s}"] for s in STAGE_KEYS]
            bkg_label = "Background (total)"
        else:
            bkg_total = [
                max(0, sc_sig[f"Stage_n_{s}"]
                    - sc_sig[f"Stage_nGold_{s}"]
                    - sc_sig[f"Stage_nSilver_{s}"]
                    - sc_sig[f"Stage_nBronze_{s}"])
                for s in STAGE_KEYS
            ]
            bkg_label = "Non-signal SVs (signal file)"

        norm_bkg = bkg_total[0] if (normalized and bkg_total[0] > 0) else 1

        if normalized:
            sig_y   = [v / norm_sig for v in sig_gold]
            bkg_y   = [v / norm_bkg for v in bkg_total]
            y_title = "Survival Fraction (relative to Seeding)"
            min_y, max_y = 1e-6, 1.2
        else:
            sig_y   = [float(v) for v in sig_gold]
            bkg_y   = [float(v) for v in bkg_total]
            y_title = "Vertex Yield"
            all_vals = sig_y + bkg_y
            nonzero  = [v for v in all_vals if v > 0]
            min_y    = max(0.5, min(nonzero) / 2) if nonzero else 0.5
            max_y    = max(nonzero) * 5 if nonzero else 10

        h = setup_axis_hist(canvas, STAGE_NAMES, y_title, min_y, max_y)
        h.Draw("AXIS")

        x     = np.arange(1, len(STAGE_NAMES) + 1, dtype=float)
        g_sig = make_tgraph(x, sig_y, _COLOR_SIG, 20, "Signal (gold)",  1)
        g_bkg = make_tgraph(x, bkg_y, COLOR_BKG,  21, bkg_label,        2)
        graphs = [g_sig, g_bkg]

        mg = ROOT.TMultiGraph()
        for g in graphs:
            mg.Add(g)
        mg.Draw("lp")
        canvas._mg = mg
        canvas._graphs = graphs

        add_legend(canvas, graphs, x1=0.54, y1=0.78, x2=0.82, y2=0.90)
        draw_grid_lines(canvas, len(STAGE_NAMES), min_y, max_y, logy=not normalized)
        draw_cms_label()
        canvas.Update()
        tdir.cd()
        canvas.Write()
        print(f"    [yield_flow_{suffix}] done")


# ── Efficiency funnel ─────────────────────────────────────────────────────────

def plot_efficiency_funnel(tdir, gf):
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    n_denom    = int(np.sum(has_tracks))
    if n_denom == 0:
        print("    [efficiency_funnel] No gen vertices with tracks — skipping")
        return

    canvas = make_canvas("efficiency_funnel", "Signal Efficiency per Stage", logy=False)
    canvas.SetGridx(0)

    n_stages   = len(STAGE_KEYS)
    eff_gold   = np.zeros(n_stages)
    eff_silver = np.zeros(n_stages)
    eff_bronze = np.zeros(n_stages)

    for i, s in enumerate(STAGE_KEYS):
        gold   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_gold_{s}"])  )[has_tracks].astype(bool)
        silver = ak.to_numpy(ak.flatten(gf[f"GenFunnel_silver_{s}"]))[has_tracks].astype(bool)
        bronze = ak.to_numpy(ak.flatten(gf[f"GenFunnel_bronze_{s}"]))[has_tracks].astype(bool)
        eff_gold[i]   = np.sum(gold)             / n_denom
        eff_silver[i] = np.sum(silver & ~gold)   / n_denom
        eff_bronze[i] = np.sum(bronze & ~silver & ~gold) / n_denom

    h_gold   = ROOT.TH1F("h_fgl", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)
    h_silver = ROOT.TH1F("h_fsi", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)
    h_bronze = ROOT.TH1F("h_fbr", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)

    for i, name in enumerate(STAGE_NAMES):
        for h in (h_gold, h_silver, h_bronze):
            h.GetXaxis().SetBinLabel(i + 1, name)
        h_gold  .SetBinContent(i + 1, eff_gold[i])
        h_silver.SetBinContent(i + 1, eff_silver[i])
        h_bronze.SetBinContent(i + 1, eff_bronze[i])

    h_gold  .SetFillColor(ROOT.kBlue + 2); h_gold  .SetLineColor(ROOT.kBlue + 2)
    h_silver.SetFillColor(ROOT.kAzure + 6);h_silver.SetLineColor(ROOT.kAzure + 6)
    h_bronze.SetFillColor(ROOT.kCyan - 7); h_bronze.SetLineColor(ROOT.kCyan - 7)

    max_y = max(eff_gold + eff_silver + eff_bronze) * 1.3
    h_gold.SetMaximum(max_y)
    h_gold.GetYaxis().SetTitle("Efficiency")
    h_gold.GetXaxis().CenterTitle(True)
    h_gold.GetYaxis().CenterTitle(True)
    h_gold.SetStats(0)
    h_gold.Draw("HIST")

    hs = ROOT.THStack("hs_funnel", "")
    hs.Add(h_bronze)
    hs.Add(h_silver)
    hs.Add(h_gold)
    hs.Draw("HIST SAME")

    leg = ROOT.TLegend(0.65, 0.72, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0)
    leg.AddEntry(h_gold,   "Gold (exact match)", "f")
    leg.AddEntry(h_silver, "Silver (partial, not gold)", "f")
    leg.AddEntry(h_bronze, "Bronze (spatial only)", "f")
    leg.Draw()
    canvas._leg = leg; canvas._hs = hs
    canvas._hgold = h_gold; canvas._hsilver = h_silver; canvas._hbronze = h_bronze

    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    print("    [efficiency_funnel] done")


# ── Efficiency vs gen variable ─────────────────────────────────────────────────

def plot_efficiency_vs_var(tdir, gf, var_key, var_label, bins, canvas_name):
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    var_vals   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{var_key}"]))[has_tracks]
    n_bins     = len(bins) - 1

    canvas = make_canvas(canvas_name, f"Gold efficiency vs {var_label}", logy=False)
    canvas.SetGridx(0)

    h_denom = ROOT.TH1F(f"h_denom_{var_key}", "", n_bins, np.array(bins, dtype=float))
    for v in var_vals:
        h_denom.Fill(v)

    graphs, hists = [], []
    for i, s in enumerate(STAGE_KEYS):
        gold_mask = ak.to_numpy(ak.flatten(gf[f"GenFunnel_gold_{s}"]))[has_tracks].astype(bool)
        h_num = ROOT.TH1F(f"h_num_{s}_{var_key}", "", n_bins, np.array(bins, dtype=float))
        for v in var_vals[gold_mask]:
            h_num.Fill(v)
        eff = ROOT.TEfficiency(h_num, h_denom)
        eff.SetStatisticOption(ROOT.TEfficiency.kFCP)
        g = eff.CreateGraph()
        g.SetTitle(STAGE_NAMES[i])
        g.SetMarkerColor(COLORS_STAGE[i]); g.SetLineColor(COLORS_STAGE[i])
        g.SetMarkerStyle(MARKERS[i]); g.SetLineWidth(2)
        graphs.append(g)
        hists.append((h_num, eff))

    h_ax = ROOT.TH1F(f"h_ax_{var_key}", f";{var_label};Gold Efficiency",
                     n_bins, np.array(bins, dtype=float))
    h_ax.SetMinimum(0); h_ax.SetMaximum(1.1)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")
    canvas._h_ax = h_ax; canvas._hists = hists; canvas._graphs = graphs

    for g in graphs:
        g.Draw("P SAME")

    add_legend(canvas, graphs, x1=0.6, y1=0.55, x2=0.88, y2=0.88)
    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    print(f"    [eff_vs_{var_key}] done")


# ── Silver-to-gold recovery ───────────────────────────────────────────────────

def plot_silver_to_gold_recovery(tdir, gf):
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    dxy        = ak.to_numpy(ak.flatten(gf["GenFunnel_dxy"])          )[has_tracks]
    gold_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_gold_merged"])  )[has_tracks].astype(bool)
    silv_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_silver_merged"]))[has_tracks].astype(bool)
    gold_clean = ak.to_numpy(ak.flatten(gf["GenFunnel_gold_cleaned"]) )[has_tracks].astype(bool)

    absorbed = silv_merge & ~gold_merge
    if int(np.sum(absorbed)) == 0:
        print("    [silver_to_gold_recovery] No absorbed vertices found — skipping")
        return

    recovered = absorbed & gold_clean
    bins      = np.linspace(0, 30, 31)
    h_denom   = ROOT.TH1F("h_sg_denom", "", len(bins) - 1, bins)
    h_numer   = ROOT.TH1F("h_sg_numer", "", len(bins) - 1, bins)
    for v, a, r in zip(dxy, absorbed, recovered):
        if a:
            h_denom.Fill(v)
            if r:
                h_numer.Fill(v)

    canvas = make_canvas("silver_to_gold_recovery",
                         "Silver-to-Gold Recovery after Cleaning", logy=False)
    canvas.SetGridx(0)

    eff = ROOT.TEfficiency(h_numer, h_denom)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)
    g = eff.CreateGraph()
    g.SetTitle("Silver#rightarrowGold recovery")
    g.SetMarkerColor(ROOT.kGreen + 2); g.SetLineColor(ROOT.kGreen + 2)
    g.SetMarkerStyle(20); g.SetLineWidth(2)

    h_ax = ROOT.TH1F("h_sg_ax", ";Gen dxy (cm);Recovery Fraction", len(bins) - 1, bins)
    h_ax.SetMinimum(0); h_ax.SetMaximum(1.1)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")
    g.Draw("P SAME")

    canvas._h_ax = h_ax; canvas._eff = eff; canvas._g = g
    canvas._hists = (h_denom, h_numer)
    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    print("    [silver_to_gold_recovery] done")


# ── Loss stage distribution ───────────────────────────────────────────────────

def plot_loss_stage_distribution(tdir, gf):
    """Bar chart: how many gen-gold vertices are lost at each stage."""
    has_tracks  = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    lost_stages = ak.to_numpy(ak.flatten(gf["GenFunnel_lostGoldAtStage"]))[has_tracks]

    # -1 = survived, 0 = never found, 1-4 = lost after that stage
    labels  = ["Never seeded", "Lost at Seeding", "Lost at Merging",
               "Lost at Cleaning", "Lost at Disambig."]
    n_cats  = len(labels)
    counts  = np.array([np.sum(lost_stages == v) for v in range(0, n_cats)], dtype=float)
    n_surv  = int(np.sum(lost_stages == -1))
    n_total = len(lost_stages)

    canvas = make_canvas("loss_stage_distribution",
                         "Gold vertex loss per stage", logy=False)
    canvas.SetGridx(0)

    h = ROOT.TH1F("h_loss", ";Loss Category;Count", n_cats, 0.5, n_cats + 0.5)
    for i, (c, lbl) in enumerate(zip(counts, labels)):
        h.SetBinContent(i + 1, c)
        h.GetXaxis().SetBinLabel(i + 1, lbl)
    h.SetFillColor(ROOT.kOrange + 1)
    h.SetLineColor(ROOT.kOrange + 1)
    h.GetYaxis().SetTitle("Gen vertex count")
    h.GetXaxis().CenterTitle(True)
    h.GetYaxis().CenterTitle(True)
    h.SetStats(0)
    h.SetMaximum(max(counts) * 1.4 if max(counts) > 0 else 1)
    h.Draw("HIST")

    # Annotation: survived count
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.035)
    latex.DrawLatex(0.18, 0.84,
                    f"Survived to filter: {n_surv}/{n_total}")
    canvas._h = h
    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    print("    [loss_stage_distribution] done")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sc_sig, sc_bkg=None):
    print("  [summary] Yield flow...")
    plot_yield_flow(tdir, sc_sig, sc_bkg)
    if gf is not None and len(gf) > 0:
        print("  [summary] Efficiency funnel...")
        plot_efficiency_funnel(tdir, gf)
        print("  [summary] Efficiency vs gen dxy...")
        plot_efficiency_vs_var(tdir, gf, "dxy", "Gen dxy (cm)", GEN_DXY_BINS, "eff_vs_dxy")
        print("  [summary] Efficiency vs gen pt...")
        plot_efficiency_vs_var(tdir, gf, "pt",  "Gen p_{T} (GeV)", GEN_PT_BINS, "eff_vs_pt")
        print("  [summary] Silver-to-gold recovery...")
        plot_silver_to_gold_recovery(tdir, gf)
        print("  [summary] Loss stage distribution...")
        plot_loss_stage_distribution(tdir, gf)
