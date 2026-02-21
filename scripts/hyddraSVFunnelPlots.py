"""
hyddraSVFunnelPlots.py

Reads the ROOT file produced by HyddraSVsDiagnosticAnalyzer and generates:

  1. Yield flow — unnormalized counts at each stage (signal gold + background total)
  2. Yield flow — gold-normalized (survival fraction relative to seeding)
  3. Efficiency funnel — fraction of gen vertices found (gold / silver / bronze) per stage
  4. Efficiency vs gen dxy — gold at each stage
  5. Efficiency vs gen pt  — gold at each stage
  6. Silver-to-gold recovery — among silver-at-merge, fraction recovered to gold after cleaning
  7. 2D cleaning variables — (compatibility, cosTheta) for signal vs background tracks

Usage:
    python3 hyddraSVFunnelPlots.py --signal signal.root --output funnel_plots.root
    python3 hyddraSVFunnelPlots.py --signal signal.root --background bkg.root --output funnel_plots.root
"""

import argparse
import sys
import numpy as np
import uproot
import awkward as ak
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)


STAGE_NAMES   = ["Seeding", "Merging", "Cleaning", "Disambiguation", "Filtering"]
STAGE_KEYS    = ["seed", "merged", "cleaned", "disambig", "filtered"]
COLORS_SIGNAL = ROOT.kBlue + 2
COLORS_BKG    = ROOT.kRed + 2
COLORS_STAGE  = [ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kRed+2, ROOT.kOrange-3, ROOT.kMagenta-2]
MARKERS       = [20, 21, 22, 23, 29]


def draw_cms_label(subtitle="Leptonic HYDDRA"):
    """Draw CMS label matching StandardPlots/src/style.py (1D plotter style)."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextFont(61)
    latex.SetTextSize(0.057)      # 0.052 * 1.1 (1D multiplier)
    latex.DrawLatex(0.16, 0.955, "CMS")
    latex.SetTextFont(52)         # italic, matching StandardPlots
    latex.SetTextSize(0.044)      # 0.04 * 1.1
    latex.DrawLatex(0.26, 0.955, subtitle)


def make_canvas(name, title, logy=True, width=800, height=600):
    c = ROOT.TCanvas(name, title, width, height)
    c.SetLeftMargin(0.16)         # margin_left(0.12) + 0.04, as 1D plotter
    c.SetRightMargin(0.10)        # margin_right(0.12) - 0.02
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.10)
    if logy:
        c.SetLogy(1)
    c.SetGridx(1)
    c.SetGridy(1)
    c.SetTicks(1, 1)
    return c


def draw_grid_lines(canvas, n_stages, min_y, max_y):
    canvas._vlines, canvas._hlines = [], []
    for i in range(1, n_stages):
        ln = ROOT.TLine(i + 0.5, min_y, i + 0.5, max_y)
        ln.SetLineStyle(3); ln.SetLineColor(ROOT.kGray + 2); ln.Draw()
        canvas._vlines.append(ln)
    import math
    lo = int(math.floor(math.log10(max(min_y, 1e-3))))
    hi = int(math.ceil(math.log10(max(max_y, 1))))
    for exp in range(lo, hi):
        y = 10**exp
        if min_y <= y <= max_y:
            ln = ROOT.TLine(0.5, y, n_stages + 0.5, y)
            ln.SetLineStyle(3); ln.SetLineColor(ROOT.kGray + 2); ln.Draw()
            canvas._hlines.append(ln)


def setup_axis_hist(canvas, stage_names, y_title, min_y, max_y, title=""):
    n = len(stage_names)
    h = ROOT.TH1F(f"h_ax_{canvas.GetName()}", title, n, 0.5, n + 0.5)
    h.SetMinimum(min_y)
    h.SetMaximum(max_y)
    h.GetXaxis().SetTitle("Algorithm Stage")
    h.GetYaxis().SetTitle(y_title)
    h.GetXaxis().CenterTitle(True)
    h.GetYaxis().CenterTitle(True)
    h.GetXaxis().SetTitleSize(0.045)
    h.GetYaxis().SetTitleSize(0.045)
    h.GetXaxis().SetLabelSize(0.04)
    h.GetYaxis().SetLabelSize(0.04)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.SetStats(0)
    for i, s in enumerate(stage_names):
        h.GetXaxis().SetBinLabel(i + 1, s)
    canvas._axis_hist = h
    return h


def make_tgraph(x, y, color, marker, label, line_style=1):
    g = ROOT.TGraph(len(x), np.array(x, dtype=float), np.array(y, dtype=float))
    g.SetTitle(label)
    g.SetMarkerColor(color); g.SetLineColor(color)
    g.SetMarkerStyle(marker); g.SetLineWidth(2); g.SetLineStyle(line_style)
    return g


def add_legend(canvas, graphs, x1=0.6, y1=0.7, x2=0.88, y2=0.88):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)
    for g in graphs:
        leg.AddEntry(g, g.GetTitle(), "lp")
    leg.Draw()
    canvas._legend = leg


# Data loading helpers
def load_gen_funnel(root_file):
    """Return genFunnel tree arrays from an open uproot file."""
    return root_file["hyddraSVsDiagAnalyzer/genFunnel"].arrays(library="ak")


def load_stage_counts(root_file):
    """Return stageCounts tree arrays (summed over all events)."""
    arr = root_file["hyddraSVsDiagAnalyzer/stageCounts"].arrays(library="ak")
    totals = {}
    for key in ["n", "nGold", "nSilver", "nBronze"]:
        for stage in STAGE_KEYS:
            branch = f"Stage_{key}_{stage}"
            totals[branch] = int(ak.sum(arr[branch]))
    return totals


def load_cleaning_tracks(root_file):
    """Return cleaningTracks tree arrays from an open uproot file."""
    return root_file["hyddraSVsDiagAnalyzer/cleaningTracks"].arrays(library="ak")


# Plot 1 & 2: Yield flow (styled after vertexYieldFlow4.py)
def plot_yield_flow(out_file, sig_counts, bkg_counts=None):
    """
    Two canvases:
      - yield_flow_unnorm : raw counts at each stage
      - yield_flow_norm   : counts normalized to seeding (survival fraction)

    Signal curve = gold-matched vertices at each stage.
    Background curve = total vertices from bkg_counts if provided, otherwise
    total non-gold vertices from the signal file itself (useful for studying
    fake-rate when running on a signal-only sample).
    """
    for normalized in [False, True]:
        suffix = "norm" if normalized else "unnorm"
        canvas = make_canvas(f"yield_flow_{suffix}",
                             f"HYDDRA Leptonic Yield Flow ({'Normalized' if normalized else 'Raw'})")

        # Signal: gold-matched vertices
        sig_gold = [sig_counts[f"Stage_nGold_{s}"] for s in STAGE_KEYS]
        norm_sig = sig_gold[0] if (normalized and sig_gold[0] > 0) else 1

        # Background: dedicated file if available, otherwise non-gold SVs from signal file
        if bkg_counts:
            bkg_total  = [bkg_counts[f"Stage_n_{s}"] for s in STAGE_KEYS]
            bkg_label  = "Background (total)"
            bkg_source = "dedicated background file"
        else:
            # Non-gold = total - (gold + silver + bronze).  Negative values are
            # clipped to 0 (can happen due to double-counting in the match flags).
            bkg_total = [
                max(0, sig_counts[f"Stage_n_{s}"]
                    - sig_counts[f"Stage_nGold_{s}"]
                    - sig_counts[f"Stage_nSilver_{s}"]
                    - sig_counts[f"Stage_nBronze_{s}"])
                for s in STAGE_KEYS
            ]
            bkg_label  = "Non-signal SVs (signal file)"
            bkg_source = "signal file non-matched vertices"

        norm_bkg = bkg_total[0] if (normalized and bkg_total[0] > 0) else 1
        print(f"    background source: {bkg_source}")

        if normalized:
            sig_y   = [v / norm_sig for v in sig_gold]
            bkg_y   = [v / norm_bkg for v in bkg_total]
            y_title = "Survival Fraction (relative to Seeding)"
            min_y, max_y = 1e-4, 2.0
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

        x = np.arange(1, len(STAGE_NAMES) + 1, dtype=float)

        g_sig = make_tgraph(x, sig_y, COLORS_SIGNAL, 20, "Signal (gold)", line_style=1)
        g_bkg = make_tgraph(x, bkg_y, COLORS_BKG,    21, bkg_label,       line_style=2)
        graphs = [g_sig, g_bkg]

        mg = ROOT.TMultiGraph()
        for g in graphs:
            mg.Add(g)
        mg.Draw("lp")
        canvas._mg = mg
        canvas._graphs = graphs

        add_legend(canvas, graphs)
        draw_grid_lines(canvas, len(STAGE_NAMES), min_y, max_y)
        draw_cms_label()
        canvas.Update()
        out_file.cd()
        canvas.Write()

        print(f"  [yield_flow_{suffix}] done")


# Plot 3: Efficiency funnel bar chart
def plot_efficiency_funnel(out_file, gf):
    """
    Fraction of gen vertices with hasTracks=True found at each stage
    broken down by gold / silver / bronze.
    One grouped bar per stage.
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    n_denom = int(np.sum(has_tracks))
    if n_denom == 0:
        print("  [efficiency_funnel] No gen vertices with tracks — skipping")
        return

    canvas = make_canvas("efficiency_funnel", "Signal Efficiency per Stage", logy=False)
    canvas.SetGridx(0)

    n_stages = len(STAGE_KEYS)
    # Stacked bar structure: gold on top, then silver (exclusive), then bronze (exclusive)
    eff_gold   = np.zeros(n_stages)
    eff_silver = np.zeros(n_stages)  # silver but not gold
    eff_bronze = np.zeros(n_stages)  # bronze but not gold or silver

    for i, s in enumerate(STAGE_KEYS):
        gold   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_gold_{s}"]))[has_tracks].astype(bool)
        silver = ak.to_numpy(ak.flatten(gf[f"GenFunnel_silver_{s}"]))[has_tracks].astype(bool)
        bronze = ak.to_numpy(ak.flatten(gf[f"GenFunnel_bronze_{s}"]))[has_tracks].astype(bool)
        eff_gold[i]   = np.sum(gold) / n_denom
        eff_silver[i] = np.sum(silver & ~gold) / n_denom
        eff_bronze[i] = np.sum(bronze & ~silver & ~gold) / n_denom

    h_gold   = ROOT.TH1F("h_fgl", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)
    h_silver = ROOT.TH1F("h_fsi", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)
    h_bronze = ROOT.TH1F("h_fbr", ";Algorithm Stage;Efficiency", n_stages, 0.5, n_stages + 0.5)

    for i, name in enumerate(STAGE_NAMES):
        h_gold  .GetXaxis().SetBinLabel(i + 1, name)
        h_silver.GetXaxis().SetBinLabel(i + 1, name)
        h_bronze.GetXaxis().SetBinLabel(i + 1, name)
        h_gold  .SetBinContent(i + 1, eff_gold[i])
        h_silver.SetBinContent(i + 1, eff_silver[i])
        h_bronze.SetBinContent(i + 1, eff_bronze[i])

    h_gold  .SetFillColor(ROOT.kBlue + 2);  h_gold  .SetLineColor(ROOT.kBlue + 2)
    h_silver.SetFillColor(ROOT.kAzure + 6); h_silver.SetLineColor(ROOT.kAzure + 6)
    h_bronze.SetFillColor(ROOT.kCyan - 7);  h_bronze.SetLineColor(ROOT.kCyan - 7)

    max_y = max(eff_gold + eff_silver + eff_bronze) * 1.3
    h_gold.SetMaximum(max_y)
    h_gold.GetYaxis().SetTitle("Efficiency")
    h_gold.GetXaxis().CenterTitle(True)
    h_gold.GetYaxis().CenterTitle(True)
    h_gold.SetStats(0)
    h_gold.Draw("HIST")

    # Draw stacked (bronze baseline, silver on top, gold on top)
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
    out_file.cd()
    canvas.Write()
    print("  [efficiency_funnel] done")


# Plots 4 & 5: Efficiency vs gen dxy / pt
def plot_efficiency_vs_var(out_file, gf, var_key, var_label, bins, canvas_name):
    """
    TEfficiency-style ratio curve for gold efficiency vs a gen-level variable,
    one curve per stage, drawn on the same canvas.
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    var_vals   = ak.to_numpy(ak.flatten(gf[f"GenFunnel_{var_key}"]))[has_tracks]
    n_bins     = len(bins) - 1

    canvas = make_canvas(f"eff_vs_{var_key}", f"Gold efficiency vs {var_label}", logy=False)
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

    # Dummy axis
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
    out_file.cd()
    canvas.Write()
    print(f"  [eff_vs_{var_key}] done")


# Plot 6: Silver-to-gold recovery
def plot_silver_to_gold_recovery(out_file, gf):
    """
    Among gen vertices that are silver-but-not-gold at the merge stage,
    show what fraction become gold after cleaning — plotted vs gen dxy.
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    dxy        = ak.to_numpy(ak.flatten(gf["GenFunnel_dxy"]))[has_tracks]
    gold_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_gold_merged"]))[has_tracks].astype(bool)
    silv_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_silver_merged"]))[has_tracks].astype(bool)
    gold_clean = ak.to_numpy(ak.flatten(gf["GenFunnel_gold_cleaned"]))[has_tracks].astype(bool)

    # Absorbed = silver at merge but NOT gold at merge
    absorbed = silv_merge & ~gold_merge
    if int(np.sum(absorbed)) == 0:
        print("  [silver_to_gold_recovery] No absorbed vertices found — skipping")
        return

    recovered = absorbed & gold_clean

    bins = np.linspace(0, 30, 31)
    h_denom = ROOT.TH1F("h_sg_denom", "", len(bins) - 1, bins)
    h_numer = ROOT.TH1F("h_sg_numer", "", len(bins) - 1, bins)
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

    canvas._h_ax = h_ax; canvas._eff = eff; canvas._g = g; canvas._hists = (h_denom, h_numer)
    draw_cms_label()
    canvas.Update()
    out_file.cd()
    canvas.Write()
    print("  [silver_to_gold_recovery] done")


# Plot 7: 2D cleaning variables
def plot_cleaning_2d(out_file, ct, max_compat=1.5, min_cos_theta=0.5):
    """
    2D histogram of (compatibility, cosTheta) for tracks in post-merge multi-track
    vertices, separated by:
      - signal tracks in signal-matched vertices (gold/silver)
      - background tracks in signal-matched vertices ("contaminating")
      - all tracks in background-only vertices
    Overlaid with the cleaning cut lines.
    """
    if ct is None or len(ak.flatten(ct["CleanTrack_compatibility"])) == 0:
        print("  [cleaning_2d] No cleaning track data — skipping")
        return

    compat    = ak.to_numpy(ak.flatten(ct["CleanTrack_compatibility"])).astype(float)
    cos_th    = ak.to_numpy(ak.flatten(ct["CleanTrack_cosTheta"])).astype(float)
    is_signal = ak.to_numpy(ak.flatten(ct["CleanTrack_isSignal"])).astype(bool)
    vtx_gold  = ak.to_numpy(ak.flatten(ct["CleanTrack_vtxIsGold"])).astype(bool)
    vtx_silver= ak.to_numpy(ak.flatten(ct["CleanTrack_vtxIsSilver"])).astype(bool)

    sig_vtx = vtx_gold | vtx_silver
    sig_track_in_sig_vtx = is_signal & sig_vtx
    bkg_track_in_sig_vtx = ~is_signal & sig_vtx
    bkg_vtx_tracks       = ~sig_vtx

    x_bins = np.linspace(0, 5, 51)
    y_bins = np.linspace(-1, 1, 51)

    def fill_h2(name, title, mask):
        ROOT.gStyle.SetPalette(ROOT.kViridis)
        h = ROOT.TH2F(name, title,
                      len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)
        h.SetDirectory(0)
        for c, ct_v in zip(compat[mask], cos_th[mask]):
            h.Fill(c, ct_v)
        h.GetXaxis().SetTitle("Track Compatibility (#sigma)")
        h.GetYaxis().SetTitle("Track cos#theta")
        h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.3); h.GetYaxis().SetTitleOffset(1.3)
        h.SetStats(0)
        return h

    h_sig   = fill_h2("h2_sig_track",   "Signal tracks in signal vtx; Compatibility; cos#theta",
                      sig_track_in_sig_vtx)
    h_cont  = fill_h2("h2_cont_track",  "Bkg tracks in signal vtx; Compatibility; cos#theta",
                      bkg_track_in_sig_vtx)
    h_bkg   = fill_h2("h2_bkg_track",   "All tracks in bkg vtx; Compatibility; cos#theta",
                      bkg_vtx_tracks)

    def draw_cut_lines(h):
        """Overlay the cleaning cut boundaries."""
        lines = []
        # Vertical line at maxCompatibility
        ln1 = ROOT.TLine(max_compat, y_bins[0], max_compat, y_bins[-1])
        ln1.SetLineColor(ROOT.kRed); ln1.SetLineWidth(2); ln1.SetLineStyle(2)
        ln1.Draw(); lines.append(ln1)
        # Horizontal line at minCleanCosTheta
        ln2 = ROOT.TLine(x_bins[0], min_cos_theta, x_bins[-1], min_cos_theta)
        ln2.SetLineColor(ROOT.kRed); ln2.SetLineWidth(2); ln2.SetLineStyle(2)
        ln2.Draw(); lines.append(ln2)
        return lines

    for h, suffix, title in [
        (h_sig,  "signal_tracks",       "Signal tracks in signal vertices"),
        (h_cont, "contaminating_tracks","Background tracks in signal vertices"),
        (h_bkg,  "bkg_tracks",          "All tracks in background vertices"),
    ]:
        c = ROOT.TCanvas(f"cleaning_2d_{suffix}", title, 800, 600)
        c.SetLeftMargin(0.16)
        c.SetRightMargin(0.18)    # extra room for z-axis colour bar
        c.SetBottomMargin(0.12)
        c.SetTopMargin(0.10)
        #h.Draw("COLZ")
        c.Update()
        lines = draw_cut_lines(h)

        leg = ROOT.TLegend(0.18, 0.75, 0.55, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0)
        # Dummy lines for legend
        l1 = ROOT.TLine(); l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.SetLineWidth(2)
        leg.AddEntry(l1, "Cleaning cuts", "l")
        leg.Draw()
        draw_cms_label()
        c.Update()

        c._h = h; c._lines = lines; c._leg = leg
        out_file.cd()
        c.Write()
        print(f"  [cleaning_2d_{suffix}] done")


# Main

def main():
    parser = argparse.ArgumentParser(description="HYDDRA diagnostic funnel plots")
    parser.add_argument("--signal",     required=True, help="Signal ROOT file")
    parser.add_argument("--background", default=None,  help="Background ROOT file (optional)")
    parser.add_argument("--output",     default="hyddraSVFunnelPlots.root",
                        help="Output ROOT file")
    parser.add_argument("--max-compat",   type=float, default=1.5,
                        help="maxCompatibility cut value (default: 1.5)")
    parser.add_argument("--min-cos-theta",type=float, default=0.5,
                        help="minCleanCosTheta cut value (default: 0.5)")
    args = parser.parse_args()

    print(f"[hyddraSVFunnelPlots] Signal:     {args.signal}")
    print(f"[hyddraSVFunnelPlots] Background: {args.background or 'none'}")
    print(f"[hyddraSVFunnelPlots] Output:     {args.output}")

    with uproot.open(args.signal) as sig_file:
        gf_sig  = load_gen_funnel(sig_file)
        sc_sig  = load_stage_counts(sig_file)
        ct_sig  = load_cleaning_tracks(sig_file)

    sc_bkg = None
    if args.background:
        with uproot.open(args.background) as bkg_file:
            sc_bkg = load_stage_counts(bkg_file)

    out_file = ROOT.TFile(args.output, "RECREATE")

    print("\n[1/7] Yield flow plots...")
    plot_yield_flow(out_file, sc_sig, sc_bkg)

    print("[2/7] Efficiency funnel bar chart...")
    plot_efficiency_funnel(out_file, gf_sig)

    print("[3/7] Efficiency vs gen dxy...")
    dxy_bins = list(np.concatenate([
        np.linspace(0, 5, 11), np.linspace(6, 20, 8), np.linspace(25, 50, 6)]))
    plot_efficiency_vs_var(out_file, gf_sig, "dxy", "Gen dxy (cm)", dxy_bins, "eff_vs_dxy")

    print("[4/7] Efficiency vs gen pt...")
    pt_bins = list(np.linspace(0, 100, 21))
    plot_efficiency_vs_var(out_file, gf_sig, "pt", "Gen p_{T} (GeV)", pt_bins, "eff_vs_pt")

    print("[5/7] Silver-to-gold recovery...")
    plot_silver_to_gold_recovery(out_file, gf_sig)

    print("[6/7] 2D cleaning variable plots...")
    plot_cleaning_2d(out_file, ct_sig, args.max_compat, args.min_cos_theta)

    out_file.Close()
    print(f"\n[hyddraSVFunnelPlots] All plots saved to {args.output}")


if __name__ == "__main__":
    main()
