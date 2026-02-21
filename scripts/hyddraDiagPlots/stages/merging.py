"""
stages/merging.py — merging-stage plots for hyddraDiagPlots.

Plots:
  - reco_merged_{cosTheta,decayAngle,pOverE,dxySignif}
  - track_absorption  (nTracks distribution at merge: gold vs silver)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import RECO_OBSERVABLES, COLOR_GOLD, COLOR_SILVER
from ..src.style   import draw_cms_label, make_canvas
from ..src.plotter import plot_reco_observable


def plot_track_absorption(tdir, gf):
    """
    Distribution of nTracks at the merge stage for gold vs silver vertices.
    Shows whether signal vertices accumulate extra tracks after merging.
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    ntrk_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_nTracks_merged"]))[has_tracks]
    gold_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_gold_merged"])   )[has_tracks].astype(bool)
    silv_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_silver_merged"]) )[has_tracks].astype(bool)

    # Only consider matched vertices (not -1 sentinels)
    matched = ntrk_merge > 0

    if int(np.sum(matched & gold_merge)) == 0 and int(np.sum(matched & silv_merge & ~gold_merge)) == 0:
        print("    [track_absorption] No merged vertices found — skipping")
        return

    n_bins = 10
    bins   = np.arange(1.5, n_bins + 2.5, 1.0)

    canvas = make_canvas("track_absorption",
                         "nTracks at merge: gold vs silver", logy=False)
    canvas.SetGridx(0)

    h_gold   = ROOT.TH1F("h_tab_gold",   "", n_bins, bins)
    h_silver = ROOT.TH1F("h_tab_silver", "", n_bins, bins)

    for v, g, s in zip(ntrk_merge, gold_merge, silv_merge):
        if v < 0:
            continue
        if g:
            h_gold.Fill(float(v))
        elif s:
            h_silver.Fill(float(v))

    def norm_and_style(h, color):
        h.SetLineColor(color); h.SetLineWidth(2)
        h.SetFillStyle(0); h.SetStats(0)
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0 / integral)

    norm_and_style(h_gold,   COLOR_GOLD)
    norm_and_style(h_silver, COLOR_SILVER)

    max_y = max(h_gold.GetMaximum(), h_silver.GetMaximum()) * 1.3

    h_ax = ROOT.TH1F("h_ax_tab",
                     ";Number of tracks at merge;Normalised to Unit Area",
                     n_bins, bins)
    h_ax.SetMinimum(0); h_ax.SetMaximum(max(max_y, 0.1))
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.04);  h_ax.GetYaxis().SetLabelSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.3); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")
    for h in (h_gold, h_silver):
        if h.Integral() > 0:
            h.Draw("HIST SAME")

    leg = ROOT.TLegend(0.62, 0.75, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.03)
    leg.AddEntry(h_gold,   "Gold (exact match)", "l")
    leg.AddEntry(h_silver, "Silver (partial)",   "l")
    leg.Draw()

    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    canvas._h_ax = h_ax; canvas._hists = [h_gold, h_silver]; canvas._leg = leg
    print("    [track_absorption] done")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg):
    print("  [merging] Reco observables...")
    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "merged", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_merged_{obs_key}] done")
    if gf is not None and len(gf) > 0:
        print("  [merging] Track absorption...")
        plot_track_absorption(tdir, gf)
