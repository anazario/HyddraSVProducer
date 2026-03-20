"""
stages/hadronic_merging.py — merging-stage plots for hadronic HYDDRA.

Plots:
  - reco_merged_{cosTheta,decayAngle,pOverE,dxySignif,mass,nTracks}
  - hadronic_track_absorption  (nTracks distribution at merge: tight vs loose signal)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import HADRONIC_RECO_OBSERVABLES, COLOR_GOLD, COLOR_SILVER
from ..src.style   import draw_cms_label, make_canvas, draw_axis_grid
from ..src.plotter import plot_reco_observable


def plot_hadronic_track_absorption(tdir, gf):
    """
    Distribution of nTracks at the merge stage for tight vs loose signal vertices.
    Shows whether signal vertices accumulate extra tracks after merging.
    Uses matchRatio thresholds: tight (>=0.5) and loose (>0, <0.5).
    """
    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    ntrk_merge = ak.to_numpy(ak.flatten(gf["GenFunnel_nTracks_merged"]))[has_tracks]
    mr_merge   = ak.to_numpy(ak.flatten(gf["GenFunnel_matchRatio_merged"]))[has_tracks]

    valid = ntrk_merge > 0  # -1 is sentinel for no match
    tight = valid & (mr_merge >= 0.5)
    loose = valid & (mr_merge > 0) & (mr_merge < 0.5)

    if int(np.sum(tight)) == 0 and int(np.sum(loose)) == 0:
        print("    [hadronic_track_absorption] No merged signal vertices — skipping")
        return

    n_bins = 20
    bins   = np.arange(1, n_bins + 2, 1)

    canvas = make_canvas("hadronic_track_absorption",
                         "Hadronic HYDDRA — nTracks at merge (signal)", logy=False)
    canvas.SetGridx(0)

    h_tight = ROOT.TH1F("h_htab_tight", "", n_bins, bins)
    h_loose = ROOT.TH1F("h_htab_loose", "", n_bins, bins)

    for v, t, l in zip(ntrk_merge, tight, loose):
        if v < 0:
            continue
        if t:
            h_tight.Fill(float(v))
        elif l:
            h_loose.Fill(float(v))

    def norm_and_style(h, color):
        h.SetLineColor(color); h.SetLineWidth(2)
        h.SetFillStyle(0); h.SetStats(0)
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0 / integral)

    norm_and_style(h_tight, COLOR_GOLD)
    norm_and_style(h_loose, COLOR_SILVER)

    max_y = max(h_tight.GetMaximum(), h_loose.GetMaximum()) * 1.3

    h_ax = ROOT.TH1F("h_ax_htab",
                     ";Number of tracks at merge;Normalised to Unit Area",
                     n_bins, bins)
    h_ax.SetMinimum(1e-4); h_ax.SetMaximum(max(max_y, 0.1))
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.04);  h_ax.GetYaxis().SetLabelSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.3); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")
    for h in (h_tight, h_loose):
        if h.Integral() > 0:
            h.Draw("HIST SAME")

    canvas._grid_lines = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.62, 0.75, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.03)
    leg.AddEntry(h_tight, "Signal (matchRatio #geq 0.5)", "l")
    leg.AddEntry(h_loose, "Signal (0 < matchRatio < 0.5)", "l")
    leg.Draw()

    draw_cms_label("Hadronic HYDDRA")
    canvas.Update()

    canvas._h_ax = h_ax; canvas._h_tight = h_tight; canvas._h_loose = h_loose
    canvas._leg  = leg
    tdir.cd()
    canvas.Write()
    print("    [hadronic_track_absorption] done")


def make_plots(tdir, gf, sv_sig, sv_bkg):
    print("  [had/merging] Reco observables...")
    for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "merged", obs_key, obs_cfg, sv_bkg, hadronic=True)
        print(f"    [reco_merged_{obs_key}] done")
    if gf is not None and len(gf) > 0:
        print("  [had/merging] Track absorption...")
        plot_hadronic_track_absorption(tdir, gf)
