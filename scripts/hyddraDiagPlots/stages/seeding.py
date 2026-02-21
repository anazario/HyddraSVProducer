"""
stages/seeding.py — seeding-stage plots for hyddraDiagPlots.

Plots:
  - reco_seed_{cosTheta,decayAngle,pOverE,dxySignif}  (via plot_reco_observable)
  - seed_track_cos_theta
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import RECO_OBSERVABLES
from ..src.style   import draw_cms_label, make_canvas, draw_axis_grid
from ..src.plotter import plot_reco_observable


def plot_seed_track_cos_theta(tdir, st):
    """
    Track cosTheta distribution for tracks inside 2-track OS seeded vertices.
    Overlays signal tracks (SeedTrack_isSignal=True) vs non-signal tracks.
    Split by whether the vertex is gold-matched.
    """
    if st is None or len(st) == 0:
        print("    [seed_track_cos_theta] No seed track data — skipping")
        return

    cos_theta  = ak.to_numpy(ak.flatten(st["SeedTrack_cosTheta"]))
    is_signal  = ak.to_numpy(ak.flatten(st["SeedTrack_isSignal"])).astype(bool)
    vtx_gold   = ak.to_numpy(ak.flatten(st["SeedTrack_vtxIsGold"])).astype(bool)

    if len(cos_theta) == 0:
        print("    [seed_track_cos_theta] Empty after flatten — skipping")
        return

    bins   = np.linspace(-1, 1, 51)
    n_bins = len(bins) - 1

    canvas = make_canvas("seed_track_cos_theta",
                         "Track cos#theta in 2-track OS seeds", logy=False)

    h_sig_gold    = ROOT.TH1F("h_stct_sig_gold",    "", n_bins, bins)
    h_sig_nongold = ROOT.TH1F("h_stct_sig_nongold", "", n_bins, bins)
    h_ns_gold     = ROOT.TH1F("h_stct_ns_gold",     "", n_bins, bins)
    h_ns_nongold  = ROOT.TH1F("h_stct_ns_nongold",  "", n_bins, bins)

    for v, sig, gold in zip(cos_theta, is_signal, vtx_gold):
        if sig  and gold:  h_sig_gold.Fill(float(v))
        if sig  and not gold: h_sig_nongold.Fill(float(v))
        if not sig and gold:  h_ns_gold.Fill(float(v))
        if not sig and not gold: h_ns_nongold.Fill(float(v))

    def norm_and_style(h, color, lstyle):
        h.SetLineColor(color); h.SetLineWidth(2); h.SetLineStyle(lstyle)
        h.SetFillStyle(0); h.SetStats(0)
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0 / integral)

    norm_and_style(h_sig_gold,    ROOT.kBlue + 2,   1)
    norm_and_style(h_sig_nongold, ROOT.kBlue + 2,   2)
    norm_and_style(h_ns_gold,     ROOT.kRed  + 2,   1)
    norm_and_style(h_ns_nongold,  ROOT.kRed  + 2,   2)

    hists = [h_sig_gold, h_sig_nongold, h_ns_gold, h_ns_nongold]
    max_y = max((h.GetMaximum() for h in hists if h.Integral() > 0), default=1.0) * 1.3

    h_ax = ROOT.TH1F("h_ax_stct", ";Track cos#theta (wrt PV);Normalised to Unit Area",
                     n_bins, bins)
    h_ax.SetMinimum(0); h_ax.SetMaximum(max_y)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.04);  h_ax.GetYaxis().SetLabelSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.3); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")
    for h in hists:
        if h.Integral() > 0:
            h.Draw("HIST SAME")

    leg = ROOT.TLegend(0.18, 0.72, 0.60, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.03)
    leg.AddEntry(h_sig_gold,    "Signal track, gold vtx",     "l")
    leg.AddEntry(h_sig_nongold, "Signal track, non-gold vtx", "l")
    leg.AddEntry(h_ns_gold,     "Non-sig track, gold vtx",    "l")
    leg.AddEntry(h_ns_nongold,  "Non-sig track, non-gold vtx","l")
    leg.Draw()

    canvas._grid_lines = draw_axis_grid(h_ax, logy=False)
    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()
    canvas._h_ax = h_ax; canvas._hists = hists; canvas._leg = leg
    print("    [seed_track_cos_theta] done")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg, st):
    print("  [seeding] Reco observables...")
    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "seed", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_seed_{obs_key}] done")
    print("  [seeding] Seed track cosTheta...")
    plot_seed_track_cos_theta(tdir, st)
