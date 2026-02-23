"""
stages/filtering.py — filtering-stage plots for hyddraDiagPlots.

Plots:
  - reco_filtered_{cosTheta,decayAngle,pOverE,dxySignif}
  - filtering_costheta_decay_{signal,nonsignal}_{before,after}
  - filter_cutflow  (sequential removal efficiency per cut, signal vs non-signal)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import RECO_OBSERVABLES, STAGE_IDX, COLOR_GOLD, COLOR_NONSIGNAL
from ..src.style   import draw_cms_label, draw_colz_grid, draw_axis_grid
from ..src.plotter import plot_reco_observable


def plot_filtering_costheta_decay_2d(tdir, sv_sig):
    """
    2D (decayAngle, cosTheta) scatter before and after filtering,
    split by signal (isGold) and non-signal.  Four canvases total.
    """
    if sv_sig is None or len(sv_sig) == 0:
        print("    [filtering_2d] No allStageVtx data — skipping")
        return

    sidx_pre  = STAGE_IDX["disambig"]
    sidx_post = STAGE_IDX["filtered"]

    stage   = ak.to_numpy(sv_sig["StageVtx_stageIdx"])
    cos_th  = ak.to_numpy(sv_sig["StageVtx_cosTheta"]).astype(float)
    decay   = ak.to_numpy(sv_sig["StageVtx_decayAngle"]).astype(float)
    is_gold = ak.to_numpy(sv_sig["StageVtx_isGold"]).astype(bool)

    x_bins = np.linspace(-1, 1, 101)
    y_bins = np.linspace(-1, 1, 101)

    def fill_h2(name, mask):
        ROOT.gStyle.SetPalette(ROOT.kViridis)
        h = ROOT.TH2F(name, name, len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)
        h.SetDirectory(0)
        for da, ct in zip(decay[mask], cos_th[mask]):
            h.Fill(da, ct)
        h.GetXaxis().SetTitle("cos#theta^{*}_{CM}")
        h.GetYaxis().SetTitle("cos#theta (wrt PV)")
        h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.3); h.GetYaxis().SetTitleOffset(1.3)
        h.SetStats(0)
        return h

    configs = [
        ("signal_before",    (stage == sidx_pre)  &  is_gold,  "Signal vertices (pre-filter)"),
        ("signal_after",     (stage == sidx_post) &  is_gold,  "Signal vertices (post-filter)"),
        ("nonsignal_before", (stage == sidx_pre)  & ~is_gold,  "Non-signal vertices (pre-filter)"),
        ("nonsignal_after",  (stage == sidx_post) & ~is_gold,  "Non-signal vertices (post-filter)"),
    ]

    for suffix, mask, title in configs:
        cname = f"filtering_costheta_decay_{suffix}"
        h = fill_h2(f"h2_filt_{suffix}", mask)

        c = ROOT.TCanvas(cname, title, 800, 600)
        c.SetLeftMargin(0.16)
        c.SetRightMargin(0.18)
        c.SetBottomMargin(0.12)
        c.SetTopMargin(0.10)
        h.Draw("COLZ")
        c.Update()
        c._grid_lines = draw_colz_grid(h)
        draw_cms_label()
        c.Update()

        pal = h.GetListOfFunctions().FindObject("palette")
        if pal:
            h.GetListOfFunctions().Remove(pal)

        c._h = h
        tdir.cd()
        c.Write()
        print(f"    [{cname}] done")


def make_cutflow_table(tdir, fv, cfg):
    """
    Sequential removal-efficiency table for the filtering cuts.
    Produces a grouped bar chart (signal vs non-signal) and prints a text table.
    """
    if fv is None or len(fv) == 0:
        print("    [filter_cutflow] No filteringVtx data — skipping")
        return

    n_tracks     = ak.to_numpy(fv["FilterVtx_nTracks"]).astype(int)
    min_ct       = ak.to_numpy(fv["FilterVtx_minTrackCosTheta"]).astype(float)
    max_abs_cm   = ak.to_numpy(fv["FilterVtx_maxAbsCosThetaCM"]).astype(float)
    max_slope    = ak.to_numpy(fv["FilterVtx_maxSlopeMetric"]).astype(float)
    charge       = ak.to_numpy(fv["FilterVtx_charge"]).astype(int)
    dxy_signif   = ak.to_numpy(fv["FilterVtx_dxySignif"]).astype(float)
    is_gold      = ak.to_numpy(fv["FilterVtx_isGold"]).astype(bool)

    # Build slope-cut mask using stored metric when slope ≈ -1, else maxAbsCosThetaCM
    slope     = float(cfg["trackCosThetaCM_Slope"])     if cfg else -1.0
    intercept = float(cfg["maxTrackCosThetaCM_Intercept"]) if cfg else 1.8
    if abs(slope + 1.0) < 1e-6:
        slope_fail = max_slope > intercept
    elif abs(slope) < 1e-6:
        slope_fail = max_abs_cm > intercept
    else:
        print(f"    [filter_cutflow] Warning: trackCosThetaCM_Slope={slope:.2f} not "
              "directly supported; skipping slope cut in table")
        slope_fail = np.zeros(len(n_tracks), dtype=bool)

    req_charge = bool(cfg["requireChargeNeutrality"]) if cfg else True

    cuts = [
        ("size #neq 2",        n_tracks != 2),
        ("trackcos#theta",     min_ct < (float(cfg["minTrackCosTheta"]) if cfg else 0.5)),
        ("|cos#theta*| slope", slope_fail),
        ("|cos#theta*| limit", max_abs_cm > (float(cfg["maxTrackCosThetaCM_Limit"]) if cfg else 0.95)),
        ("charge",             (charge != 0) if req_charge else np.zeros(len(n_tracks), dtype=bool)),
        ("dxy significance",   dxy_signif <= (float(cfg["minDxySignificance"]) if cfg else 25.0)),
    ]

    surviving = np.ones(len(n_tracks), dtype=bool)
    cut_labels, sig_fracs, nonsig_fracs = [], [], []

    header = f"  {'Cut':<22} | {'Signal':>18} | {'Non-signal':>18}"
    print("  [filter_cutflow] Sequential removal per cut:")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for label, fail_mask in cuts:
        n_sig_in    = int(np.sum( surviving &  is_gold))
        n_nonsig_in = int(np.sum( surviving & ~is_gold))
        fails       = surviving & fail_mask
        n_sig_fail    = int(np.sum(fails &  is_gold))
        n_nonsig_fail = int(np.sum(fails & ~is_gold))

        sig_frac    = n_sig_fail    / n_sig_in    if n_sig_in    > 0 else 0.0
        nonsig_frac = n_nonsig_fail / n_nonsig_in if n_nonsig_in > 0 else 0.0

        print(f"  {label:<22} | "
              f"{n_sig_fail:>5}/{n_sig_in:<5} ({sig_frac*100:>5.1f}%) | "
              f"{n_nonsig_fail:>5}/{n_nonsig_in:<5} ({nonsig_frac*100:>5.1f}%)")

        cut_labels.append(label)
        sig_fracs.append(sig_frac)
        nonsig_fracs.append(nonsig_frac)
        surviving &= ~fail_mask

    # ── ROOT bar chart ────────────────────────────────────────────────────────
    n_bins = len(cut_labels)
    h_sig    = ROOT.TH1F("h_cutflow_sig",    "", n_bins, 0.5, n_bins + 0.5)
    h_nonsig = ROOT.TH1F("h_cutflow_nonsig", "", n_bins, 0.5, n_bins + 0.5)

    for i, (lbl, sf, nf) in enumerate(zip(cut_labels, sig_fracs, nonsig_fracs), 1):
        h_sig   .GetXaxis().SetBinLabel(i, lbl)
        h_nonsig.GetXaxis().SetBinLabel(i, lbl)
        h_sig   .SetBinContent(i, sf)
        h_nonsig.SetBinContent(i, nf)

    h_sig.SetFillColor(COLOR_GOLD);     h_sig.SetLineColor(COLOR_GOLD)
    h_sig.SetBarWidth(0.4);             h_sig.SetBarOffset(0.05)
    h_nonsig.SetFillColor(COLOR_NONSIGNAL); h_nonsig.SetLineColor(COLOR_NONSIGNAL)
    h_nonsig.SetFillStyle(3004)
    h_nonsig.SetBarWidth(0.4);          h_nonsig.SetBarOffset(0.55)

    max_y = max(max(sig_fracs, default=0.0), max(nonsig_fracs, default=0.0))
    h_ax = ROOT.TH1F("h_ax_cutflow", ";Filter cut;Fraction removed",
                     n_bins, 0.5, n_bins + 0.5)
    h_ax.SetMinimum(0.0)
    h_ax.SetMaximum(max_y * 1.4 if max_y > 0 else 1.0)
    for i, lbl in enumerate(cut_labels, 1):
        h_ax.GetXaxis().SetBinLabel(i, lbl)
    h_ax.GetXaxis().SetLabelSize(0.048)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.04); h_ax.GetYaxis().SetTitleSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.4); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)

    c = ROOT.TCanvas("filter_cutflow", "Filter cutflow", 950, 600)
    c.SetLeftMargin(0.12); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.16); c.SetTopMargin(0.10)
    h_ax.Draw("AXIS")
    h_sig.Draw("BAR HIST SAME")
    h_nonsig.Draw("BAR HIST SAME")
    c._grid_lines = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.60, 0.76, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.035)
    leg.AddEntry(h_sig,    "Signal (gold)", "f")
    leg.AddEntry(h_nonsig, "Non-signal",    "f")
    leg.Draw()
    draw_cms_label()
    c.Update()

    c._h_ax = h_ax; c._h_sig = h_sig; c._h_nonsig = h_nonsig; c._leg = leg
    tdir.cd()
    c.Write()
    print("    [filter_cutflow] done")


# ── Orchestrator ───────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg, fv=None, cfg=None):
    print("  [filtering] Reco observables...")
    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "filtered", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_filtered_{obs_key}] done")
    print("  [filtering] 2D cosTheta vs decay angle...")
    plot_filtering_costheta_decay_2d(tdir, sv_sig)
    print("  [filtering] Sequential cut-flow table...")
    make_cutflow_table(tdir, fv, cfg)
