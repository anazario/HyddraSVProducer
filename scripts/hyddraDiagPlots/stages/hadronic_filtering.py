"""
stages/hadronic_filtering.py — filtering-stage plots for hadronic HYDDRA.

Hadronic filtering applies only a dxy significance cut (after rejecting
vertices with dxyError <= 0).  This stage is therefore much simpler than
the leptonic equivalent.

Plots:
  - reco_filtered_{cosTheta,decayAngle,pOverE,dxySignif,mass,nTracks}
  - hadronic_dxySignif_filter  (pre/post-filter dxySignif, signal vs non-signal)
  - hadronic_filter_cutflow    (two-bar chart: dxyError and dxySignif cuts)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import HADRONIC_RECO_OBSERVABLES, STAGE_IDX, COLOR_GOLD, COLOR_NONSIGNAL, COLOR_BKG
from ..src.style   import draw_cms_label, draw_axis_grid
from ..src.plotter import had_signal_mask, plot_reco_observable


def plot_dxySignif_pre_post(tdir, sv_sig, sv_bkg, cfg):
    """
    1D dxySignif: pre-filter (disambig, solid) vs post-filter (filtered, dashed)
    for signal and non-signal.  A vertical line marks the applied cut threshold.
    """
    tag = "hadronic_dxySignif_filter"

    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping")
        return

    sidx_pre  = STAGE_IDX["disambig"]
    sidx_post = STAGE_IDX["filtered"]
    bins = np.linspace(0, 150, 76)

    stage     = ak.to_numpy(sv_sig["StageVtx_stageIdx"]).astype(int)
    is_signal = had_signal_mask(sv_sig)
    vals      = ak.to_numpy(sv_sig["StageVtx_dxySignif"]).astype(float)

    # Mask out sentinel values (-1)
    valid = vals >= 0

    def make_th1(name, mask, color, lstyle):
        h = ROOT.TH1F(name, name, len(bins)-1, bins)
        h.SetDirectory(0)
        combined = mask & valid
        for v in vals[combined]:
            h.Fill(v)
        h.SetLineColor(color); h.SetLineWidth(2); h.SetLineStyle(lstyle)
        h.SetFillStyle(0); h.SetStats(0)
        h.GetXaxis().SetTitle("dxy Significance")
        h.GetYaxis().SetTitle("Vertices")
        h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.2); h.GetYaxis().SetTitleOffset(1.3)
        return h

    h_sig_pre    = make_th1(f"h_{tag}_sig_pre",    (stage == sidx_pre)  &  is_signal, COLOR_GOLD,      1)
    h_sig_post   = make_th1(f"h_{tag}_sig_post",   (stage == sidx_post) &  is_signal, COLOR_GOLD,      2)
    h_nonsig_pre = make_th1(f"h_{tag}_nonsig_pre", (stage == sidx_pre)  & ~is_signal, COLOR_NONSIGNAL, 1)
    h_nonsig_post= make_th1(f"h_{tag}_nonsig_post",(stage == sidx_post) & ~is_signal, COLOR_NONSIGNAL, 2)

    hists  = [h_sig_pre, h_sig_post, h_nonsig_pre, h_nonsig_post]
    labels = ["Signal (pre-filter)", "Signal (post-filter)",
              "Non-signal (pre-filter)", "Non-signal (post-filter)"]

    if sv_bkg is not None and len(sv_bkg) > 0:
        bkg_stage = ak.to_numpy(sv_bkg["StageVtx_stageIdx"]).astype(int)
        bkg_vals  = ak.to_numpy(sv_bkg["StageVtx_dxySignif"]).astype(float)
        bkg_valid = bkg_vals >= 0
        h_bkg = ROOT.TH1F(f"h_{tag}_bkg", "", len(bins)-1, bins)
        h_bkg.SetDirectory(0)
        for v in bkg_vals[(bkg_stage == sidx_post) & bkg_valid]:
            h_bkg.Fill(v)
        h_bkg.SetLineColor(COLOR_BKG); h_bkg.SetLineWidth(2); h_bkg.SetLineStyle(1)
        h_bkg.SetFillStyle(0); h_bkg.SetStats(0)
        hists.append(h_bkg)
        labels.append("Background (post-filter)")

    max_y = max((h.GetMaximum() for h in hists), default=1.0)

    c = ROOT.TCanvas(tag, "dxy Significance (hadronic filter)", 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)
    c.SetLogy(True)

    for i, h in enumerate(hists):
        h.SetMaximum(max_y * 5)
        h.SetMinimum(0.5)
        h.Draw("HIST" if i == 0 else "HIST SAME")

    c._grid_lines = draw_axis_grid(hists[0], logy=True)

    # Draw cut threshold line
    cut_val = float(cfg["minDxySignificance"]) if cfg else 40.0
    ln = ROOT.TLine(cut_val, 0.5, cut_val, max_y * 5)
    ln.SetLineColor(ROOT.kRed); ln.SetLineWidth(2); ln.SetLineStyle(2)
    ln.Draw()

    leg = ROOT.TLegend(0.50, 0.68, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.033)
    for h, lbl in zip(hists, labels):
        leg.AddEntry(h, lbl, "l")
    l1 = ROOT.TLine()
    l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.SetLineWidth(2)
    leg.AddEntry(l1, f"Cut: dxySignif > {cut_val:.0f}", "l")
    leg.Draw()
    draw_cms_label("Hadronic HYDDRA")
    c.Update()

    c._hists = hists; c._leg = leg; c._cut_line = ln
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


def make_filter_cutflow(tdir, sv_sig, cfg):
    """
    Two-bar chart showing independent cut fractions at the disambiguation stage:
      1. dxyError <= 0  (vertex fit quality failure)
      2. dxySignif <= threshold  (displacement significance cut)
    """
    tag = "hadronic_filter_cutflow"

    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping")
        return

    sidx_disambig = STAGE_IDX["disambig"]
    mask = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx_disambig
    if not np.any(mask):
        print(f"    [{tag}] No disambiguated vertices — skipping")
        return

    dxy_signif = ak.to_numpy(sv_sig["StageVtx_dxySignif" ][mask]).astype(float)
    is_signal  = had_signal_mask(sv_sig)[mask]

    min_dxy_sig = float(cfg["minDxySignificance"]) if cfg else 40.0

    # dxySignif == -1 is the sentinel for dxyError <= 0
    cuts = [
        ("dxyError #leq 0",                   dxy_signif <  0),
        (f"dxySignif #leq {min_dxy_sig:.0f}", (dxy_signif >= 0) & (dxy_signif <= min_dxy_sig)),
    ]

    n_gold_total   = int(np.sum( is_signal))
    n_nonsig_total = int(np.sum(~is_signal))

    cut_labels, sig_fracs, nonsig_fracs = [], [], []
    for label, fail_mask in cuts:
        n_sig_fail    = int(np.sum(fail_mask &  is_signal))
        n_nonsig_fail = int(np.sum(fail_mask & ~is_signal))
        sig_fracs.append(n_sig_fail    / n_gold_total   if n_gold_total   > 0 else 0.0)
        nonsig_fracs.append(n_nonsig_fail / n_nonsig_total if n_nonsig_total > 0 else 0.0)
        cut_labels.append(label)

    n_bins = len(cut_labels)
    h_sig    = ROOT.TH1F(f"h_{tag}_sig",    "", n_bins, 0.5, n_bins + 0.5)
    h_nonsig = ROOT.TH1F(f"h_{tag}_nonsig", "", n_bins, 0.5, n_bins + 0.5)

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
    h_ax = ROOT.TH1F(f"h_ax_{tag}", ";Hadronic filter cut;Fraction of disambig SVs removed",
                     n_bins, 0.5, n_bins + 0.5)
    h_ax.SetMinimum(0.0)
    h_ax.SetMaximum(max_y * 1.5 if max_y > 0 else 1.0)
    for i, lbl in enumerate(cut_labels, 1):
        h_ax.GetXaxis().SetBinLabel(i, lbl)
    h_ax.GetXaxis().SetLabelSize(0.05)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.04); h_ax.GetYaxis().SetTitleSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.4); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)

    c = ROOT.TCanvas(tag, "Hadronic filter cutflow", 700, 500)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.14); c.SetTopMargin(0.10)
    h_ax.Draw("AXIS")
    h_sig.Draw("BAR HIST SAME")
    h_nonsig.Draw("BAR HIST SAME")
    c._grid_lines = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.55, 0.76, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.035)
    leg.AddEntry(h_sig,    "Signal (matchRatio > 0)", "f")
    leg.AddEntry(h_nonsig, "Non-signal",    "f")
    leg.Draw()
    draw_cms_label("Hadronic HYDDRA")
    c.Update()

    c._h_ax = h_ax; c._h_sig = h_sig; c._h_nonsig = h_nonsig; c._leg = leg
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


# ── Orchestrator ───────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg, cfg=None):
    print("  [had/filtering] Reco observables...")
    for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "filtered", obs_key, obs_cfg, sv_bkg, hadronic=True)
        print(f"    [reco_filtered_{obs_key}] done")
    print("  [had/filtering] dxySignif pre/post filter...")
    plot_dxySignif_pre_post(tdir, sv_sig, sv_bkg, cfg)
    print("  [had/filtering] Filter cutflow...")
    make_filter_cutflow(tdir, sv_sig, cfg)
