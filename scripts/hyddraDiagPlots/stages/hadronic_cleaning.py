"""
stages/hadronic_cleaning.py — cleaning-stage plots for hadronic HYDDRA.

Unlike leptonic cleaning (per-track removal), hadronic cleaning is purely
vertex-level: vertices failing any of the following cuts are dropped entirely.
  - nTracks < minSize
  - normChi2 > maxNormChi2
  - cosTheta < minCosTheta
  - |decayAngle| > maxDecayAngle
  - mass < minMass
  - p/E < minPOverE

Plots:
  - reco_cleaned_{cosTheta,decayAngle,pOverE,dxySignif,mass,nTracks}
  - hadronic_cleaning_cutflow  (independent per-cut removal bar chart)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import HADRONIC_RECO_OBSERVABLES, STAGE_IDX, COLOR_GOLD, COLOR_NONSIGNAL
from ..src.style   import draw_cms_label, draw_axis_grid
from ..src.plotter import had_signal_mask, plot_reco_observable


def make_cleaning_cutflow(tdir, sv_sig, cfg):
    """
    Independent removal-efficiency table for the hadronic cleaning cuts.
    Each cut is evaluated against ALL post-merge vertices so that every
    fraction shares the same denominator and cuts are directly comparable.
    Produces a grouped bar chart (signal gold vs non-signal).
    """
    if sv_sig is None or len(sv_sig) == 0:
        print("    [had_cleaning_cutflow] No allStageVtx data — skipping")
        return
    if "StageVtx_nTracks" not in sv_sig.fields:
        print("    [had_cleaning_cutflow] StageVtx_nTracks not in ntuple — skipping "
              "(regenerate with hadronic analyzer)")
        return

    sidx_merged = STAGE_IDX["merged"]
    mask = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx_merged
    if not np.any(mask):
        print("    [had_cleaning_cutflow] No merged vertices found — skipping")
        return

    n_tracks   = ak.to_numpy(sv_sig["StageVtx_nTracks"    ][mask]).astype(int)
    norm_chi2  = ak.to_numpy(sv_sig["StageVtx_normChi2"   ][mask]).astype(float)
    cos_theta  = ak.to_numpy(sv_sig["StageVtx_cosTheta"   ][mask]).astype(float)
    decay_angle= ak.to_numpy(sv_sig["StageVtx_decayAngle" ][mask]).astype(float)
    mass       = ak.to_numpy(sv_sig["StageVtx_mass"       ][mask]).astype(float)
    p_over_e   = ak.to_numpy(sv_sig["StageVtx_pOverE"     ][mask]).astype(float)
    is_signal  = had_signal_mask(sv_sig)[mask]

    min_size      = int  (cfg["minSize"])       if cfg else 5
    max_norm_chi2 = float(cfg["maxNormChi2"])   if cfg else 5.0
    min_cos       = float(cfg["minCosTheta"])   if cfg else 0.0
    max_da        = float(cfg["maxDecayAngle"]) if cfg else 0.9
    min_mass      = float(cfg["minMass"])       if cfg else 2.0
    min_poe       = float(cfg["minPOverE"])     if cfg else 0.6

    cuts = [
        (f"nTracks < {min_size}",              n_tracks    <  min_size),
        (f"chi2 > {max_norm_chi2:.1f}",        norm_chi2   >  max_norm_chi2),
        (f"cos#theta < {min_cos:.1f}",         cos_theta   <  min_cos),
        (f"|decay angle| > {max_da:.2f}",      np.fabs(decay_angle) > max_da),
        (f"mass < {min_mass:.1f} GeV",         mass        <  min_mass),
        (f"p/E < {min_poe:.2f}",               p_over_e    <  min_poe),
    ]

    n_gold_total   = int(np.sum( is_signal))
    n_nonsig_total = int(np.sum(~is_signal))

    cut_labels, sig_fracs, nonsig_fracs = [], [], []

    col_w = 28
    header = (f"  {'Cut':<26} | {'Signal lost':{col_w}} | "
              f"{'Non-sig reduction':{col_w}} | {'Non-sig removed':>10}")
    sep = "  " + "-" * (len(header) - 2)
    lines = [
        "  [had_cleaning_cutflow] Per-cut removal (denominator = all merged SVs):",
        f"  {'Denominator':<26} | "
        f"{'':>5}{n_gold_total:<5}{'':>18}| "
        f"{'':>5}{n_nonsig_total:<5}{'':>18}|",
        header, sep,
    ]

    for label, fail_mask in cuts:
        n_sig_fail    = int(np.sum(fail_mask &  is_signal))
        n_nonsig_fail = int(np.sum(fail_mask & ~is_signal))
        sig_frac    = n_sig_fail    / n_gold_total   if n_gold_total   > 0 else 0.0
        nonsig_frac = n_nonsig_fail / n_nonsig_total if n_nonsig_total > 0 else 0.0
        lines.append(
            f"  {label:<26} | "
            f"{n_sig_fail:>5}/{n_gold_total:<5} ({sig_frac*100:>5.1f}%) | "
            f"{n_nonsig_fail:>5}/{n_nonsig_total:<5} ({nonsig_frac*100:>5.1f}%) | "
            f"{n_nonsig_fail:>10}"
        )
        cut_labels.append(label)
        sig_fracs.append(sig_frac)
        nonsig_fracs.append(nonsig_frac)

    table_str = "\n".join(lines)

    # ── ROOT bar chart ────────────────────────────────────────────────────────
    n_bins = len(cut_labels)
    h_sig    = ROOT.TH1F("h_had_cutflow_sig",    "", n_bins, 0.5, n_bins + 0.5)
    h_nonsig = ROOT.TH1F("h_had_cutflow_nonsig", "", n_bins, 0.5, n_bins + 0.5)

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
    h_ax = ROOT.TH1F("h_ax_had_cutflow",
                     ";Hadronic cleaning cut;Fraction of merged SVs removed",
                     n_bins, 0.5, n_bins + 0.5)
    h_ax.SetMinimum(0.0)
    h_ax.SetMaximum(max_y * 1.4 if max_y > 0 else 1.0)
    for i, lbl in enumerate(cut_labels, 1):
        h_ax.GetXaxis().SetBinLabel(i, lbl)
    h_ax.GetXaxis().SetLabelSize(0.042)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.04); h_ax.GetYaxis().SetTitleSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.5); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)

    c = ROOT.TCanvas("hadronic_cleaning_cutflow", "Hadronic cleaning cutflow", 1000, 600)
    c.SetLeftMargin(0.12); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.18); c.SetTopMargin(0.10)
    h_ax.Draw("AXIS")
    h_sig.Draw("BAR HIST SAME")
    h_nonsig.Draw("BAR HIST SAME")
    c._grid_lines = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.60, 0.76, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.035)
    leg.AddEntry(h_sig,    "Signal (matchRatio > 0)", "f")
    leg.AddEntry(h_nonsig, "Non-signal",    "f")
    leg.Draw()
    draw_cms_label("Hadronic HYDDRA")
    c.Update()

    c._h_ax = h_ax; c._h_sig = h_sig; c._h_nonsig = h_nonsig; c._leg = leg
    tdir.cd()
    c.Write()
    return table_str


# ── Orchestrator ───────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg, cfg=None):
    for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "cleaned", obs_key, obs_cfg, sv_bkg, hadronic=True)
    return make_cleaning_cutflow(tdir, sv_sig, cfg)
