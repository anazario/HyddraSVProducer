"""
stages/filtering.py — filtering-stage plots for hyddraDiagPlots.

Plots:
  - reco_filtered_{cosTheta,decayAngle,pOverE,dxySignif}
  - filtering_mintrackcos_decay_{signal,nonsignal}_{before,after}
  - reco_filtered_minTrackCosTheta  (1D, filtering stage only)
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
    2D (trackCosThetaCM, trackCosTheta) scatter before and after filtering,
    split by signal (isGold) and non-signal.  Four canvases total.
    One point per track (two points per 2-track vertex), using the exact
    per-track variables evaluated in LeptonicHYDDRA::filteringImpl so that
    the slope cut boundary |cos theta*_CM| = slope * cos theta + intercept
    maps directly onto the plotted axes.
    """
    if sv_sig is None or len(sv_sig) == 0:
        print("    [filtering_2d] No allStageVtx data — skipping")
        return

    if "StageVtx_trackCosTheta" not in sv_sig.fields:
        print("    [filtering_2d] StageVtx_trackCosTheta not in ntuple — skipping")
        return

    sidx_pre  = STAGE_IDX["disambig"]
    sidx_post = STAGE_IDX["filtered"]

    # Jagged arrays: one inner list per vertex, one element per track
    trk_ct_jagged   = sv_sig["StageVtx_trackCosTheta"]
    trk_ctcm_jagged = sv_sig["StageVtx_trackCosThetaCM"]

    # Broadcast vertex-level scalars to per-track level, then flatten
    stage   = ak.to_numpy(ak.flatten(
                  ak.broadcast_arrays(sv_sig["StageVtx_stageIdx"], trk_ct_jagged)[0]))
    is_gold = ak.to_numpy(ak.flatten(
                  ak.broadcast_arrays(sv_sig["StageVtx_isGold"], trk_ct_jagged)[0])).astype(bool)
    trk_ct   = ak.to_numpy(ak.flatten(trk_ct_jagged)).astype(float)
    trk_ctcm = ak.to_numpy(ak.flatten(trk_ctcm_jagged)).astype(float)

    x_bins = np.linspace(-1, 1, 101)
    y_bins = np.linspace(-1, 1, 101)

    def fill_h2(name, mask):
        ROOT.gStyle.SetPalette(ROOT.kViridis)
        h = ROOT.TH2F(name, name, len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)
        h.SetDirectory(0)
        for da, ct in zip(trk_ctcm[mask], trk_ct[mask]):
            h.Fill(da, ct)
        h.GetXaxis().SetTitle("cos#theta^{*}_{CM}")
        h.GetYaxis().SetTitle("track cos#theta")
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
        cname = f"filtering_mintrackcos_decay_{suffix}"
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


def _plot_filtering_1d(tdir, sv_sig, sv_bkg, branch, x_title, cname, tag):
    """
    Generic helper: 1D per-track distribution of a jagged StageVtx branch,
    showing pre-filter (disambig, solid) and post-filter (filtered, dashed)
    for signal and non-signal.  Background (post-filter only) added if sv_bkg
    is provided and contains the branch.
    """
    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping")
        return
    if branch not in sv_sig.fields:
        print(f"    [{tag}] {branch} not in ntuple — skipping")
        return

    from ..src.config import COLOR_GOLD, COLOR_NONSIGNAL, COLOR_BKG

    sidx_pre  = STAGE_IDX["disambig"]
    sidx_post = STAGE_IDX["filtered"]
    bins = np.linspace(-1, 1, 51)

    jag     = sv_sig[branch]
    stage   = ak.to_numpy(ak.flatten(ak.broadcast_arrays(sv_sig["StageVtx_stageIdx"], jag)[0]))
    is_gold = ak.to_numpy(ak.flatten(ak.broadcast_arrays(sv_sig["StageVtx_isGold"],   jag)[0])).astype(bool)
    vals    = ak.to_numpy(ak.flatten(jag)).astype(float)

    def make_th1(name, data, mask, color, lstyle):
        h = ROOT.TH1F(name, name, len(bins)-1, bins)
        h.SetDirectory(0)
        for v in data[mask]:
            h.Fill(v)
        h.SetLineColor(color); h.SetLineWidth(2); h.SetLineStyle(lstyle)
        h.SetFillStyle(0); h.SetStats(0)
        h.GetXaxis().SetTitle(x_title)
        h.GetYaxis().SetTitle("Tracks")
        h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.2); h.GetYaxis().SetTitleOffset(1.3)
        return h

    h_sig_pre    = make_th1(f"h_{tag}_sig_pre",    vals, (stage == sidx_pre)  &  is_gold, COLOR_GOLD,      1)
    h_sig_post   = make_th1(f"h_{tag}_sig_post",   vals, (stage == sidx_post) &  is_gold, COLOR_GOLD,      2)
    h_nonsig_pre = make_th1(f"h_{tag}_nonsig_pre", vals, (stage == sidx_pre)  & ~is_gold, COLOR_NONSIGNAL, 1)
    h_nonsig_post= make_th1(f"h_{tag}_nonsig_post",vals, (stage == sidx_post) & ~is_gold, COLOR_NONSIGNAL, 2)

    hists  = [h_sig_pre, h_sig_post, h_nonsig_pre, h_nonsig_post]
    labels = ["Signal (pre-filter)", "Signal (post-filter)",
              "Non-signal (pre-filter)", "Non-signal (post-filter)"]

    if sv_bkg is not None and len(sv_bkg) > 0 and branch in sv_bkg.fields:
        bkg_jag   = sv_bkg[branch]
        bkg_stage = ak.to_numpy(ak.flatten(
                        ak.broadcast_arrays(sv_bkg["StageVtx_stageIdx"], bkg_jag)[0]))
        bkg_vals  = ak.to_numpy(ak.flatten(bkg_jag)).astype(float)
        h_bkg = make_th1(f"h_{tag}_bkg", bkg_vals, bkg_stage == sidx_post, COLOR_BKG, 1)
        hists.append(h_bkg)
        labels.append("Background (post-filter)")

    max_y = max((h.GetMaximum() for h in hists), default=1.0)

    c = ROOT.TCanvas(cname, x_title, 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)
    c.SetLogy(True)

    for i, h in enumerate(hists):
        h.SetMaximum(max_y * 5)
        h.SetMinimum(0.5)
        h.Draw("HIST" if i == 0 else "HIST SAME")

    c._grid_lines = draw_axis_grid(hists[0], logy=True)

    leg = ROOT.TLegend(0.15, 0.70, 0.55, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.033)
    for h, lbl in zip(hists, labels):
        leg.AddEntry(h, lbl, "l")
    leg.Draw()
    draw_cms_label()
    c.Update()

    c._hists = hists
    c._leg   = leg
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


def plot_filtering_mintrackcos_1d(tdir, sv_sig, sv_bkg=None):
    """
    1D per-track distribution of track cos theta (wrt PV) at the filtering
    stage, with pre-filter (solid) and post-filter (dashed) overlaid.
    """
    _plot_filtering_1d(tdir, sv_sig, sv_bkg,
                       branch="StageVtx_trackCosTheta",
                       x_title="track cos#theta",
                       cname="reco_filtering_trackCosTheta",
                       tag="reco_filtering_trackCosTheta")


def plot_filtering_trackcosThetaCM_1d(tdir, sv_sig, sv_bkg=None):
    """
    1D per-track distribution of track decay angle in the CM frame at the
    filtering stage, with pre-filter (solid) and post-filter (dashed) overlaid.
    """
    _plot_filtering_1d(tdir, sv_sig, sv_bkg,
                       branch="StageVtx_trackCosThetaCM",
                       x_title="cos#theta^{*}_{CM}",
                       cname="reco_filtering_trackCosThetaCM",
                       tag="reco_filtering_trackCosThetaCM")


def make_cutflow_table(tdir, sv_sig, cfg):
    """
    Independent removal-efficiency table for the filtering cuts.
    Each cut is evaluated against ALL disambiguated vertices so that every
    fraction shares the same denominator and cuts are directly comparable.
    Produces a grouped bar chart (signal vs non-signal) and prints a text table.
    The bar chart y-axis is "fraction of disambig SVs removed" — same denominator
    for all bars, so heights are directly comparable across cuts.
    """
    if sv_sig is None or len(sv_sig) == 0:
        print("    [filter_cutflow] No allStageVtx data — skipping")
        return

    # Slice to the disambiguation stage
    sidx_disambig = STAGE_IDX["disambig"]
    mask = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx_disambig
    if not np.any(mask):
        print("    [filter_cutflow] No disambiguated vertices found — skipping")
        return

    trk_ct_jag   = sv_sig["StageVtx_trackCosTheta"][mask]
    trk_ctcm_jag = sv_sig["StageVtx_trackCosThetaCM"][mask]

    n_tracks   = ak.to_numpy(ak.num(trk_ct_jag)).astype(int)
    min_ct     = ak.to_numpy(ak.min(trk_ct_jag,             axis=1)).astype(float)
    max_abs_cm = ak.to_numpy(ak.max(ak.abs(trk_ctcm_jag),   axis=1)).astype(float)
    max_slope  = ak.to_numpy(ak.max(ak.abs(trk_ctcm_jag) + trk_ct_jag, axis=1)).astype(float)

    vtx_cos    = ak.to_numpy(sv_sig["StageVtx_cosTheta"  ][mask]).astype(float)
    decay_angle= ak.to_numpy(sv_sig["StageVtx_decayAngle"][mask]).astype(float)
    dxy_signif = ak.to_numpy(sv_sig["StageVtx_dxySignif" ][mask]).astype(float)
    charge     = ak.to_numpy(sv_sig["StageVtx_charge"    ][mask]).astype(int)
    sv_mass    = ak.to_numpy(sv_sig["StageVtx_mass"      ][mask]).astype(float)
    sv_beta    = ak.to_numpy(sv_sig["StageVtx_pOverE"    ][mask]).astype(float)
    is_gold    = ak.to_numpy(sv_sig["StageVtx_isGold"    ][mask]).astype(bool)

    # Build slope-cut mask using pre-computed metric when slope ≈ -1
    slope     = float(cfg["trackCosThetaCM_Slope"])        if cfg else -1.0
    intercept = float(cfg["maxTrackCosThetaCM_Intercept"]) if cfg else 1.8
    if abs(slope + 1.0) < 1e-6:
        slope_fail = max_slope > intercept
    elif abs(slope) < 1e-6:
        slope_fail = max_abs_cm > intercept
    else:
        print(f"    [filter_cutflow] Warning: trackCosThetaCM_Slope={slope:.2f} not "
              "directly supported; skipping slope cut in table")
        slope_fail = np.zeros(len(n_tracks), dtype=bool)

    req_charge      = bool(cfg["requireChargeNeutrality"]) if cfg else True
    min_vtx_cos     = float(cfg["minVtxCosTheta"])   if cfg and "minVtxCosTheta"    in cfg else -1.0
    max_vtx_cos     = float(cfg["maxVtxCosTheta"])   if cfg and "maxVtxCosTheta"    in cfg else  1.0
    use_abs_vtx_cos = bool(cfg["useAbsVtxCosTheta"]) if cfg and "useAbsVtxCosTheta" in cfg else False
    val_cos = np.fabs(vtx_cos) if use_abs_vtx_cos else vtx_cos
    vtx_cos_fail = (val_cos < min_vtx_cos) | (val_cos > max_vtx_cos)

    apply_da_filter = bool(cfg["applyVtxDecayAngleFiltering"]) if cfg and "applyVtxDecayAngleFiltering" in cfg else False
    max_vtx_da      = float(cfg["maxVtxDecayAngle"])   if cfg and "maxVtxDecayAngle"    in cfg else 1.0
    use_abs_vtx_da  = bool(cfg["useAbsVtxDecayAngle"]) if cfg and "useAbsVtxDecayAngle" in cfg else False
    if apply_da_filter:
        val_da = np.fabs(decay_angle) if use_abs_vtx_da else decay_angle
        vtx_da_fail = val_da > max_vtx_da
    else:
        vtx_da_fail = np.zeros(len(n_tracks), dtype=bool)

    min_mass_filter = float(cfg["minMassFilter"]) if cfg and "minMassFilter" in cfg else 0.0
    min_beta_filter = float(cfg["minBetaFilter"]) if cfg and "minBetaFilter" in cfg else 0.0

    cuts = [
        ("size #neq 2",        n_tracks != 2),
        ("trackcos#theta",     min_ct < (float(cfg["minTrackCosTheta"]) if cfg else 0.5)),
        ("|cos#theta*| slope", slope_fail),
        ("|cos#theta*| limit", max_abs_cm > (float(cfg["maxTrackCosThetaCM_Limit"]) if cfg else 0.95)),
        ("charge",             (charge != 0) if req_charge else np.zeros(len(n_tracks), dtype=bool)),
        ("vtxcos#theta",       vtx_cos_fail),
        ("vtx decay angle",    vtx_da_fail),
        ("dxy significance",   dxy_signif <= (float(cfg["minDxySignificance"]) if cfg else 25.0)),
        ("mass filter",        sv_mass < min_mass_filter),
        ("beta filter",        sv_beta < min_beta_filter),
    ]

    # Fixed denominators: all disambiguated SVs, split by gold / non-signal.
    # Every cut fraction is relative to these totals so cuts are comparable.
    n_gold_total   = int(np.sum( is_gold))
    n_nonsig_total = int(np.sum(~is_gold))

    cut_labels, sig_fracs, nonsig_fracs = [], [], []

    col_w = 24
    header = (f"  {'Cut':<22} | {'Gold lost':{col_w}} | "
              f"{'Non-sig reduction':{col_w}} | {'Non-sig removed':>10}")
    sep    = "  " + "-" * (len(header) - 2)
    print("  [filter_cutflow] Per-cut removal (denominator = all disambig SVs):")
    print(f"  {'Denominator':<22} | "
          f"{'':>5}{n_gold_total:<5}{'':>14}| "
          f"{'':>5}{n_nonsig_total:<5}{'':>14}|")
    print(header)
    print(sep)

    for label, fail_mask in cuts:
        n_sig_fail    = int(np.sum(fail_mask &  is_gold))
        n_nonsig_fail = int(np.sum(fail_mask & ~is_gold))

        sig_frac    = n_sig_fail    / n_gold_total   if n_gold_total   > 0 else 0.0
        nonsig_frac = n_nonsig_fail / n_nonsig_total if n_nonsig_total > 0 else 0.0

        print(f"  {label:<22} | "
              f"{n_sig_fail:>5}/{n_gold_total:<5} ({sig_frac*100:>5.1f}%) | "
              f"{n_nonsig_fail:>5}/{n_nonsig_total:<5} ({nonsig_frac*100:>5.1f}%) | "
              f"{n_nonsig_fail:>10}")

        cut_labels.append(label)
        sig_fracs.append(sig_frac)
        nonsig_fracs.append(nonsig_frac)

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
    h_ax = ROOT.TH1F("h_ax_cutflow", ";Filter cut;Fraction of disambig SVs removed",
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

def make_plots(tdir, gf, sv_sig, sv_bkg, cfg=None):
    print("  [filtering] Reco observables...")
    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "filtered", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_filtered_{obs_key}] done")
    print("  [filtering] 1D track cos theta (pre/post filter)...")
    plot_filtering_mintrackcos_1d(tdir, sv_sig, sv_bkg)
    print("  [filtering] 1D track cos theta CM (pre/post filter)...")
    plot_filtering_trackcosThetaCM_1d(tdir, sv_sig, sv_bkg)
    print("  [filtering] 2D track cos theta vs decay angle...")
    plot_filtering_costheta_decay_2d(tdir, sv_sig)
    print("  [filtering] Sequential cut-flow table...")
    make_cutflow_table(tdir, sv_sig, cfg)
