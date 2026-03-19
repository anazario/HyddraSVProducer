"""
stages/hadronic_id.py — ID-stage plots for hadronic HYDDRA.

The hadronic ID cut requires nTracks >= 3 AND mass / nTracks > 1.

Plots:
  - reco_id_{cosTheta,decayAngle,pOverE,dxySignif,mass,nTracks}
  - hadronic_id_massOverNTracks  (mass/nTracks pre/post ID, signal vs non-signal)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import HADRONIC_RECO_OBSERVABLES, STAGE_IDX, COLOR_GOLD, COLOR_NONSIGNAL, COLOR_BKG, GEN_DXY_BINS
from ..src.style   import draw_cms_label, draw_axis_grid
from ..src.plotter import (plot_reco_observable, plot_costheta_zoom,
                            plot_dxy_1d, plot_2d_dxy_vs_x,
                            geom_cut_lines_z, geom_cut_lines_eta)


def plot_mass_over_ntracks(tdir, sv_sig, sv_bkg=None):
    """
    1D mass/nTracks distribution comparing pre-ID (filtered stage) vs post-ID
    for signal and non-signal.  A dashed vertical line marks the ID cut at 1.0.
    """
    tag = "hadronic_id_massOverNTracks"

    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping")
        return

    if "StageVtx_nTracks" not in sv_sig.fields:
        print(f"    [{tag}] StageVtx_nTracks not in ntuple — skipping")
        return

    sidx_pre  = STAGE_IDX["filtered"]
    sidx_post = STAGE_IDX["id"]
    bins = np.linspace(0, 10, 51)

    stage     = ak.to_numpy(sv_sig["StageVtx_stageIdx"]).astype(int)
    mass      = ak.to_numpy(sv_sig["StageVtx_mass"    ]).astype(float)
    ntracks   = ak.to_numpy(sv_sig["StageVtx_nTracks" ]).astype(int)
    is_signal = ak.to_numpy(sv_sig["StageVtx_matchRatio"]) > 0

    valid = ntracks > 0
    mass_over_n = np.where(valid, mass / np.where(valid, ntracks, 1).astype(float), -1.0)

    def make_th1(name, mask, color, lstyle):
        h = ROOT.TH1F(name, name, len(bins) - 1, bins)
        h.SetDirectory(0)
        combined = mask & valid
        for v in mass_over_n[combined]:
            h.Fill(v)
        h.SetLineColor(color); h.SetLineWidth(2); h.SetLineStyle(lstyle)
        h.SetFillStyle(0); h.SetStats(0)
        h.GetXaxis().SetTitle("Mass / n_{tracks} (GeV)")
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
    labels = ["Signal (pre-ID)", "Signal (post-ID)",
              "Non-signal (pre-ID)", "Non-signal (post-ID)"]

    if sv_bkg is not None and len(sv_bkg) > 0 and "StageVtx_nTracks" in sv_bkg.fields:
        bkg_stage   = ak.to_numpy(sv_bkg["StageVtx_stageIdx"]).astype(int)
        bkg_mass    = ak.to_numpy(sv_bkg["StageVtx_mass"    ]).astype(float)
        bkg_ntracks = ak.to_numpy(sv_bkg["StageVtx_nTracks" ]).astype(int)
        bkg_valid   = bkg_ntracks > 0
        bkg_mon     = np.where(bkg_valid, bkg_mass / np.where(bkg_valid, bkg_ntracks, 1).astype(float), -1.0)
        h_bkg = ROOT.TH1F(f"h_{tag}_bkg", f"h_{tag}_bkg", len(bins) - 1, bins)
        h_bkg.SetDirectory(0)
        for v in bkg_mon[(bkg_stage == sidx_post) & bkg_valid]:
            h_bkg.Fill(v)
        h_bkg.SetLineColor(COLOR_BKG); h_bkg.SetLineWidth(2); h_bkg.SetLineStyle(1)
        h_bkg.SetFillStyle(0); h_bkg.SetStats(0)
        h_bkg.GetXaxis().SetTitle("Mass / n_{tracks} (GeV)")
        h_bkg.GetYaxis().SetTitle("Vertices")
        hists.append(h_bkg)
        labels.append("Background (post-ID)")

    max_y = max((h.GetMaximum() for h in hists), default=1.0)

    c = ROOT.TCanvas(tag, "Mass / nTracks at ID stage", 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)
    c.SetLogy(True)

    for i, h in enumerate(hists):
        h.SetMaximum(max_y * 5)
        h.SetMinimum(0.5)
        h.Draw("HIST" if i == 0 else "HIST SAME")

    line = ROOT.TLine(1.0, 0.5, 1.0, max_y * 5)
    line.SetLineColor(ROOT.kBlack); line.SetLineWidth(2); line.SetLineStyle(2)
    line.Draw()

    c._grid_lines = draw_axis_grid(hists[0], logy=True)

    leg = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.033)
    for h, lbl in zip(hists, labels):
        leg.AddEntry(h, lbl, "l")
    leg.Draw()
    draw_cms_label()
    c.Update()

    c._hists = hists; c._leg = leg; c._line = line
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


def plot_eff_vs_dxy(tdir, gf):
    """
    Signal efficiency vs gen dxy at the hadronic ID stage.
    Two curves: tight (matchRatio >= 0.5) and loose (matchRatio > 0).
    """
    tag = "hadronic_id_eff_vs_dxy"

    if gf is None or len(gf) == 0:
        print(f"    [{tag}] No gen funnel data — skipping")
        return
    if "GenFunnel_matchRatio_id" not in gf.fields:
        print(f"    [{tag}] GenFunnel_matchRatio_id not in ntuple — skipping")
        return

    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
    mr_id      = ak.to_numpy(ak.flatten(gf["GenFunnel_matchRatio_id"]))[has_tracks]
    dxy        = ak.to_numpy(ak.flatten(gf["GenFunnel_dxy"          ]))[has_tracks]

    bins   = np.array(GEN_DXY_BINS, dtype=float)
    n_bins = len(bins) - 1

    h_denom = ROOT.TH1F(f"h_denom_{tag}", "", n_bins, bins)
    h_denom.SetDirectory(0)
    for v in dxy:
        h_denom.Fill(v)

    sig_configs = [
        ("tight", mr_id >= 0.5, ROOT.kRed+1,   20, 1),
        ("loose", mr_id > 0,    ROOT.kAzure+6,  24, 2),
    ]

    graphs, hists = [], []
    for label, sig_mask, color, marker, lstyle in sig_configs:
        h_num = ROOT.TH1F(f"h_num_{tag}_{label}", "", n_bins, bins)
        h_num.SetDirectory(0)
        for v in dxy[sig_mask]:
            h_num.Fill(v)

        eff = ROOT.TEfficiency(h_num, h_denom)
        eff.SetStatisticOption(ROOT.TEfficiency.kFCP)
        g = eff.CreateGraph()
        leg_title = (f"Signal tight (matchRatio #geq 0.5)"
                     if label == "tight" else "Signal loose (matchRatio > 0)")
        g.SetTitle(leg_title)
        g.SetMarkerColor(color); g.SetLineColor(color)
        g.SetMarkerStyle(marker); g.SetLineWidth(2); g.SetLineStyle(lstyle)
        graphs.append(g)
        hists.append((h_num, eff))

    c = ROOT.TCanvas(tag, "Hadronic ID signal efficiency vs gen dxy", 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)

    h_ax = ROOT.TH1F(f"h_ax_{tag}", ";Gen dxy (cm);Signal Efficiency at ID Stage", n_bins, bins)
    h_ax.SetMinimum(0.0); h_ax.SetMaximum(1.2)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.04);  h_ax.GetYaxis().SetLabelSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.2); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")

    for g in graphs:
        g.Draw("P SAME")

    c._grid_lines = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.40, 0.20, 0.88, 0.35)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.035)
    for g in graphs:
        leg.AddEntry(g, g.GetTitle(), "lp")
    leg.Draw()
    draw_cms_label("Hadronic HYDDRA")
    c.Update()

    c._h_ax = h_ax; c._h_denom = h_denom; c._hists = hists; c._graphs = graphs; c._leg = leg
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


# ── Orchestrator ───────────────────────────────────────────────────────────────

def _had_id_sig_mask(sv_sig):
    return ak.to_numpy(sv_sig["StageVtx_matchRatio"]).astype(float) > 0


def _had_id_pos(sv_sig, sidx):
    mask     = ak.to_numpy(sv_sig["StageVtx_stageIdx"]).astype(int) == sidx
    x        = ak.to_numpy(sv_sig["StageVtx_x"]).astype(float)[mask]
    y        = ak.to_numpy(sv_sig["StageVtx_y"]).astype(float)[mask]
    z        = ak.to_numpy(sv_sig["StageVtx_z"]).astype(float)[mask]
    dxy      = np.sqrt(x**2 + y**2)
    eta_orig = np.arcsinh(z / np.where(dxy > 0, dxy, 1e-6))
    is_sig   = _had_id_sig_mask(sv_sig)[mask]
    return z, dxy, eta_orig, is_sig


def make_plots(tdir, gf, sv_sig, sv_bkg=None, cfg=None):
    _cms = "Hadronic HYDDRA"
    sidx   = STAGE_IDX["id"]
    is_sig = _had_id_sig_mask(sv_sig) if sv_sig is not None and len(sv_sig) > 0 else None

    for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "id", obs_key, obs_cfg, sv_bkg, hadronic=True)
    plot_mass_over_ntracks(tdir, sv_sig, sv_bkg)

    if is_sig is not None:
        plot_costheta_zoom(tdir, sv_sig, sv_bkg, sidx, is_sig,
                           "hadronic_id_cosTheta_zoom",
                           sig_label="Signal (matchRatio > 0)", cms_label=_cms)
        plot_dxy_1d(tdir, sv_sig, sv_bkg, sidx, is_sig,
                    "hadronic_id_dxy",
                    sig_label="Signal (matchRatio > 0)", cms_label=_cms)
        if "StageVtx_z" in sv_sig.fields:
            z, dxy, eta, is_sig_pos = _had_id_pos(sv_sig, sidx)
            plot_2d_dxy_vs_x(tdir, z,   dxy, is_sig_pos, "Vertex z (cm)",
                             np.linspace(-30, 30, 61), np.linspace(0, 80, 41),
                             "hadronic_id_dxy_vs_z",
                             sig_label="Signal (matchRatio > 0)",
                             cut_lines_fn=geom_cut_lines_z, cms_label=_cms)
            plot_2d_dxy_vs_x(tdir, eta, dxy, is_sig_pos, "#eta_{origin}",
                             np.linspace(-5, 5, 51), np.linspace(0, 80, 41),
                             "hadronic_id_dxy_vs_etaorigin",
                             sig_label="Signal (matchRatio > 0)",
                             cut_lines_fn=geom_cut_lines_eta, cms_label=_cms)

    plot_eff_vs_dxy(tdir, gf)
