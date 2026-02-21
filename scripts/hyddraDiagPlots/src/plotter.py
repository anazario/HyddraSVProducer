"""
plotter.py — reusable plot helpers for hyddraDiagPlots.
"""
import numpy as np
import awkward as ak
import ROOT

from .config import (
    STAGE_NAMES, STAGE_IDX,
    COLOR_GOLD, COLOR_SILVER, COLOR_BRONZE, COLOR_NONSIGNAL, COLOR_BKG,
)
from .style import make_canvas, draw_cms_label, draw_axis_grid


# ── Axis / graph helpers ──────────────────────────────────────────────────────

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
    g.SetMarkerColor(color)
    g.SetLineColor(color)
    g.SetMarkerStyle(marker)
    g.SetLineWidth(2)
    g.SetLineStyle(line_style)
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


# ── Observable plot ───────────────────────────────────────────────────────────

def _fill_hist(h, vals, mask=None):
    """Fill a TH1F from a numpy array."""
    arr = vals if mask is None else vals[mask]
    for v in arr:
        h.Fill(float(v))


def _style_hist(h, color, line_style=1):
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetLineStyle(line_style)
    h.SetFillStyle(0)
    h.SetStats(0)


def plot_reco_observable(tdir, gf, sv_sig, stage_key, obs_key, obs_cfg, sv_bkg=None):
    """
    Four normalised shape distributions on one canvas:
      - Gold         : genFunnel vertices gold-matched at this stage
      - Silver excl. : silver but not gold
      - Bronze excl. : bronze but not silver
      - Non-signal   : sv_sig allStageVtx where stageIdx==stage and isGold=False
      - Background   : sv_bkg allStageVtx at this stage (dashed, if provided)

    All normalised to unit area.  -1 sentinels are skipped.
    """
    bins   = np.array(obs_cfg["bins"], dtype=float)
    n_bins = len(bins) - 1
    logy   = obs_cfg["log_y"]
    label  = obs_cfg["label"]
    sidx   = STAGE_IDX[stage_key]

    # genFunnel branch name
    gf_branch_map = {
        "cosTheta":   f"GenFunnel_cosTheta_{stage_key}",
        "decayAngle": f"GenFunnel_decayAngle_{stage_key}",
        "pOverE":     f"GenFunnel_pOverE_{stage_key}",
        "dxySignif":  f"GenFunnel_dxySignif_{stage_key}",
        "mass":       f"GenFunnel_mass_{stage_key}",
    }
    sv_branch_map = {
        "cosTheta":   "StageVtx_cosTheta",
        "decayAngle": "StageVtx_decayAngle",
        "pOverE":     "StageVtx_pOverE",
        "dxySignif":  "StageVtx_dxySignif",
        "mass":       "StageVtx_mass",
    }

    cname  = f"reco_{stage_key}_{obs_key}"
    canvas = make_canvas(cname, f"{label} — {stage_key}", logy=logy)

    def make_h(tag):
        return ROOT.TH1F(f"h_{cname}_{tag}", "", n_bins, bins)

    h_gold   = make_h("gold")
    h_silver = make_h("silver")
    h_bronze = make_h("bronze")
    h_nonsig = make_h("nonsig")
    h_bkg    = make_h("bkg") if sv_bkg is not None else None

    # Fill from genFunnel (signal categories)
    if gf is not None and len(gf) > 0:
        has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
        vals_gf    = ak.to_numpy(ak.flatten(gf[gf_branch_map[obs_key]]))[has_tracks]
        gold       = ak.to_numpy(ak.flatten(gf[f"GenFunnel_gold_{stage_key}"]))[has_tracks].astype(bool)
        silver     = ak.to_numpy(ak.flatten(gf[f"GenFunnel_silver_{stage_key}"]))[has_tracks].astype(bool)
        bronze     = ak.to_numpy(ak.flatten(gf[f"GenFunnel_bronze_{stage_key}"]))[has_tracks].astype(bool)
        _fill_hist(h_gold,   vals_gf, gold)
        _fill_hist(h_silver, vals_gf, silver & ~gold)
        _fill_hist(h_bronze, vals_gf, bronze & ~silver & ~gold)

    # Fill non-signal from signal file's allStageVtx
    if sv_sig is not None and len(sv_sig) > 0:
        sv_stage  = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx
        sv_gold   = ak.to_numpy(sv_sig["StageVtx_isGold"]).astype(bool)
        sv_vals   = ak.to_numpy(sv_sig[sv_branch_map[obs_key]])
        _fill_hist(h_nonsig, sv_vals, sv_stage & ~sv_gold)

    # Fill background from dedicated file's allStageVtx
    if sv_bkg is not None and len(sv_bkg) > 0:
        sv_stage_b = ak.to_numpy(sv_bkg["StageVtx_stageIdx"]) == sidx
        sv_vals_b  = ak.to_numpy(sv_bkg[sv_branch_map[obs_key]])
        _fill_hist(h_bkg, sv_vals_b, sv_stage_b)

    # Normalise to unit area
    for h, col, lstyle in [
        (h_gold,   COLOR_GOLD,      1),
        (h_silver, COLOR_SILVER,    1),
        (h_bronze, COLOR_BRONZE,    1),
        (h_nonsig, COLOR_NONSIGNAL, 1),
    ]:
        _style_hist(h, col, lstyle)
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0 / integral)
    if h_bkg is not None:
        _style_hist(h_bkg, COLOR_BKG, 2)
        integral = h_bkg.Integral()
        if integral > 0:
            h_bkg.Scale(1.0 / integral)

    active_hists = [(h_gold,   "Gold"),
                    (h_silver, "Silver (excl.)"),
                    (h_bronze, "Bronze (excl.)"),
                    (h_nonsig, "Non-signal")]
    if h_bkg is not None:
        active_hists.append((h_bkg, "Background"))

    max_y = max((h.GetMaximum() for h, _ in active_hists if h.Integral() > 0), default=1.0)
    y_max = max_y * (5 if logy else 1.3)
    y_min = 1e-4

    h_ax = ROOT.TH1F(f"h_ax_{cname}", f";{label};Normalised to Unit Area", n_bins, bins)
    h_ax.SetMinimum(y_min)
    h_ax.SetMaximum(y_max)
    h_ax.GetXaxis().CenterTitle(True)
    h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045)
    h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.04)
    h_ax.GetYaxis().SetLabelSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.3)
    h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")

    for h, _ in active_hists:
        if h.Integral() > 0:
            h.Draw("HIST SAME")

    canvas._grid_lines = draw_axis_grid(h_ax, logy=logy)

    leg = ROOT.TLegend(0.62, 0.68, 0.88, 0.88)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)
    for h, lbl in active_hists:
        if h.Integral() > 0:
            leg.AddEntry(h, lbl, "l")
    leg.Draw()

    draw_cms_label()
    canvas.Update()
    tdir.cd()
    canvas.Write()

    canvas._h_ax    = h_ax
    canvas._hists   = [h for h, _ in active_hists]
    canvas._leg     = leg
