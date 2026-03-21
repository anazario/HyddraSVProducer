"""
plotter.py — reusable plot helpers for hyddraDiagPlots.
"""
import numpy as np
import awkward as ak
import ROOT

from .config import (
    STAGE_NAMES, STAGE_KEYS, STAGE_IDX,
    COLOR_GOLD, COLOR_SILVER, COLOR_BRONZE, COLOR_NONSIGNAL, COLOR_BKG,
    COLORS_STAGE, MARKERS, GEN_DXY_BINS, HAD_MIN3D_CUT,
)
from .style import make_canvas, draw_cms_label, draw_axis_grid


# ── Hadronic signal masks ─────────────────────────────────────────────────────

def had_signal_mask(data, key=None, loose=True):
    """
    Boolean signal mask for hadronic data.  Works for both data sources:
      key=None  → allStageVtx (StageVtx_matchRatio / StageVtx_min3D, flat arrays)
      key=<str> → GenFunnel   (GenFunnel_matchRatio_{key} / GenFunnel_min3D_{key}, jagged)

    loose=True  (default): matchRatio > 0
    loose=False           : matchRatio >= 0.5

    Gates on min3D < HAD_MIN3D_CUT when the branch exists; falls back to
    matchRatio threshold only for old ntuples.
    """
    jagged      = key is not None
    mr_field    = f"GenFunnel_matchRatio_{key}" if jagged else "StageVtx_matchRatio"
    min3d_field = f"GenFunnel_min3D_{key}"      if jagged else "StageVtx_min3D"

    def _flat(arr):
        return ak.to_numpy(ak.flatten(arr)).astype(float) if jagged else ak.to_numpy(arr).astype(float)

    mr   = _flat(data[mr_field])
    mask = (mr > 0) if loose else (mr >= 0.5)
    if min3d_field in data.fields:
        min3d = _flat(data[min3d_field])
        mask  = mask & (min3d >= 0) & (min3d < HAD_MIN3D_CUT)
    return mask


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


def plot_reco_observable(tdir, gf, sv_sig, stage_key, obs_key, obs_cfg, sv_bkg=None,
                         hadronic=False):
    """
    Normalised shape distributions on one canvas.

    Leptonic (hadronic=False):
      - Gold         : genFunnel vertices gold-matched at this stage
      - Silver excl. : silver but not gold
      - Bronze excl. : bronze but not silver
      - Non-signal   : sv_sig allStageVtx where stageIdx==stage and isGold=False

    Hadronic (hadronic=True):
      - Signal tight (matchRatio >= 0.5)
      - Signal loose excl. (0 < matchRatio < 0.5)
      - Non-signal   : sv_sig allStageVtx where matchRatio <= 0

    All normalised to unit area.  -1 sentinels are skipped.
    """
    bins   = np.array(obs_cfg["bins"], dtype=float)
    n_bins = len(bins) - 1
    logy   = obs_cfg["log_y"]
    label  = obs_cfg["label"]
    sidx   = STAGE_IDX[stage_key]

    gf_branch_map = {
        "cosTheta":   f"GenFunnel_cosTheta_{stage_key}",
        "decayAngle": f"GenFunnel_decayAngle_{stage_key}",
        "pOverE":     f"GenFunnel_pOverE_{stage_key}",
        "dxySignif":  f"GenFunnel_dxySignif_{stage_key}",
        "mass":       f"GenFunnel_mass_{stage_key}",
        "nTracks":    f"GenFunnel_nTracks_{stage_key}",
    }
    sv_branch_map = {
        "cosTheta":   "StageVtx_cosTheta",
        "decayAngle": "StageVtx_decayAngle",
        "pOverE":     "StageVtx_pOverE",
        "dxySignif":  "StageVtx_dxySignif",
        "mass":       "StageVtx_mass",
        "nTracks":    "StageVtx_nTracks",
    }

    cname  = f"reco_{stage_key}_{obs_key}"
    canvas = make_canvas(cname, f"{label} — {stage_key}", logy=logy)

    def make_h(tag):
        return ROOT.TH1F(f"h_{cname}_{tag}", "", n_bins, bins)

    h_bkg = make_h("bkg") if sv_bkg is not None else None

    if hadronic:
        h_tight  = make_h("tight")
        h_loose  = make_h("loose")
        h_nonsig = make_h("nonsig")

        if gf is not None and len(gf) > 0:
            has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_isHadronic"])).astype(bool)
            vals_gf    = ak.to_numpy(ak.flatten(gf[gf_branch_map[obs_key]]))[has_tracks]
            mr         = ak.to_numpy(ak.flatten(gf[f"GenFunnel_matchRatio_{stage_key}"]))[has_tracks]
            valid      = mr >= 0  # -1 is sentinel for no match
            _fill_hist(h_tight, vals_gf, valid & (mr >= 0.5))
            _fill_hist(h_loose, vals_gf, valid & (mr > 0) & (mr < 0.5))

        if sv_sig is not None and len(sv_sig) > 0:
            sv_stage = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx
            sv_vals  = ak.to_numpy(sv_sig[sv_branch_map[obs_key]])
            _fill_hist(h_nonsig, sv_vals, sv_stage & ~had_signal_mask(sv_sig))

        if sv_bkg is not None and len(sv_bkg) > 0:
            sv_stage_b = ak.to_numpy(sv_bkg["StageVtx_stageIdx"]) == sidx
            sv_vals_b  = ak.to_numpy(sv_bkg[sv_branch_map[obs_key]])
            _fill_hist(h_bkg, sv_vals_b, sv_stage_b)

        for h, col, lstyle in [
            (h_tight,  COLOR_GOLD,      1),
            (h_loose,  COLOR_SILVER,    1),
            (h_nonsig, COLOR_NONSIGNAL, 1),
        ]:
            _style_hist(h, col, lstyle)
            integral = h.Integral()
            if integral > 0:
                h.Scale(1.0 / integral)

        active_hists = [
            (h_tight,  "Signal (matchRatio #geq 0.5)"),
            (h_loose,  "Signal (0 < matchRatio < 0.5)"),
            (h_nonsig, "Non-signal"),
        ]

    else:
        h_gold   = make_h("gold")
        h_silver = make_h("silver")
        h_bronze = make_h("bronze")
        h_nonsig = make_h("nonsig")

        if gf is not None and len(gf) > 0:
            has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
            vals_gf    = ak.to_numpy(ak.flatten(gf[gf_branch_map[obs_key]]))[has_tracks]
            gold       = ak.to_numpy(ak.flatten(gf[f"GenFunnel_gold_{stage_key}"]))[has_tracks].astype(bool)
            silver     = ak.to_numpy(ak.flatten(gf[f"GenFunnel_silver_{stage_key}"]))[has_tracks].astype(bool)
            bronze     = ak.to_numpy(ak.flatten(gf[f"GenFunnel_bronze_{stage_key}"]))[has_tracks].astype(bool)
            _fill_hist(h_gold,   vals_gf, gold)
            _fill_hist(h_silver, vals_gf, silver & ~gold)
            _fill_hist(h_bronze, vals_gf, bronze & ~silver & ~gold)

        if sv_sig is not None and len(sv_sig) > 0:
            sv_stage  = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx
            sv_gold   = ak.to_numpy(sv_sig["StageVtx_isGold"]).astype(bool)
            sv_vals   = ak.to_numpy(sv_sig[sv_branch_map[obs_key]])
            _fill_hist(h_nonsig, sv_vals, sv_stage & ~sv_gold)

        if sv_bkg is not None and len(sv_bkg) > 0:
            sv_stage_b = ak.to_numpy(sv_bkg["StageVtx_stageIdx"]) == sidx
            sv_vals_b  = ak.to_numpy(sv_bkg[sv_branch_map[obs_key]])
            _fill_hist(h_bkg, sv_vals_b, sv_stage_b)

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

        active_hists = [
            (h_gold,   "Gold"),
            (h_silver, "Silver (excl.)"),
            (h_bronze, "Bronze (excl.)"),
            (h_nonsig, "Non-signal"),
        ]

    if h_bkg is not None:
        _style_hist(h_bkg, COLOR_BKG, 2)
        integral = h_bkg.Integral()
        if integral > 0:
            h_bkg.Scale(1.0 / integral)
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

    draw_cms_label("Hadronic HYDDRA" if hadronic else "Leptonic HYDDRA")
    canvas.Update()
    tdir.cd()
    canvas.Write()

    canvas._h_ax    = h_ax
    canvas._hists   = [h for h, _ in active_hists]
    canvas._leg     = leg


# ── Fakes vs reco dxy (shared, called from leptonic and hadronic ID stages) ──

def plot_fakes_vs_dxy(tdir, sv_sig, is_signal, cms_label="", tag="fakes_vs_dxy"):
    """
    Fake fraction (non-signal SVs / all SVs) vs reconstructed dxy, one curve
    per algorithm stage.  Complementary to the per-stage eff_vs_dxy plots.

    is_signal : full-length boolean array aligned with sv_sig rows.
                True = signal (gold for leptonic, matchRatio>0 for hadronic).
    """
    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping"); return
    if "StageVtx_x" not in sv_sig.fields:
        print(f"    [{tag}] StageVtx_x not in ntuple — skipping"); return

    bins   = np.array(GEN_DXY_BINS, dtype=float)
    n_bins = len(bins) - 1

    stage = ak.to_numpy(sv_sig["StageVtx_stageIdx"]).astype(int)
    dxy   = np.sqrt(ak.to_numpy(sv_sig["StageVtx_x"]).astype(float)**2 +
                    ak.to_numpy(sv_sig["StageVtx_y"]).astype(float)**2)

    graphs, hists = [], []
    for i, s in enumerate(STAGE_KEYS):
        mask = stage == i
        h_denom = ROOT.TH1F(f"h_denom_{tag}_{s}", "", n_bins, bins)
        h_numer = ROOT.TH1F(f"h_numer_{tag}_{s}", "", n_bins, bins)
        h_denom.SetDirectory(0); h_numer.SetDirectory(0)
        for v in dxy[mask]:
            h_denom.Fill(v)
        for v in dxy[mask & ~is_signal]:
            h_numer.Fill(v)
        if h_denom.Integral() == 0:
            hists.append((h_denom, h_numer, None))
            continue
        eff = ROOT.TEfficiency(h_numer, h_denom)
        eff.SetStatisticOption(ROOT.TEfficiency.kFCP)
        g = eff.CreateGraph()
        g.SetTitle(STAGE_NAMES[i])
        g.SetMarkerColor(COLORS_STAGE[i]); g.SetLineColor(COLORS_STAGE[i])
        g.SetMarkerStyle(MARKERS[i]); g.SetLineWidth(2)
        graphs.append(g)
        hists.append((h_denom, h_numer, eff))

    c = ROOT.TCanvas(tag, tag, 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)

    h_ax = ROOT.TH1F(f"h_ax_{tag}", ";Reco dxy (cm);Fake fraction at ID stage",
                     n_bins, bins)
    h_ax.SetMinimum(0.0); h_ax.SetMaximum(1.1)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.GetXaxis().SetLabelSize(0.04);  h_ax.GetYaxis().SetLabelSize(0.04)
    h_ax.GetXaxis().SetTitleOffset(1.2); h_ax.GetYaxis().SetTitleOffset(1.3)
    h_ax.SetStats(0)
    h_ax.Draw("AXIS")

    for g in graphs:
        g.Draw("P SAME")

    c._grid_lines = draw_axis_grid(h_ax, logy=False)
    add_legend(c, graphs, x1=0.55, y1=0.55, x2=0.88, y2=0.88)
    draw_cms_label(cms_label)
    c.Update()

    c._h_ax = h_ax; c._hists = hists; c._graphs = graphs
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


# ── Shared geometry-cut line drawers ─────────────────────────────────────────

def geom_cut_lines_z(c, x_bins, y_bins):
    """Draw KUCMSSkimmer baseline cut boundaries on a dxy-vs-z canvas."""
    lines = []
    for zval in [-15.0, 15.0]:
        ln = ROOT.TLine(zval, 0.0, zval, 15.0)
        ln.SetLineColor(ROOT.kRed); ln.SetLineWidth(2); ln.SetLineStyle(2)
        ln.Draw(); lines.append(ln)
    ln2 = ROOT.TLine(x_bins[0], 15.0, x_bins[-1], 15.0)
    ln2.SetLineColor(ROOT.kRed); ln2.SetLineWidth(2); ln2.SetLineStyle(2)
    ln2.Draw(); lines.append(ln2)
    ln3 = ROOT.TLine(x_bins[0], 2.0, x_bins[-1], 2.0)
    ln3.SetLineColor(ROOT.kOrange+1); ln3.SetLineWidth(2); ln3.SetLineStyle(3)
    ln3.Draw(); lines.append(ln3)
    c._cut_lines = lines


def geom_cut_lines_eta(c, x_bins, y_bins):
    """Draw KUCMSSkimmer baseline cut boundaries on a dxy-vs-eta_origin canvas."""
    lines = []
    for etaval in [-1.0, 1.0]:
        ln = ROOT.TLine(etaval, 15.0, etaval, y_bins[-1])
        ln.SetLineColor(ROOT.kRed); ln.SetLineWidth(2); ln.SetLineStyle(2)
        ln.Draw(); lines.append(ln)
    ln2 = ROOT.TLine(x_bins[0], 15.0, x_bins[-1], 15.0)
    ln2.SetLineColor(ROOT.kRed); ln2.SetLineWidth(2); ln2.SetLineStyle(2)
    ln2.Draw(); lines.append(ln2)
    ln3 = ROOT.TLine(x_bins[0], 2.0, x_bins[-1], 2.0)
    ln3.SetLineColor(ROOT.kOrange+1); ln3.SetLineWidth(2); ln3.SetLineStyle(3)
    ln3.Draw(); lines.append(ln3)
    c._cut_lines = lines


# ── Shared ID-stage position plots ───────────────────────────────────────────

def _make_th1_line(name, vals, mask, bins, color, lstyle, x_title):
    h = ROOT.TH1F(name, name, len(bins) - 1, bins)
    h.SetDirectory(0)
    for v in vals[mask]:
        h.Fill(float(v))
    h.SetLineColor(color); h.SetLineWidth(2); h.SetLineStyle(lstyle)
    h.SetFillStyle(0); h.SetStats(0)
    h.GetXaxis().SetTitle(x_title)
    h.GetYaxis().SetTitle("Vertices")
    h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
    h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
    h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
    h.GetXaxis().SetTitleOffset(1.2); h.GetYaxis().SetTitleOffset(1.3)
    return h


def plot_costheta_zoom(tdir, sv_sig, sv_bkg, sidx, is_sig, tag,
                       sig_label="Signal", cms_label="", cut_val=0.995):
    """
    cosTheta distribution zoomed into [0.90, 1.002] with fine binning.
    is_sig: full-length bool array aligned with sv_sig rows.
    """
    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping"); return

    bins    = np.linspace(0.90, 1.002, 52)
    stage   = ak.to_numpy(sv_sig["StageVtx_stageIdx"]).astype(int)
    ct      = ak.to_numpy(sv_sig["StageVtx_cosTheta"]).astype(float)
    mask_id = stage == sidx

    h_sig    = _make_th1_line(f"h_{tag}_sig",    ct, mask_id &  is_sig, bins, COLOR_GOLD,      1, "cos#theta (wrt PV)")
    h_nonsig = _make_th1_line(f"h_{tag}_nonsig", ct, mask_id & ~is_sig, bins, COLOR_NONSIGNAL, 1, "cos#theta (wrt PV)")
    hists  = [h_sig, h_nonsig]
    labels = [sig_label, "Non-signal"]

    if sv_bkg is not None and len(sv_bkg) > 0:
        bkg_stage = ak.to_numpy(sv_bkg["StageVtx_stageIdx"]).astype(int)
        bkg_ct    = ak.to_numpy(sv_bkg["StageVtx_cosTheta"]).astype(float)
        h_bkg = _make_th1_line(f"h_{tag}_bkg", bkg_ct, bkg_stage == sidx, bins, COLOR_BKG, 1, "cos#theta (wrt PV)")
        hists.append(h_bkg); labels.append("Background")

    max_y = max((h.GetMaximum() for h in hists), default=1.0)
    c = ROOT.TCanvas(tag, tag, 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10); c.SetLogy(True)
    for i, h in enumerate(hists):
        h.SetMaximum(max_y * 5); h.SetMinimum(0.5)
        h.Draw("HIST" if i == 0 else "HIST SAME")

    line = ROOT.TLine(cut_val, 0.5, cut_val, max_y * 5)
    line.SetLineColor(ROOT.kBlack); line.SetLineWidth(2); line.SetLineStyle(2)
    line.Draw()
    c._grid_lines = draw_axis_grid(hists[0], logy=True)

    leg = ROOT.TLegend(0.15, 0.70, 0.55, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.033)
    for h, lbl in zip(hists, labels):
        leg.AddEntry(h, lbl, "l")
    l1 = ROOT.TLine()
    l1.SetLineColor(ROOT.kBlack); l1.SetLineStyle(2); l1.SetLineWidth(2)
    leg.AddEntry(l1, f"Skimmer cut: cosTheta < {cut_val}", "l")
    leg.Draw()
    draw_cms_label(cms_label); c.Update()

    c._hists = hists; c._leg = leg; c._line = line
    tdir.cd(); c.Write()
    print(f"    [{tag}] done")


def plot_dxy_1d(tdir, sv_sig, sv_bkg, sidx, is_sig, tag,
                sig_label="Signal", cms_label=""):
    """
    1D dxy = sqrt(x^2+y^2) distribution at a given stage, split by signal/non-signal.
    is_sig: full-length bool array aligned with sv_sig rows.
    """
    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping"); return
    if "StageVtx_x" not in sv_sig.fields:
        print(f"    [{tag}] StageVtx_x not in ntuple — skipping"); return

    bins    = np.linspace(0, 80, 41)
    stage   = ak.to_numpy(sv_sig["StageVtx_stageIdx"]).astype(int)
    dxy     = np.sqrt(ak.to_numpy(sv_sig["StageVtx_x"]).astype(float)**2 +
                      ak.to_numpy(sv_sig["StageVtx_y"]).astype(float)**2)
    mask_id = stage == sidx

    h_sig    = _make_th1_line(f"h_{tag}_sig",    dxy, mask_id &  is_sig, bins, COLOR_GOLD,      1, "dxy (cm)")
    h_nonsig = _make_th1_line(f"h_{tag}_nonsig", dxy, mask_id & ~is_sig, bins, COLOR_NONSIGNAL, 1, "dxy (cm)")
    hists  = [h_sig, h_nonsig]
    labels = [sig_label, "Non-signal"]

    if sv_bkg is not None and len(sv_bkg) > 0 and "StageVtx_x" in sv_bkg.fields:
        bkg_stage = ak.to_numpy(sv_bkg["StageVtx_stageIdx"]).astype(int)
        bkg_dxy   = np.sqrt(ak.to_numpy(sv_bkg["StageVtx_x"]).astype(float)**2 +
                            ak.to_numpy(sv_bkg["StageVtx_y"]).astype(float)**2)
        h_bkg = _make_th1_line(f"h_{tag}_bkg", bkg_dxy, bkg_stage == sidx, bins, COLOR_BKG, 1, "dxy (cm)")
        hists.append(h_bkg); labels.append("Background")

    max_y = max((h.GetMaximum() for h in hists), default=1.0)
    c = ROOT.TCanvas(tag, tag, 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10); c.SetLogy(True)
    for i, h in enumerate(hists):
        h.SetMaximum(max_y * 5); h.SetMinimum(0.5)
        h.Draw("HIST" if i == 0 else "HIST SAME")
    c._grid_lines = draw_axis_grid(hists[0], logy=True)

    leg = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.035)
    for h, lbl in zip(hists, labels):
        leg.AddEntry(h, lbl, "l")
    leg.Draw()
    draw_cms_label(cms_label); c.Update()

    c._hists = hists; c._leg = leg
    tdir.cd(); c.Write()
    print(f"    [{tag}] done")


def plot_2d_dxy_vs_x(tdir, x_vals, dxy_vals, is_sig, x_title, x_bins, y_bins,
                     canvas_prefix, sig_label="Signal",
                     cut_lines_fn=None, cms_label=""):
    """
    Two COLZ TH2F canvases (signal / non-signal) for dxy (y) vs x_vals (x).
    All arrays must already be filtered to the desired stage.
    cut_lines_fn(c, x_bins, y_bins) is called after Draw to overlay cut lines.
    """
    x_bins_a = np.array(x_bins, dtype=float)
    y_bins_a = np.array(y_bins, dtype=float)
    ROOT.gStyle.SetPalette(ROOT.kViridis)

    for cname, mask, title in [
        (f"{canvas_prefix}_signal",    is_sig,  sig_label),
        (f"{canvas_prefix}_nonsignal", ~is_sig, "Non-signal"),
    ]:
        h = ROOT.TH2F(cname, cname,
                      len(x_bins_a) - 1, x_bins_a,
                      len(y_bins_a) - 1, y_bins_a)
        h.SetDirectory(0)
        for xv, yv in zip(x_vals[mask], dxy_vals[mask]):
            h.Fill(float(xv), float(yv))
        h.GetXaxis().SetTitle(x_title)
        h.GetYaxis().SetTitle("dxy (cm)")
        h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.3); h.GetYaxis().SetTitleOffset(1.3)
        h.SetStats(0)

        c = ROOT.TCanvas(cname, f"{title} — {x_title} vs dxy", 800, 640)
        c.SetLeftMargin(0.14); c.SetRightMargin(0.18)
        c.SetBottomMargin(0.12); c.SetTopMargin(0.10)
        h.Draw("COLZ"); c.Update()

        if cut_lines_fn:
            cut_lines_fn(c, x_bins_a, y_bins_a)

        draw_cms_label(cms_label); c.Update()
        pal = h.GetListOfFunctions().FindObject("palette")
        if pal:
            h.GetListOfFunctions().Remove(pal)

        c._h = h
        tdir.cd(); c.Write()
        print(f"    [{cname}] done")
