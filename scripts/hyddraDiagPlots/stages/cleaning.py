"""
stages/cleaning.py — cleaning-stage plots for hyddraDiagPlots.

Plots:
  - reco_cleaned_{cosTheta,decayAngle,pOverE,dxySignif}
  - cleaning_2d_{signal_tracks,contaminating_tracks,bkg_tracks}
  - cleaning_track_{compatibility,cos_theta}  (1D signal vs contaminating vs bkg)
  - cleaning_compat_cosTheta_gt_{wp}  (compatibility at cosTheta WPs: -0.5, 0.0, 0.5, 0.8)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import RECO_OBSERVABLES
from ..src.style   import draw_cms_label, make_canvas, draw_axis_grid
from ..src.plotter import plot_reco_observable


def plot_cleaning_2d(tdir, ct, max_compat=1.5, min_cos_theta=0.5,
                     use_diagonal_cut=False, clean_cut_slope=0.0):
    """
    2D (compatibility, cosTheta) for tracks in post-merge multi-track vertices.
    Three canvases: signal tracks, contaminating tracks, background-only tracks.
    """
    if ct is None or len(ak.flatten(ct["CleanTrack_compatibility"])) == 0:
        print("    [cleaning_2d] No cleaning track data — skipping")
        return

    compat     = ak.to_numpy(ak.flatten(ct["CleanTrack_compatibility"])).astype(float)
    cos_th     = ak.to_numpy(ak.flatten(ct["CleanTrack_cosTheta"])).astype(float)
    is_signal  = ak.to_numpy(ak.flatten(ct["CleanTrack_isSignal"])).astype(bool)
    vtx_gold   = ak.to_numpy(ak.flatten(ct["CleanTrack_vtxIsGold"])).astype(bool)
    vtx_silver = ak.to_numpy(ak.flatten(ct["CleanTrack_vtxIsSilver"])).astype(bool)

    sig_vtx              = vtx_gold | vtx_silver
    sig_track_in_sig_vtx = is_signal & sig_vtx
    bkg_track_in_sig_vtx = ~is_signal & sig_vtx
    bkg_vtx_tracks       = ~sig_vtx

    x_bins = np.linspace(0, 5, 101)
    y_bins = np.linspace(-1, 1, 101)

    def fill_h2(name, title, mask):
        ROOT.gStyle.SetPalette(ROOT.kViridis)
        h = ROOT.TH2F(name, title, len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)
        h.SetDirectory(0)
        for c, ct_v in zip(compat[mask], cos_th[mask]):
            h.Fill(c, ct_v)
        h.GetXaxis().SetTitle("Track Compatibility (#sigma)")
        h.GetYaxis().SetTitle("Track cos#theta")
        h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.3); h.GetYaxis().SetTitleOffset(1.3)
        h.SetStats(0)
        return h

    h_sig  = fill_h2("h2_sig_track",  "Signal tracks in signal vtx",
                     sig_track_in_sig_vtx)
    h_cont = fill_h2("h2_cont_track", "Bkg tracks in signal vtx",
                     bkg_track_in_sig_vtx)
    h_bkg  = fill_h2("h2_bkg_track",  "All tracks in bkg vtx",
                     bkg_vtx_tracks)

    def draw_cut_lines():
        lines = []
        if use_diagonal_cut:
            ln = ROOT.TLine(x_bins[0],
                            clean_cut_slope * x_bins[0] + min_cos_theta,
                            x_bins[-1],
                            clean_cut_slope * x_bins[-1] + min_cos_theta)
            ln.SetLineColor(ROOT.kRed); ln.SetLineWidth(2); ln.SetLineStyle(2)
            ln.Draw(); lines.append(ln)
        else:
            ln1 = ROOT.TLine(max_compat, y_bins[0], max_compat, y_bins[-1])
            ln1.SetLineColor(ROOT.kRed); ln1.SetLineWidth(2); ln1.SetLineStyle(2)
            ln1.Draw(); lines.append(ln1)
            ln2 = ROOT.TLine(x_bins[0], min_cos_theta, x_bins[-1], min_cos_theta)
            ln2.SetLineColor(ROOT.kRed); ln2.SetLineWidth(2); ln2.SetLineStyle(2)
            ln2.Draw(); lines.append(ln2)
        return lines

    for h, suffix, title in [
        (h_sig,  "signal_tracks",        "Signal tracks in signal vertices"),
        (h_cont, "contaminating_tracks", "Background tracks in signal vertices"),
        (h_bkg,  "bkg_tracks",           "All tracks in background vertices"),
    ]:
        c = ROOT.TCanvas(f"cleaning_2d_{suffix}", title, 800, 600)
        c.SetLeftMargin(0.16)
        c.SetRightMargin(0.18)
        c.SetBottomMargin(0.12)
        c.SetTopMargin(0.10)
        h.Draw("COLZ")
        c.Update()
        lines = draw_cut_lines()

        leg = ROOT.TLegend(0.18, 0.75, 0.55, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0)
        l1 = ROOT.TLine()
        l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.SetLineWidth(2)
        leg.AddEntry(l1, "Cleaning cuts", "l")
        leg.Draw()
        draw_cms_label()
        c.Update()

        pal = h.GetListOfFunctions().FindObject("palette")
        if pal:
            h.GetListOfFunctions().Remove(pal)

        c._h = h; c._lines = lines; c._leg = leg
        tdir.cd()
        c.Write()
        print(f"    [cleaning_2d_{suffix}] done")


def plot_cleaning_track_distributions(tdir, ct, cos_theta_min=None):
    """
    1D normalised distributions of compatibility and cosTheta for:
      - signal tracks in signal-matched vertices
      - contaminating tracks (non-signal in signal-matched vertices)
      - all tracks in background-only vertices

    If cos_theta_min is given, only tracks with cosTheta > cos_theta_min are
    included and only the compatibility plot is produced (the cosTheta plot is
    skipped since we are cutting on it).  A TLatex label showing the applied
    cut is drawn on the canvas.
    """
    if ct is None or len(ak.flatten(ct["CleanTrack_compatibility"])) == 0:
        print("    [cleaning_track_1d] No cleaning track data — skipping")
        return

    compat     = ak.to_numpy(ak.flatten(ct["CleanTrack_compatibility"])).astype(float)
    cos_th     = ak.to_numpy(ak.flatten(ct["CleanTrack_cosTheta"])).astype(float)
    is_signal  = ak.to_numpy(ak.flatten(ct["CleanTrack_isSignal"])).astype(bool)
    vtx_gold   = ak.to_numpy(ak.flatten(ct["CleanTrack_vtxIsGold"])).astype(bool)
    vtx_silver = ak.to_numpy(ak.flatten(ct["CleanTrack_vtxIsSilver"])).astype(bool)

    sig_vtx  = vtx_gold | vtx_silver
    m_sig    = is_signal  & sig_vtx
    m_cont   = ~is_signal & sig_vtx
    m_bkg    = ~sig_vtx

    # Keep uncut masks so WP plots can normalize to the full population
    m_sig_base, m_cont_base, m_bkg_base = m_sig, m_cont, m_bkg

    if cos_theta_min is not None:
        cos_cut = cos_th > cos_theta_min
        m_sig   = m_sig  & cos_cut
        m_cont  = m_cont & cos_cut
        m_bkg   = m_bkg  & cos_cut

    cat_colors = [(m_sig,  ROOT.kBlue + 2, 1, "Signal track in sig. vtx",  m_sig_base),
                  (m_cont, ROOT.kOrange+ 1, 1, "Contam. track in sig. vtx", m_cont_base),
                  (m_bkg,  ROOT.kRed  + 2, 2, "Track in bkg-only vtx",     m_bkg_base)]

    wp_str = (f"_cosTheta_gt_{f'{cos_theta_min:.1f}'.replace('-','m').replace('.','p')}"
              if cos_theta_min is not None else "")

    var_list = [(compat, np.linspace(0, 5, 51), "Track Compatibility (#sigma)",
                 f"cleaning_compat_1d{wp_str}")]
    if cos_theta_min is None:
        var_list.append((cos_th, np.linspace(-1, 1, 51), "Track cos#theta",
                         "cleaning_costheta_1d"))

    y_label = ("Fraction of tracks (no cut)" if cos_theta_min is not None
               else "Normalised to Unit Area")

    for var_vals, bins, x_label, cname in var_list:
        n_bins = len(bins) - 1
        canvas = make_canvas(cname, x_label, logy=False)
        hists  = []
        for mask, col, lstyle, lbl, base_mask in cat_colors:
            h = ROOT.TH1F(f"h_{cname}_{lbl[:4]}", "", n_bins, bins)
            for v in var_vals[mask]:
                h.Fill(float(v))
            h.SetLineColor(col); h.SetLineWidth(2); h.SetLineStyle(lstyle)
            h.SetFillStyle(0); h.SetStats(0)
            denom = int(np.sum(base_mask)) if cos_theta_min is not None else h.Integral()
            if denom > 0:
                h.Scale(1.0 / denom)
            hists.append((h, lbl))

        max_y = max((h.GetMaximum() for h, _ in hists if h.Integral() > 0), default=1.0)
        h_ax = ROOT.TH1F(f"h_ax_{cname}",
                         f";{x_label};{y_label}", n_bins, bins)
        h_ax.SetMinimum(1e-4); h_ax.SetMaximum(max_y * 1.3)
        h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
        h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
        h_ax.GetXaxis().SetLabelSize(0.04);  h_ax.GetYaxis().SetLabelSize(0.04)
        h_ax.GetXaxis().SetTitleOffset(1.3); h_ax.GetYaxis().SetTitleOffset(1.3)
        h_ax.SetStats(0)
        h_ax.Draw("AXIS")
        for h, _ in hists:
            if h.Integral() > 0:
                h.Draw("HIST SAME")

        leg = ROOT.TLegend(0.18, 0.72, 0.72, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.03)
        for h, lbl in hists:
            if h.Integral() > 0:
                leg.AddEntry(h, lbl, "l")
        leg.Draw()
        canvas._grid_lines = draw_axis_grid(h_ax, logy=False)

        if cos_theta_min is not None:
            lat = ROOT.TLatex()
            lat.SetNDC(); lat.SetTextFont(42)
            lat.SetTextSize(0.038); lat.SetTextAlign(11)
            lat.DrawLatex(0.18, 0.67, f"cos#theta > {cos_theta_min}")
            canvas._lat = lat

        draw_cms_label()
        canvas.Update()
        tdir.cd()
        canvas.Write()
        canvas._h_ax = h_ax; canvas._hists = hists; canvas._leg = leg
        print(f"    [{cname}] done")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg, ct, max_compat=1.5, min_cos_theta=0.5,
               use_diagonal_cut=False, clean_cut_slope=0.0):
    print("  [cleaning] Reco observables...")
    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "cleaned", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_cleaned_{obs_key}] done")
    print("  [cleaning] 2D cleaning variables...")
    plot_cleaning_2d(tdir, ct, max_compat, min_cos_theta, use_diagonal_cut, clean_cut_slope)
    print("  [cleaning] 1D cleaning track distributions...")
    plot_cleaning_track_distributions(tdir, ct)
    print("  [cleaning] Compatibility vs cosTheta working points...")
    for wp in [-0.5, 0.0, 0.5, 0.8]:
        plot_cleaning_track_distributions(tdir, ct, cos_theta_min=wp)
