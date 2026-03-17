"""
stages/hadronic_disambiguation.py — disambiguation-stage plots for hadronic HYDDRA.

Hadronic disambiguation keeps the largest vertex (most tracks) when two vertices
share a track.  The key plot here is the 2D (nTracks, nSignalTracks) scatter
at both pre- and post-disambiguation, with an optional third version restricted
to vertices that share a track with at least one other vertex ("non-unique").

Plots:
  - reco_disambig_{cosTheta,decayAngle,pOverE,dxySignif,mass,nTracks}
  - disambiguation_purity_{pre,post,shared}
      2D: nTracks (x) vs nSignalTracks (y) at the given stage subset.
      Two working-point lines are drawn:
        * y = 1       (loose: >=1 signal track)
        * y = 0.5 * x (tight: >=50% purity diagonal)
      Background vertices (nSignalTracks=0) and signal-contaminated vertices
      are shown on the same 2D canvas using separate colour palettes drawn
      as a hatched 1D profile overlaid on the 2D scatter.
      In practice we use a TH2 for all vertices and overlay a TH2 for the
      signal-contaminated population drawn with a contrasting palette.
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import HADRONIC_RECO_OBSERVABLES, STAGE_IDX
from ..src.style   import draw_cms_label, draw_colz_grid
from ..src.plotter import plot_reco_observable


def _plot_purity_2d(tdir, sv_sig, stage_key, title_suffix, shared_only=False):
    """
    2D histogram: x = StageVtx_nTracks, y = StageVtx_nSignalTracks.
    Restricted to the given stage.  If shared_only=True, further restricted to
    vertices with StageVtx_hasSharedTrack == True (non-unique, pre-disambig).

    Two diagonal lines mark the working points:
      y = 1          (loose: at least 1 signal track)
      y = 0.5 * x   (tight: >= 50% of tracks are signal)
    """
    tag = f"disambiguation_purity_{title_suffix}"

    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping")
        return
    if "StageVtx_nTracks" not in sv_sig.fields:
        print(f"    [{tag}] StageVtx_nTracks not in ntuple — skipping "
              "(regenerate with hadronic analyzer)")
        return

    sidx    = STAGE_IDX[stage_key]
    st_mask = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx

    if shared_only:
        if "StageVtx_hasSharedTrack" not in sv_sig.fields:
            print(f"    [{tag}] StageVtx_hasSharedTrack not available — skipping")
            return
        shared_mask = ak.to_numpy(sv_sig["StageVtx_hasSharedTrack"]).astype(bool)
        base_mask = st_mask & shared_mask
    else:
        base_mask = st_mask

    if not np.any(base_mask):
        print(f"    [{tag}] No vertices in selection — skipping")
        return

    n_trk = ak.to_numpy(sv_sig["StageVtx_nTracks"     ][base_mask]).astype(int)
    n_sig = ak.to_numpy(sv_sig["StageVtx_nSignalTracks"][base_mask]).astype(int)

    # Clip nSignalTracks to [0, nTracks] to guard against any sentinel values
    n_sig = np.clip(n_sig, 0, n_trk)

    n_max = max(int(n_trk.max()) if len(n_trk) > 0 else 30, 10)
    x_bins = np.arange(0.5, n_max + 1.5, 1.0)
    y_bins = np.arange(-0.5, n_max + 1.5, 1.0)

    ROOT.gStyle.SetPalette(ROOT.kViridis)
    h_all = ROOT.TH2F(f"h2_{tag}_all", "", len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)
    h_sig_cont = ROOT.TH2F(f"h2_{tag}_sig", "", len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)

    for nt, ns in zip(n_trk, n_sig):
        h_all.Fill(float(nt), float(ns))

    # Signal-contaminated = nSignalTracks >= 1
    sig_cont_mask = n_sig >= 1
    for nt, ns in zip(n_trk[sig_cont_mask], n_sig[sig_cont_mask]):
        h_sig_cont.Fill(float(nt), float(ns))

    def _style_h2(h, xtitle, ytitle):
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.GetXaxis().CenterTitle(True); h.GetYaxis().CenterTitle(True)
        h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.04);  h.GetYaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.3); h.GetYaxis().SetTitleOffset(1.3)
        h.SetStats(0)

    _style_h2(h_all,      "Number of tracks",         "Number of signal tracks")
    _style_h2(h_sig_cont, "Number of tracks (signal-contaminated)", "Number of signal tracks")

    # ── Draw the two canvases: all vertices and signal-contaminated only ──────
    for h, suffix, canvas_title in [
        (h_all,      "all",          f"All vertices ({title_suffix})"),
        (h_sig_cont, "sig_cont",     f"Signal-contaminated vertices ({title_suffix})"),
    ]:
        cname = f"{tag}_{suffix}"
        c = ROOT.TCanvas(cname, canvas_title, 800, 700)
        c.SetLeftMargin(0.14); c.SetRightMargin(0.18)
        c.SetBottomMargin(0.12); c.SetTopMargin(0.10)

        h.Draw("COLZ")
        c.Update()
        c._grid_lines = draw_colz_grid(h)

        # Draw working-point lines
        lines = []
        x_max = float(x_bins[-1])

        # Loose: y = 1 (horizontal)
        ln_loose = ROOT.TLine(x_bins[0], 1.0, x_max, 1.0)
        ln_loose.SetLineColor(ROOT.kRed); ln_loose.SetLineWidth(2); ln_loose.SetLineStyle(2)
        ln_loose.Draw(); lines.append(ln_loose)

        # Tight: y = 0.5 * x (diagonal from origin)
        y_at_xmax = 0.5 * x_max
        y_lim = float(y_bins[-1])
        if y_at_xmax > y_lim:
            x_at_ylim = y_lim / 0.5
            ln_tight = ROOT.TLine(x_bins[0], 0.0, x_at_ylim, y_lim)
        else:
            ln_tight = ROOT.TLine(x_bins[0], 0.0, x_max, y_at_xmax)
        ln_tight.SetLineColor(ROOT.kOrange + 7); ln_tight.SetLineWidth(2); ln_tight.SetLineStyle(2)
        ln_tight.Draw(); lines.append(ln_tight)

        leg = ROOT.TLegend(0.15, 0.75, 0.58, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.03)
        leg.AddEntry(ln_loose, "Loose: #geq 1 signal track",     "l")
        leg.AddEntry(ln_tight, "Tight: #geq 50% signal tracks",  "l")
        leg.Draw()

        pal = h.GetListOfFunctions().FindObject("palette")
        if pal:
            h.GetListOfFunctions().Remove(pal)

        draw_cms_label("Hadronic HYDDRA")
        c.Update()
        c._h = h; c._lines = lines; c._leg = leg
        tdir.cd()
        c.Write()
        print(f"    [{cname}] done")


# ── Orchestrator ───────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg):
    print("  [had/disambiguation] Reco observables...")
    for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "disambig", obs_key, obs_cfg, sv_bkg, hadronic=True)
        print(f"    [reco_disambig_{obs_key}] done")

    print("  [had/disambiguation] Purity 2D: pre-disambiguation (all)...")
    _plot_purity_2d(tdir, sv_sig, "cleaned",  "pre",    shared_only=False)

    print("  [had/disambiguation] Purity 2D: pre-disambiguation (non-unique only)...")
    _plot_purity_2d(tdir, sv_sig, "cleaned",  "shared", shared_only=True)

    print("  [had/disambiguation] Purity 2D: post-disambiguation...")
    _plot_purity_2d(tdir, sv_sig, "disambig", "post",   shared_only=False)
