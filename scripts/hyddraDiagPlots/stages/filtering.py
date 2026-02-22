"""
stages/filtering.py — filtering-stage plots for hyddraDiagPlots.

Plots:
  - reco_filtered_{cosTheta,decayAngle,pOverE,dxySignif}
  - filtering_costheta_decay_{signal,nonsignal}_{before,after}
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import RECO_OBSERVABLES, STAGE_IDX
from ..src.style   import draw_cms_label
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
        draw_cms_label()
        c.Update()

        pal = h.GetListOfFunctions().FindObject("palette")
        if pal:
            h.GetListOfFunctions().Remove(pal)

        c._h = h
        tdir.cd()
        c.Write()
        print(f"    [{cname}] done")


# ── Orchestrator ───────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg):
    print("  [filtering] Reco observables...")
    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "filtered", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_filtered_{obs_key}] done")
    print("  [filtering] 2D cosTheta vs decay angle...")
    plot_filtering_costheta_decay_2d(tdir, sv_sig)
