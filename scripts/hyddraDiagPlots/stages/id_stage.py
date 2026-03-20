"""
stages/id.py — ID-stage plots for leptonic HYDDRA.

The leptonic ID cut requires tracksSize == 2 and (passLooseMuonID OR passLooseElectronID).

Plots:
  - reco_id_{cosTheta,decayAngle,pOverE,dxySignif,mass}
  - id_flavor_breakdown  (muon-ID vs electron-ID, signal vs non-signal)
"""
import numpy as np
import awkward as ak
import ROOT

from ..src.config  import RECO_OBSERVABLES, STAGE_IDX, COLOR_GOLD, COLOR_NONSIGNAL, COLOR_BKG, GEN_DXY_BINS
from ..src.style   import draw_cms_label, draw_axis_grid
from ..src.plotter import (plot_reco_observable, plot_costheta_zoom,
                            plot_dxy_1d, plot_2d_dxy_vs_x,
                            plot_fakes_vs_dxy,
                            geom_cut_lines_z, geom_cut_lines_eta)


def plot_id_flavor_breakdown(tdir, sv_sig, sv_bkg=None):
    """
    Bar chart showing the number of ID-passing vertices split by muon-ID vs
    electron-ID, for signal (gold) and non-signal separately.
    """
    tag = "id_flavor_breakdown"

    if sv_sig is None or len(sv_sig) == 0:
        print(f"    [{tag}] No allStageVtx data — skipping")
        return

    if "StageVtx_passLooseMuonID" not in sv_sig.fields:
        print(f"    [{tag}] StageVtx_passLooseMuonID not in ntuple — skipping")
        return

    sidx_id = STAGE_IDX["id"]
    mask_id = ak.to_numpy(sv_sig["StageVtx_stageIdx"]) == sidx_id

    if not np.any(mask_id):
        print(f"    [{tag}] No ID-stage vertices found — skipping")
        return

    is_gold = ak.to_numpy(sv_sig["StageVtx_isGold"          ])[mask_id].astype(bool)
    is_muon = ak.to_numpy(sv_sig["StageVtx_passLooseMuonID"    ])[mask_id].astype(bool)
    is_elec = ak.to_numpy(sv_sig["StageVtx_passLooseElectronID"])[mask_id].astype(bool)

    n_sig_muon = int(np.sum( is_gold &  is_muon))
    n_sig_elec = int(np.sum( is_gold &  is_elec))
    n_ns_muon  = int(np.sum(~is_gold &  is_muon))
    n_ns_elec  = int(np.sum(~is_gold &  is_elec))

    h_sig    = ROOT.TH1F("h_id_flavor_sig",    "", 2, 0.5, 2.5)
    h_nonsig = ROOT.TH1F("h_id_flavor_nonsig", "", 2, 0.5, 2.5)
    for i, (ns, nn) in enumerate([(n_sig_muon, n_ns_muon), (n_sig_elec, n_ns_elec)], 1):
        h_sig   .SetBinContent(i, ns)
        h_nonsig.SetBinContent(i, nn)
        h_sig   .GetXaxis().SetBinLabel(i, ["Muon ID", "Electron ID"][i - 1])
        h_nonsig.GetXaxis().SetBinLabel(i, ["Muon ID", "Electron ID"][i - 1])

    h_sig.SetFillColor(COLOR_GOLD);         h_sig.SetLineColor(COLOR_GOLD)
    h_sig.SetBarWidth(0.35);                h_sig.SetBarOffset(0.10)
    h_nonsig.SetFillColor(COLOR_NONSIGNAL); h_nonsig.SetLineColor(COLOR_NONSIGNAL)
    h_nonsig.SetFillStyle(3004)
    h_nonsig.SetBarWidth(0.35);             h_nonsig.SetBarOffset(0.55)
    h_sig.SetStats(0);  h_nonsig.SetStats(0)

    max_y = max(h_sig.GetMaximum(), h_nonsig.GetMaximum())
    h_ax = ROOT.TH1F("h_ax_id_flavor", ";Lepton ID type;Vertices at ID stage", 2, 0.5, 2.5)
    h_ax.SetMinimum(0)
    h_ax.SetMaximum(max_y * 1.5 if max_y > 0 else 10)
    for i, lbl in enumerate(["Muon ID", "Electron ID"], 1):
        h_ax.GetXaxis().SetBinLabel(i, lbl)
    h_ax.GetXaxis().SetLabelSize(0.06)
    h_ax.GetXaxis().CenterTitle(True); h_ax.GetYaxis().CenterTitle(True)
    h_ax.GetXaxis().SetTitleSize(0.045); h_ax.GetYaxis().SetTitleSize(0.045)
    h_ax.SetStats(0)

    c = ROOT.TCanvas(tag, "ID flavor breakdown", 700, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.14); c.SetTopMargin(0.10)
    h_ax.Draw("AXIS")
    h_sig.Draw("BAR HIST SAME")
    h_nonsig.Draw("BAR HIST SAME")
    c._grid_lines = draw_axis_grid(h_ax, logy=False)

    leg = ROOT.TLegend(0.55, 0.76, 0.88, 0.88)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.035)
    leg.AddEntry(h_sig,    "Signal (gold)", "f")
    leg.AddEntry(h_nonsig, "Non-signal",    "f")
    leg.Draw()
    draw_cms_label()
    c.Update()

    c._h_ax = h_ax; c._h_sig = h_sig; c._h_nonsig = h_nonsig; c._leg = leg
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


def _id_sig_mask(sv_sig):
    return ak.to_numpy(sv_sig["StageVtx_isGold"]).astype(bool)


def _id_pos(sv_sig, sidx):
    mask     = ak.to_numpy(sv_sig["StageVtx_stageIdx"]).astype(int) == sidx
    x        = ak.to_numpy(sv_sig["StageVtx_x"]).astype(float)[mask]
    y        = ak.to_numpy(sv_sig["StageVtx_y"]).astype(float)[mask]
    z        = ak.to_numpy(sv_sig["StageVtx_z"]).astype(float)[mask]
    dxy      = np.sqrt(x**2 + y**2)
    eta_orig = np.arcsinh(z / np.where(dxy > 0, dxy, 1e-6))
    is_sig   = _id_sig_mask(sv_sig)[mask]
    return z, dxy, eta_orig, is_sig


def plot_eff_vs_dxy(tdir, gf):
    """
    Gold efficiency vs gen dxy at the ID stage, split by muon and electron flavor.
    Denominator: gen vertices with hasTracks by flavor.
    Numerator: those also gold at the ID stage.
    """
    tag = "id_eff_vs_dxy"

    if gf is None or len(gf) == 0:
        print(f"    [{tag}] No gen funnel data — skipping")
        return
    if "GenFunnel_gold_id" not in gf.fields:
        print(f"    [{tag}] GenFunnel_gold_id not in ntuple — skipping")
        return

    has_tracks = ak.to_numpy(ak.flatten(gf["GenFunnel_hasTracks"])).astype(bool)
    is_muon    = ak.to_numpy(ak.flatten(gf["GenFunnel_isMuon"    ]))[has_tracks].astype(bool)
    is_elec    = ak.to_numpy(ak.flatten(gf["GenFunnel_isElectron"]))[has_tracks].astype(bool)
    gold_id    = ak.to_numpy(ak.flatten(gf["GenFunnel_gold_id"   ]))[has_tracks].astype(bool)
    dxy        = ak.to_numpy(ak.flatten(gf["GenFunnel_dxy"       ]))[has_tracks]

    bins   = np.array(GEN_DXY_BINS, dtype=float)
    n_bins = len(bins) - 1

    flavor_configs = [
        ("Muon",     is_muon, ROOT.kRed+1,   20),
        ("Electron", is_elec, ROOT.kBlue+2,  24),
    ]

    graphs, hists = [], []
    for flavor, flavor_mask, color, marker in flavor_configs:
        dxy_f  = dxy[flavor_mask]
        gold_f = gold_id[flavor_mask]

        h_denom = ROOT.TH1F(f"h_denom_{tag}_{flavor}", "", n_bins, bins)
        h_denom.SetDirectory(0)
        for v in dxy_f:
            h_denom.Fill(v)

        h_num = ROOT.TH1F(f"h_num_{tag}_{flavor}", "", n_bins, bins)
        h_num.SetDirectory(0)
        for v in dxy_f[gold_f]:
            h_num.Fill(v)

        eff = ROOT.TEfficiency(h_num, h_denom)
        eff.SetStatisticOption(ROOT.TEfficiency.kFCP)
        g = eff.CreateGraph()
        g.SetTitle(flavor)
        g.SetMarkerColor(color); g.SetLineColor(color)
        g.SetMarkerStyle(marker); g.SetLineWidth(2)
        graphs.append(g)
        hists.append((h_denom, h_num, eff))

    c = ROOT.TCanvas(tag, "ID gold efficiency vs gen dxy", 800, 600)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06)
    c.SetBottomMargin(0.12); c.SetTopMargin(0.10)

    h_ax = ROOT.TH1F(f"h_ax_{tag}", ";Gen dxy (cm);Gold Efficiency at ID Stage", n_bins, bins)
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

    leg = ROOT.TLegend(0.55, 0.20, 0.88, 0.38)
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.035)
    for g in graphs:
        leg.AddEntry(g, g.GetTitle(), "lp")
    leg.Draw()
    draw_cms_label()
    c.Update()

    c._h_ax = h_ax; c._hists = hists; c._graphs = graphs; c._leg = leg
    tdir.cd()
    c.Write()
    print(f"    [{tag}] done")


# ── Orchestrator ───────────────────────────────────────────────────────────────

def make_plots(tdir, gf, sv_sig, sv_bkg=None):
    sidx   = STAGE_IDX["id"]
    is_sig = _id_sig_mask(sv_sig) if sv_sig is not None and len(sv_sig) > 0 else None

    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "id", obs_key, obs_cfg, sv_bkg)
    plot_id_flavor_breakdown(tdir, sv_sig)

    if is_sig is not None:
        plot_costheta_zoom(tdir, sv_sig, sv_bkg, sidx, is_sig,
                           "id_cosTheta_zoom", sig_label="Signal (gold)")
        plot_dxy_1d(tdir, sv_sig, sv_bkg, sidx, is_sig,
                    "id_dxy", sig_label="Signal (gold)")
        if "StageVtx_z" in sv_sig.fields:
            z, dxy, eta, is_sig_pos = _id_pos(sv_sig, sidx)
            plot_2d_dxy_vs_x(tdir, z,   dxy, is_sig_pos, "Vertex z (cm)",
                             np.linspace(-30, 30, 61), np.linspace(0, 80, 41),
                             "id_dxy_vs_z", sig_label="Signal (gold)",
                             cut_lines_fn=geom_cut_lines_z)
            plot_2d_dxy_vs_x(tdir, eta, dxy, is_sig_pos, "#eta_{origin}",
                             np.linspace(-5, 5, 51), np.linspace(0, 80, 41),
                             "id_dxy_vs_etaorigin", sig_label="Signal (gold)",
                             cut_lines_fn=geom_cut_lines_eta)

    plot_eff_vs_dxy(tdir, gf)
    if is_sig is not None:
        plot_fakes_vs_dxy(tdir, sv_sig, is_sig, tag="id_fakes_vs_dxy")
