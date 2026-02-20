#!/usr/bin/env python3
"""
Gen-match threshold study script

Reads output from processLLPNanoSV.py and produces diagnostic plots of
per-leg deltaR and relPtDiff distributions, split by gold match status.
Helps choose optimal thresholds for gold matching.

Usage:
    python genMatchThresholdStudy.py -i llpNanoSV.root
    python genMatchThresholdStudy.py -i llpNanoSV.root --delta-r 0.2 --rel-pt-diff 0.3
    python genMatchThresholdStudy.py -i "output_*.root" -o combined_study.root
"""

import argparse
import glob
import os
import numpy as np
import uproot
import ROOT


def load_data(filename, tree_path=None):
    """Load SV data from ROOT file produced by processLLPNanoSV.py"""
    file = uproot.open(filename)

    if tree_path:
        tree = file[tree_path]
    else:
        possible_paths = ['hyddraSVAnalyzer/tree', 'llpNanoSVAnalyzer/tree', 'tree']
        tree = None
        for path in possible_paths:
            if path in file:
                tree = file[path]
                print(f"Found tree at: {path}")
                break
        if tree is None:
            print(f"Available keys: {file.keys()}")
            raise KeyError("Could not find tree. Try specifying --tree-path")

    branches = [
        'HyddraSV_deltaR1', 'HyddraSV_deltaR2',
        'HyddraSV_relPtDiff1', 'HyddraSV_relPtDiff2',
        'HyddraSV_isGold',
    ]

    available = tree.keys()
    to_load = [b for b in branches if b in available]
    missing = [b for b in branches if b not in available]
    if missing:
        raise KeyError(f"Missing required branches: {missing}")

    data = tree.arrays(to_load, library='np')
    return data


def flatten_legs(data):
    """Flatten per-SV leg branches into combined arrays.

    Returns (deltaR, relPtDiff, isGold) where each entry is one leg.
    Legs with deltaR < 0 (no gen match found) are excluded.
    isGold is duplicated per-SV so both legs share the same gold status.
    """
    dr1_list, dr2_list = [], []
    rpt1_list, rpt2_list = [], []
    gold1_list, gold2_list = [], []

    n_events = len(data['HyddraSV_deltaR1'])
    for ei in range(n_events):
        dr1 = data['HyddraSV_deltaR1'][ei]
        dr2 = data['HyddraSV_deltaR2'][ei]
        rpt1 = data['HyddraSV_relPtDiff1'][ei]
        rpt2 = data['HyddraSV_relPtDiff2'][ei]
        gold = data['HyddraSV_isGold'][ei]

        for i in range(len(dr1)):
            dr1_list.append(dr1[i])
            dr2_list.append(dr2[i])
            rpt1_list.append(rpt1[i])
            rpt2_list.append(rpt2[i])
            gold1_list.append(gold[i])
            gold2_list.append(gold[i])

    # Combine leg1 and leg2
    all_dr = np.concatenate([np.array(dr1_list), np.array(dr2_list)])
    all_rpt = np.concatenate([np.array(rpt1_list), np.array(rpt2_list)])
    all_gold = np.concatenate([np.array(gold1_list), np.array(gold2_list)]).astype(bool)

    # Filter out legs with no gen match (deltaR < 0)
    valid = all_dr >= 0
    print(f"Total legs: {len(all_dr)}, with gen match: {np.sum(valid)}, "
          f"gold: {np.sum(all_gold[valid])}, non-gold: {np.sum(~all_gold[valid])}")

    return all_dr[valid], all_rpt[valid], all_gold[valid]


def make_scatter_all(deltaR, relPtDiff, isGold, dr_cut, relpt_cut):
    """2D scatter: deltaR vs relPtDiff, colored by gold status, with threshold lines."""
    c = ROOT.TCanvas('c_scatter_all', 'deltaR vs relPtDiff (all legs)', 800, 600)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.08)
    c.SetBottomMargin(0.12)

    # Non-gold (red) drawn first, then gold (green) on top
    ng = ~isGold
    g_nongold = ROOT.TGraph(int(np.sum(ng)))
    for i, (dr, rpt) in enumerate(zip(deltaR[ng], relPtDiff[ng])):
        g_nongold.SetPoint(i, float(dr), float(rpt))
    g_nongold.SetMarkerStyle(6)
    g_nongold.SetMarkerColor(ROOT.kRed)
    g_nongold.SetTitle(';#DeltaR(reco, gen);|p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen}')

    g_gold = ROOT.TGraph(int(np.sum(isGold)))
    for i, (dr, rpt) in enumerate(zip(deltaR[isGold], relPtDiff[isGold])):
        g_gold.SetPoint(i, float(dr), float(rpt))
    g_gold.SetMarkerStyle(6)
    g_gold.SetMarkerColor(ROOT.kGreen + 2)

    # Set axis ranges
    dr_max = min(float(np.percentile(deltaR, 99.5)), 5.0)
    rpt_max = min(float(np.percentile(relPtDiff, 99.5)), 5.0)

    g_nongold.GetXaxis().SetRangeUser(0, dr_max)
    g_nongold.GetYaxis().SetRangeUser(0, rpt_max)
    g_nongold.GetXaxis().SetTitleSize(0.05)
    g_nongold.GetXaxis().SetLabelSize(0.04)
    g_nongold.GetYaxis().SetTitleSize(0.05)
    g_nongold.GetYaxis().SetLabelSize(0.04)

    g_nongold.Draw('AP')
    g_gold.Draw('P SAME')

    # Threshold lines
    line_dr = ROOT.TLine(dr_cut, 0, dr_cut, rpt_max)
    line_dr.SetLineColor(ROOT.kBlue)
    line_dr.SetLineStyle(2)
    line_dr.SetLineWidth(2)
    line_dr.Draw()

    line_rpt = ROOT.TLine(0, relpt_cut, dr_max, relpt_cut)
    line_rpt.SetLineColor(ROOT.kBlue)
    line_rpt.SetLineStyle(2)
    line_rpt.SetLineWidth(2)
    line_rpt.Draw()

    legend = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(g_gold, 'Gold-matched', 'p')
    legend.AddEntry(g_nongold, 'Non-gold', 'p')
    legend.AddEntry(line_dr, f'Thresholds (#DeltaR={dr_cut}, relPt={relpt_cut})', 'l')
    legend.Draw()

    c.Update()
    return c, [g_nongold, g_gold, line_dr, line_rpt, legend]


def make_scatter_gold(deltaR, relPtDiff, isGold, dr_cut, relpt_cut):
    """2D scatter: deltaR vs relPtDiff, gold-matched legs only (zoomed)."""
    c = ROOT.TCanvas('c_scatter_gold', 'deltaR vs relPtDiff (gold only)', 800, 600)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.08)
    c.SetBottomMargin(0.12)

    gold_dr = deltaR[isGold]
    gold_rpt = relPtDiff[isGold]

    g = ROOT.TGraph(len(gold_dr))
    for i, (dr, rpt) in enumerate(zip(gold_dr, gold_rpt)):
        g.SetPoint(i, float(dr), float(rpt))
    g.SetMarkerStyle(6)
    g.SetMarkerColor(ROOT.kGreen + 2)
    g.SetTitle(';#DeltaR(reco, gen);|p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen}')

    # Zoomed range based on gold distribution
    if len(gold_dr) > 0:
        dr_max = min(float(np.percentile(gold_dr, 99.5)) * 1.2, 2.0)
        rpt_max = min(float(np.percentile(gold_rpt, 99.5)) * 1.2, 2.0)
    else:
        dr_max = 1.0
        rpt_max = 1.0

    g.GetXaxis().SetRangeUser(0, dr_max)
    g.GetYaxis().SetRangeUser(0, rpt_max)
    g.GetXaxis().SetTitleSize(0.05)
    g.GetXaxis().SetLabelSize(0.04)
    g.GetYaxis().SetTitleSize(0.05)
    g.GetYaxis().SetLabelSize(0.04)

    g.Draw('AP')

    line_dr = ROOT.TLine(dr_cut, 0, dr_cut, rpt_max)
    line_dr.SetLineColor(ROOT.kBlue)
    line_dr.SetLineStyle(2)
    line_dr.SetLineWidth(2)
    line_dr.Draw()

    line_rpt = ROOT.TLine(0, relpt_cut, dr_max, relpt_cut)
    line_rpt.SetLineColor(ROOT.kBlue)
    line_rpt.SetLineStyle(2)
    line_rpt.SetLineWidth(2)
    line_rpt.Draw()

    legend = ROOT.TLegend(0.55, 0.78, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(g, 'Gold-matched legs', 'p')
    legend.AddEntry(line_dr, f'Thresholds (#DeltaR={dr_cut}, relPt={relpt_cut})', 'l')
    legend.Draw()

    c.Update()
    return c, [g, line_dr, line_rpt, legend]


def make_1d_deltaR(deltaR, isGold):
    """1D deltaR distribution, gold vs non-gold, log y-axis."""
    c = ROOT.TCanvas('c_deltaR', 'deltaR distribution', 800, 600)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetLogy()

    nbins = 100
    dr_max = min(float(np.percentile(deltaR, 99.5)), 5.0)

    h_gold = ROOT.TH1F('h_deltaR_gold', ';#DeltaR(reco, gen);Legs / bin', nbins, 0, dr_max)
    h_nongold = ROOT.TH1F('h_deltaR_nongold', ';#DeltaR(reco, gen);Legs / bin', nbins, 0, dr_max)

    for dr, g in zip(deltaR, isGold):
        if g:
            h_gold.Fill(float(dr))
        else:
            h_nongold.Fill(float(dr))

    h_gold.SetLineColor(ROOT.kGreen + 2)
    h_gold.SetFillColor(ROOT.kGreen + 2)
    h_gold.SetFillStyle(3004)
    h_gold.SetLineWidth(2)

    h_nongold.SetLineColor(ROOT.kRed)
    h_nongold.SetFillColor(ROOT.kRed)
    h_nongold.SetFillStyle(3005)
    h_nongold.SetLineWidth(2)

    max_val = max(h_gold.GetMaximum(), h_nongold.GetMaximum(), 1)
    h_nongold.SetMinimum(0.5)
    h_nongold.SetMaximum(max_val * 5)
    h_nongold.SetStats(0)
    h_nongold.GetXaxis().SetTitleSize(0.05)
    h_nongold.GetXaxis().SetLabelSize(0.04)
    h_nongold.GetYaxis().SetTitleSize(0.05)
    h_nongold.GetYaxis().SetLabelSize(0.04)

    h_nongold.Draw('HIST')
    h_gold.Draw('HIST SAME')

    legend = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_gold, 'Gold-matched', 'f')
    legend.AddEntry(h_nongold, 'Non-gold', 'f')
    legend.Draw()

    c.Update()
    return c, [h_gold, h_nongold, legend]


def make_1d_relPtDiff(relPtDiff, isGold):
    """1D relPtDiff distribution, gold vs non-gold, log y-axis."""
    c = ROOT.TCanvas('c_relPtDiff', 'relPtDiff distribution', 800, 600)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetLogy()

    nbins = 100
    rpt_max = min(float(np.percentile(relPtDiff, 99.5)), 5.0)

    h_gold = ROOT.TH1F('h_relPtDiff_gold',
                        ';|p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen};Legs / bin',
                        nbins, 0, rpt_max)
    h_nongold = ROOT.TH1F('h_relPtDiff_nongold',
                           ';|p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen};Legs / bin',
                           nbins, 0, rpt_max)

    for rpt, g in zip(relPtDiff, isGold):
        if g:
            h_gold.Fill(float(rpt))
        else:
            h_nongold.Fill(float(rpt))

    h_gold.SetLineColor(ROOT.kGreen + 2)
    h_gold.SetFillColor(ROOT.kGreen + 2)
    h_gold.SetFillStyle(3004)
    h_gold.SetLineWidth(2)

    h_nongold.SetLineColor(ROOT.kRed)
    h_nongold.SetFillColor(ROOT.kRed)
    h_nongold.SetFillStyle(3005)
    h_nongold.SetLineWidth(2)

    max_val = max(h_gold.GetMaximum(), h_nongold.GetMaximum(), 1)
    h_nongold.SetMinimum(0.5)
    h_nongold.SetMaximum(max_val * 5)
    h_nongold.SetStats(0)
    h_nongold.GetXaxis().SetTitleSize(0.05)
    h_nongold.GetXaxis().SetLabelSize(0.04)
    h_nongold.GetYaxis().SetTitleSize(0.05)
    h_nongold.GetYaxis().SetLabelSize(0.04)

    h_nongold.Draw('HIST')
    h_gold.Draw('HIST SAME')

    legend = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_gold, 'Gold-matched', 'f')
    legend.AddEntry(h_nongold, 'Non-gold', 'f')
    legend.Draw()

    c.Update()
    return c, [h_gold, h_nongold, legend]


def make_2d_relPtDiff_vs_deltaR(deltaR, relPtDiff, isGold):
    """TH2D of relPtDiff vs deltaR: all, gold-only, and non-gold-only.

    Note: relPtDiff as stored is |reco-gen|/gen (non-negative).
    Returns list of (canvas, objects) tuples for 3 plots.
    """
    dr_max = min(float(np.percentile(deltaR, 99.5)), 5.0)
    rpt_max = min(float(np.percentile(relPtDiff, 99.5)), 5.0)
    nbins_dr = 100
    nbins_rpt = 100

    results = []

    configs = [
        ('all', 'All legs', np.ones(len(deltaR), dtype=bool)),
        ('gold', 'Gold-matched legs', isGold),
        ('nongold', 'Non-gold legs', ~isGold),
    ]

    canvases = []
    histos = []

    for tag, title, mask in configs:
        c = ROOT.TCanvas(f'c_2d_{tag}', f'relPtDiff vs deltaR ({title})', 800, 600)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.14)
        c.SetBottomMargin(0.12)
        c.SetLogz()

        h = ROOT.TH2D(f'h2_relPtDiff_vs_deltaR_{tag}',
                       f'{title};#DeltaR(reco, gen);'
                       f'|p_{{T}}^{{reco}} - p_{{T}}^{{gen}}| / p_{{T}}^{{gen}};Legs',
                       nbins_dr, 0, dr_max, nbins_rpt, 0, rpt_max)
        h.SetDirectory(0)
        h.SetStats(0)
        h.SetMinimum(0.5)
        h.GetXaxis().SetTitleSize(0.05)
        h.GetXaxis().SetLabelSize(0.04)
        h.GetYaxis().SetTitleSize(0.05)
        h.GetYaxis().SetLabelSize(0.04)
        h.GetZaxis().SetTitleSize(0.05)
        h.GetZaxis().SetLabelSize(0.04)

        for dr, rpt in zip(deltaR[mask], relPtDiff[mask]):
            h.Fill(float(dr), float(rpt))

        h.Draw('COLZ')

        c.Update()
        canvases.append(c)
        histos.append(h)

    return [(c, [h]) for c, h in zip(canvases, histos)]


def main():
    parser = argparse.ArgumentParser(
        description='Gen-match threshold study: visualize deltaR and relPtDiff distributions')
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='Input ROOT file(s) or glob pattern(s) (from processLLPNanoSV.py)')
    parser.add_argument('-o', '--output', default='gen_match_threshold_study.root',
                        help='Output ROOT file name')
    parser.add_argument('--tree-path', default=None,
                        help='Path to tree in ROOT file (default: auto-detect)')
    parser.add_argument('--delta-r', type=float, default=0.3,
                        help='deltaR threshold line to draw on plots (default: 0.3)')
    parser.add_argument('--rel-pt-diff', type=float, default=0.5,
                        help='relPtDiff threshold line to draw on plots (default: 0.5)')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # Expand glob patterns and collect unique files
    input_files = []
    for pattern in args.input:
        expanded = sorted(glob.glob(pattern))
        if expanded:
            input_files.extend(expanded)
        else:
            # Not a glob, treat as literal path
            input_files.append(pattern)
    # Deduplicate while preserving order
    seen = set()
    unique_files = []
    for f in input_files:
        fabs = os.path.abspath(f)
        if fabs not in seen:
            seen.add(fabs)
            unique_files.append(f)
    input_files = unique_files

    if not input_files:
        print("ERROR: No input files found.")
        return

    print(f"Processing {len(input_files)} file(s)...")

    # Load and flatten all files
    all_deltaR, all_relPtDiff, all_isGold = [], [], []
    for fname in input_files:
        print(f"  Loading {fname}...")
        data = load_data(fname, args.tree_path)
        dr, rpt, gold = flatten_legs(data)
        all_deltaR.append(dr)
        all_relPtDiff.append(rpt)
        all_isGold.append(gold)

    deltaR = np.concatenate(all_deltaR)
    relPtDiff = np.concatenate(all_relPtDiff)
    isGold = np.concatenate(all_isGold)

    print(f"\nCombined: {len(deltaR)} legs, {np.sum(isGold)} gold, {np.sum(~isGold)} non-gold")

    if len(deltaR) == 0:
        print("ERROR: No valid legs found. Check your input file.")
        return

    # Create plots
    all_canvases = []
    all_objects = []

    def save(canvas, objects):
        all_canvases.append(canvas)
        all_objects.extend(objects)

    print("\nCreating plots...")

    c, o = make_scatter_all(deltaR, relPtDiff, isGold, args.delta_r, args.rel_pt_diff)
    save(c, o)

    c, o = make_scatter_gold(deltaR, relPtDiff, isGold, args.delta_r, args.rel_pt_diff)
    save(c, o)

    c, o = make_1d_deltaR(deltaR, isGold)
    save(c, o)

    c, o = make_1d_relPtDiff(relPtDiff, isGold)
    save(c, o)

    for c, o in make_2d_relPtDiff_vs_deltaR(deltaR, relPtDiff, isGold):
        save(c, o)

    # Save to ROOT file
    output_file = ROOT.TFile(args.output, 'RECREATE')
    for canvas in all_canvases:
        canvas.Write()
    for obj in all_objects:
        if hasattr(obj, 'Write'):
            obj.Write()
    output_file.Close()

    print(f"\nDone! Output saved to {args.output}")
    print(f"Created {len(all_canvases)} plots")


if __name__ == "__main__":
    main()
