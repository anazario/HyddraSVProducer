#!/usr/bin/env python3
"""
SV Efficiency Analysis Script

Analyzes secondary vertex efficiency from HyddraSVAnalyzer output.
- Leptonic SVs: Uses isGold matching (exclusive to leptonic decays)
- Hadronic SVs: Uses matchRatio thresholds (>=50% and >0%)
Uses uproot for reading and ROOT/cmsstyle for plotting.

Usage:
    python svEfficiencyAnalysis.py -i hyddraSV_ntuple.root
    python svEfficiencyAnalysis.py -i hyddraSV_ntuple.root --dxy-min 1e-5 --dxy-max 1.0
"""

import argparse
import uproot
import numpy as np
import ROOT
import cmsstyle as CMS
from array import array


def load_data(filename, tree_path=None):
    """Load SV data from ROOT file produced by HyddraSVAnalyzer"""
    file = uproot.open(filename)

    # Try common tree paths
    if tree_path:
        tree = file[tree_path]
    else:
        possible_paths = ['hyddraSVAnalyzer/tree', 'tree']
        tree = None
        for path in possible_paths:
            if path in file:
                tree = file[path]
                print(f"Found tree at: {path}")
                break

        if tree is None:
            print(f"Available keys: {file.keys()}")
            raise KeyError("Could not find tree. Try specifying --tree-path")

    # Reco SV branches
    reco_branches = [
        'HyddraSV_nVertices', 'HyddraSV_isLeptonic', 'HyddraSV_nTracks',
        'HyddraSV_x', 'HyddraSV_y', 'HyddraSV_z', 'HyddraSV_dxy',
        'HyddraSV_mass', 'HyddraSV_pt', 'HyddraSV_eta', 'HyddraSV_phi',
        'HyddraSV_chi2', 'HyddraSV_normalizedChi2',
        'HyddraSV_passLooseMuonID', 'HyddraSV_passLooseElectronID', 'HyddraSV_passLooseID',
        'HyddraSV_genVertexIndex', 'HyddraSV_nearestGenVertexIndex',
        'HyddraSV_isBronze', 'HyddraSV_isSilver', 'HyddraSV_isGold',
        'HyddraSV_matchRatio', 'HyddraSV_min3D',
    ]

    # Gen vertex branches
    gen_branches = [
        'HyddraGenVertex_nTotal', 'HyddraGenVertex_nElectron', 'HyddraGenVertex_nMuon', 'HyddraGenVertex_nHadronic',
        'HyddraGenVertex_dxy', 'HyddraGenVertex_pt', 'HyddraGenVertex_eta', 'HyddraGenVertex_phi',
        'HyddraGenVertex_x', 'HyddraGenVertex_y', 'HyddraGenVertex_z',
        'HyddraGenVertex_isElectron', 'HyddraGenVertex_isMuon', 'HyddraGenVertex_isHadronic',
        'HyddraGenVertex_passSelection', 'HyddraGenVertex_passSelectionAndCuts',
        'HyddraGenVertex_nTracks', 'HyddraGenVertex_mass',
    ]

    all_branches = reco_branches + gen_branches

    # Filter to only branches that exist
    available_branches = tree.keys()
    branches_to_load = [b for b in all_branches if b in available_branches]

    data = tree.arrays(branches_to_load, library='np')
    return data


def calculate_efficiency(values, pass_mask, bins):
    """Calculate efficiency in bins with binomial errors"""
    bin_centers = []
    efficiencies = []
    errors = []
    n_total_list = []
    n_pass_list = []

    for i in range(len(bins) - 1):
        bin_low = bins[i]
        bin_high = bins[i + 1]
        bin_center = (bin_low + bin_high) / 2.0

        in_bin = (values >= bin_low) & (values < bin_high)
        n_total = np.sum(in_bin)
        n_pass = np.sum(pass_mask[in_bin])

        if n_total > 0:
            efficiency = n_pass / n_total
            error = np.sqrt(efficiency * (1 - efficiency) / n_total)
        else:
            efficiency = 0
            error = 0

        bin_centers.append(bin_center)
        efficiencies.append(efficiency)
        errors.append(error)
        n_total_list.append(n_total)
        n_pass_list.append(n_pass)

    return (np.array(bin_centers), np.array(efficiencies), np.array(errors),
            np.array(n_total_list), np.array(n_pass_list))


def create_tgraph(x, y, yerr, name, title):
    """Create ROOT TGraphErrors"""
    n_points = len(x)
    graph = ROOT.TGraphErrors(n_points)
    graph.SetName(name)
    graph.SetTitle(title)

    for i in range(n_points):
        graph.SetPoint(i, float(x[i]), float(y[i]))
        graph.SetPointError(i, 0.0, float(yerr[i]))

    return graph


def draw_cms_labels():
    """Draw standard CMS labels on plots"""
    cms_label = ROOT.TLatex()
    cms_label.SetNDC()
    cms_label.SetTextSize(0.065)
    cms_label.SetTextFont(61)
    cms_label.DrawLatex(0.122, 0.945, "CMS")
    cms_label.SetTextFont(52)
    cms_label.SetTextSize(0.05)
    cms_label.DrawLatex(0.23, 0.945, "Simulation Preliminary")
    return cms_label


def analyze_sv_efficiency(data, dxy_min=1e-5, dxy_max=1.0, dxy_nbins=21):
    """Analyze SV efficiency from gen vertex perspective.

    Leptonic SVs: Uses isGold matching (exclusive to leptonic decays)
    Hadronic SVs: Uses matchRatio thresholds (>=50% and >0%)
    """

    # Flatten jagged arrays for gen vertices, tracking matches from reco
    all_gen_dxy = []
    all_gen_pt = []
    all_gen_eta = []
    all_gen_isElectron = []
    all_gen_isMuon = []
    all_gen_isHadronic = []
    all_gen_hasGoldMatch = []          # For leptonic (isGold)
    all_gen_hasMatchRatio50 = []       # For hadronic (matchRatio >= 0.5)
    all_gen_hasMatchRatioAny = []      # For hadronic (matchRatio > 0)

    n_events = len(data['HyddraGenVertex_dxy'])

    for event_idx in range(n_events):
        gen_dxy = data['HyddraGenVertex_dxy'][event_idx]
        gen_pt = data['HyddraGenVertex_pt'][event_idx]
        gen_eta = data['HyddraGenVertex_eta'][event_idx]
        gen_isElectron = data['HyddraGenVertex_isElectron'][event_idx]
        gen_isMuon = data['HyddraGenVertex_isMuon'][event_idx]
        gen_isHadronic = data['HyddraGenVertex_isHadronic'][event_idx]

        # Get reco vertex info for this event
        reco_genVertexIndex = data['HyddraSV_genVertexIndex'][event_idx]
        reco_nearestGenVertexIndex = data['HyddraSV_nearestGenVertexIndex'][event_idx]
        reco_isGold = data['HyddraSV_isGold'][event_idx]
        reco_matchRatio = data['HyddraSV_matchRatio'][event_idx]

        # Build sets of gen vertex indices for different matching criteria
        gold_matched_gen_indices = set()
        matchRatio50_gen_indices = set()
        matchRatioAny_gen_indices = set()

        for reco_idx in range(len(reco_genVertexIndex)):
            # Gold matching uses genVertexIndex (for leptonic)
            gen_idx = reco_genVertexIndex[reco_idx]
            if gen_idx >= 0 and reco_isGold[reco_idx]:
                gold_matched_gen_indices.add(gen_idx)

            # matchRatio matching uses nearestGenVertexIndex (for hadronic)
            nearest_gen_idx = reco_nearestGenVertexIndex[reco_idx]
            if nearest_gen_idx >= 0:
                # matchRatio >= 50%
                if reco_matchRatio[reco_idx] >= 0.5:
                    matchRatio50_gen_indices.add(nearest_gen_idx)
                # matchRatio > 0%
                if reco_matchRatio[reco_idx] > 0:
                    matchRatioAny_gen_indices.add(nearest_gen_idx)

        for i in range(len(gen_dxy)):
            all_gen_dxy.append(gen_dxy[i])
            all_gen_pt.append(gen_pt[i])
            all_gen_eta.append(gen_eta[i])
            all_gen_isElectron.append(gen_isElectron[i])
            all_gen_isMuon.append(gen_isMuon[i])
            all_gen_isHadronic.append(gen_isHadronic[i])
            all_gen_hasGoldMatch.append(i in gold_matched_gen_indices)
            all_gen_hasMatchRatio50.append(i in matchRatio50_gen_indices)
            all_gen_hasMatchRatioAny.append(i in matchRatioAny_gen_indices)

    all_gen_dxy = np.array(all_gen_dxy)
    all_gen_pt = np.array(all_gen_pt)
    all_gen_eta = np.array(all_gen_eta)
    all_gen_isElectron = np.array(all_gen_isElectron)
    all_gen_isMuon = np.array(all_gen_isMuon)
    all_gen_isHadronic = np.array(all_gen_isHadronic)
    all_gen_hasGoldMatch = np.array(all_gen_hasGoldMatch)
    all_gen_hasMatchRatio50 = np.array(all_gen_hasMatchRatio50)
    all_gen_hasMatchRatioAny = np.array(all_gen_hasMatchRatioAny)

    # Define masks
    is_leptonic = all_gen_isElectron | all_gen_isMuon

    # Define bins
    dxy_bins = np.logspace(np.log10(dxy_min), np.log10(dxy_max), dxy_nbins)
    pt_bins = np.array([0, 5, 10, 15, 20, 30, 50, 100, 200])
    eta_bins = np.linspace(-2.5, 2.5, 21)

    results = {
        'all': {},
        'leptonic': {},
        'electron': {},
        'muon': {},
        'hadronic_50': {},      # matchRatio >= 50%
        'hadronic_any': {},     # matchRatio > 0%
    }

    # Calculate efficiencies for leptonic categories using isGold matching
    leptonic_categories = [
        ('leptonic', is_leptonic),
        ('electron', all_gen_isElectron),
        ('muon', all_gen_isMuon),
    ]

    for label, mask in leptonic_categories:
        if np.sum(mask) == 0:
            continue

        # Efficiency vs dxy (isGold matching)
        centers, eff, err, n_total, n_pass = calculate_efficiency(
            all_gen_dxy[mask], all_gen_hasGoldMatch[mask], dxy_bins
        )
        results[label]['dxy'] = {
            'centers': centers, 'efficiency': eff, 'error': err,
            'n_total': n_total, 'n_pass': n_pass
        }

        # Efficiency vs pt (isGold matching)
        centers, eff, err, n_total, n_pass = calculate_efficiency(
            all_gen_pt[mask], all_gen_hasGoldMatch[mask], pt_bins
        )
        results[label]['pt'] = {
            'centers': centers, 'efficiency': eff, 'error': err,
            'n_total': n_total, 'n_pass': n_pass
        }

        # Efficiency vs eta (isGold matching)
        centers, eff, err, n_total, n_pass = calculate_efficiency(
            all_gen_eta[mask], all_gen_hasGoldMatch[mask], eta_bins
        )
        results[label]['eta'] = {
            'centers': centers, 'efficiency': eff, 'error': err,
            'n_total': n_total, 'n_pass': n_pass
        }

    # Calculate efficiencies for hadronic using matchRatio thresholds
    hadronic_mask = all_gen_isHadronic.astype(bool)
    if np.sum(hadronic_mask) > 0:
        # Hadronic with matchRatio >= 50%
        for var_name, var_data, bins in [
            ('dxy', all_gen_dxy, dxy_bins),
            ('pt', all_gen_pt, pt_bins),
            ('eta', all_gen_eta, eta_bins),
        ]:
            centers, eff, err, n_total, n_pass = calculate_efficiency(
                var_data[hadronic_mask], all_gen_hasMatchRatio50[hadronic_mask], bins
            )
            results['hadronic_50'][var_name] = {
                'centers': centers, 'efficiency': eff, 'error': err,
                'n_total': n_total, 'n_pass': n_pass
            }

        # Hadronic with matchRatio > 0%
        for var_name, var_data, bins in [
            ('dxy', all_gen_dxy, dxy_bins),
            ('pt', all_gen_pt, pt_bins),
            ('eta', all_gen_eta, eta_bins),
        ]:
            centers, eff, err, n_total, n_pass = calculate_efficiency(
                var_data[hadronic_mask], all_gen_hasMatchRatioAny[hadronic_mask], bins
            )
            results['hadronic_any'][var_name] = {
                'centers': centers, 'efficiency': eff, 'error': err,
                'n_total': n_total, 'n_pass': n_pass
            }

    # Summary statistics
    n_total_gen = len(all_gen_dxy)
    n_leptonic = np.sum(is_leptonic)
    n_hadronic = np.sum(all_gen_isHadronic)
    n_electrons = np.sum(all_gen_isElectron)
    n_muons = np.sum(all_gen_isMuon)

    # Leptonic uses isGold matching
    n_lep_gold = np.sum(all_gen_hasGoldMatch[is_leptonic])
    n_ele_gold = np.sum(all_gen_hasGoldMatch[all_gen_isElectron])
    n_mu_gold = np.sum(all_gen_hasGoldMatch[all_gen_isMuon])

    # Hadronic uses matchRatio thresholds
    n_had_ratio50 = np.sum(all_gen_hasMatchRatio50[all_gen_isHadronic])
    n_had_ratioAny = np.sum(all_gen_hasMatchRatioAny[all_gen_isHadronic])

    print(f"\n=== SV Efficiency Summary ===")
    print(f"Total gen vertices: {n_total_gen}")
    print(f"\n  Leptonic: {n_leptonic} (electrons: {n_electrons}, muons: {n_muons})")
    print(f"    [isGold matching]")
    print(f"    Gold matched: {n_lep_gold} ({100*n_lep_gold/max(1,n_leptonic):.1f}%)")
    print(f"    - Electrons: {n_ele_gold} ({100*n_ele_gold/max(1,n_electrons):.1f}%)")
    print(f"    - Muons: {n_mu_gold} ({100*n_mu_gold/max(1,n_muons):.1f}%)")
    print(f"\n  Hadronic: {n_hadronic}")
    print(f"    [matchRatio thresholds]")
    print(f"    matchRatio >= 50%: {n_had_ratio50} ({100*n_had_ratio50/max(1,n_hadronic):.1f}%)")
    print(f"    matchRatio > 0%:   {n_had_ratioAny} ({100*n_had_ratioAny/max(1,n_hadronic):.1f}%)")

    return results, {
        'all_gen_dxy': all_gen_dxy,
        'all_gen_pt': all_gen_pt,
        'all_gen_eta': all_gen_eta,
        'is_leptonic': is_leptonic,
        'is_hadronic': all_gen_isHadronic,
        'hasGoldMatch': all_gen_hasGoldMatch,
        'hasMatchRatio50': all_gen_hasMatchRatio50,
        'hasMatchRatioAny': all_gen_hasMatchRatioAny,
    }


def analyze_reco_sv_properties(data):
    """Analyze reconstructed SV properties"""

    all_dxy = []
    all_nTracks = []
    all_mass = []
    all_isLeptonic = []
    all_isGold = []
    all_isSilver = []
    all_isBronze = []
    all_matchRatio = []
    all_min3D = []

    n_events = len(data['HyddraSV_dxy'])

    for event_idx in range(n_events):
        dxy = data['HyddraSV_dxy'][event_idx]
        nTracks = data['HyddraSV_nTracks'][event_idx]
        mass = data['HyddraSV_mass'][event_idx]
        isLeptonic = data['HyddraSV_isLeptonic'][event_idx]
        isGold = data['HyddraSV_isGold'][event_idx]
        isSilver = data['HyddraSV_isSilver'][event_idx]
        isBronze = data['HyddraSV_isBronze'][event_idx]
        matchRatio = data['HyddraSV_matchRatio'][event_idx]
        min3D = data['HyddraSV_min3D'][event_idx]

        for i in range(len(dxy)):
            all_dxy.append(dxy[i])
            all_nTracks.append(nTracks[i])
            all_mass.append(mass[i])
            all_isLeptonic.append(isLeptonic[i])
            all_isGold.append(isGold[i])
            all_isSilver.append(isSilver[i])
            all_isBronze.append(isBronze[i])
            all_matchRatio.append(matchRatio[i])
            all_min3D.append(min3D[i])

    return {
        'dxy': np.array(all_dxy),
        'nTracks': np.array(all_nTracks),
        'mass': np.array(all_mass),
        'isLeptonic': np.array(all_isLeptonic),
        'isGold': np.array(all_isGold),
        'isSilver': np.array(all_isSilver),
        'isBronze': np.array(all_isBronze),
        'matchRatio': np.array(all_matchRatio),
        'min3D': np.array(all_min3D),
    }


def create_efficiency_vs_dxy_plot(results, output_name, dxy_min=1e-5, dxy_max=1.0, dxy_nbins=21):
    """Create efficiency vs gen dxy plot.

    Leptonic: uses isGold matching
    Hadronic: uses matchRatio thresholds (>=50% and >0%)
    """

    COLOR_LEPTONIC = ROOT.kOrange + 6
    COLOR_HADRONIC_50 = ROOT.kAzure + 1
    COLOR_HADRONIC_ANY = ROOT.kAzure - 4
    COLOR_ELECTRON = ROOT.kRed + 1
    COLOR_MUON = ROOT.kGreen + 2

    x_label = 'd_{xy}^{gen} [cm]'
    y_label = 'SV Reconstruction Efficiency'

    canvas = CMS.cmsCanvas(output_name, dxy_min, dxy_max, 0.0001, 1.5, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetLogx()
    canvas.SetLogy()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)

    dxy_bins = np.logspace(np.log10(dxy_min), np.log10(dxy_max), dxy_nbins)
    h_dummy = ROOT.TH1F(f'h_dummy_{output_name}', f';{x_label};{y_label}',
                        len(dxy_bins) - 1, dxy_bins)
    h_dummy.SetStats(0)
    h_dummy.SetMinimum(0.0001)
    h_dummy.SetMaximum(1.5)
    h_dummy.GetXaxis().SetMoreLogLabels()
    h_dummy.GetXaxis().SetNoExponent()
    h_dummy.GetXaxis().SetTitleSize(0.05)
    h_dummy.GetXaxis().SetLabelSize(0.05)
    h_dummy.GetXaxis().SetTitleOffset(1.25)
    h_dummy.GetXaxis().CenterTitle(True)
    h_dummy.GetYaxis().SetTitleSize(0.05)
    h_dummy.GetYaxis().SetLabelSize(0.05)
    h_dummy.GetYaxis().SetTitleOffset(1.25)
    h_dummy.GetYaxis().CenterTitle(True)
    h_dummy.Draw()

    graphs = []

    # Leptonic categories (isGold matching)
    plot_configs_leptonic = [
        ('leptonic', 'Leptonic (isGold)', COLOR_LEPTONIC, 20),
        ('electron', 'Electrons (isGold)', COLOR_ELECTRON, 22),
        ('muon', 'Muons (isGold)', COLOR_MUON, 23),
    ]

    for key, label, color, marker in plot_configs_leptonic:
        if key in results and 'dxy' in results[key]:
            r = results[key]['dxy']
            g = create_tgraph(r['centers'], r['efficiency'], r['error'],
                              f'efficiency_vs_dxy_{key}', label)
            g.SetMarkerStyle(marker)
            g.SetMarkerSize(1.2)
            g.SetMarkerColor(color)
            g.SetLineColor(color)
            g.Draw('P SAME')
            graphs.append((label, g, color))

    # Hadronic categories (matchRatio thresholds)
    plot_configs_hadronic = [
        ('hadronic_50', 'Hadronic (ratio #geq 50%)', COLOR_HADRONIC_50, 21),
        ('hadronic_any', 'Hadronic (ratio > 0%)', COLOR_HADRONIC_ANY, 25),
    ]

    for key, label, color, marker in plot_configs_hadronic:
        if key in results and 'dxy' in results[key]:
            r = results[key]['dxy']
            g = create_tgraph(r['centers'], r['efficiency'], r['error'],
                              f'efficiency_vs_dxy_{key}', label)
            g.SetMarkerStyle(marker)
            g.SetMarkerSize(1.2)
            g.SetMarkerColor(color)
            g.SetLineColor(color)
            g.Draw('P SAME')
            graphs.append((label, g, color))

    n_entries = len(graphs)
    legend_height = 0.045 * n_entries
    legend = ROOT.TLegend(0.45, 0.91 - legend_height, 0.88, 0.91)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)
    for label, g, _ in graphs:
        legend.AddEntry(g, label, 'p')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [g for _, g, _ in graphs] + [h_dummy, legend, cms_label]
    return canvas, all_objects


def create_efficiency_vs_pt_plot(results, output_name):
    """Create efficiency vs gen pt plot.

    Leptonic: uses isGold matching
    Hadronic: uses matchRatio thresholds (>=50% and >0%)
    """

    COLOR_LEPTONIC = ROOT.kOrange + 6
    COLOR_HADRONIC_50 = ROOT.kAzure + 1
    COLOR_HADRONIC_ANY = ROOT.kAzure - 4

    x_label = 'p_{T}^{gen} [GeV]'
    y_label = 'SV Reconstruction Efficiency'

    canvas = CMS.cmsCanvas(output_name, 0, 200, 0.0, 1.2, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)

    h_dummy = ROOT.TH1F(f'h_dummy_{output_name}', f';{x_label};{y_label}', 100, 0, 200)
    h_dummy.SetStats(0)
    h_dummy.SetMinimum(0.0)
    h_dummy.SetMaximum(1.2)
    h_dummy.GetXaxis().SetTitleSize(0.05)
    h_dummy.GetXaxis().SetLabelSize(0.05)
    h_dummy.GetXaxis().SetTitleOffset(1.25)
    h_dummy.GetXaxis().CenterTitle(True)
    h_dummy.GetYaxis().SetTitleSize(0.05)
    h_dummy.GetYaxis().SetLabelSize(0.05)
    h_dummy.GetYaxis().SetTitleOffset(1.25)
    h_dummy.GetYaxis().CenterTitle(True)
    h_dummy.Draw()

    graphs = []

    plot_configs = [
        ('Leptonic (isGold)', 'leptonic', COLOR_LEPTONIC, 20),
        ('Hadronic (ratio #geq 50%)', 'hadronic_50', COLOR_HADRONIC_50, 21),
        ('Hadronic (ratio > 0%)', 'hadronic_any', COLOR_HADRONIC_ANY, 25),
    ]

    for label, key, color, marker in plot_configs:
        if key in results and 'pt' in results[key]:
            r = results[key]['pt']
            g = create_tgraph(r['centers'], r['efficiency'], r['error'],
                              f'efficiency_vs_pt_{key}', label)
            g.SetMarkerStyle(marker)
            g.SetMarkerSize(1.2)
            g.SetMarkerColor(color)
            g.SetLineColor(color)
            g.Draw('P SAME')
            graphs.append((label, g, color))

    n_entries = len(graphs)
    legend_height = 0.045 * n_entries
    legend = ROOT.TLegend(0.45, 0.20, 0.88, 0.20 + legend_height)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)
    for label, g, _ in graphs:
        legend.AddEntry(g, label, 'p')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [g for _, g, _ in graphs] + [h_dummy, legend, cms_label]
    return canvas, all_objects


def create_matching_quality_plot(reco_props, output_name, dxy_min=1e-5, dxy_max=1.0):
    """Create plot showing gold/silver/bronze matching quality vs dxy"""

    dxy = reco_props['dxy']
    isGold = reco_props['isGold']
    isSilver = reco_props['isSilver']
    isBronze = reco_props['isBronze']
    isLeptonic = reco_props['isLeptonic']

    # Filter to leptonic only (where gold/silver/bronze makes sense)
    lep_mask = isLeptonic.astype(bool)

    x_label = 'd_{xy}^{reco} [cm]'
    y_label = 'Fraction'

    canvas = CMS.cmsCanvas(output_name, dxy_min, dxy_max, 0.0, 1.2, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetLogx()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)

    dxy_bins = np.logspace(np.log10(dxy_min), np.log10(dxy_max), 21)

    h_dummy = ROOT.TH1F(f'h_dummy_{output_name}', f';{x_label};{y_label}',
                        len(dxy_bins) - 1, dxy_bins)
    h_dummy.SetStats(0)
    h_dummy.SetMinimum(0.0)
    h_dummy.SetMaximum(1.2)
    h_dummy.GetXaxis().SetMoreLogLabels()
    h_dummy.GetXaxis().SetNoExponent()
    h_dummy.GetXaxis().SetTitleSize(0.05)
    h_dummy.GetXaxis().SetLabelSize(0.05)
    h_dummy.GetXaxis().SetTitleOffset(1.25)
    h_dummy.GetXaxis().CenterTitle(True)
    h_dummy.GetYaxis().SetTitleSize(0.05)
    h_dummy.GetYaxis().SetLabelSize(0.05)
    h_dummy.GetYaxis().SetTitleOffset(1.25)
    h_dummy.GetYaxis().CenterTitle(True)
    h_dummy.Draw()

    graphs = []

    # Calculate fraction in each quality category
    for label, mask, color, marker in [
        ('Gold', isGold[lep_mask], ROOT.kGreen + 2, 20),
        ('Silver', isSilver[lep_mask], ROOT.kBlue + 1, 21),
        ('Bronze', isBronze[lep_mask], ROOT.kOrange + 1, 22),
    ]:
        centers, frac, err, _, _ = calculate_efficiency(
            dxy[lep_mask], mask.astype(bool), dxy_bins
        )
        g = create_tgraph(centers, frac, err, f'quality_{label}', label)
        g.SetMarkerStyle(marker)
        g.SetMarkerSize(1.2)
        g.SetMarkerColor(color)
        g.SetLineColor(color)
        g.Draw('P SAME')
        graphs.append((label, g, color))

    legend = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for label, g, _ in graphs:
        legend.AddEntry(g, label, 'p')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [g for _, g, _ in graphs] + [h_dummy, legend, cms_label]
    return canvas, all_objects


def create_matchRatio_distribution_plot(reco_props, output_name):
    """Create matchRatio distribution plot"""

    matchRatio = reco_props['matchRatio']
    isLeptonic = reco_props['isLeptonic'].astype(bool)

    # Filter valid matchRatio
    valid = matchRatio >= 0

    x_label = 'Match Ratio'
    y_label = 'Entries'

    canvas = CMS.cmsCanvas(output_name, 0, 1, 0.1, 1e6, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetLogy()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)

    h_lep = ROOT.TH1F('h_matchRatio_leptonic', f';{x_label};{y_label}', 50, 0, 1)
    h_had = ROOT.TH1F('h_matchRatio_hadronic', f';{x_label};{y_label}', 50, 0, 1)

    for mr, is_lep in zip(matchRatio[valid], isLeptonic[valid]):
        if is_lep:
            h_lep.Fill(mr)
        else:
            h_had.Fill(mr)

    h_lep.SetStats(0)
    h_lep.SetLineColor(ROOT.kOrange + 6)
    h_lep.SetFillColor(ROOT.kOrange + 6)
    h_lep.SetFillStyle(3004)
    h_lep.SetLineWidth(2)

    h_had.SetLineColor(ROOT.kAzure + 1)
    h_had.SetFillColor(ROOT.kAzure + 1)
    h_had.SetFillStyle(3005)
    h_had.SetLineWidth(2)

    max_val = max(h_lep.GetMaximum(), h_had.GetMaximum(), 1)
    h_lep.SetMinimum(0.1)
    h_lep.SetMaximum(max_val * 5)
    h_lep.GetXaxis().SetTitleSize(0.05)
    h_lep.GetXaxis().SetLabelSize(0.05)
    h_lep.GetXaxis().SetTitleOffset(1.25)
    h_lep.GetXaxis().CenterTitle(True)
    h_lep.GetYaxis().SetTitleSize(0.05)
    h_lep.GetYaxis().SetLabelSize(0.05)
    h_lep.GetYaxis().SetTitleOffset(1.25)
    h_lep.GetYaxis().CenterTitle(True)

    h_lep.Draw('HIST')
    h_had.Draw('HIST SAME')

    legend = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_lep, 'Leptonic', 'f')
    legend.AddEntry(h_had, 'Hadronic', 'f')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [h_lep, h_had, legend, cms_label]
    return canvas, all_objects


def create_min3D_distribution_plot(reco_props, output_name):
    """Create min3D (distance to nearest gen vertex) distribution"""

    min3D = reco_props['min3D']
    isLeptonic = reco_props['isLeptonic'].astype(bool)

    # Filter valid min3D
    valid = min3D >= 0

    x_label = 'Min 3D Distance to Gen Vertex [cm]'
    y_label = 'Entries'

    canvas = CMS.cmsCanvas(output_name, 0.001, 100, 0.1, 1e6, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetLogx()
    canvas.SetLogy()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)

    dist_bins = np.logspace(-3, 2, 51)

    h_lep = ROOT.TH1F('h_min3D_leptonic', f';{x_label};{y_label}', len(dist_bins) - 1, dist_bins)
    h_had = ROOT.TH1F('h_min3D_hadronic', f';{x_label};{y_label}', len(dist_bins) - 1, dist_bins)

    for d, is_lep in zip(min3D[valid], isLeptonic[valid]):
        if is_lep:
            h_lep.Fill(d)
        else:
            h_had.Fill(d)

    h_lep.SetStats(0)
    h_lep.SetLineColor(ROOT.kOrange + 6)
    h_lep.SetFillColor(ROOT.kOrange + 6)
    h_lep.SetFillStyle(3004)
    h_lep.SetLineWidth(2)

    h_had.SetLineColor(ROOT.kAzure + 1)
    h_had.SetFillColor(ROOT.kAzure + 1)
    h_had.SetFillStyle(3005)
    h_had.SetLineWidth(2)

    max_val = max(h_lep.GetMaximum(), h_had.GetMaximum(), 1)
    h_lep.SetMinimum(0.1)
    h_lep.SetMaximum(max_val * 5)
    h_lep.GetXaxis().SetMoreLogLabels()
    h_lep.GetXaxis().SetNoExponent()
    h_lep.GetXaxis().SetTitleSize(0.05)
    h_lep.GetXaxis().SetLabelSize(0.05)
    h_lep.GetXaxis().SetTitleOffset(1.25)
    h_lep.GetXaxis().CenterTitle(True)
    h_lep.GetYaxis().SetTitleSize(0.05)
    h_lep.GetYaxis().SetLabelSize(0.05)
    h_lep.GetYaxis().SetTitleOffset(1.25)
    h_lep.GetYaxis().CenterTitle(True)

    h_lep.Draw('HIST')
    h_had.Draw('HIST SAME')

    # Vertical line at 0.05 cm (spatial matching cut)
    line = ROOT.TLine(0.05, 0.1, 0.05, max_val * 5)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()

    legend = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_lep, 'Leptonic', 'f')
    legend.AddEntry(h_had, 'Hadronic', 'f')
    legend.AddEntry(line, '0.05 cm cut', 'l')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [h_lep, h_had, line, legend, cms_label]
    return canvas, all_objects


def main():
    parser = argparse.ArgumentParser(description='SV efficiency analysis from HyddraSVAnalyzer output')
    parser.add_argument('-i', '--input', required=True,
                        help='Input ROOT file (from HyddraSVAnalyzer)')
    parser.add_argument('-o', '--output', default='sv_efficiency_analysis.root',
                        help='Output ROOT file name')
    parser.add_argument('--tree-path', default=None,
                        help='Path to tree in ROOT file (default: auto-detect)')
    parser.add_argument('--dxy-min', type=float, default=1e-5,
                        help='Minimum dxy for efficiency plot in cm (default: 1e-5)')
    parser.add_argument('--dxy-max', type=float, default=1.0,
                        help='Maximum dxy for efficiency plot in cm (default: 1.0)')
    parser.add_argument('--dxy-nbins', type=int, default=21,
                        help='Number of dxy bins (default: 21)')
    args = parser.parse_args()

    # Set CMS style
    CMS.SetExtraText("Simulation")
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # Load data
    print(f"Loading data from {args.input}...")
    data = load_data(args.input, args.tree_path)

    # Analyze efficiency
    results, summary = analyze_sv_efficiency(data, args.dxy_min, args.dxy_max, args.dxy_nbins)

    # Analyze reco SV properties
    reco_props = analyze_reco_sv_properties(data)

    # Collect all canvases and objects
    all_canvases = []
    all_objects = []

    def save(canvas, objects):
        all_canvases.append(canvas)
        all_objects.extend(objects)

    # Create plots
    print("\nCreating plots...")

    # Efficiency vs dxy (isGold matching)
    c, o = create_efficiency_vs_dxy_plot(results, 'sv_efficiency_vs_dxy',
                                          args.dxy_min, args.dxy_max, args.dxy_nbins)
    save(c, o)

    # Efficiency vs pt
    c, o = create_efficiency_vs_pt_plot(results, 'sv_efficiency_vs_pt')
    save(c, o)

    # Matching quality (gold/silver/bronze)
    c, o = create_matching_quality_plot(reco_props, 'sv_matching_quality',
                                         args.dxy_min, args.dxy_max)
    save(c, o)

    # Match ratio distribution
    c, o = create_matchRatio_distribution_plot(reco_props, 'sv_matchRatio_distribution')
    save(c, o)

    # Min3D distribution
    c, o = create_min3D_distribution_plot(reco_props, 'sv_min3D_distribution')
    save(c, o)

    # Save to ROOT file
    output_file = ROOT.TFile(args.output, 'RECREATE')

    for canvas in all_canvases:
        canvas.Write()
    for obj in all_objects:
        obj.Write()

    output_file.Close()

    print(f"\nAnalysis complete! Output saved to {args.output}")
    print(f"Created {len(all_canvases)} plots")


if __name__ == "__main__":
    main()
