#!/usr/bin/env python3
"""
Track Efficiency Analysis Script

Analyzes single-track signal efficiency from TrackAnalyzer output.
Uses uproot for reading and ROOT/cmsstyle for plotting.

Usage:
    python trackEfficiencyAnalysis.py -i track_ntuple.root
    python trackEfficiencyAnalysis.py -i track_ntuple.root -o my_output.root
    python trackEfficiencyAnalysis.py -i track_ntuple.root --deltaR-cut 0.05
"""

import argparse
import uproot
import numpy as np
import ROOT
import cmsstyle as CMS
from array import array


def load_data(filename, tree_path=None):
    """Load track data from ROOT file produced by TrackAnalyzer"""
    file = uproot.open(filename)

    # Try common tree paths
    if tree_path:
        tree = file[tree_path]
    else:
        possible_paths = ['trackAnalyzer/tree', 'tree']
        tree = None
        for path in possible_paths:
            if path in file:
                tree = file[path]
                print(f"Found tree at: {path}")
                break

        if tree is None:
            print(f"Available keys: {file.keys()}")
            raise KeyError(f"Could not find tree. Try specifying --tree-path")

    # Track branches
    track_branches = [
        'Track_nTotal', 'Track_pt', 'Track_eta', 'Track_phi',
        'Track_dxy', 'Track_dxySignificance', 'Track_sip2D', 'Track_sip3D',
        'Track_nValidHits', 'Track_nValidPixelHits', 'Track_normalizedChi2',
        'Track_isSignal', 'Track_isSignalElectron', 'Track_isSignalMuon',
        'Track_genMatchDeltaR', 'Track_genMatchPt', 'Track_genMatchDxy',
        'Track_genMatchPdgId', 'Track_genMatchRelPtDiff',
    ]

    # Signal gen branches
    gen_branches = [
        'SignalGen_nTotal', 'SignalGen_nElectron', 'SignalGen_nMuon',
        'SignalGen_pt', 'SignalGen_eta', 'SignalGen_phi', 'SignalGen_dxy',
        'SignalGen_pdgId', 'SignalGen_isMatched', 'SignalGen_matchedTrackPt',
        'SignalGen_matchedDeltaR',
    ]

    all_branches = track_branches + gen_branches

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


def analyze_track_efficiency(data, dxy_min=1e-5, dxy_max=1.0, dxy_nbins=21):
    """Analyze single-track efficiency from gen particle perspective"""

    # Flatten jagged arrays for gen particles
    all_gen_dxy = []
    all_gen_pt = []
    all_gen_eta = []
    all_gen_pdgId = []
    all_gen_isMatched = []
    all_gen_matchedDeltaR = []

    n_events = len(data['SignalGen_dxy'])

    for event_idx in range(n_events):
        gen_dxy = data['SignalGen_dxy'][event_idx]
        gen_pt = data['SignalGen_pt'][event_idx]
        gen_eta = data['SignalGen_eta'][event_idx]
        gen_pdgId = data['SignalGen_pdgId'][event_idx]
        gen_isMatched = data['SignalGen_isMatched'][event_idx]
        gen_matchedDeltaR = data['SignalGen_matchedDeltaR'][event_idx]

        for i in range(len(gen_dxy)):
            all_gen_dxy.append(gen_dxy[i])
            all_gen_pt.append(gen_pt[i])
            all_gen_eta.append(gen_eta[i])
            all_gen_pdgId.append(gen_pdgId[i])
            all_gen_isMatched.append(gen_isMatched[i])
            all_gen_matchedDeltaR.append(gen_matchedDeltaR[i])

    all_gen_dxy = np.array(all_gen_dxy)
    all_gen_pt = np.array(all_gen_pt)
    all_gen_eta = np.array(all_gen_eta)
    all_gen_pdgId = np.array(all_gen_pdgId)
    all_gen_isMatched = np.array(all_gen_isMatched)
    all_gen_matchedDeltaR = np.array(all_gen_matchedDeltaR)

    # Separate electrons and muons
    is_electron = np.abs(all_gen_pdgId) == 11
    is_muon = np.abs(all_gen_pdgId) == 13

    # Define bins (configurable dxy range)
    dxy_bins = np.logspace(np.log10(dxy_min), np.log10(dxy_max), dxy_nbins)
    pt_bins = np.array([0, 5, 10, 15, 20, 30, 50, 100, 200])
    eta_bins = np.linspace(-2.5, 2.5, 21)

    results = {
        'all': {},
        'electron': {},
        'muon': {},
    }

    # Calculate efficiencies
    for label, mask in [('all', np.ones(len(all_gen_dxy), dtype=bool)),
                        ('electron', is_electron),
                        ('muon', is_muon)]:

        if np.sum(mask) == 0:
            continue

        # Efficiency vs dxy
        centers, eff, err, n_total, n_pass = calculate_efficiency(
            all_gen_dxy[mask], all_gen_isMatched[mask], dxy_bins
        )
        results[label]['dxy'] = {
            'centers': centers, 'efficiency': eff, 'error': err,
            'n_total': n_total, 'n_pass': n_pass
        }

        # Efficiency vs pt
        centers, eff, err, n_total, n_pass = calculate_efficiency(
            all_gen_pt[mask], all_gen_isMatched[mask], pt_bins
        )
        results[label]['pt'] = {
            'centers': centers, 'efficiency': eff, 'error': err,
            'n_total': n_total, 'n_pass': n_pass
        }

        # Efficiency vs eta
        centers, eff, err, n_total, n_pass = calculate_efficiency(
            all_gen_eta[mask], all_gen_isMatched[mask], eta_bins
        )
        results[label]['eta'] = {
            'centers': centers, 'efficiency': eff, 'error': err,
            'n_total': n_total, 'n_pass': n_pass
        }

    # Summary statistics
    n_total_gen = len(all_gen_dxy)
    n_matched = np.sum(all_gen_isMatched)
    n_electrons = np.sum(is_electron)
    n_muons = np.sum(is_muon)
    n_electrons_matched = np.sum(all_gen_isMatched[is_electron])
    n_muons_matched = np.sum(all_gen_isMatched[is_muon])

    print(f"\n=== Track Efficiency Summary ===")
    print(f"Total signal gen particles: {n_total_gen}")
    print(f"  Electrons: {n_electrons} ({n_electrons_matched} matched, {100*n_electrons_matched/max(1,n_electrons):.1f}%)")
    print(f"  Muons: {n_muons} ({n_muons_matched} matched, {100*n_muons_matched/max(1,n_muons):.1f}%)")
    print(f"Overall efficiency: {100*n_matched/max(1,n_total_gen):.1f}%")

    return results, {
        'n_total': n_total_gen,
        'n_matched': n_matched,
        'n_electrons': n_electrons,
        'n_muons': n_muons,
        'n_electrons_matched': n_electrons_matched,
        'n_muons_matched': n_muons_matched,
        'all_gen_dxy': all_gen_dxy,
        'all_gen_matchedDeltaR': all_gen_matchedDeltaR,
        'all_gen_isMatched': all_gen_isMatched,
    }


def analyze_deltaR_distribution(data):
    """Analyze deltaR distribution for matched and unmatched tracks"""

    all_deltaR = []
    all_isMatched = []
    all_gen_dxy = []

    n_events = len(data['SignalGen_dxy'])

    for event_idx in range(n_events):
        gen_matchedDeltaR = data['SignalGen_matchedDeltaR'][event_idx]
        gen_isMatched = data['SignalGen_isMatched'][event_idx]
        gen_dxy = data['SignalGen_dxy'][event_idx]

        for i in range(len(gen_matchedDeltaR)):
            if gen_matchedDeltaR[i] >= 0:  # Has a match (even if not passing cut)
                all_deltaR.append(gen_matchedDeltaR[i])
                all_isMatched.append(gen_isMatched[i])
                all_gen_dxy.append(gen_dxy[i])

    return np.array(all_deltaR), np.array(all_isMatched), np.array(all_gen_dxy)


def create_efficiency_vs_dxy_plot(results, output_name, dxy_min=1e-5, dxy_max=1.0, dxy_nbins=21):
    """Create efficiency vs gen dxy plot"""

    COLOR_ALL = ROOT.kBlack
    COLOR_ELECTRON = ROOT.kOrange + 6
    COLOR_MUON = ROOT.kAzure + 1

    x_label = 'd_{xy}^{gen} [cm]'
    y_label = 'Single-Track Efficiency'

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

    # All particles
    if 'all' in results and 'dxy' in results['all']:
        r = results['all']['dxy']
        g = create_tgraph(r['centers'], r['efficiency'], r['error'],
                          'efficiency_vs_dxy_all', 'All Signal')
        g.SetMarkerStyle(20)
        g.SetMarkerSize(1.2)
        g.SetMarkerColor(COLOR_ALL)
        g.SetLineColor(COLOR_ALL)
        g.Draw('P SAME')
        graphs.append(('All', g, COLOR_ALL))

    # Electrons
    if 'electron' in results and 'dxy' in results['electron']:
        r = results['electron']['dxy']
        g = create_tgraph(r['centers'], r['efficiency'], r['error'],
                          'efficiency_vs_dxy_electron', 'Electrons')
        g.SetMarkerStyle(22)
        g.SetMarkerSize(1.2)
        g.SetMarkerColor(COLOR_ELECTRON)
        g.SetLineColor(COLOR_ELECTRON)
        g.Draw('P SAME')
        graphs.append(('Electrons', g, COLOR_ELECTRON))

    # Muons
    if 'muon' in results and 'dxy' in results['muon']:
        r = results['muon']['dxy']
        g = create_tgraph(r['centers'], r['efficiency'], r['error'],
                          'efficiency_vs_dxy_muon', 'Muons')
        g.SetMarkerStyle(23)
        g.SetMarkerSize(1.2)
        g.SetMarkerColor(COLOR_MUON)
        g.SetLineColor(COLOR_MUON)
        g.Draw('P SAME')
        graphs.append(('Muons', g, COLOR_MUON))

    # Legend
    n_entries = len(graphs)
    legend_height = 0.05 * n_entries
    legend = ROOT.TLegend(0.55, 0.91 - legend_height, 0.88, 0.91)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for label, g, _ in graphs:
        legend.AddEntry(g, label, 'p')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [g for _, g, _ in graphs] + [h_dummy, legend, cms_label]
    return canvas, all_objects


def create_efficiency_vs_pt_plot(results, output_name):
    """Create efficiency vs gen pt plot"""

    COLOR_ALL = ROOT.kBlack
    COLOR_ELECTRON = ROOT.kOrange + 6
    COLOR_MUON = ROOT.kAzure + 1

    x_label = 'p_{T}^{gen} [GeV]'
    y_label = 'Single-Track Efficiency'

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

    for label, key, color, marker in [('All', 'all', COLOR_ALL, 20),
                                       ('Electrons', 'electron', COLOR_ELECTRON, 22),
                                       ('Muons', 'muon', COLOR_MUON, 23)]:
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
    legend_height = 0.05 * n_entries
    legend = ROOT.TLegend(0.55, 0.25, 0.88, 0.25 + legend_height)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for label, g, _ in graphs:
        legend.AddEntry(g, label, 'p')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [g for _, g, _ in graphs] + [h_dummy, legend, cms_label]
    return canvas, all_objects


def create_efficiency_vs_eta_plot(results, output_name):
    """Create efficiency vs gen eta plot"""

    COLOR_ALL = ROOT.kBlack
    COLOR_ELECTRON = ROOT.kOrange + 6
    COLOR_MUON = ROOT.kAzure + 1

    x_label = '#eta^{gen}'
    y_label = 'Single-Track Efficiency'

    canvas = CMS.cmsCanvas(output_name, -2.5, 2.5, 0.0, 1.2, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)

    h_dummy = ROOT.TH1F(f'h_dummy_{output_name}', f';{x_label};{y_label}', 100, -2.5, 2.5)
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

    for label, key, color, marker in [('All', 'all', COLOR_ALL, 20),
                                       ('Electrons', 'electron', COLOR_ELECTRON, 22),
                                       ('Muons', 'muon', COLOR_MUON, 23)]:
        if key in results and 'eta' in results[key]:
            r = results[key]['eta']
            g = create_tgraph(r['centers'], r['efficiency'], r['error'],
                              f'efficiency_vs_eta_{key}', label)
            g.SetMarkerStyle(marker)
            g.SetMarkerSize(1.2)
            g.SetMarkerColor(color)
            g.SetLineColor(color)
            g.Draw('P SAME')
            graphs.append((label, g, color))

    n_entries = len(graphs)
    legend_height = 0.05 * n_entries
    legend = ROOT.TLegend(0.55, 0.25, 0.88, 0.25 + legend_height)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for label, g, _ in graphs:
        legend.AddEntry(g, label, 'p')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [g for _, g, _ in graphs] + [h_dummy, legend, cms_label]
    return canvas, all_objects


def create_deltaR_distribution_plot(deltaR_values, is_matched, output_name):
    """Create deltaR distribution plot"""

    x_label = '#DeltaR(track, gen)'
    y_label = 'Entries'

    canvas = CMS.cmsCanvas(output_name, 0, 0.5, 0.1, 1e6, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetLogy()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)

    # Create histograms
    h_all = ROOT.TH1F('h_deltaR_all', f';{x_label};{y_label}', 100, 0, 0.5)
    h_matched = ROOT.TH1F('h_deltaR_matched', f';{x_label};{y_label}', 100, 0, 0.5)
    h_unmatched = ROOT.TH1F('h_deltaR_unmatched', f';{x_label};{y_label}', 100, 0, 0.5)

    for i, dr in enumerate(deltaR_values):
        h_all.Fill(dr)
        if is_matched[i]:
            h_matched.Fill(dr)
        else:
            h_unmatched.Fill(dr)

    h_all.SetStats(0)
    h_all.SetLineColor(ROOT.kBlack)
    h_all.SetLineWidth(2)

    h_matched.SetLineColor(ROOT.kGreen + 2)
    h_matched.SetFillColor(ROOT.kGreen + 2)
    h_matched.SetFillStyle(3004)
    h_matched.SetLineWidth(2)

    h_unmatched.SetLineColor(ROOT.kRed + 1)
    h_unmatched.SetFillColor(ROOT.kRed + 1)
    h_unmatched.SetFillStyle(3005)
    h_unmatched.SetLineWidth(2)

    max_val = max(h_all.GetMaximum(), 1)
    h_all.SetMinimum(0.1)
    h_all.SetMaximum(max_val * 5)
    h_all.GetXaxis().SetTitleSize(0.05)
    h_all.GetXaxis().SetLabelSize(0.05)
    h_all.GetXaxis().SetTitleOffset(1.25)
    h_all.GetXaxis().CenterTitle(True)
    h_all.GetYaxis().SetTitleSize(0.05)
    h_all.GetYaxis().SetLabelSize(0.05)
    h_all.GetYaxis().SetTitleOffset(1.25)
    h_all.GetYaxis().CenterTitle(True)

    h_all.Draw('HIST')
    h_matched.Draw('HIST SAME')
    h_unmatched.Draw('HIST SAME')

    # Vertical line at deltaR cut (0.02)
    line = ROOT.TLine(0.02, 0.1, 0.02, max_val * 5)
    line.SetLineColor(ROOT.kBlue)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()

    legend = ROOT.TLegend(0.55, 0.70, 0.88, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_all, 'All matches', 'l')
    legend.AddEntry(h_matched, 'Passing #DeltaR cut', 'f')
    legend.AddEntry(h_unmatched, 'Failing #DeltaR cut', 'f')
    legend.AddEntry(line, '#DeltaR = 0.02 cut', 'l')
    legend.Draw()

    cms_label = draw_cms_labels()

    all_objects = [h_all, h_matched, h_unmatched, line, legend, cms_label]
    return canvas, all_objects


def create_2d_efficiency_plot(summary_data, output_name, dxy_min=1e-5, dxy_max=1.0):
    """Create 2D efficiency plot: dxy vs deltaR"""

    gen_dxy = summary_data['all_gen_dxy']
    gen_deltaR = summary_data['all_gen_matchedDeltaR']
    is_matched = summary_data['all_gen_isMatched']

    # Only include entries with valid deltaR
    valid_mask = gen_deltaR >= 0
    gen_dxy = gen_dxy[valid_mask]
    gen_deltaR = gen_deltaR[valid_mask]
    is_matched = is_matched[valid_mask]

    x_label = 'd_{xy}^{gen} [cm]'
    y_label = '#DeltaR(track, gen)'

    canvas = CMS.cmsCanvas(output_name, dxy_min, dxy_max, 0, 0.3, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetLogx()
    canvas.SetLogz()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.15)

    dxy_bins = np.logspace(np.log10(dxy_min), np.log10(dxy_max), 31)
    deltaR_bins = np.linspace(0, 0.3, 31)

    h2d = ROOT.TH2F('h2d_dxy_deltaR', f';{x_label};{y_label}',
                    len(dxy_bins) - 1, dxy_bins,
                    len(deltaR_bins) - 1, deltaR_bins)
    h2d.SetStats(0)

    for i in range(len(gen_dxy)):
        h2d.Fill(gen_dxy[i], gen_deltaR[i])

    h2d.GetXaxis().SetMoreLogLabels()
    h2d.GetXaxis().SetNoExponent()
    h2d.GetXaxis().SetTitleSize(0.05)
    h2d.GetXaxis().SetLabelSize(0.05)
    h2d.GetXaxis().SetTitleOffset(1.25)
    h2d.GetXaxis().CenterTitle(True)
    h2d.GetYaxis().SetTitleSize(0.05)
    h2d.GetYaxis().SetLabelSize(0.05)
    h2d.GetYaxis().SetTitleOffset(1.25)
    h2d.GetYaxis().CenterTitle(True)
    h2d.GetZaxis().SetTitleSize(0.05)
    h2d.GetZaxis().SetLabelSize(0.04)
    h2d.Draw('COLZ')

    # Horizontal line at deltaR cut
    line = ROOT.TLine(dxy_min, 0.02, dxy_max, 0.02)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()

    cms_label = draw_cms_labels()

    all_objects = [h2d, line, cms_label]
    return canvas, all_objects


def create_deltaR_vs_relPt_plot(data, output_name, deltaR_max=0.1):
    """Create 2D plot: deltaR vs relative pT difference"""

    # Extract track-level gen match info
    all_deltaR = []
    all_relPtDiff = []

    n_events = len(data['Track_genMatchDeltaR'])

    for event_idx in range(n_events):
        deltaR = data['Track_genMatchDeltaR'][event_idx]
        relPtDiff = data['Track_genMatchRelPtDiff'][event_idx]

        for i in range(len(deltaR)):
            # Only include valid matches (deltaR >= 0)
            if deltaR[i] >= 0:
                all_deltaR.append(deltaR[i])
                all_relPtDiff.append(relPtDiff[i])

    all_deltaR = np.array(all_deltaR)
    all_relPtDiff = np.array(all_relPtDiff)

    print(f"2D deltaR vs relPt: {len(all_deltaR)} matched tracks")

    x_label = '(p_{T}^{gen} - p_{T}^{reco}) / p_{T}^{gen}'
    y_label = '#DeltaR(track, gen)'

    canvas = CMS.cmsCanvas(output_name, -1, 1, 0, deltaR_max, x_label, y_label,
                           square=False, extraSpace=0.01, iPos=0)
    canvas.SetCanvasSize(800, 600)
    canvas.SetLogz()
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.15)

    relPt_bins = np.linspace(-1, 1, 51)
    deltaR_bins = np.linspace(0, deltaR_max, 51)

    h2d = ROOT.TH2F('h2d_deltaR_vs_relPt', f';{x_label};{y_label}',
                    len(relPt_bins) - 1, relPt_bins,
                    len(deltaR_bins) - 1, deltaR_bins)
    h2d.SetStats(0)

    for i in range(len(all_deltaR)):
        h2d.Fill(all_relPtDiff[i], all_deltaR[i])

    h2d.GetXaxis().SetTitleSize(0.05)
    h2d.GetXaxis().SetLabelSize(0.05)
    h2d.GetXaxis().SetTitleOffset(1.25)
    h2d.GetXaxis().CenterTitle(True)
    h2d.GetYaxis().SetTitleSize(0.05)
    h2d.GetYaxis().SetLabelSize(0.05)
    h2d.GetYaxis().SetTitleOffset(1.25)
    h2d.GetYaxis().CenterTitle(True)
    h2d.GetZaxis().SetTitleSize(0.05)
    h2d.GetZaxis().SetLabelSize(0.04)
    h2d.Draw('COLZ')

    # Horizontal line at deltaR cut (0.02)
    line_deltaR = ROOT.TLine(-1, 0.02, 1, 0.02)
    line_deltaR.SetLineColor(ROOT.kRed)
    line_deltaR.SetLineStyle(2)
    line_deltaR.SetLineWidth(2)
    line_deltaR.Draw()

    # Vertical line at relPt = 0
    line_relPt = ROOT.TLine(0, 0, 0, deltaR_max)
    line_relPt.SetLineColor(ROOT.kBlack)
    line_relPt.SetLineStyle(2)
    line_relPt.SetLineWidth(1)
    line_relPt.Draw()

    cms_label = draw_cms_labels()

    all_objects = [h2d, line_deltaR, line_relPt, cms_label]
    return canvas, all_objects


def main():
    parser = argparse.ArgumentParser(description='Track efficiency analysis from TrackAnalyzer output')
    parser.add_argument('-i', '--input', required=True,
                        help='Input ROOT file (from TrackAnalyzer)')
    parser.add_argument('-o', '--output', default='track_efficiency_analysis.root',
                        help='Output ROOT file name')
    parser.add_argument('--tree-path', default=None,
                        help='Path to tree in ROOT file (default: auto-detect)')
    parser.add_argument('--dxy-min', type=float, default=1e-5,
                        help='Minimum dxy for efficiency plot in cm (default: 1e-5)')
    parser.add_argument('--dxy-max', type=float, default=1.0,
                        help='Maximum dxy for efficiency plot in cm (default: 1.0)')
    parser.add_argument('--dxy-nbins', type=int, default=21,
                        help='Number of dxy bins (default: 21)')
    parser.add_argument('--deltaR-max', type=float, default=0.1,
                        help='Maximum deltaR for 2D relPt plot (default: 0.1)')
    args = parser.parse_args()

    # Set CMS style
    CMS.SetExtraText("Simulation")
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # Load data
    print(f"Loading data from {args.input}...")
    data = load_data(args.input, args.tree_path)

    # Analyze efficiency
    results, summary = analyze_track_efficiency(data, args.dxy_min, args.dxy_max, args.dxy_nbins)

    # Analyze deltaR distribution
    deltaR_values, is_matched, gen_dxy = analyze_deltaR_distribution(data)

    # Collect all canvases and objects
    all_canvases = []
    all_objects = []

    def save(canvas, objects):
        all_canvases.append(canvas)
        all_objects.extend(objects)

    # Create plots
    print("\nCreating plots...")

    # Efficiency vs dxy
    c, o = create_efficiency_vs_dxy_plot(results, 'track_efficiency_vs_dxy',
                                          args.dxy_min, args.dxy_max, args.dxy_nbins)
    save(c, o)

    # Efficiency vs pt
    c, o = create_efficiency_vs_pt_plot(results, 'track_efficiency_vs_pt')
    save(c, o)

    # Efficiency vs eta
    c, o = create_efficiency_vs_eta_plot(results, 'track_efficiency_vs_eta')
    save(c, o)

    # DeltaR distribution
    if len(deltaR_values) > 0:
        c, o = create_deltaR_distribution_plot(deltaR_values, is_matched, 'deltaR_distribution')
        save(c, o)

    # 2D dxy vs deltaR
    if len(summary['all_gen_dxy']) > 0:
        c, o = create_2d_efficiency_plot(summary, 'dxy_vs_deltaR_2d',
                                          args.dxy_min, args.dxy_max)
        save(c, o)

    # 2D deltaR vs relative pT difference
    if 'Track_genMatchDeltaR' in data and 'Track_genMatchRelPtDiff' in data:
        c, o = create_deltaR_vs_relPt_plot(data, 'deltaR_vs_relPt_2d',
                                           deltaR_max=args.deltaR_max)
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
