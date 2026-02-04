#!/usr/bin/env python3
"""
Track Signal vs Background Comparison Script

Compares Track_isSignal vs !Track_isSignal distributions for key variables.
Designed for large files - uses chunked reading and basic ROOT styling (no cmsstyle).

Usage:
    python trackSignalComparison.py -i track_ntuple.root
    python trackSignalComparison.py -i track_ntuple.root -o output.root --chunk-size 100000
"""

import argparse
import uproot
import numpy as np
import ROOT
from array import array


def setup_root_style():
    """Setup basic ROOT style without cmsstyle"""
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.08)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetTitleXOffset(1.2)
    ROOT.gStyle.SetTitleYOffset(1.4)
    ROOT.gStyle.SetHistLineWidth(2)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)


def load_data_chunked(filename, tree_path=None, chunk_size=500000):
    """Load track data from ROOT file using chunked reading for large files"""

    file = uproot.open(filename)

    # Find tree
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
            raise KeyError("Could not find tree. Try specifying --tree-path")

    # Branches we need
    branches = [
        'Track_pt', 'Track_ptError', 'Track_eta',
        'Track_normalizedChi2', 'Track_sip2D', 'Track_isSignal'
    ]

    # Check which branches exist
    available = tree.keys()
    branches = [b for b in branches if b in available]

    if 'Track_isSignal' not in branches:
        raise KeyError("Track_isSignal branch not found - was this run with hasGenInfo=True?")

    print(f"Loading branches: {branches}")
    print(f"Total entries: {tree.num_entries}")
    print(f"Chunk size: {chunk_size}")

    # Accumulate data in chunks
    all_data = {b: [] for b in branches}

    n_chunks = 0
    for chunk in tree.iterate(branches, step_size=chunk_size, library='np'):
        n_chunks += 1
        for branch in branches:
            # Flatten jagged arrays (each event has multiple tracks)
            arr = chunk[branch]
            if hasattr(arr, 'flatten'):
                # Already flat numpy array from iterate
                flat = np.concatenate([np.asarray(x) for x in arr])
            else:
                flat = np.concatenate([np.asarray(x) for x in arr])
            all_data[branch].append(flat)

        if n_chunks % 10 == 0:
            print(f"  Processed {n_chunks} chunks...")

    # Concatenate all chunks
    print(f"Concatenating {n_chunks} chunks...")
    for branch in branches:
        all_data[branch] = np.concatenate(all_data[branch])

    n_tracks = len(all_data['Track_pt'])
    print(f"Total tracks loaded: {n_tracks:,}")

    return all_data


def create_comparison_histogram(sig_data, bkg_data, name, title, nbins, xmin, xmax,
                                 x_label, cut_value=None, cut_direction=None):
    """Create signal vs background comparison histogram"""

    # Create histograms
    h_sig = ROOT.TH1F(f'{name}_signal', title, nbins, xmin, xmax)
    h_bkg = ROOT.TH1F(f'{name}_background', title, nbins, xmin, xmax)

    # Fill histograms
    for val in sig_data:
        h_sig.Fill(val)
    for val in bkg_data:
        h_bkg.Fill(val)

    # Normalize to unit area
    if h_sig.Integral() > 0:
        h_sig.Scale(1.0 / h_sig.Integral())
    if h_bkg.Integral() > 0:
        h_bkg.Scale(1.0 / h_bkg.Integral())

    # Style
    h_sig.SetLineColor(ROOT.kRed + 1)
    h_sig.SetFillColor(ROOT.kRed + 1)
    h_sig.SetFillStyle(3004)
    h_sig.SetLineWidth(2)

    h_bkg.SetLineColor(ROOT.kBlue + 1)
    h_bkg.SetFillColor(ROOT.kBlue + 1)
    h_bkg.SetFillStyle(3005)
    h_bkg.SetLineWidth(2)

    # Axis labels
    h_sig.GetXaxis().SetTitle(x_label)
    h_sig.GetYaxis().SetTitle('Normalized')
    h_bkg.GetXaxis().SetTitle(x_label)
    h_bkg.GetYaxis().SetTitle('Normalized')

    # Calculate cut efficiencies if cut is specified
    cut_info = None
    if cut_value is not None and cut_direction is not None:
        n_sig = len(sig_data)
        n_bkg = len(bkg_data)

        if cut_direction == '<':
            sig_pass = np.sum(sig_data < cut_value)
            bkg_pass = np.sum(bkg_data < cut_value)
        elif cut_direction == '>':
            sig_pass = np.sum(sig_data > cut_value)
            bkg_pass = np.sum(bkg_data > cut_value)
        elif cut_direction == 'abs<':
            sig_pass = np.sum(np.abs(sig_data) < cut_value)
            bkg_pass = np.sum(np.abs(bkg_data) < cut_value)

        sig_eff = sig_pass / n_sig if n_sig > 0 else 0
        bkg_eff = bkg_pass / n_bkg if n_bkg > 0 else 0
        bkg_rej = 1 - bkg_eff

        cut_info = {
            'cut_value': cut_value,
            'cut_direction': cut_direction,
            'sig_eff': sig_eff,
            'bkg_eff': bkg_eff,
            'bkg_rej': bkg_rej,
            'n_sig': n_sig,
            'n_bkg': n_bkg,
            'sig_pass': sig_pass,
            'bkg_pass': bkg_pass,
        }

    return h_sig, h_bkg, cut_info


def draw_comparison_canvas(h_sig, h_bkg, name, cut_info=None, logy=False):
    """Draw comparison canvas with legend and cut line"""

    canvas = ROOT.TCanvas(name, name, 800, 600)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.08)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)

    if logy:
        canvas.SetLogy()

    # Get max for y-axis
    max_val = max(h_sig.GetMaximum(), h_bkg.GetMaximum())
    h_sig.SetMaximum(max_val * 1.4)
    h_sig.SetMinimum(0.0001 if logy else 0)

    h_sig.Draw('HIST')
    h_bkg.Draw('HIST SAME')

    # Legend
    legend = ROOT.TLegend(0.55, 0.75, 0.88, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(h_sig, f'Signal ({h_sig.GetEntries():.0f})', 'f')
    legend.AddEntry(h_bkg, f'Background ({h_bkg.GetEntries():.0f})', 'f')
    legend.Draw()

    # Cut line and efficiency text
    line = None
    cut_text = None
    if cut_info is not None:
        # Draw cut line
        cut_val = cut_info['cut_value']
        if cut_info['cut_direction'] == 'abs<':
            # Draw two lines for absolute value cut
            line = ROOT.TLine(cut_val, 0, cut_val, max_val * 1.2)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            line2 = ROOT.TLine(-cut_val, 0, -cut_val, max_val * 1.2)
            line2.SetLineColor(ROOT.kBlack)
            line2.SetLineStyle(2)
            line2.SetLineWidth(2)
            line2.Draw()
        else:
            line = ROOT.TLine(cut_val, 0, cut_val, max_val * 1.2)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()

        # Efficiency text
        cut_text = ROOT.TLatex()
        cut_text.SetNDC()
        cut_text.SetTextSize(0.032)
        cut_text.SetTextFont(42)
        cut_text.DrawLatex(0.15, 0.85, f"Cut: {cut_info['cut_direction']} {cut_val}")
        cut_text.DrawLatex(0.15, 0.80, f"Signal eff: {cut_info['sig_eff']*100:.1f}%")
        cut_text.DrawLatex(0.15, 0.75, f"Bkg rejection: {cut_info['bkg_rej']*100:.1f}%")

    # CMS label
    cms_label = ROOT.TLatex()
    cms_label.SetNDC()
    cms_label.SetTextSize(0.05)
    cms_label.SetTextFont(61)
    cms_label.DrawLatex(0.12, 0.93, "CMS")
    cms_label.SetTextFont(52)
    cms_label.SetTextSize(0.04)
    cms_label.DrawLatex(0.22, 0.93, "Simulation")

    return canvas, [h_sig, h_bkg, legend, line, cut_text, cms_label]


def main():
    parser = argparse.ArgumentParser(description='Track signal vs background comparison')
    parser.add_argument('-i', '--input', required=True,
                        help='Input ROOT file (from TrackAnalyzer)')
    parser.add_argument('-o', '--output', default='track_signal_comparison.root',
                        help='Output ROOT file name')
    parser.add_argument('--tree-path', default=None,
                        help='Path to tree in ROOT file (default: auto-detect)')
    parser.add_argument('--chunk-size', type=int, default=500000,
                        help='Chunk size for reading large files (default: 500000)')
    args = parser.parse_args()

    setup_root_style()

    # Load data
    print(f"\nLoading data from {args.input}...")
    data = load_data_chunked(args.input, args.tree_path, args.chunk_size)

    # Separate signal and background
    is_signal = data['Track_isSignal'].astype(bool)
    is_bkg = ~is_signal

    n_sig = np.sum(is_signal)
    n_bkg = np.sum(is_bkg)
    print(f"\nSignal tracks: {n_sig:,}")
    print(f"Background tracks: {n_bkg:,}")
    print(f"Signal fraction: {100*n_sig/(n_sig+n_bkg):.2f}%")

    # Calculate derived variables
    pt = data['Track_pt']
    pt_res = data['Track_ptError'] / np.maximum(data['Track_pt'], 1e-10)
    norm_chi2 = data['Track_normalizedChi2']
    sip2d = data['Track_sip2D']
    eta = data['Track_eta']

    # Define plots: (variable, name, title, nbins, xmin, xmax, xlabel, cut_val, cut_dir, logy)
    plots = [
        (pt, 'pt', 'Track p_{T}', 100, 0, 50, 'p_{T} [GeV]', 5, '<', True),
        (pt_res, 'ptResolution', 'Track p_{T} Resolution', 100, 0, 0.5, '#sigma_{p_{T}}/p_{T}', 0.1, '<', True),
        (norm_chi2, 'normalizedChi2', 'Track #chi^{2}/ndof', 100, 0, 20, '#chi^{2}/ndof', 10, '<', True),
        (sip2d, 'sip2D', 'Track SIP_{2D}', 100, -50, 50, 'SIP_{2D}', 10, 'abs<', True),
        (eta, 'eta', 'Track #eta', 100, -3, 3, '#eta', None, None, False),
    ]

    # Create output file
    output_file = ROOT.TFile(args.output, 'RECREATE')

    all_canvases = []
    all_objects = []

    print("\n=== Cut Efficiency Summary ===")
    print(f"{'Variable':<20} {'Cut':<15} {'Sig Eff':<12} {'Bkg Rej':<12}")
    print("-" * 60)

    for var, name, title, nbins, xmin, xmax, xlabel, cut_val, cut_dir, logy in plots:
        # Get signal and background data
        sig_data = var[is_signal]
        bkg_data = var[is_bkg]

        # Create histograms
        h_sig, h_bkg, cut_info = create_comparison_histogram(
            sig_data, bkg_data, name, title, nbins, xmin, xmax, xlabel, cut_val, cut_dir
        )

        # Draw canvas
        canvas, objects = draw_comparison_canvas(h_sig, h_bkg, f'c_{name}', cut_info, logy)

        # Print summary
        if cut_info:
            cut_str = f"{cut_dir} {cut_val}"
            print(f"{name:<20} {cut_str:<15} {cut_info['sig_eff']*100:>6.1f}%      {cut_info['bkg_rej']*100:>6.1f}%")
        else:
            print(f"{name:<20} {'N/A':<15} {'N/A':<12} {'N/A':<12}")

        # Save to file
        canvas.Write()
        h_sig.Write()
        h_bkg.Write()

        all_canvases.append(canvas)
        all_objects.extend(objects)

    print("-" * 60)

    # Create combined efficiency plot
    print("\nCreating ROC-style summary...")

    # Summary canvas with all distributions
    c_summary = ROOT.TCanvas('c_summary', 'Summary', 1200, 800)
    c_summary.Divide(3, 2)

    for i, (var, name, title, nbins, xmin, xmax, xlabel, cut_val, cut_dir, logy) in enumerate(plots):
        c_summary.cd(i + 1)
        if logy:
            ROOT.gPad.SetLogy()

        sig_data = var[is_signal]
        bkg_data = var[is_bkg]

        h_sig, h_bkg, _ = create_comparison_histogram(
            sig_data, bkg_data, f'{name}_sum', title, nbins, xmin, xmax, xlabel
        )

        max_val = max(h_sig.GetMaximum(), h_bkg.GetMaximum())
        h_sig.SetMaximum(max_val * 1.5)
        h_sig.SetTitle(title)
        h_sig.Draw('HIST')
        h_bkg.Draw('HIST SAME')

        leg = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(h_sig, 'Signal', 'f')
        leg.AddEntry(h_bkg, 'Background', 'f')
        leg.Draw()

    c_summary.Write()

    output_file.Close()

    print(f"\nOutput saved to: {args.output}")
    print(f"Created {len(all_canvases) + 1} plots")


if __name__ == "__main__":
    main()
