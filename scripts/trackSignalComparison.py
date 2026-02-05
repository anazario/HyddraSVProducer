#!/usr/bin/env python3
"""
Track Signal vs Background Comparison Script

Compares Track_isSignal vs !Track_isSignal distributions for key variables.
Zooms into relevant regions to help identify optimal cut values.
Designed for large files - uses chunked reading and basic ROOT styling (no cmsstyle).

Usage:
    python trackSignalComparison.py -i track_ntuple.root
    python trackSignalComparison.py -i track_ntuple.root -o output.root --chunk-size 100000

    # Custom cuts
    python trackSignalComparison.py -i track_ntuple.root --pt-cut 10 --sip2d-cut 3

    # Invert sip2D cut (keep displaced tracks)
    python trackSignalComparison.py -i track_ntuple.root --sip2d-cut 3 --invert-sip2d
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


def create_comparison_histogram(sig_data, bkg_data, name, title, nbins, xmin, xmax, x_label):
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

    return h_sig, h_bkg


def draw_comparison_canvas(h_sig, h_bkg, name, logy=False):
    """Draw comparison canvas with legend"""

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

    # CMS label
    cms_label = ROOT.TLatex()
    cms_label.SetNDC()
    cms_label.SetTextSize(0.05)
    cms_label.SetTextFont(61)
    cms_label.DrawLatex(0.12, 0.93, "CMS")
    cms_label.SetTextFont(52)
    cms_label.SetTextSize(0.04)
    cms_label.DrawLatex(0.22, 0.93, "Simulation")

    return canvas, [h_sig, h_bkg, legend, cms_label]


def draw_comparison_canvas_with_cut(h_sig, h_bkg, name, cut_val, cut_dir, logy=False):
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
    h_sig.SetMaximum(max_val * 1.5)
    h_sig.SetMinimum(0.0001 if logy else 0)

    h_sig.Draw('HIST')
    h_bkg.Draw('HIST SAME')

    # Draw cut line(s)
    lines = []
    if cut_val is not None:
        ymin = 0.0001 if logy else 0
        ymax = max_val * 1.3

        if cut_dir in ['abs<', 'abs>=']:
            # Two lines for absolute value cut
            line1 = ROOT.TLine(cut_val, ymin, cut_val, ymax)
            line1.SetLineColor(ROOT.kBlack)
            line1.SetLineStyle(2)
            line1.SetLineWidth(2)
            line1.Draw()
            lines.append(line1)

            line2 = ROOT.TLine(-cut_val, ymin, -cut_val, ymax)
            line2.SetLineColor(ROOT.kBlack)
            line2.SetLineStyle(2)
            line2.SetLineWidth(2)
            line2.Draw()
            lines.append(line2)
        else:
            line = ROOT.TLine(cut_val, ymin, cut_val, ymax)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            lines.append(line)

    # Legend
    legend = ROOT.TLegend(0.55, 0.70, 0.88, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(h_sig, f'Signal ({h_sig.GetEntries():.0f})', 'f')
    legend.AddEntry(h_bkg, f'Background ({h_bkg.GetEntries():.0f})', 'f')
    if lines:
        legend.AddEntry(lines[0], 'Cut', 'l')
    legend.Draw()

    # Cut label
    cut_label = None
    if cut_val is not None:
        cut_label = ROOT.TLatex()
        cut_label.SetNDC()
        cut_label.SetTextSize(0.035)
        cut_label.SetTextFont(42)
        if cut_dir == '>':
            cut_label.DrawLatex(0.15, 0.85, f"Cut: > {cut_val}")
        elif cut_dir == '<':
            cut_label.DrawLatex(0.15, 0.85, f"Cut: < {cut_val}")
        elif cut_dir == 'abs<':
            cut_label.DrawLatex(0.15, 0.85, f"Cut: |x| < {cut_val}")
        elif cut_dir == 'abs>=':
            cut_label.DrawLatex(0.15, 0.85, f"Cut: |x| #geq {cut_val}")

    # CMS label
    cms_label = ROOT.TLatex()
    cms_label.SetNDC()
    cms_label.SetTextSize(0.05)
    cms_label.SetTextFont(61)
    cms_label.DrawLatex(0.12, 0.93, "CMS")
    cms_label.SetTextFont(52)
    cms_label.SetTextSize(0.04)
    cms_label.DrawLatex(0.22, 0.93, "Simulation (N-1)")

    return canvas, [h_sig, h_bkg, legend, cms_label, cut_label] + lines


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

    # Configurable cuts
    parser.add_argument('--pt-cut', type=float, default=5.0,
                        help='Minimum track pT cut in GeV (default: 5.0)')
    parser.add_argument('--pt-res-cut', type=float, default=0.08,
                        help='Maximum pT resolution cut (default: 0.08)')
    parser.add_argument('--chi2-cut', type=float, default=5.0,
                        help='Maximum normalized chi2 cut (default: 5.0)')
    parser.add_argument('--sip2d-cut', type=float, default=5.0,
                        help='SIP2D cut value (default: 5.0)')
    parser.add_argument('--invert-sip2d', action='store_true',
                        help='Invert sip2D cut: keep |sip2D| >= cut instead of < cut (for displaced tracks)')

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

    # Determine sip2D cut direction based on --invert-sip2d flag
    # Default (invert=False): keep |sip2D| < cut (prompt-like tracks)
    # Inverted (invert=True): keep |sip2D| >= cut (displaced tracks)
    sip2d_cut_dir = 'abs>=' if args.invert_sip2d else 'abs<'

    # Print cut configuration
    print(f"\n--- Cut Configuration ---")
    print(f"  pT > {args.pt_cut} GeV")
    print(f"  pT resolution < {args.pt_res_cut}")
    print(f"  normalized chi2 < {args.chi2_cut}")
    if args.invert_sip2d:
        print(f"  |sip2D| >= {args.sip2d_cut} (inverted - selecting displaced)")
    else:
        print(f"  |sip2D| < {args.sip2d_cut} (default - selecting prompt-like)")

    # Define plots: (variable, name, title, nbins, xmin, xmax, xlabel, logy, cut_val, cut_dir)
    # cut_dir: '>' means keep if var > cut, '<' means keep if var < cut,
    #          'abs<' means keep if |var| < cut, 'abs>=' means keep if |var| >= cut
    plots = [
        (pt, 'pt', 'Track p_{T}', 100, 0, 50, 'p_{T} [GeV]', True, args.pt_cut, '>'),
        (pt_res, 'ptResolution', 'Track p_{T} Resolution', 100, 0, 0.5, '#sigma_{p_{T}}/p_{T}', True, args.pt_res_cut, '<'),
        (norm_chi2, 'normalizedChi2', 'Track #chi^{2}/ndof', 100, 0, 20, '#chi^{2}/ndof', True, args.chi2_cut, '<'),
        (sip2d, 'sip2D', 'Track SIP_{2D}', 100, -50, 50, 'SIP_{2D}', False, args.sip2d_cut, sip2d_cut_dir),
        (eta, 'eta', 'Track #eta', 100, -3, 3, '#eta', False, None, None),
    ]

    # Build cut masks for each variable
    def apply_cut(var_data, cut_val, cut_dir):
        if cut_dir == '>':
            return var_data > cut_val
        elif cut_dir == '<':
            return var_data < cut_val
        elif cut_dir == 'abs<':
            return np.abs(var_data) < cut_val
        elif cut_dir == 'abs>=':
            return np.abs(var_data) >= cut_val
        return np.ones(len(var_data), dtype=bool)

    # Pre-compute all cut masks
    cut_masks = {}
    for var, name, title, nbins, xmin, xmax, xlabel, logy, cut_val, cut_dir in plots:
        if cut_val is not None:
            cut_masks[name] = {
                'signal': apply_cut(var[is_signal], cut_val, cut_dir),
                'background': apply_cut(var[is_bkg], cut_val, cut_dir),
            }

    # Create output file
    output_file = ROOT.TFile(args.output, 'RECREATE')

    all_canvases = []
    all_objects = []

    print("\nCreating N-1 plots (each variable with all OTHER cuts applied)...")

    for var, name, title, nbins, xmin, xmax, xlabel, logy, cut_val, cut_dir in plots:
        # Build N-1 mask: apply all cuts EXCEPT the current variable
        n_minus_1_sig = np.ones(n_sig, dtype=bool)
        n_minus_1_bkg = np.ones(n_bkg, dtype=bool)

        for other_name, masks in cut_masks.items():
            if other_name != name:  # Skip the current variable's cut
                n_minus_1_sig &= masks['signal']
                n_minus_1_bkg &= masks['background']

        # Get signal and background data with N-1 selection
        sig_data = var[is_signal][n_minus_1_sig]
        bkg_data = var[is_bkg][n_minus_1_bkg]

        # Create histograms
        h_sig, h_bkg = create_comparison_histogram(
            sig_data, bkg_data, f'{name}_Nminus1', f'{title} (N-1)', nbins, xmin, xmax, xlabel
        )

        # Draw canvas with cut line
        canvas, objects = draw_comparison_canvas_with_cut(
            h_sig, h_bkg, f'c_{name}_Nminus1', cut_val, cut_dir, logy
        )

        n_sig_nminus1 = len(sig_data)
        n_bkg_nminus1 = len(bkg_data)
        print(f"  {name}: {n_sig_nminus1:,} sig, {n_bkg_nminus1:,} bkg after other cuts")

        # Save to file
        canvas.Write()
        h_sig.Write()
        h_bkg.Write()

        all_canvases.append(canvas)
        all_objects.extend(objects)

    # =========================================================================
    # Cut Efficiency Analysis
    # =========================================================================
    print("\n" + "=" * 70)
    print("  CUT EFFICIENCY ANALYSIS")
    print("=" * 70)

    print(f"\nTotal tracks: {n_sig + n_bkg:,}")
    print(f"  Signal:     {n_sig:,} ({100*n_sig/(n_sig+n_bkg):.2f}%)")
    print(f"  Background: {n_bkg:,} ({100*n_bkg/(n_sig+n_bkg):.2f}%)")

    # Individual cut efficiencies
    print("\n--- Individual Cut Efficiencies ---")
    print(f"{'Cut':<25} {'Sig Eff':<12} {'Bkg Eff':<12} {'Bkg Rej':<12} {'Sig Yield':<12} {'Bkg Yield':<12}")
    print("-" * 85)

    eff_cut_masks_sig = {}
    eff_cut_masks_bkg = {}

    for var, name, title, nbins, xmin, xmax, xlabel, logy, cut_val, cut_dir in plots:
        if cut_val is None:
            continue

        sig_data = var[is_signal]
        bkg_data = var[is_bkg]

        if cut_dir == '>':
            sig_pass = sig_data > cut_val
            bkg_pass = bkg_data > cut_val
            cut_str = f"{name} > {cut_val}"
        elif cut_dir == '<':
            sig_pass = sig_data < cut_val
            bkg_pass = bkg_data < cut_val
            cut_str = f"{name} < {cut_val}"
        elif cut_dir == 'abs<':
            sig_pass = np.abs(sig_data) < cut_val
            bkg_pass = np.abs(bkg_data) < cut_val
            cut_str = f"|{name}| < {cut_val}"
        elif cut_dir == 'abs>=':
            sig_pass = np.abs(sig_data) >= cut_val
            bkg_pass = np.abs(bkg_data) >= cut_val
            cut_str = f"|{name}| >= {cut_val}"

        eff_cut_masks_sig[name] = sig_pass
        eff_cut_masks_bkg[name] = bkg_pass

        n_sig_pass = np.sum(sig_pass)
        n_bkg_pass = np.sum(bkg_pass)
        sig_eff = n_sig_pass / n_sig if n_sig > 0 else 0
        bkg_eff = n_bkg_pass / n_bkg if n_bkg > 0 else 0
        bkg_rej = 1 - bkg_eff

        print(f"{cut_str:<25} {sig_eff*100:>6.2f}%      {bkg_eff*100:>6.2f}%      {bkg_rej*100:>6.2f}%      {n_sig_pass:>10,}  {n_bkg_pass:>10,}")

    # Combined cut efficiency
    print("\n--- Combined Cut Efficiency (All Cuts Applied) ---")

    # Build combined mask
    combined_sig_mask = np.ones(n_sig, dtype=bool)
    combined_bkg_mask = np.ones(n_bkg, dtype=bool)

    for name in eff_cut_masks_sig:
        combined_sig_mask &= eff_cut_masks_sig[name]
        combined_bkg_mask &= eff_cut_masks_bkg[name]

    n_sig_combined = np.sum(combined_sig_mask)
    n_bkg_combined = np.sum(combined_bkg_mask)
    sig_eff_combined = n_sig_combined / n_sig if n_sig > 0 else 0
    bkg_eff_combined = n_bkg_combined / n_bkg if n_bkg > 0 else 0
    bkg_rej_combined = 1 - bkg_eff_combined

    print(f"\nCuts applied:")
    for var, name, title, nbins, xmin, xmax, xlabel, logy, cut_val, cut_dir in plots:
        if cut_val is None:
            continue
        if cut_dir == '>':
            print(f"  - {name} > {cut_val}")
        elif cut_dir == '<':
            print(f"  - {name} < {cut_val}")
        elif cut_dir == 'abs<':
            print(f"  - |{name}| < {cut_val}")
        elif cut_dir == 'abs>=':
            print(f"  - |{name}| >= {cut_val}")

    print(f"\n{'Metric':<25} {'Signal':<20} {'Background':<20}")
    print("-" * 65)
    print(f"{'Initial yield':<25} {n_sig:>15,}      {n_bkg:>15,}")
    print(f"{'Final yield':<25} {n_sig_combined:>15,}      {n_bkg_combined:>15,}")
    print(f"{'Efficiency':<25} {sig_eff_combined*100:>14.2f}%      {bkg_eff_combined*100:>14.2f}%")
    print(f"{'Rejection':<25} {(1-sig_eff_combined)*100:>14.2f}%      {bkg_rej_combined*100:>14.2f}%")

    # Signal purity after cuts
    total_after_cuts = n_sig_combined + n_bkg_combined
    purity = n_sig_combined / total_after_cuts if total_after_cuts > 0 else 0
    print(f"\n{'Signal purity after cuts':<25} {purity*100:>14.2f}%")
    print(f"{'S/sqrt(B)':<25} {n_sig_combined/np.sqrt(max(1,n_bkg_combined)):>14.2f}")

    print("\n" + "=" * 70)

    # Summary canvas with all distributions (N-1 selection applied)
    print("\nCreating summary canvas...")

    c_summary = ROOT.TCanvas('c_summary', 'Summary (N-1 plots)', 1200, 800)
    c_summary.Divide(3, 2)

    # Keep references to all objects so they don't get garbage collected
    summary_objects = []

    for i, (var, name, title, nbins, xmin, xmax, xlabel, logy, cut_val, cut_dir) in enumerate(plots):
        c_summary.cd(i + 1)
        if logy:
            ROOT.gPad.SetLogy()

        # Apply N-1 selection
        n_minus_1_sig = np.ones(n_sig, dtype=bool)
        n_minus_1_bkg = np.ones(n_bkg, dtype=bool)
        for other_name, masks in cut_masks.items():
            if other_name != name:
                n_minus_1_sig &= masks['signal']
                n_minus_1_bkg &= masks['background']

        sig_data = var[is_signal][n_minus_1_sig]
        bkg_data = var[is_bkg][n_minus_1_bkg]

        h_sig, h_bkg = create_comparison_histogram(
            sig_data, bkg_data, f'{name}_sum', title, nbins, xmin, xmax, xlabel
        )
        summary_objects.extend([h_sig, h_bkg])

        max_val = max(h_sig.GetMaximum(), h_bkg.GetMaximum())
        h_sig.SetMaximum(max_val * 1.5)
        h_sig.SetMinimum(0.0001 if logy else 0)
        h_sig.SetTitle(title)
        h_sig.Draw('HIST')
        h_bkg.Draw('HIST SAME')

        # Draw cut line
        if cut_val is not None:
            ymax = max_val * 1.3
            ymin = 0.0001 if logy else 0
            line = ROOT.TLine(cut_val, ymin, cut_val, ymax)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            summary_objects.append(line)
            if cut_dir in ['abs<', 'abs>=']:
                line2 = ROOT.TLine(-cut_val, ymin, -cut_val, ymax)
                line2.SetLineColor(ROOT.kBlack)
                line2.SetLineStyle(2)
                line2.SetLineWidth(2)
                line2.Draw()
                summary_objects.append(line2)

        leg = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(h_sig, 'Signal', 'f')
        leg.AddEntry(h_bkg, 'Background', 'f')
        leg.Draw()
        summary_objects.append(leg)

    c_summary.Write()

    output_file.Close()

    print(f"\nOutput saved to: {args.output}")
    print(f"Created {len(all_canvases) + 1} plots")


if __name__ == "__main__":
    main()
