#!/usr/bin/env python3
"""
hyddraEXOElectronEfficiency.py — Electron track collection efficiency comparison

Reads the 18 HyddraSVsEXOAnalyzer output files (2 mass groups × 3 ctau × 3
collections: GED, low-pT, merged) and writes a single ROOT file containing
efficiency canvases organized by comparison type.

Canvas organization:
  collection_comparison/{sample}/  — GED vs low-pT vs merged per sample
  ctau_scan/{collection}/{mass}/   — ctau 1/10/100 mm per collection and mass
  mass_comparison/{collection}/{ctau}/ — soft vs hard electrons per collection and ctau

Each canvas shows absolute reconstruction efficiency (isGold / all gen vertices).
For collection_comparison, solid markers = inclusive, open markers = isolated.
For ctau_scan and mass_comparison, inclusive only.

Usage:
    python scripts/hyddraEXOElectronEfficiency.py -i iDMroot/*.root -o eff.root
    python scripts/hyddraEXOElectronEfficiency.py --input-dir iDMroot/ -o eff.root
"""

import argparse
import glob
import os
import re
import numpy as np
import ROOT
import uproot
import cmsstyle as CMS

from efficiency_utils import calculate_efficiency, create_tgraph, draw_cms_labels
from hyddraEXOEfficiency import _make_canvas, _make_dummy

# ── Collection metadata ────────────────────────────────────────────────────────

COLLECTIONS = ['gedElectronTracks', 'lowPtElectronTracks', 'mergedElectronTracks']

COLL_SHORT = {
    'gedElectronTracks':    'ged',
    'lowPtElectronTracks':  'lowPt',
    'mergedElectronTracks': 'merged',
}
COLL_LABEL = {
    'gedElectronTracks':    'GED',
    'lowPtElectronTracks':  'Low-p_{T}',
    'mergedElectronTracks': 'GED+Low-p_{T}',
}
COLL_COLOR = {
    'gedElectronTracks':    ROOT.kAzure + 1,
    'lowPtElectronTracks':  ROOT.kOrange + 1,
    'mergedElectronTracks': ROOT.kTeal + 4,
}
COLL_SOLID = {'gedElectronTracks': 20, 'lowPtElectronTracks': 21, 'mergedElectronTracks': 22}
COLL_OPEN  = {'gedElectronTracks': 24, 'lowPtElectronTracks': 25, 'mergedElectronTracks': 26}

CTAU_COLORS  = {1: ROOT.kBlue + 1, 10: ROOT.kGreen + 2, 100: ROOT.kRed - 4}
CTAU_MARKERS = {1: 20,             10: 21,               100: 22}
MASS_COLORS  = [ROOT.kAzure + 1, ROOT.kOrange + 1]
MASS_MARKERS = [20, 21]

# ── Branches ──────────────────────────────────────────────────────────────────

BRANCHES = [
    'HyddraSV_genVtxIdx',
    'HyddraSV_isGold',
    'HyddraSV_dxy',
    'HyddraSV_passesIsolation',
    'HyddraGenSV_dxy',
    'HyddraGenSV_pt',
    'HyddraGenSV_trk1Pt',
    'HyddraGenSV_trk2Pt',
    'HyddraGenSV_passSelection',
    # Event-level branches for AN-style skim
    'Event_genMET',
    'nRecoElectrons',
]

# AN gen-level filter thresholds (Table 7 of AN2023-091)
AN_GEN_MET_CUT      = 80.   # GeV  (gen p_T^miss > 80 GeV gen filter)
AN_MIN_RECO_ELECTRONS = 2   # reco electrons required (AN footnote 20)

TREE_PATHS = ['Events', 'hyddraEXOAnalyzer/Events', 'tree']

# ── Axis definitions ──────────────────────────────────────────────────────────

X_BINS = {
    'gen_dxy':    np.logspace(-1,  2,    21),   # 0.1 – 100 cm
    'gen_pt':     np.logspace( 0,  np.log10(300.), 21),   # 1 – 300 GeV
    'gen_minpt':  np.logspace(-1,  np.log10(50.),  21),   # 0.1 – 50 GeV
}
X_LABEL = {
    'gen_dxy':    'd_{xy}^{gen} [cm]',
    'gen_pt':     'p_{T}^{gen} (di-e system) [GeV]',
    'gen_minpt':  'p_{T}^{gen} (softer electron) [GeV]',
}
X_NAMES = {
    'gen_dxy':    'eff_vs_dxy',
    'gen_pt':     'eff_vs_pt',
    'gen_minpt':  'eff_vs_mintrk_pt',
}

# ── Filename parsing ──────────────────────────────────────────────────────────

_FNAME_RE = re.compile(
    r'Mchi-([\dp]+)_dMchi-([\dp]+)_ctau-(\d+)'
    r'.*?_(gedElectronTracks|lowPtElectronTracks|mergedElectronTracks)(?:_[^.]+)?\.root$'
)


def parse_filename(path):
    m = _FNAME_RE.search(os.path.basename(path))
    if not m:
        return None
    mchi, dmchi, ctau_s, coll = m.group(1), m.group(2), m.group(3), m.group(4)
    ctau = int(ctau_s)
    mass_key  = f'Mchi{mchi}_dMchi{dmchi}'
    mchi_d    = mchi.replace('p', '.')
    dmchi_d   = dmchi.replace('p', '.')
    return {
        'path':         path,
        'collection':   coll,
        'mass_key':     mass_key,
        'ctau':         ctau,
        'mass_label':   f'M_{{#chi}}={mchi_d} GeV, #Deltam={dmchi_d} GeV',
        'ctau_label':   f'c#tau = {ctau} mm',
        'sample_label': (f'M_{{#chi}}={mchi_d} GeV, '
                         f'#Deltam={dmchi_d} GeV, c#tau={ctau} mm'),
    }

# ── Data loading ──────────────────────────────────────────────────────────────

def load_and_analyze(path, an_skim=False):
    """Return efficiency arrays for all gen vertices in the file.

    If *an_skim* is True, apply the two event-level skims used in the AN
    (AN2023-091, Sec. 7) before counting gen vertices in the denominator:
      1. Gen-level MET > 80 GeV  (mirrors the generator-level p_T^miss filter
         applied to the official signal samples, Table 7).
      2. nRecoElectrons >= 2  (mirrors the "two or more reco electrons"
         requirement described in footnote 20 of Sec. 7.2).
    These raise the apparent efficiency by removing events that are
    unlikely to be triggered on or to have any reconstruction activity.
    """
    f = uproot.open(path)
    tree = None
    for tp in TREE_PATHS:
        if tp in f:
            tree = f[tp]
            break
    if tree is None:
        raise KeyError(f"Cannot find tree in {path}. Keys: {list(f.keys())}")

    available = set(tree.keys())
    to_load   = [b for b in BRANCHES if b in available]
    missing   = [b for b in BRANCHES if b not in available]
    if missing:
        print(f"  WARNING: missing branches: {missing}")

    data     = tree.arrays(to_load, library='np')
    n_events = len(data['HyddraGenSV_dxy'])

    gen_dxy, gen_pt, gen_minpt = [], [], []
    gen_pass_sel = []
    gen_pass_gold_inc = []
    gen_pass_gold_iso = []

    for ev in range(n_events):
        # ── AN-style event skim ─────────────────────────────────────────────
        if an_skim:
            gen_met   = float(data['Event_genMET'][ev])   if 'Event_genMET'   in data else -1.
            n_reco_e  = int  (data['nRecoElectrons'][ev]) if 'nRecoElectrons' in data else -1
            if gen_met < AN_GEN_MET_CUT or n_reco_e < AN_MIN_RECO_ELECTRONS:
                continue

        sv_gvIdx  = data['HyddraSV_genVtxIdx'][ev]
        sv_isGold = data['HyddraSV_isGold'][ev]
        sv_iso    = data['HyddraSV_passesIsolation'][ev]

        gold_inc, gold_iso = set(), set()
        for ri in range(len(sv_gvIdx)):
            gi = int(sv_gvIdx[ri])
            if gi >= 0 and sv_isGold[ri]:
                gold_inc.add(gi)
                if sv_iso[ri]:
                    gold_iso.add(gi)

        gv_dxy  = data['HyddraGenSV_dxy'][ev]
        gv_pt   = data['HyddraGenSV_pt'][ev]
        gv_t1pt = data['HyddraGenSV_trk1Pt'][ev]
        gv_t2pt = data['HyddraGenSV_trk2Pt'][ev]
        gv_psel = data['HyddraGenSV_passSelection'][ev]

        for i in range(len(gv_dxy)):
            gen_dxy.append(float(gv_dxy[i]))
            gen_pt.append(float(gv_pt[i]))
            gen_minpt.append(min(float(gv_t1pt[i]), float(gv_t2pt[i])))
            gen_pass_sel.append(bool(gv_psel[i]))
            gen_pass_gold_inc.append(i in gold_inc)
            gen_pass_gold_iso.append(i in gold_iso)

    return {
        'gen_dxy':           np.array(gen_dxy),
        'gen_pt':            np.array(gen_pt),
        'gen_minpt':         np.array(gen_minpt),
        'gen_pass_sel':      np.array(gen_pass_sel),
        'gen_pass_gold_inc': np.array(gen_pass_gold_inc),
        'gen_pass_gold_iso': np.array(gen_pass_gold_iso),
        'n_events':          n_events,
    }

# ── Canvas helpers ────────────────────────────────────────────────────────────

# Keep all intermediate ROOT objects alive until the output file is closed.
_KEEP_ALIVE = []


def _draw_graph(result, x_key, pass_key, bins, gname, label, color, marker):
    x    = result[x_key]
    mask = result[pass_key]
    if len(x) == 0:
        return None
    c, e, err = calculate_efficiency(x, mask, bins)
    g = create_tgraph(c, e, err, gname, label)
    g.SetMarkerStyle(marker)
    g.SetMarkerSize(1.2)
    g.SetMarkerColor(color)
    g.SetLineColor(color)
    g.Draw('P SAME')
    _KEEP_ALIVE.append(g)
    return g


def _make_legend(graphs, x1=0.55, y1=0.20, header=None):
    lh  = 0.045 * max(len(graphs), 1)
    leg = ROOT.TLegend(x1, y1, 0.87, y1 + lh)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if header:
        leg.SetHeader(header)
    for g in graphs:
        leg.AddEntry(g, g.GetTitle(), 'p')
    leg.Draw()
    _KEEP_ALIVE.append(leg)
    return leg


def _new_canvas(name, x_key):
    bins = X_BINS[x_key]
    xl   = X_LABEL[x_key]
    yl   = 'Reconstruction Efficiency'
    xmin, xmax = float(bins[0]), float(bins[-1])
    ymin, ymax = 1e-4, 1.
    c = _make_canvas(name, xl, yl, xmin, xmax, y_min=ymin, y_max=ymax,
                     logx=True, logy=True)
    h = _make_dummy(name, xl, yl, xmin, xmax, ymin, ymax, len(bins) - 1, bins)
    _KEEP_ALIVE.append(h)
    return c

# ── Per-group canvas builders ─────────────────────────────────────────────────

def make_collection_comparison(results_by_coll, x_key, canvas_name, header=None):
    """One canvas: GED / low-pT / merged overlaid; solid=inclusive, open=isolated."""
    bins = X_BINS[x_key]
    c    = _new_canvas(canvas_name, x_key)
    graphs = []
    for coll in COLLECTIONS:
        if coll not in results_by_coll:
            continue
        res    = results_by_coll[coll]
        color  = COLL_COLOR[coll]
        label  = COLL_LABEL[coll]
        for pass_key, marker, suffix in [
            ('gen_pass_gold_inc', COLL_SOLID[coll], 'inc'),
            ('gen_pass_gold_iso', COLL_OPEN[coll],  'iso'),
        ]:
            gname = f'{canvas_name}_{COLL_SHORT[coll]}_{suffix}'
            g = _draw_graph(res, x_key, pass_key, bins, gname,
                            f'{label} ({suffix})', color, marker)
            if g:
                graphs.append(g)

    _make_legend(graphs, header=header)
    lbl = draw_cms_labels()
    _KEEP_ALIVE.append(lbl)
    return c


def make_ctau_scan(results_by_ctau, x_key, canvas_name, header=None):
    """One canvas: ctau = 1 / 10 / 100 mm, inclusive efficiency."""
    bins = X_BINS[x_key]
    c    = _new_canvas(canvas_name, x_key)
    graphs = []
    for ctau in sorted(results_by_ctau):
        res    = results_by_ctau[ctau]
        color  = CTAU_COLORS.get(ctau, ROOT.kBlack)
        marker = CTAU_MARKERS.get(ctau, 20)
        gname  = f'{canvas_name}_ctau{ctau}'
        g = _draw_graph(res, x_key, 'gen_pass_gold_inc', bins, gname,
                        f'c#tau = {ctau} mm', color, marker)
        if g:
            graphs.append(g)

    _make_legend(graphs, header=header)
    lbl = draw_cms_labels()
    _KEEP_ALIVE.append(lbl)
    return c


def make_mass_comparison(results_by_mass, mass_keys_ordered, mass_labels,
                         x_key, canvas_name, header=None):
    """One canvas: soft vs hard mass point, inclusive efficiency."""
    bins = X_BINS[x_key]
    c    = _new_canvas(canvas_name, x_key)
    graphs = []
    for idx, mk in enumerate(mass_keys_ordered):
        if mk not in results_by_mass:
            continue
        res    = results_by_mass[mk]
        color  = MASS_COLORS[idx % len(MASS_COLORS)]
        marker = MASS_MARKERS[idx % len(MASS_MARKERS)]
        gname  = f'{canvas_name}_{mk}'
        g = _draw_graph(res, x_key, 'gen_pass_gold_inc', bins, gname,
                        mass_labels[mk], color, marker)
        if g:
            graphs.append(g)

    _make_legend(graphs, header=header)
    lbl = draw_cms_labels()
    _KEEP_ALIVE.append(lbl)
    return c


def save_canvas(tdir, canvas):
    tdir.cd()
    canvas.Write()
    ROOT.gROOT.cd()

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('-i', '--input', nargs='+',
                     help='Input ROOT files (explicit list)')
    grp.add_argument('--input-dir', metavar='DIR',
                     help='Directory to search for matching *.root files')
    parser.add_argument('-o', '--output',
                        default='hyddraEXO_electron_efficiency.root',
                        help='Output ROOT file')
    parser.add_argument('--an-skim', action='store_true',
                        help='Apply AN-style event skims to the efficiency '
                             'denominator: gen-MET > 80 GeV and '
                             'nRecoElectrons >= 2 (see AN2023-091 Sec. 7 '
                             'footnote 20 and Table 7). Requires files '
                             'produced with the updated analyzer that stores '
                             'Event_genMET and nRecoElectrons branches.')
    args = parser.parse_args()

    paths = (args.input if args.input
             else sorted(glob.glob(os.path.join(args.input_dir, '*.root'))))

    file_meta = {}
    for p in paths:
        meta = parse_filename(p)
        if meta is None:
            print(f"  WARNING: skipping unrecognized filename: {os.path.basename(p)}")
            continue
        file_meta[p] = meta

    if not file_meta:
        print("No recognized input files found.")
        return

    # ── Load all data ────────────────────────────────────────────────────────
    # Nested dict: data[mass_key][ctau][collection] = result
    print(f"\nLoading {len(file_meta)} files...")
    data = {}
    mass_labels = {}
    for path, meta in sorted(file_meta.items()):
        mk, ct, coll = meta['mass_key'], meta['ctau'], meta['collection']
        print(f"  {mk}  ctau={ct:3d}  {COLL_SHORT[coll]:7s}  {os.path.basename(path)}")
        result = load_and_analyze(path, an_skim=args.an_skim)
        data.setdefault(mk, {}).setdefault(ct, {})[coll] = result
        mass_labels[mk] = meta['mass_label']

    mass_keys = sorted(data.keys())
    ctau_vals = sorted({ct for mk in data for ct in data[mk]})

    CMS.SetExtraText("Simulation Preliminary")
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    out = ROOT.TFile(args.output, 'RECREATE')

    # ════════════════════════════════════════════════════════════════════════
    # 1. collection_comparison: per sample, 3 collections overlaid
    # ════════════════════════════════════════════════════════════════════════
    d_coll = out.mkdir('collection_comparison')

    for mk in mass_keys:
        for ct in sorted(data[mk]):
            sample_key = f'{mk}_ctau{ct}'
            sample_label = next(
                meta['sample_label']
                for p, meta in file_meta.items()
                if meta['mass_key'] == mk and meta['ctau'] == ct
            )
            d_samp = d_coll.mkdir(sample_key)
            results_by_coll = data[mk][ct]

            for xk in ['gen_dxy', 'gen_pt', 'gen_minpt']:
                cname = f'coll_{sample_key}_{X_NAMES[xk]}'
                c = make_collection_comparison(results_by_coll, xk, cname,
                                               header=sample_label)
                save_canvas(d_samp, c)

    # ════════════════════════════════════════════════════════════════════════
    # 2. ctau_scan: per collection and mass, 3 ctau curves
    # ════════════════════════════════════════════════════════════════════════
    d_ctau_dir = out.mkdir('ctau_scan')

    for coll in COLLECTIONS:
        cs = COLL_SHORT[coll]
        d_c = d_ctau_dir.mkdir(cs)

        for mk in mass_keys:
            d_m = d_c.mkdir(mk)
            results_by_ctau = {
                ct: data[mk][ct][coll]
                for ct in data[mk]
                if coll in data[mk][ct]
            }
            if not results_by_ctau:
                continue
            header = f'{COLL_LABEL[coll]}, {mass_labels[mk]}'

            for xk in ['gen_dxy', 'gen_pt', 'gen_minpt']:
                cname = f'ctau_{cs}_{mk}_{X_NAMES[xk]}'
                c = make_ctau_scan(results_by_ctau, xk, cname, header=header)
                save_canvas(d_m, c)

    # ════════════════════════════════════════════════════════════════════════
    # 3. mass_comparison: per collection and ctau, 2 mass curves
    # ════════════════════════════════════════════════════════════════════════
    d_mass_dir = out.mkdir('mass_comparison')

    for coll in COLLECTIONS:
        cs = COLL_SHORT[coll]
        d_c = d_mass_dir.mkdir(cs)

        for ct in ctau_vals:
            d_ct = d_c.mkdir(f'ctau{ct}')
            results_by_mass = {
                mk: data[mk][ct][coll]
                for mk in mass_keys
                if ct in data[mk] and coll in data[mk][ct]
            }
            if not results_by_mass:
                continue
            header = f'{COLL_LABEL[coll]}, c#tau = {ct} mm'

            for xk in ['gen_dxy', 'gen_pt', 'gen_minpt']:
                cname = f'mass_{cs}_ctau{ct}_{X_NAMES[xk]}'
                c = make_mass_comparison(results_by_mass, mass_keys, mass_labels,
                                         xk, cname, header=header)
                save_canvas(d_ct, c)

    out.Close()
    _KEEP_ALIVE.clear()
    print(f"\nCanvases saved to {args.output}")


if __name__ == '__main__':
    main()
