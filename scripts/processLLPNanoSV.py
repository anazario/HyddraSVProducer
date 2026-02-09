#!/usr/bin/env python3
"""
Process LLPNanoAOD files to produce HyddraSV-compatible ntuples.

Reads NanoAOD flat ROOT files directly (no CMSSW required) and writes
output trees compatible with muonSVEfficiencyComparison.py and
svEfficiencyAnalysis.py.

Supports PatMuonVertex and PatDSAMuonVertex collections.
Uses multiprocessing for parallel file processing.

Usage:
    python processLLPNanoSV.py -i file1.root file2.root --collections PatMuonVertex
    python processLLPNanoSV.py -i file1.root file2.root --collections PatMuonVertex PatDSAMuonVertex -j 8
    python processLLPNanoSV.py --input-file-list files.txt --collections PatDSAMuonVertex -o output.root
"""

import argparse
import os
import sys
import time
import numpy as np
import uproot
import awkward as ak
from multiprocessing import Pool

MUON_MASS = 0.10566  # GeV


# ============================================================================
# Gen vertex reconstruction
# ============================================================================

def build_gen_vertices(evt, mother_pdg_id):
    """Build gen muon vertices from GenPart arrays for one event.

    Finds status-1 muons whose mother chain leads to the signal particle,
    groups them by signal mother, and computes vertex kinematics.
    """
    pdg_id = evt['GenPart_pdgId']
    status = evt['GenPart_status']
    mother_idx = evt['GenPart_genPartIdxMother']
    n = len(pdg_id)

    mother_to_muons = {}
    for i in range(n):
        if abs(int(pdg_id[i])) != 13 or int(status[i]) != 1:
            continue
        mid = int(mother_idx[i])
        if mid < 0 or mid >= n:
            continue
        # Walk up chain to find signal mother
        current = mid
        found = False
        for _ in range(10):
            if current < 0 or current >= n:
                break
            if int(pdg_id[current]) == mother_pdg_id:
                found = True
                break
            if abs(int(pdg_id[current])) == 13:
                current = int(mother_idx[current])
                continue
            break
        if found:
            mother_to_muons.setdefault(current, []).append(i)

    vertices = []
    for mother, muon_indices in mother_to_muons.items():
        if len(muon_indices) < 2:
            continue
        first = muon_indices[0]
        gv_x = float(evt['GenPart_vx'][first])
        gv_y = float(evt['GenPart_vy'][first])
        gv_z = float(evt['GenPart_vz'][first])

        # 4-vector sum
        px_tot = py_tot = pz_tot = e_tot = 0.
        for idx in muon_indices:
            pt_i = float(evt['GenPart_pt'][idx])
            eta_i = float(evt['GenPart_eta'][idx])
            phi_i = float(evt['GenPart_phi'][idx])
            m_i = float(evt['GenPart_mass'][idx])
            px_tot += pt_i * np.cos(phi_i)
            py_tot += pt_i * np.sin(phi_i)
            pz_tot += pt_i * np.sinh(eta_i)
            e_tot += np.sqrt((pt_i * np.cosh(eta_i))**2 + m_i**2)

        p_tot = np.sqrt(px_tot**2 + py_tot**2 + pz_tot**2)
        vertices.append({
            'x': gv_x, 'y': gv_y, 'z': gv_z,
            'dxy': np.sqrt(gv_x**2 + gv_y**2),
            'pt': np.sqrt(px_tot**2 + py_tot**2),
            'eta': np.arctanh(pz_tot / p_tot) if p_tot > 0 else 0.,
            'phi': np.arctan2(py_tot, px_tot),
            'mass': np.sqrt(max(0., e_tot**2 - p_tot**2)),
            'muon_indices': muon_indices,
        })
    return vertices


# ============================================================================
# DSA-to-gen matching
# ============================================================================

def match_dsa_to_gen(dsa_eta, dsa_phi, dsa_pt, evt, gen_vertices, dr_cut, relpt_cut):
    """Find best deltaR match between a DSA muon and signal gen muons.

    Returns (gen_vertex_index, deltaR, relPtDiff). (-1, 999, 999) if no match.
    """
    gen_eta = evt['GenPart_eta']
    gen_phi = evt['GenPart_phi']
    gen_pt = evt['GenPart_pt']

    best_gv = -1
    best_dr = dr_cut
    best_relpt = 999.
    for gv_idx, gv in enumerate(gen_vertices):
        for mu_idx in gv['muon_indices']:
            deta = dsa_eta - float(gen_eta[mu_idx])
            dphi = float(dsa_phi) - float(gen_phi[mu_idx])
            dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
            dr = np.sqrt(deta**2 + dphi**2)
            relpt = abs(float(dsa_pt) - float(gen_pt[mu_idx])) / max(float(gen_pt[mu_idx]), 0.001)
            if relpt > relpt_cut:
                continue
            if dr < best_dr:
                best_dr = dr
                best_gv = gv_idx
                best_relpt = relpt
    return best_gv, best_dr, best_relpt


# ============================================================================
# Gold matching
# ============================================================================

def find_gold_match(iv, evt, vtx_pre, gen_vertices, dr_cut, relpt_cut):
    """Check if both vertex legs match to the same gen vertex. Returns index or -1."""
    if not gen_vertices:
        return -1

    def trace_leg(orig_idx_f, is_dsa_f):
        is_dsa = float(is_dsa_f) > 0.5
        orig_idx = int(orig_idx_f)
        if is_dsa:
            if orig_idx < 0 or orig_idx >= len(evt['DSAMuon_pt']):
                return -1
            gv_idx, _, _ = match_dsa_to_gen(
                float(evt['DSAMuon_eta'][orig_idx]),
                float(evt['DSAMuon_phi'][orig_idx]),
                float(evt['DSAMuon_pt'][orig_idx]),
                evt, gen_vertices, dr_cut, relpt_cut)
            return gv_idx
        else:
            if orig_idx < 0 or orig_idx >= len(evt['Muon_genPartIdx']):
                return -1
            gp_idx = int(evt['Muon_genPartIdx'][orig_idx])
            if gp_idx < 0 or gp_idx >= len(evt['GenPart_pdgId']):
                return -1
            if abs(int(evt['GenPart_pdgId'][gp_idx])) != 13:
                return -1
            for gv_idx, gv in enumerate(gen_vertices):
                if gp_idx in gv['muon_indices']:
                    return gv_idx
            return -1

    vtx = {k.replace(vtx_pre, ''): evt[k] for k in evt if k.startswith(vtx_pre)}
    gv1 = trace_leg(vtx['originalMuonIdx1'][iv], vtx['isDSAMuon1'][iv])
    gv2 = trace_leg(vtx['originalMuonIdx2'][iv], vtx['isDSAMuon2'][iv])
    if gv1 >= 0 and gv1 == gv2:
        return gv1
    return -1


# ============================================================================
# Per-file processing (runs in worker process)
# ============================================================================

def process_file(args):
    """Process one NanoAOD file for one collection. Returns output dict."""
    filename, collection, mother_pdg_id, dr_cut, relpt_cut = args

    vtx_pre = collection + '_'
    ref_pre = collection + 'RefittedTracks_'

    try:
        f = uproot.open(filename)
    except Exception as e:
        print(f'  WARNING: cannot open {filename}: {e}', file=sys.stderr)
        return None

    if 'Events' not in f:
        print(f'  WARNING: no Events tree in {filename}', file=sys.stderr)
        return None

    tree = f['Events']
    available = set(tree.keys())

    # Branches to read
    need = [
        'GenPart_pdgId', 'GenPart_genPartIdxMother', 'GenPart_status',
        'GenPart_vx', 'GenPart_vy', 'GenPart_vz',
        'GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'GenPart_mass',
        'Muon_genPartIdx', 'Muon_pt', 'Muon_eta', 'Muon_phi',
        'DSAMuon_pt', 'DSAMuon_eta', 'DSAMuon_phi',
    ]
    vtx_cols = ['isValid', 'vx', 'vy', 'vz', 'vxy',
                'chi2', 'normChi2', 'ndof',
                'originalMuonIdx1', 'originalMuonIdx2',
                'isDSAMuon1', 'isDSAMuon2']
    ref_cols = ['px', 'py', 'pz', 'idx']
    need += [vtx_pre + c for c in vtx_cols]
    need += [ref_pre + c for c in ref_cols]

    to_read = [b for b in need if b in available]
    missing = [b for b in need if b not in available]
    if missing:
        print(f'  WARNING [{os.path.basename(filename)}]: missing branches: {missing}',
              file=sys.stderr)

    data = tree.arrays(to_read, library='np')
    n_events = len(data['GenPart_pdgId'])

    # Output accumulators
    scalar_keys = ['nVertices', 'genVertex_nTotal', 'genVertex_nMuon',
                   'genVertex_nElectron', 'genVertex_nHadronic']
    vec_float_keys = ['x', 'y', 'z', 'dxy', 'mass', 'pt', 'eta', 'phi',
                      'chi2', 'normalizedChi2', 'ndof', 'min3D']
    vec_int_keys = ['genVertexIndex', 'nearestGenVertexIndex']
    vec_bool_keys = ['isLeptonic', 'isBronze', 'isSilver', 'isGold']
    gen_float_keys = ['dxy', 'x', 'y', 'z', 'pt', 'eta', 'phi', 'mass']
    gen_bool_keys = ['isMuon', 'isElectron', 'isHadronic', 'passSelection']

    out = {k: [] for k in scalar_keys}
    out['nTracks'] = []
    for k in vec_float_keys + vec_int_keys + vec_bool_keys:
        out['sv_' + k] = []
    for k in gen_float_keys + gen_bool_keys:
        out['gv_' + k] = []

    dsa_dr_vals = []
    dsa_relpt_vals = []

    for ei in range(n_events):
        # Per-event branch data
        evt = {b: data[b][ei] for b in to_read}

        gen_vertices = build_gen_vertices(evt, mother_pdg_id)

        # ---- Gen vertex output ----
        n_gv = len(gen_vertices)
        out['genVertex_nTotal'].append(n_gv)
        out['genVertex_nMuon'].append(n_gv)
        out['genVertex_nElectron'].append(0)
        out['genVertex_nHadronic'].append(0)

        for k in gen_float_keys:
            out['gv_' + k].append(np.array([gv[k] for gv in gen_vertices], dtype=np.float32))
        out['gv_isMuon'].append(np.ones(n_gv, dtype=np.bool_))
        out['gv_isElectron'].append(np.zeros(n_gv, dtype=np.bool_))
        out['gv_isHadronic'].append(np.zeros(n_gv, dtype=np.bool_))
        out['gv_passSelection'].append(np.ones(n_gv, dtype=np.bool_))

        # ---- Reco vertices ----
        is_valid_key = vtx_pre + 'isValid'
        n_vtx = len(evt.get(is_valid_key, []))

        ev = {k: [] for k in vec_float_keys + vec_int_keys + vec_bool_keys + ['nTracks']}

        n_valid = 0
        for iv in range(n_vtx):
            if float(evt[is_valid_key][iv]) < 0.5:
                continue
            n_valid += 1

            ev['isLeptonic'].append(True)
            ev['nTracks'].append(2)
            ev['x'].append(float(evt[vtx_pre + 'vx'][iv]))
            ev['y'].append(float(evt[vtx_pre + 'vy'][iv]))
            ev['z'].append(float(evt[vtx_pre + 'vz'][iv]))
            ev['dxy'].append(float(evt[vtx_pre + 'vxy'][iv]))
            ev['chi2'].append(float(evt[vtx_pre + 'chi2'][iv]))
            ev['normalizedChi2'].append(float(evt[vtx_pre + 'normChi2'][iv]))
            ev['ndof'].append(float(evt[vtx_pre + 'ndof'][iv]))

            # Kinematics from refitted tracks
            idx_key = ref_pre + 'idx'
            px_tot = py_tot = pz_tot = e_tot = 0.
            found = False
            if idx_key in evt:
                for it in range(len(evt[idx_key])):
                    if int(evt[idx_key][it]) == iv:
                        px = float(evt[ref_pre + 'px'][it])
                        py = float(evt[ref_pre + 'py'][it])
                        pz = float(evt[ref_pre + 'pz'][it])
                        e_tot += np.sqrt(px**2 + py**2 + pz**2 + MUON_MASS**2)
                        px_tot += px; py_tot += py; pz_tot += pz
                        found = True

            if found:
                p_tot = np.sqrt(px_tot**2 + py_tot**2 + pz_tot**2)
                ev['mass'].append(np.sqrt(max(0., e_tot**2 - p_tot**2)))
                ev['pt'].append(np.sqrt(px_tot**2 + py_tot**2))
                ev['eta'].append(np.arctanh(pz_tot / p_tot) if p_tot > 0 else 0.)
                ev['phi'].append(np.arctan2(py_tot, px_tot))
            else:
                ev['mass'].append(0.); ev['pt'].append(0.)
                ev['eta'].append(0.); ev['phi'].append(0.)

            # Gold matching
            gold_idx = find_gold_match(iv, evt, vtx_pre, gen_vertices, dr_cut, relpt_cut)
            is_gold = gold_idx >= 0
            ev['genVertexIndex'].append(gold_idx)
            ev['isGold'].append(is_gold)
            ev['isSilver'].append(is_gold)
            ev['isBronze'].append(is_gold)

            # DSA histogram (no cuts)
            for oidx_key, dsa_key in [('originalMuonIdx1', 'isDSAMuon1'),
                                       ('originalMuonIdx2', 'isDSAMuon2')]:
                if float(evt[vtx_pre + dsa_key][iv]) > 0.5:
                    leg = int(evt[vtx_pre + oidx_key][iv])
                    if 0 <= leg < len(evt.get('DSAMuon_pt', [])):
                        _, dr, rpt = match_dsa_to_gen(
                            float(evt['DSAMuon_eta'][leg]),
                            float(evt['DSAMuon_phi'][leg]),
                            float(evt['DSAMuon_pt'][leg]),
                            evt, gen_vertices, 999., 999.)
                        if dr < 900.:
                            dsa_dr_vals.append(dr)
                            dsa_relpt_vals.append(rpt)

            # Spatial matching
            min_dist = -1.
            nearest_idx = -1
            vx_r = float(evt[vtx_pre + 'vx'][iv])
            vy_r = float(evt[vtx_pre + 'vy'][iv])
            vz_r = float(evt[vtx_pre + 'vz'][iv])
            for gv_idx, gv in enumerate(gen_vertices):
                d = np.sqrt((vx_r - gv['x'])**2 + (vy_r - gv['y'])**2 + (vz_r - gv['z'])**2)
                if nearest_idx < 0 or d < min_dist:
                    min_dist = d
                    nearest_idx = gv_idx
            ev['nearestGenVertexIndex'].append(nearest_idx)
            ev['min3D'].append(min_dist)

        out['nVertices'].append(n_valid)
        out['nTracks'].append(np.array(ev['nTracks'], dtype=np.uint32))
        for k in vec_float_keys:
            out['sv_' + k].append(np.array(ev[k], dtype=np.float32))
        for k in vec_int_keys:
            out['sv_' + k].append(np.array(ev[k], dtype=np.int32))
        for k in vec_bool_keys:
            out['sv_' + k].append(np.array(ev[k], dtype=np.bool_))

    return out, dsa_dr_vals, dsa_relpt_vals, n_events


# ============================================================================
# Output writing
# ============================================================================

def write_output(output_file, all_results, dsa_dr, dsa_relpt):
    """Concatenate results from all files and write output ROOT file."""

    # Concatenate
    combined = None
    for result in all_results:
        if result is None:
            continue
        out, dr_vals, relpt_vals, _ = result
        dsa_dr.extend(dr_vals)
        dsa_relpt.extend(relpt_vals)
        if combined is None:
            combined = out
        else:
            for key in combined:
                combined[key].extend(out[key])

    if combined is None:
        print('ERROR: No data to write.', file=sys.stderr)
        return

    # Build tree dict with HyddraSV/HyddraGenVertex branch names
    tree_data = {}

    # Scalar branches
    tree_data['HyddraSV_nVertices'] = np.array(combined['nVertices'], dtype=np.uint32)
    tree_data['HyddraGenVertex_nTotal'] = np.array(combined['genVertex_nTotal'], dtype=np.uint32)
    tree_data['HyddraGenVertex_nMuon'] = np.array(combined['genVertex_nMuon'], dtype=np.uint32)
    tree_data['HyddraGenVertex_nElectron'] = np.array(combined['genVertex_nElectron'], dtype=np.uint32)
    tree_data['HyddraGenVertex_nHadronic'] = np.array(combined['genVertex_nHadronic'], dtype=np.uint32)

    # Jagged reco branches
    sv_renames = {
        'x': 'x', 'y': 'y', 'z': 'z', 'dxy': 'dxy',
        'mass': 'mass', 'pt': 'pt', 'eta': 'eta', 'phi': 'phi',
        'chi2': 'chi2', 'normalizedChi2': 'normalizedChi2', 'ndof': 'ndof',
        'min3D': 'min3D', 'genVertexIndex': 'genVertexIndex',
        'nearestGenVertexIndex': 'nearestGenVertexIndex',
        'isLeptonic': 'isLeptonic', 'isBronze': 'isBronze',
        'isSilver': 'isSilver', 'isGold': 'isGold',
    }
    for internal, branch in sv_renames.items():
        tree_data['HyddraSV_' + branch] = ak.Array(combined['sv_' + internal])
    tree_data['HyddraSV_nTracks'] = ak.Array(combined['nTracks'])

    # Jagged gen branches
    gv_float = ['dxy', 'x', 'y', 'z', 'pt', 'eta', 'phi', 'mass']
    gv_bool = ['isMuon', 'isElectron', 'isHadronic', 'passSelection']
    for k in gv_float + gv_bool:
        tree_data['HyddraGenVertex_' + k] = ak.Array(combined['gv_' + k])

    with uproot.recreate(output_file) as f:
        f['llpNanoSVAnalyzer/tree'] = tree_data

        # DSA matching histograms
        if dsa_dr:
            f['h_dsaGenDeltaR'] = np.histogram(
                dsa_dr, bins=150, range=(0., 3.))
            f['h_dsaGenRelPtDiff'] = np.histogram(
                dsa_relpt, bins=100, range=(0., 5.))

    n_events = len(combined['nVertices'])
    print(f'  Wrote {n_events} events to {output_file}')


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Process LLPNanoAOD files to HyddraSV-compatible ntuples')
    parser.add_argument('-i', '--input', nargs='+', default=[],
                        help='Input NanoAOD ROOT files')
    parser.add_argument('--input-file-list', default='',
                        help='Text file with input file paths (one per line)')
    parser.add_argument('--collections', nargs='+',
                        default=['PatMuonVertex'],
                        help='Collections to process (default: PatMuonVertex)')
    parser.add_argument('-o', '--output', default='llpNanoSV',
                        help='Output file basename (default: llpNanoSV)')
    parser.add_argument('-j', '--workers', type=int, default=4,
                        help='Number of parallel workers (default: 4)')
    parser.add_argument('--mother-pdg-id', type=int, default=54,
                        help='Signal mother PDG ID (default: 54)')
    parser.add_argument('--dsa-delta-r', type=float, default=0.3,
                        help='DSA-gen deltaR matching threshold (default: 0.3)')
    parser.add_argument('--dsa-rel-pt-diff', type=float, default=0.5,
                        help='DSA-gen relative pT diff threshold (default: 0.5)')
    args = parser.parse_args()

    # Build input file list
    input_files = list(args.input)
    if args.input_file_list:
        if not os.path.exists(args.input_file_list):
            print(f'ERROR: file list not found: {args.input_file_list}', file=sys.stderr)
            sys.exit(1)
        with open(args.input_file_list) as flist:
            for line in flist:
                line = line.strip()
                if line and not line.startswith('#'):
                    input_files.append(line)

    if not input_files:
        parser.error('No input files specified. Use -i or --input-file-list.')

    print(f'Input files: {len(input_files)}')
    print(f'Collections: {", ".join(args.collections)}')
    print(f'Workers:     {args.workers}')
    print()

    for collection in args.collections:
        if len(args.collections) > 1:
            out_file = f'{args.output}_{collection}.root'
        else:
            out_file = args.output if args.output.endswith('.root') else args.output + '.root'

        print(f'Processing {collection} ...')
        t0 = time.time()

        worker_args = [(fn, collection, args.mother_pdg_id,
                        args.dsa_delta_r, args.dsa_rel_pt_diff)
                       for fn in input_files]

        dsa_dr = []
        dsa_relpt = []

        if args.workers > 1 and len(input_files) > 1:
            with Pool(args.workers) as pool:
                results = []
                for i, result in enumerate(pool.imap_unordered(process_file, worker_args)):
                    if result is not None:
                        results.append(result)
                    print(f'\r  Files processed: {i + 1}/{len(input_files)}', end='', flush=True)
                print()
        else:
            results = []
            for i, wa in enumerate(worker_args):
                result = process_file(wa)
                if result is not None:
                    results.append(result)
                print(f'\r  Files processed: {i + 1}/{len(input_files)}', end='', flush=True)
            print()

        total_events = sum(r[3] for r in results)
        elapsed = time.time() - t0
        print(f'  Total events: {total_events}  ({elapsed:.1f}s)')

        write_output(out_file, results, dsa_dr, dsa_relpt)
        print()


if __name__ == '__main__':
    main()
