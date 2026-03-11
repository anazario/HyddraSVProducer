"""
loader.py — uproot-based data loading utilities for hyddraDiagPlots.
"""
import awkward as ak
from .config import STAGE_KEYS

_BASE          = "hyddraSVsDiagAnalyzer"
_HADRONIC_BASE = "hyddraSVsHadronicDiagAnalyzer"


def load_gen_funnel(root_file, base=_BASE):
    """Return genFunnel tree arrays (vectors per event, per gen vertex)."""
    return root_file[f"{base}/genFunnel"].arrays(library="ak")


def load_stage_counts(root_file, base=_BASE):
    """Return stageCounts summed over all events as a flat dict."""
    arr = root_file[f"{base}/stageCounts"].arrays(library="ak")
    totals = {}
    for key in ["n", "nGold", "nSilver", "nBronze"]:
        for stage in STAGE_KEYS:
            branch = f"Stage_{key}_{stage}"
            totals[branch] = int(ak.sum(arr[branch]))
    return totals


def load_cleaning_tracks(root_file, base=_BASE):
    """Return cleaningTracks tree arrays (vectors per event)."""
    return root_file[f"{base}/cleaningTracks"].arrays(library="ak")


def load_all_stage_vtx(root_file, base=_BASE):
    """Return allStageVtx tree arrays (one row per vertex per stage, flat)."""
    return root_file[f"{base}/allStageVtx"].arrays(library="ak")


def load_seed_tracks(root_file, base=_BASE):
    """Return seedTracks tree arrays (vectors per event, tracks in 2-track OS seeds)."""
    return root_file[f"{base}/seedTracks"].arrays(library="ak")


def load_leptonic_config(root_file):
    """Return the stored leptonic HYDDRA config as a flat dict, or None if absent."""
    try:
        arr = root_file[f"{_BASE}/leptonicConfig"].arrays(library="np")
        return {k: v[0] for k, v in arr.items()}
    except KeyError:
        return None


def load_hadronic_config(root_file):
    """Return the stored hadronic HYDDRA config as a flat dict, or None if absent."""
    try:
        arr = root_file[f"{_HADRONIC_BASE}/hadronicConfig"].arrays(library="np")
        return {k: v[0] for k, v in arr.items()}
    except KeyError:
        return None


def has_leptonic_data(root_file):
    """Return True if the file contains leptonic diagnostic trees."""
    try:
        root_file[f"{_BASE}/leptonicConfig"]
        return True
    except KeyError:
        return False


def has_hadronic_data(root_file):
    """Return True if the file contains hadronic diagnostic trees."""
    try:
        root_file[f"{_HADRONIC_BASE}/hadronicConfig"]
        return True
    except KeyError:
        return False
