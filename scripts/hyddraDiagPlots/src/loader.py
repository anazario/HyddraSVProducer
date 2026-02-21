"""
loader.py — uproot-based data loading utilities for hyddraDiagPlots.
"""
import awkward as ak
from .config import STAGE_KEYS

_BASE = "hyddraSVsDiagAnalyzer"


def load_gen_funnel(root_file):
    """Return genFunnel tree arrays (vectors per event, per gen vertex)."""
    return root_file[f"{_BASE}/genFunnel"].arrays(library="ak")


def load_stage_counts(root_file):
    """Return stageCounts summed over all events as a flat dict."""
    arr = root_file[f"{_BASE}/stageCounts"].arrays(library="ak")
    totals = {}
    for key in ["n", "nGold", "nSilver", "nBronze"]:
        for stage in STAGE_KEYS:
            branch = f"Stage_{key}_{stage}"
            totals[branch] = int(ak.sum(arr[branch]))
    return totals


def load_cleaning_tracks(root_file):
    """Return cleaningTracks tree arrays (vectors per event)."""
    return root_file[f"{_BASE}/cleaningTracks"].arrays(library="ak")


def load_all_stage_vtx(root_file):
    """Return allStageVtx tree arrays (one row per vertex per stage, flat)."""
    return root_file[f"{_BASE}/allStageVtx"].arrays(library="ak")


def load_seed_tracks(root_file):
    """Return seedTracks tree arrays (vectors per event, tracks in 2-track OS seeds)."""
    return root_file[f"{_BASE}/seedTracks"].arrays(library="ak")
