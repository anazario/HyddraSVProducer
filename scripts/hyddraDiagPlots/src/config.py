"""
config.py — shared constants for the hyddraDiagPlots package.
"""
import numpy as np
import ROOT

# ── Hadronic signal definition ────────────────────────────────────────────────
# An allStageVtx row is considered signal if matchRatio > 0 AND the reco vertex
# is within HAD_MIN3D_CUT cm of the nearest hadronic gen vertex in 3D.
# Without the position cut, momentum-only track matching lets completely unrelated
# SVs accumulate matchRatio > 0 just by having a track pointing the same way.
HAD_MIN3D_CUT = 0.05  # cm

STAGE_NAMES = ["Seeding", "Merging", "Cleaning", "Disambiguation", "Filtering", "ID"]
STAGE_KEYS  = ["seed", "merged", "cleaned", "disambig", "filtered", "id"]
STAGE_DIRS  = ["seeding", "merging", "cleaning", "disambiguation", "filtering", "id", "summary"]
STAGE_IDX   = {k: i for i, k in enumerate(STAGE_KEYS)}

# ── Colours ──────────────────────────────────────────────────────────────────
COLOR_GOLD      = ROOT.kOrange + 1   # gold-ish
COLOR_SILVER    = ROOT.kAzure  + 6   # light blue
COLOR_BRONZE    = ROOT.kRed    - 3
COLOR_NONSIGNAL = ROOT.kGray   + 2   # non-signal from signal file
COLOR_BKG       = ROOT.kRed    + 2   # dedicated background file
COLORS_STAGE    = [ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kRed+2, ROOT.kOrange-3, ROOT.kMagenta-2, ROOT.kCyan+2]
MARKERS         = [20, 21, 22, 23, 29, 33]

# ── Reco observables shown per stage ─────────────────────────────────────────
RECO_OBSERVABLES = {
    "cosTheta":   {"label": "cos#theta (wrt PV)",       "bins": list(np.linspace(-1,   1,  51)), "log_y": True},
    "decayAngle": {"label": "cos#theta* (decay angle)", "bins": list(np.linspace(-1,   1,  51)), "log_y": True},
    "pOverE":     {"label": "p/E",                      "bins": list(np.linspace( 0,   1,  51)), "log_y": True},
    "dxySignif":  {"label": "dxy Significance",         "bins": list(np.linspace( 0, 150,  76)), "log_y": True},
    "mass":       {"label": "Invariant mass (GeV)",     "bins": list(np.linspace( 0, 100,  51)), "log_y": True},
}

# ── Hadronic-specific reco observables (adds nTracks for post-seeding stages) ─
HADRONIC_RECO_OBSERVABLES = {
    "cosTheta":   {"label": "cos#theta (wrt PV)",       "bins": list(np.linspace(-1,   1,  51)), "log_y": True},
    "decayAngle": {"label": "cos#theta* (decay angle)", "bins": list(np.linspace(-1,   1,  51)), "log_y": True},
    "pOverE":     {"label": "p/E",                      "bins": list(np.linspace( 0,   1,  51)), "log_y": True},
    "dxySignif":  {"label": "dxy Significance",         "bins": list(np.linspace( 0, 150,  76)), "log_y": True},
    "mass":       {"label": "Invariant mass (GeV)",     "bins": list(np.linspace( 0, 100,  51)), "log_y": True},
    "nTracks":    {"label": "Number of tracks",         "bins": list(range(2, 33)),               "log_y": True},
}

# ── Gen-level binning for efficiency plots ────────────────────────────────────
GEN_DXY_BINS = list(np.concatenate([
    np.linspace(0,  5, 11),
    np.linspace(6, 20,  8),
    np.linspace(25, 50, 6)
]))
GEN_PT_BINS = list(np.linspace(0, 100, 21))
