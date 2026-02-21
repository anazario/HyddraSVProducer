"""
stages/filtering.py — filtering-stage plots for hyddraDiagPlots.

Plots:
  - reco_filtered_{cosTheta,decayAngle,pOverE,dxySignif}
"""
from ..src.config  import RECO_OBSERVABLES
from ..src.plotter import plot_reco_observable


def make_plots(tdir, gf, sv_sig, sv_bkg):
    print("  [filtering] Reco observables...")
    for obs_key, obs_cfg in RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "filtered", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_filtered_{obs_key}] done")
