"""
stages/hadronic_merging.py — merging-stage plots for hadronic HYDDRA.

Plots:
  - reco_merged_{cosTheta,decayAngle,pOverE,dxySignif,mass,nTracks}
  - track_absorption  (nTracks at merge: gold vs silver, reused from leptonic)
"""
from ..src.config  import HADRONIC_RECO_OBSERVABLES
from ..src.plotter import plot_reco_observable
from .merging      import plot_track_absorption


def make_plots(tdir, gf, sv_sig, sv_bkg):
    print("  [had/merging] Reco observables...")
    for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
        plot_reco_observable(tdir, gf, sv_sig, "merged", obs_key, obs_cfg, sv_bkg)
        print(f"    [reco_merged_{obs_key}] done")
    if gf is not None and len(gf) > 0:
        print("  [had/merging] Track absorption...")
        plot_track_absorption(tdir, gf)
