"""
stages/hadronic_seeding.py — seeding-stage plots for hadronic HYDDRA.

Plots:
  - reco_seed_{cosTheta,decayAngle,pOverE,dxySignif,mass}
    (nTracks omitted: all seeds are 2-track by construction)
"""
from ..src.config  import HADRONIC_RECO_OBSERVABLES
from ..src.plotter import plot_reco_observable


def make_plots(tdir, gf, sv_sig, sv_bkg):
    print("  [had/seeding] Reco observables...")
    for obs_key, obs_cfg in HADRONIC_RECO_OBSERVABLES.items():
        if obs_key == "nTracks":
            continue  # all seeds are 2-track; plot is uninformative
        plot_reco_observable(tdir, gf, sv_sig, "seed", obs_key, obs_cfg, sv_bkg, hadronic=True)
        print(f"    [reco_seed_{obs_key}] done")
