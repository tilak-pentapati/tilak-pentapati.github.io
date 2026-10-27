"""Entry point for the descriptive and dynamics analysis.

    /Users/jay/anaconda3/envs/envs563/bin/python -I -W ignore \
        projects/python/energy_inequality/analysis/run_dynamics.py [--only esda,markov] [--fresh]

Steps (each writes tidy tables to results/dynamics/tables, figures to
results/dynamics/figures and small JSON summaries to results/dynamics):

    esda         a_01_esda.py         summary stats, KDE, Fisher-Jenks maps, Moran's I, LISA
    inequality   a_02_inequality.py   Theil/MLD decompositions, Gini, spatial Gini, dispersion vs mean, shocks
    markov       a_03_markov.py       classic/spatial/LISA Markov, homogeneity tests, rose diagrams
    mobility     a_04_mobility.py     Van Kerm decomposition, Kendall/spatial tau, LIMA
    sensitivity  a_05_sensitivity.py  excluding outliers, gas comparison, MSOA scale, permutation x2

Heavy intermediate results (spatial weights, LISA for every year with 9,999
permutations, local tau) are cached in results/dynamics/cache and reused;
``--fresh`` deletes the cache first. ``markov`` reads the LISA cache written
by ``esda``. A manifest of every output is written to results/dynamics/MANIFEST.json.
Random seeds are fixed (a_common.SEED), so reruns reproduce the tables exactly.
"""
from __future__ import annotations

import argparse
import importlib
import json
import platform
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a_common import CACHE, RESULTS, log  # noqa: E402

STEPS = {"esda": "a_01_esda", "inequality": "a_02_inequality", "markov": "a_03_markov",
         "mobility": "a_04_mobility", "sensitivity": "a_05_sensitivity"}


def manifest() -> None:
    import esda
    import giddy
    import inequality
    import libpysal
    import mapclassify
    import numpy
    import pandas

    files = sorted(str(p.relative_to(RESULTS)) for p in RESULTS.rglob("*")
                   if p.is_file() and "cache" not in p.parts and p.name != "MANIFEST.json")
    meta = {"generated": time.strftime("%Y-%m-%d %H:%M"), "python": platform.python_version(),
            "packages": {m.__name__: getattr(m, "__version__", "?") for m in
                         (esda, giddy, inequality, libpysal, mapclassify, numpy, pandas)},
            "files": files}
    (RESULTS / "MANIFEST.json").write_text(json.dumps(meta, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=",".join(STEPS), help="comma-separated subset of " + ", ".join(STEPS))
    ap.add_argument("--fresh", action="store_true", help="delete cached weights/LISA before running")
    a = ap.parse_args()
    if a.fresh and CACHE.exists():
        shutil.rmtree(CACHE)
        CACHE.mkdir(parents=True)
    for name in [s.strip() for s in a.only.split(",") if s.strip()]:
        t0 = time.time()
        log(f"=== {name} ===")
        importlib.import_module(STEPS[name]).run()
        log(f"=== {name} finished in {time.time() - t0:.0f}s ===")
    manifest()


if __name__ == "__main__":
    main()
