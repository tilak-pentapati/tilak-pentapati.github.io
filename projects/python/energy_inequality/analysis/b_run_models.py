"""Entry point for the convergence and explanatory-modelling half of the
energy-inequality analysis (the `run_models.py` step; prefixed `b_` to keep
it apart from the parallel ESDA/dynamics scripts).

    /Users/jay/anaconda3/envs/envs563/bin/python -I -W ignore \
        projects/python/energy_inequality/analysis/b_run_models.py [--force] [--only 1,2,3,4]

Sections
  1  b_01_clubs.py         Phillips-Sul club convergence, join counts, club profiles
  2  b_02_cross_section.py OLS / spatial econometric suite (LSOA, MSOA), GWR + MGWR (MSOA)
  3  b_03_change.py        Change models (2010-12 -> 2022-24; 2019 -> 2022), beta-convergence, LAD solar
  4  b_04_spacetime.py     Two-way FE event study and basis-function trend surface

Outputs go to data/data_python/energy_inequality/results/models/{tables,figures,cache}.
Heavy fits (clubs, spatial models, GWR/MGWR, trend surface) are pickled in
cache/ and reused unless --force is given, so the page only reads results.
Section 3 reads the VIF-screened covariate list cached by section 2.
Approximate run time from an empty cache on a 10-core laptop: ~6-7 h,
dominated by the two MGWR fits at n = 7,264 (b_mgwr_fast.py: ~40-50 min of
backfitting plus ~2.5-3 h of exact inference each); with a warm cache, ~4 min.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import b_common as bc


def main():
    force = "--force" in sys.argv
    only = None
    if "--only" in sys.argv:
        only = {int(s) for s in sys.argv[sys.argv.index("--only") + 1].split(",")}
    import b_01_clubs
    import b_02_cross_section
    import b_03_change
    import b_04_spacetime

    steps = {1: b_01_clubs.main, 2: b_02_cross_section.main, 3: b_03_change.main, 4: b_04_spacetime.main}
    for k, fn in steps.items():
        if only and k not in only:
            continue
        t0 = time.time()
        bc.log(f"=== section {k}: {fn.__module__} ===")
        fn(force=force)
        bc.log(f"=== section {k} done in {(time.time() - t0) / 60:.1f} min ===")


if __name__ == "__main__":
    main()
