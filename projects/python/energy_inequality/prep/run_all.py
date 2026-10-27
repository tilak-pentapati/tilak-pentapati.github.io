"""Rebuild every processed output from raw/ in order.

    /Users/jay/anaconda3/envs/envs563/bin/python -I projects/python/energy_inequality/prep/run_all.py

Pass --fetch to download any raw file missing from raw/ (otherwise the run is offline).
"""
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = ["00_fetch_raw.py", "01_geography.py", "02_energy_panels.py", "03_covariates.py", "04_audit.py"]


def main() -> None:
    for step in STEPS:
        args = [sys.executable, "-I", "-B", str(HERE / step)]
        if step == "00_fetch_raw.py" and "--fetch" not in sys.argv:
            args.append("--offline")
        t = time.time()
        print(f"\n=== {step}", flush=True)
        subprocess.run(args, check=True)
        print(f"=== {step} done in {time.time() - t:.1f}s", flush=True)


if __name__ == "__main__":
    main()
