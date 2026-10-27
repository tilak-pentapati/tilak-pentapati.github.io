"""Shared paths and readers for the energy-inequality data prep.

Every script in this folder is run with the project env in isolated mode:

    /Users/jay/anaconda3/envs/envs563/bin/python -I projects/python/energy_inequality/prep/<script>.py

`-I` stops Python adding the script folder to sys.path, so each script adds
it explicitly before importing this module.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
DATA = REPO / "data" / "data_python" / "energy_inequality"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
PROCESSED = DATA / "processed"
AUDIT = PROCESSED / "audit"

for _d in (INTERIM, PROCESSED, AUDIT):
    _d.mkdir(parents=True, exist_ok=True)

# Projected CRS of every ONS boundary used here (British National Grid).
CRS = "EPSG:27700"

# Newest DESNZ subnational release (published 18 Dec 2025).
ELEC_LSOA = RAW / "desnz" / "LSOA_domestic_elec_2010-2024.xlsx"
GAS_LSOA = RAW / "desnz" / "LSOA_domestic_gas_2010-2024.xlsx"
ELEC_LA = RAW / "desnz" / "Subnational_electricity_consumption_statistics_2005-2024.xlsx"
GAS_LA = RAW / "desnz" / "Subnational_gas_consumption_statistics_2005-2024.xlsx"
# Superseded releases kept only for revision comparisons.
ELEC_LSOA_OLD21 = RAW / "superseded_local" / "LSOA_domestic_elec_2010-21 (1).xlsx"
ELEC_LSOA_OLD22 = RAW / "superseded_local" / "LSOA_domestic_elec_2010-22 (1).xlsx"

YEARS = list(range(2010, 2025))
# Years the newest release publishes on 2011 LSOAs (all later years are on 2021 LSOAs).
YEARS_ON_LSOA11 = list(range(2010, 2015))

GEO = RAW / "geography"


def log(msg: str) -> None:
    print(msg, flush=True)


def is_ew(code: pd.Series) -> pd.Series:
    """True for England (E...) and Wales (W...) GSS codes; drops Scotland (S...)."""
    return code.astype(str).str[0].isin(["E", "W"])


def read_lookup(name: str) -> pd.DataFrame:
    """Read an ONS lookup CSV (they carry a UTF-8 BOM)."""
    return pd.read_csv(GEO / name, encoding="utf-8-sig", dtype=str)


def read_desnz_lsoa(path: Path, fuel: str) -> pd.DataFrame:
    """Read every year sheet of a DESNZ LSOA domestic consumption workbook.

    Returns a long frame: year, lad_code_src, msoa_code_src, lsoa_code_src,
    meters, total_kwh, mean_kwh, median_kwh. The 'Unallocated' row of each
    year is kept with lsoa_code_src == 'Unallocated'.
    """
    sheets = pd.read_excel(path, sheet_name=None, header=None, engine="calamine")
    out = []
    for name, d in sheets.items():
        if not re.fullmatch(r"\d{4}", str(name).strip()):
            continue
        first = d[0].astype(str)
        hdr = first.index[first.str.startswith("Local authority code")]
        if len(hdr) != 1:
            raise ValueError(f"{path.name}:{name}: header row not found")
        body = d.iloc[hdr[0] + 1 :, :10].copy()
        body.columns = [
            "lad_code_src", "lad_name_src", "msoa_code_src", "msoa_name_src",
            "lsoa_code_src", "lsoa_name_src", "meters", "total_kwh", "mean_kwh", "median_kwh",
        ]
        body = body[body["lsoa_code_src"].notna()]
        for c in ["meters", "total_kwh", "mean_kwh", "median_kwh"]:
            body[c] = pd.to_numeric(body[c], errors="coerce")
        body["year"] = int(name)
        out.append(body)
    df = pd.concat(out, ignore_index=True)
    df["fuel"] = fuel
    for c in ["lad_code_src", "msoa_code_src", "lsoa_code_src"]:
        df[c] = df[c].astype(str).str.strip()
    return df.sort_values(["year", "lsoa_code_src"]).reset_index(drop=True)


def cached_desnz(path: Path, fuel: str, tag: str) -> pd.DataFrame:
    """read_desnz_lsoa with a parquet cache in interim/ keyed on file size+mtime."""
    st = path.stat()
    cache = INTERIM / f"desnz_{tag}_{st.st_size}_{int(st.st_mtime)}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    for old in INTERIM.glob(f"desnz_{tag}_*.parquet"):
        old.unlink()
    df = read_desnz_lsoa(path, fuel)
    df.to_parquet(cache, index=False)
    return df


def read_census(table: str, level: str = "lsoa") -> pd.DataFrame:
    """Read one Nomis Census 2021 bulk table at a geography level.

    Columns are renamed to the category label with the table prefix and the
    '; measures: Value' suffix removed. Indexed by GSS code.
    """
    zpath = RAW / "census2021" / f"census2021-{table}.zip"
    with zipfile.ZipFile(zpath) as z:
        with z.open(f"census2021-{table}-{level}.csv") as fh:
            d = pd.read_csv(fh)
    d = d.rename(columns={"geography code": "code"}).set_index("code")
    d = d.drop(columns=["date", "geography"])
    d.columns = [re.sub(r"; measures: Value$", "", c).split(": ", 1)[1].strip() for c in d.columns]
    return d


def safe_div(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
    r[~np.isfinite(r)] = np.nan
    return r


def md_table(df: pd.DataFrame, floatfmt: str = "{:,.1f}") -> str:
    """Tiny markdown table writer (avoids a tabulate dependency)."""
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c, v in zip(cols, r.values):
            if c == "year" and pd.notna(v):
                cells.append(str(int(v)))
            elif isinstance(v, (float, np.floating)):
                if pd.isna(v):
                    cells.append("")
                elif float(v).is_integer() and abs(v) < 1e15 and "e" not in floatfmt:
                    cells.append(f"{int(v):,}")
                else:
                    cells.append(floatfmt.format(v))
            elif isinstance(v, (int, np.integer)):
                cells.append(f"{v:,}")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)
