"""DESNZ domestic electricity and gas panels on 2021 LSOAs (England & Wales), 2010-2024.

Outputs (processed/)
- elec_lsoa_panel.parquet, gas_lsoa_panel.parquet  balanced long panels, 35,672 LSOA21 x 15 years
- elec_lad_panel.parquet                            LAD (May 2026) panel, mean = total kWh / total meters
- audit/*.csv                                       numbers quoted in DATA_AUDIT.md

Harmonisation
- 2015-2024: taken as published by the newest release (already on 2021 LSOAs).
- 2010-2014: published on 2011 LSOAs; moved to 2021 LSOAs with the weights from
  01_geography.py. Meters and total kWh are apportioned (share = 1 for unchanged
  and merged LSOAs, household-weighted for splits / complex changes) and summed;
  mean is recomputed as total / meters. Median is kept only where it is exact,
  i.e. the LSOA21 is unchanged (U) from a single LSOA11.

Outlier rule (applied identically to every year and both fuels; components kept as flags)
- flag_small_n      meters < 50: too few meters for a stable mean.
- flag_low_mean     mean < 0.25 x the E&W median LSOA mean of that year.
- flag_extreme_mean |robust z| > 5, where z is the within-year robust z-score of
                    log(mean) (median / 1.4826*MAD).
- flag_spike        one-year spike: log(mean) differs from BOTH neighbouring years by
                    more than log(1.75) in the same direction.
- flag_outlier      any of the above. Nothing is dropped; analyses decide what to exclude.
Informational (not in flag_outlier): flag_missing, flag_meter_jump (|dlog meters| > log 1.5),
flag_partial_source (2010-14 only: some LSOA11 sources unpublished).
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # keep __pycache__ out of the repo
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from _common import (  # noqa: E402
    AUDIT, ELEC_LA, ELEC_LSOA, GAS_LA, GAS_LSOA, PROCESSED, YEARS, YEARS_ON_LSOA11,
    cached_desnz, is_ew, log, safe_div,
)

MIN_METERS = 50
LOW_MEAN_FRAC = 0.25
ROBUST_Z = 5.0
SPIKE = np.log(1.75)
METER_JUMP = np.log(1.5)


def convert_lsoa11(d11: pd.DataFrame, w: pd.DataFrame) -> pd.DataFrame:
    """Apportion 2011-LSOA rows to 2021 LSOAs."""
    x = d11.merge(w[["lsoa11cd", "lsoa21cd", "share", "chgind"]], left_on="lsoa_code_src", right_on="lsoa11cd",
                  how="left", validate="m:m")
    unmatched = x["lsoa21cd"].isna()
    if unmatched.any():
        raise ValueError(f"{unmatched.sum()} LSOA11 rows not in lookup, e.g. {x.loc[unmatched, 'lsoa_code_src'].iloc[0]}")
    x["src_missing"] = x["meters"].isna()
    x["m"] = x["meters"] * x["share"]
    x["t"] = x["total_kwh"] * x["share"]
    g = x.groupby(["year", "lsoa21cd"])
    out = g.agg(meters=("m", lambda s: s.sum(min_count=1)), total_kwh=("t", lambda s: s.sum(min_count=1)),
                n_src_lsoa11=("lsoa11cd", "nunique"), n_src_missing=("src_missing", "sum")).reset_index()
    # Median is exact only for a 1:1 unchanged LSOA (share 1, single source).
    exact = x[(x["chgind"] == "U") & (x["share"] == 1.0)][["year", "lsoa21cd", "median_kwh"]]
    exact = exact[exact.groupby(["year", "lsoa21cd"])["lsoa21cd"].transform("size") == 1]
    out = out.merge(exact, on=["year", "lsoa21cd"], how="left")
    out["median_exact"] = out["median_kwh"].notna()
    out["geog_src"] = "LSOA11 converted"
    return out


def build_panel(fuel: str, path: Path, lsoa: pd.DataFrame, w: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = cached_desnz(path, fuel, f"{fuel}2024")
    unalloc = raw[raw["lsoa_code_src"].str.contains("nallocated")][["year", "meters", "total_kwh"]]
    gb = raw[~raw["lsoa_code_src"].str.contains("nallocated")]
    gb_tot = gb.groupby("year")[["meters", "total_kwh"]].sum()
    ua = unalloc.groupby("year")[["meters", "total_kwh"]].sum()
    unalloc_tbl = pd.DataFrame({
        "gb_lsoa_meters": gb_tot["meters"], "gb_unallocated_meters": ua["meters"],
        "gb_unallocated_meter_share": ua["meters"] / (gb_tot["meters"] + ua["meters"]),
        "gb_unallocated_kwh_share": ua["total_kwh"] / (gb_tot["total_kwh"] + ua["total_kwh"]),
    })

    ew = gb[is_ew(gb["lsoa_code_src"])].copy()
    d21 = ew[~ew["year"].isin(YEARS_ON_LSOA11)].rename(columns={"lsoa_code_src": "lsoa21cd"})
    bad = set(d21["lsoa21cd"]) - set(lsoa["lsoa21cd"])
    if bad:
        raise ValueError(f"{fuel}: {len(bad)} published LSOA21 codes not in boundaries, e.g. {sorted(bad)[:3]}")
    lad_src = d21.drop_duplicates("lsoa21cd").set_index("lsoa21cd")["lad_code_src"]
    lad_mismatch = (lad_src.reindex(lsoa["lsoa21cd"]).values != lsoa["lad_cd"].values).sum()
    log(f"{fuel}: LSOA rows whose DESNZ LA code differs from LAD May 2026: {lad_mismatch}")
    d21 = d21[["year", "lsoa21cd", "meters", "total_kwh", "median_kwh"]].assign(
        n_src_lsoa11=np.nan, n_src_missing=np.nan, median_exact=lambda d: d["median_kwh"].notna(),
        geog_src="LSOA21 published")
    d11 = convert_lsoa11(ew[ew["year"].isin(YEARS_ON_LSOA11)], w)
    d = pd.concat([d11, d21], ignore_index=True)

    grid = pd.MultiIndex.from_product([lsoa["lsoa21cd"], YEARS], names=["lsoa21cd", "year"]).to_frame(index=False)
    p = grid.merge(d, on=["lsoa21cd", "year"], how="left", validate="1:1")
    p = p.merge(lsoa[["lsoa21cd", "lad_cd", "households_2021"]], on="lsoa21cd", how="left")
    chg = w.drop_duplicates("lsoa21cd").set_index("lsoa21cd")["lsoa21_change"]
    p["lsoa21_change"] = p["lsoa21cd"].map(chg)
    default_src = pd.Series(np.where(p["year"].isin(YEARS_ON_LSOA11), "LSOA11 converted", "LSOA21 published"), index=p.index)
    p["geog_src"] = p["geog_src"].fillna(default_src)
    p["mean_kwh"] = safe_div(p["total_kwh"], p["meters"])
    p["median_exact"] = p["median_exact"].fillna(False).astype(bool)
    p.loc[~p["median_exact"], "median_kwh"] = np.nan
    p["meters_per_hh2021"] = safe_div(p["meters"], p["households_2021"])
    p = add_flags(p)
    p["fuel"] = fuel
    cols = ["lsoa21cd", "lad_cd", "year", "fuel", "meters", "total_kwh", "mean_kwh", "median_kwh", "median_exact",
            "geog_src", "lsoa21_change", "n_src_lsoa11", "n_src_missing", "meters_per_hh2021", "robust_z",
            "flag_missing", "flag_partial_source", "flag_small_n", "flag_low_mean", "flag_extreme_mean",
            "flag_spike", "flag_meter_jump", "flag_outlier", "outlier_reason"]
    p = p[cols].sort_values(["lsoa21cd", "year"]).reset_index(drop=True)
    return p, unalloc_tbl


def add_flags(p: pd.DataFrame) -> pd.DataFrame:
    p = p.sort_values(["lsoa21cd", "year"]).copy()
    p["flag_missing"] = p["mean_kwh"].isna()
    p["flag_partial_source"] = p["n_src_missing"].fillna(0).gt(0) & ~p["flag_missing"]
    lm = np.log(p["mean_kwh"])
    yr = p.groupby("year")
    med = lm.groupby(p["year"]).transform("median")
    mad = (lm - med).abs().groupby(p["year"]).transform("median") * 1.4826
    p["robust_z"] = (lm - med) / mad
    p["flag_small_n"] = p["meters"] < MIN_METERS
    p["flag_low_mean"] = p["mean_kwh"] < LOW_MEAN_FRAC * yr["mean_kwh"].transform("median")
    p["flag_extreme_mean"] = p["robust_z"].abs() > ROBUST_Z
    g = lm.groupby(p["lsoa21cd"])
    d_prev = lm - g.shift(1)
    d_next = lm - g.shift(-1)
    p["flag_spike"] = ((d_prev > SPIKE) & (d_next > SPIKE)) | ((d_prev < -SPIKE) & (d_next < -SPIKE))
    dm = np.log(p["meters"]) - np.log(p["meters"]).groupby(p["lsoa21cd"]).shift(1)
    p["flag_meter_jump"] = dm.abs() > METER_JUMP
    comps = ["flag_small_n", "flag_low_mean", "flag_extreme_mean", "flag_spike"]
    p["flag_outlier"] = p[comps].any(axis=1)
    reason = pd.Series("", index=p.index)
    for c in comps:
        reason = reason + np.where(p[c], c.removeprefix("flag_") + ";", "")
    p["outlier_reason"] = reason.str.rstrip(";")
    return p


def read_la_table(path: Path) -> pd.DataFrame:
    """Domestic meters and kWh per LA (and E&W / unallocated rows) from a DESNZ subnational workbook."""
    sheets = pd.read_excel(path, sheet_name=None, header=None, engine="calamine")
    out = []
    for name, d in sheets.items():
        if not str(name).strip().isdigit() or int(name) not in YEARS:
            continue
        h = d.index[d[0].astype(str).str.strip() == "Code"][0]
        hdr = [" ".join(str(v).split()) for v in d.iloc[h]]
        def col(prefix: str) -> int:
            hits = [i for i, c in enumerate(hdr) if c.startswith(prefix) and (c.endswith(": Domestic") or c.endswith("All Domestic"))]
            assert len(hits) == 1, (path.name, name, prefix, hdr)
            return hits[0]
        b = d.iloc[h + 1:, [0, 1, 2, col("Number of meters"), col("Total consumption")]].copy()
        b.columns = ["code", "area", "la", "meters_k", "gwh"]
        b["year"] = int(name)
        out.append(b)
    t = pd.concat(out, ignore_index=True)
    t["code"] = t["code"].astype(str).str.strip()
    t["meters"] = pd.to_numeric(t["meters_k"], errors="coerce") * 1000
    t["total_kwh"] = pd.to_numeric(t["gwh"], errors="coerce") * 1e6
    return t[["year", "code", "area", "la", "meters", "total_kwh"]]


def ew_reconciliation(panel: pd.DataFrame, la: pd.DataFrame) -> pd.DataFrame:
    """Compare sum of E&W LSOA rows with the LA-table E&W total (incl. unallocated where published)."""
    s = panel.groupby("year")[["meters", "total_kwh"]].sum()
    rows = []
    for y in YEARS:
        t = la[la["year"] == y]
        k04 = t[t["code"] == "K04000001"]
        if len(k04):
            ew_m, ew_k, basis = k04["meters"].iloc[0], k04["total_kwh"].iloc[0], "LA table E&W incl. unallocated (K04000001)"
        else:
            sub = t[t["area"].astype(str).str.strip().isin(["England", "Wales"]) & (t["la"] == "All local authorities")]
            ew_m, ew_k, basis = sub["meters"].sum(), sub["total_kwh"].sum(), "LA table England + Wales (GB unallocated row not split by country)"
        ua = t[t["code"] == "Unallocated England and Wales"]
        rows.append({
            "year": y, "lsoa_sum_meters": s.loc[y, "meters"], "la_ew_meters": ew_m,
            "lsoa_sum_kwh": s.loc[y, "total_kwh"], "la_ew_kwh": ew_k,
            "lsoa_vs_la_meter_gap_share": 1 - s.loc[y, "meters"] / ew_m,
            "lsoa_vs_la_kwh_gap_share": 1 - s.loc[y, "total_kwh"] / ew_k,
            "la_table_ew_unallocated_meter_share": (ua["meters"].iloc[0] / ew_m) if len(ua) else np.nan,
            "basis": basis,
        })
    return pd.DataFrame(rows)


def lad_panel(p: pd.DataFrame, la: pd.DataFrame, lad_names: pd.Series, recon: pd.DataFrame) -> pd.DataFrame:
    q = p.assign(m_ok=p["meters"].where(~p["flag_outlier"]), t_ok=p["total_kwh"].where(~p["flag_outlier"]))
    g = q.groupby(["lad_cd", "year"])
    lad = g.agg(meters=("meters", "sum"), total_kwh=("total_kwh", "sum"), n_lsoa=("lsoa21cd", "size"),
                n_lsoa_missing=("flag_missing", "sum"), n_lsoa_outlier=("flag_outlier", "sum"),
                meters_ex_outliers=("m_ok", "sum"), total_kwh_ex_outliers=("t_ok", "sum")).reset_index()
    lad["mean_kwh"] = safe_div(lad["total_kwh"], lad["meters"])
    lad["mean_kwh_ex_outliers"] = safe_div(lad["total_kwh_ex_outliers"], lad["meters_ex_outliers"])
    lad["lad_nm"] = lad["lad_cd"].map(lad_names)
    t = la[la["code"].isin(set(lad["lad_cd"]))][["year", "code", "meters", "total_kwh"]]
    t = t.rename(columns={"code": "lad_cd", "meters": "la_table_meters", "total_kwh": "la_table_total_kwh"})
    lad = lad.merge(t, on=["lad_cd", "year"], how="left", validate="1:1")
    lad["la_table_mean_kwh"] = safe_div(lad["la_table_total_kwh"], lad["la_table_meters"])
    lad["lsoa_share_of_la_meters"] = safe_div(lad["meters"], lad["la_table_meters"])
    lad["lsoa_share_of_la_kwh"] = safe_div(lad["total_kwh"], lad["la_table_total_kwh"])
    lad = lad.merge(recon[["year", "lsoa_vs_la_meter_gap_share"]].rename(
        columns={"lsoa_vs_la_meter_gap_share": "ew_unallocated_meter_share"}), on="year", how="left")
    cols = ["lad_cd", "lad_nm", "year", "meters", "total_kwh", "mean_kwh", "n_lsoa", "n_lsoa_missing", "n_lsoa_outlier",
            "meters_ex_outliers", "total_kwh_ex_outliers", "mean_kwh_ex_outliers", "la_table_meters",
            "la_table_total_kwh", "la_table_mean_kwh", "lsoa_share_of_la_meters", "lsoa_share_of_la_kwh",
            "ew_unallocated_meter_share"]
    return lad[cols].sort_values(["lad_cd", "year"]).reset_index(drop=True)


def coverage_table(p: pd.DataFrame) -> pd.DataFrame:
    g = p.groupby("year")
    t = pd.DataFrame({
        "lsoa21_rows": g.size(), "published_or_converted": g["flag_missing"].apply(lambda s: (~s).sum()),
        "missing": g["flag_missing"].sum(), "partial_source": g["flag_partial_source"].sum(),
        "small_n": g["flag_small_n"].sum(), "low_mean": g["flag_low_mean"].sum(),
        "extreme_mean": g["flag_extreme_mean"].sum(), "spike": g["flag_spike"].sum(),
        "outlier_any": g["flag_outlier"].sum(), "meter_jump": g["flag_meter_jump"].sum(),
        "median_exact": g["median_exact"].sum(),
        "meters": g["meters"].sum(), "total_kwh": g["total_kwh"].sum(),
        "mean_lsoa_mean_kwh": g["mean_kwh"].mean(), "median_lsoa_mean_kwh": g["mean_kwh"].median(),
    })
    t["ew_kwh_per_meter"] = t["total_kwh"] / t["meters"]
    return t.reset_index()


def main() -> None:
    lsoa = pd.read_parquet(PROCESSED / "lsoa21_ew.parquet", columns=["lsoa21cd", "lad_cd", "lad_nm", "households_2021"])
    w = pd.read_parquet(PROCESSED / "lsoa11_lsoa21_weights.parquet")
    lad_names = lsoa.drop_duplicates("lad_cd").set_index("lad_cd")["lad_nm"]

    for fuel, path, la_path in [("elec", ELEC_LSOA, ELEC_LA), ("gas", GAS_LSOA, GAS_LA)]:
        p, unalloc = build_panel(fuel, path, lsoa, w)
        p.to_parquet(PROCESSED / f"{fuel}_lsoa_panel.parquet", index=False)
        la = read_la_table(la_path)
        recon = ew_reconciliation(p, la)
        recon = recon.merge(unalloc.reset_index(), on="year", how="left")
        recon.to_csv(AUDIT / f"{fuel}_unallocated_by_year.csv", index=False)
        cov = coverage_table(p)
        cov.to_csv(AUDIT / f"{fuel}_coverage_by_year.csv", index=False)
        log(f"{fuel}_lsoa_panel: {len(p):,} rows; outliers {int(p['flag_outlier'].sum()):,}; missing {int(p['flag_missing'].sum()):,}")
        log(cov[["year", "missing", "outlier_any", "ew_kwh_per_meter"]].to_string(index=False))
        if fuel == "elec":
            lad = lad_panel(p, la, lad_names, recon)
            lad.to_parquet(PROCESSED / "elec_lad_panel.parquet", index=False)
            log(f"elec_lad_panel: {len(lad):,} rows, {lad['lad_cd'].nunique()} LADs; "
                f"LA-table matched rows {int(lad['la_table_meters'].notna().sum())}")


if __name__ == "__main__":
    main()
