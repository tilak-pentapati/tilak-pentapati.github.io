"""Audit the processed outputs and (re)write data/data_python/energy_inequality/DATA_AUDIT.md.

Includes revision comparisons against the two superseded DESNZ releases
(2010-21 and 2010-22) and a validation of the
2011 -> 2021 LSOA conversion using the 2010-21 release, which published 2015-2021
on 2011 LSOAs.
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # keep __pycache__ out of the repo
sys.path.insert(0, str(Path(__file__).resolve().parent))
import geopandas as gpd  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from _common import (  # noqa: E402
    AUDIT, DATA, ELEC_LSOA_OLD21, ELEC_LSOA_OLD22, PROCESSED, YEARS, cached_desnz, is_ew, log, md_table, safe_div,
)

import importlib  # noqa: E402

energy = importlib.import_module("02_energy_panels")


def compare(a: pd.DataFrame, b: pd.DataFrame, key: str, label: str) -> dict:
    """Compare mean kWh/meter of two LSOA tables for one year on a shared key."""
    m = a[[key, "meters", "total_kwh"]].merge(b[[key, "meters", "total_kwh"]], on=key, suffixes=("_a", "_b"))
    ma, mb = safe_div(m["total_kwh_a"], m["meters_a"]), safe_div(m["total_kwh_b"], m["meters_b"])
    ok = np.isfinite(ma) & np.isfinite(mb)
    pct = (mb[ok] - ma[ok]) / ma[ok]
    return {
        "comparison": label, "matched_lsoas": int(ok.sum()),
        "corr_mean": float(np.corrcoef(ma[ok], mb[ok])[0, 1]),
        "median_abs_pct_diff": float(np.median(np.abs(pct)) * 100),
        "p95_abs_pct_diff": float(np.quantile(np.abs(pct), 0.95) * 100),
        "n_abs_diff_gt_10pct": int((np.abs(pct) > 0.10).sum()),
        "ew_kwh_per_meter_a": float(m["total_kwh_a"].sum() / m["meters_a"].sum()),
        "ew_kwh_per_meter_b": float(m["total_kwh_b"].sum() / m["meters_b"].sum()),
    }


def revisions(new_raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    old21 = cached_desnz(ELEC_LSOA_OLD21, "elec", "elecold21")
    old22 = cached_desnz(ELEC_LSOA_OLD22, "elec", "elecold22")
    ew = lambda d: d[is_ew(d["lsoa_code_src"])]  # noqa: E731
    new_raw, old21, old22 = ew(new_raw), ew(old21), ew(old22)
    rows = []
    for y in YEARS:
        n = new_raw[new_raw["year"] == y]
        o22 = old22[old22["year"] == y]
        if len(o22):
            rows.append({"year": y, **compare(o22, n, "lsoa_code_src", "2010-22 release (a) vs 2010-24 release (b)")})
        o21 = old21[old21["year"] == y]
        if len(o21) and y <= 2014:
            rows.append({"year": y, **compare(o21, n, "lsoa_code_src", "2010-21 release (a) vs 2010-24 release (b)")})
    rev = pd.DataFrame(rows)

    # Validation of the 2011 -> 2021 conversion: convert the 2010-21 release's 2015-2021 (2011 LSOAs,
    # old processing system) and compare with the newest release's 2021-LSOA figures by change class.
    w = pd.read_parquet(PROCESSED / "lsoa11_lsoa21_weights.parquet")
    cls = w.drop_duplicates("lsoa21cd").set_index("lsoa21cd")["lsoa21_change"]
    val = []
    for y in [2015, 2018, 2021]:
        conv = energy.convert_lsoa11(old21[old21["year"] == y], w)
        pub = new_raw[new_raw["year"] == y].rename(columns={"lsoa_code_src": "lsoa21cd"})
        for c in ["U", "S", "M", "X"]:
            keep = set(cls[cls == c].index)
            r = compare(conv[conv["lsoa21cd"].isin(keep)], pub[pub["lsoa21cd"].isin(keep)], "lsoa21cd",
                        f"2010-21 release converted (a) vs 2010-24 published (b), class {c}")
            val.append({"year": y, "lsoa21_change": c, **r})
    val = pd.DataFrame(val)

    o18 = old21[old21["year"] == 2018]
    m18 = safe_div(o18["total_kwh"], o18["meters"])
    i = int(np.nanargmin(m18))
    code = o18["lsoa_code_src"].iloc[i]
    new18 = new_raw[(new_raw["year"] == 2018)]
    notes = {
        "old2018_min_code": code, "old2018_min_mean": float(m18[i]), "old2018_min_meters": float(o18["meters"].iloc[i]),
        "new2018_min_mean": float(np.nanmin(safe_div(new18["total_kwh"], new18["meters"]))),
        "lsoa21_codes_not_in_2011": int(len(set(new_raw.loc[new_raw["year"] == 2022, "lsoa_code_src"])
                                            - set(w["lsoa11cd"]))),
        "lsoa11_codes_without_2022_row": int(len(set(w["lsoa11cd"])
                                                 - set(new_raw.loc[new_raw["year"] == 2022, "lsoa_code_src"]))),
    }
    return rev, val, notes


def main() -> None:
    elec = pd.read_parquet(PROCESSED / "elec_lsoa_panel.parquet")
    gas = pd.read_parquet(PROCESSED / "gas_lsoa_panel.parquet")
    lad = pd.read_parquet(PROCESSED / "elec_lad_panel.parquet")
    cov = pd.read_parquet(PROCESSED / "covariates_lsoa21.parquet")
    pop = pd.read_parquet(PROCESSED / "population_lsoa_panel.parquet")
    lsoa = gpd.read_parquet(PROCESSED / "lsoa21_ew.parquet")
    lad_g = gpd.read_parquet(PROCESSED / "lad21_ew.parquet")
    w = pd.read_parquet(PROCESSED / "lsoa11_lsoa21_weights.parquet")

    new_raw = cached_desnz(energy.ELEC_LSOA, "elec", "elec2024")
    rev, val, notes = revisions(new_raw)
    rev.to_csv(AUDIT / "elec_release_revisions.csv", index=False)
    val.to_csv(AUDIT / "elec_conversion_validation.csv", index=False)

    cov_e = pd.read_csv(AUDIT / "elec_coverage_by_year.csv")
    cov_g = pd.read_csv(AUDIT / "gas_coverage_by_year.csv")
    ua_e = pd.read_csv(AUDIT / "elec_unallocated_by_year.csv")
    ua_g = pd.read_csv(AUDIT / "gas_unallocated_by_year.csv")
    chg = pd.read_csv(AUDIT / "lsoa11_lsoa21_change_summary.csv")

    # Conservation check for 2010-14: E&W meters before and after conversion.
    ew_raw = new_raw[is_ew(new_raw["lsoa_code_src"]) & new_raw["year"].le(2014)]
    cons = pd.DataFrame({
        "meters_on_2011_lsoas": ew_raw.groupby("year")["meters"].sum(),
        "meters_on_2021_lsoas": elec[elec["year"].le(2014)].groupby("year")["meters"].sum(),
        "kwh_on_2011_lsoas": ew_raw.groupby("year")["total_kwh"].sum(),
        "kwh_on_2021_lsoas": elec[elec["year"].le(2014)].groupby("year")["total_kwh"].sum(),
    }).reset_index()
    cons["max_rel_diff"] = np.maximum((cons.meters_on_2021_lsoas / cons.meters_on_2011_lsoas - 1).abs(),
                                      (cons.kwh_on_2021_lsoas / cons.kwh_on_2011_lsoas - 1).abs())

    conv_rows = elec[elec["year"].le(2014)].drop_duplicates("lsoa21cd")["lsoa21_change"].value_counts()
    lsoa_any_outlier = elec.groupby("lsoa21cd")["flag_outlier"].any().sum()
    lsoa_out_years = elec.groupby("lsoa21cd")["flag_outlier"].sum()

    def fmt_cov(c: pd.DataFrame) -> pd.DataFrame:
        t = c[["year", "lsoa21_rows", "published_or_converted", "missing", "partial_source", "small_n", "low_mean",
               "extreme_mean", "spike", "outlier_any", "meter_jump", "median_exact"]].copy()
        t.insert(1, "geography", np.where(t["year"] <= 2014, "LSOA11→21", "LSOA21"))
        return t.astype({k: int for k in t.columns if k not in ("geography",)})

    sanity = pd.DataFrame({
        "year": cov_e["year"],
        "elec E&W kWh/meter (LSOA sum)": cov_e["ew_kwh_per_meter"],
        "elec DESNZ LA-table E&W kWh/meter": ua_e["la_ew_kwh"] / ua_e["la_ew_meters"],
        "elec median LSOA mean": cov_e["median_lsoa_mean_kwh"],
        "elec meters (m)": cov_e["meters"] / 1e6,
        "gas E&W kWh/meter (LSOA sum)": cov_g["ew_kwh_per_meter"],
        "gas median LSOA mean": cov_g["median_lsoa_mean_kwh"],
        "gas meters (m)": cov_g["meters"] / 1e6,
        "population (m)": pop.groupby("year")["pop_total"].sum().reindex(cov_e["year"]).values / 1e6,
    })

    ua_tbl = pd.DataFrame({
        "year": ua_e["year"],
        "elec: LSOA-sum shortfall vs LA-table E&W meters %": ua_e["lsoa_vs_la_meter_gap_share"] * 100,
        "elec: E&W unallocated row, % of E&W meters": ua_e["la_table_ew_unallocated_meter_share"] * 100,
        "elec: GB LSOA-file unallocated meters": ua_e["gb_unallocated_meters"].astype(int),
        "elec: GB unallocated % meters": ua_e["gb_unallocated_meter_share"] * 100,
        "gas: LSOA-sum shortfall vs LA-table E&W meters %": ua_g["lsoa_vs_la_meter_gap_share"] * 100,
        "gas: GB unallocated % meters": ua_g["gb_unallocated_meter_share"] * 100,
    })

    lad_cov = lad.groupby("year").agg(lads=("lad_cd", "size"), la_table_matched=("la_table_meters", lambda s: s.notna().sum()),
                                      min_lsoa_share_meters=("lsoa_share_of_la_meters", "min"),
                                      median_lsoa_share_meters=("lsoa_share_of_la_meters", "median")).reset_index()

    rev_tbl = rev[["year", "comparison", "matched_lsoas", "corr_mean", "median_abs_pct_diff", "p95_abs_pct_diff",
                   "n_abs_diff_gt_10pct", "ew_kwh_per_meter_a", "ew_kwh_per_meter_b"]].copy()
    rev_tbl["comparison"] = rev_tbl["comparison"].str.replace(" release (a) vs 2010-24 release (b)", "", regex=False)
    val_tbl = val[["year", "lsoa21_change", "matched_lsoas", "corr_mean", "median_abs_pct_diff", "p95_abs_pct_diff",
                   "n_abs_diff_gt_10pct", "ew_kwh_per_meter_a", "ew_kwh_per_meter_b"]]
    vs, vu = val[val["lsoa21_change"] == "S"], val[val["lsoa21_change"] == "U"]
    ua_max_e = max(ua_e["gb_unallocated_meter_share"].max(), ua_e["lsoa_vs_la_meter_gap_share"].max()) * 100
    ua_max_g = max(ua_g["gb_unallocated_meter_share"].max(), ua_g["lsoa_vs_la_meter_gap_share"].max()) * 100

    miss = cov.isna().sum()
    miss = miss[miss > 0]
    files = []
    for name in ["lsoa21_ew.parquet", "lad21_ew.parquet", "elec_lsoa_panel.parquet", "gas_lsoa_panel.parquet",
                 "elec_lad_panel.parquet", "covariates_lsoa21.parquet", "population_lsoa_panel.parquet",
                 "lsoa11_lsoa21_weights.parquet"]:
        p = PROCESSED / name
        d = gpd.read_parquet(p) if name in ("lsoa21_ew.parquet", "lad21_ew.parquet") else pd.read_parquet(p)
        files.append({"file": f"processed/{name}", "rows": len(d), "cols": d.shape[1], "MB": p.stat().st_size / 1e6})
    files = pd.DataFrame(files)

    f4 = "{:,.4f}"
    md = f"""# Energy inequality project: data audit

Generated by `projects/python/energy_inequality/prep/04_audit.py` (re-run `run_all.py` to refresh).
Study area: England & Wales, 2021 LSOAs ({len(lsoa):,} LSOAs, {lsoa['msoa21cd'].nunique():,} MSOAs,
{len(lad_g)} LADs in the May 2026 vintage). Scotland is dropped. Years: **2010-2024** (15 years) from the
newest DESNZ release (*LSOA domestic electricity / gas 2010-2024*, published 18 December 2025; next update due
December 2026). Sources and licences: see `SOURCES.md`.

## Output files

{md_table(files, "{:,.1f}")}

- `elec_lsoa_panel` / `gas_lsoa_panel`: balanced long panels (every LSOA21 x every year; rows with no published
  value have `flag_missing = True`). Columns: `lsoa21cd, lad_cd, year, fuel, meters, total_kwh, mean_kwh,
  median_kwh, median_exact, geog_src, lsoa21_change, n_src_lsoa11, n_src_missing, meters_per_hh2021, robust_z`,
  the `flag_*` columns and `outlier_reason`.
- `elec_lad_panel`: LAD (May 2026) x year. `mean_kwh = total_kwh / meters` from LSOA sums (never an average of
  LSOA means). Also carries the DESNZ LA-table domestic figures (`la_table_*`, 2015-2024 for all LADs; 2012-14
  only where the LA code is unchanged) and the share of LA-table meters/kWh that the LSOA sums capture.
- `covariates_lsoa21`: Census 2021 shares, RUC 2021, IMD 2025 / WIMD 2025, IMD 2019 / WIMD 2019 converted to
  LSOA21, EPC (MSOA-level) measures. `population_lsoa_panel`: ONS mid-year estimates 2011-2024.
- `lsoa11_lsoa21_weights`: the apportionment weights used for 2010-14 (kept so the conversion is inspectable).
- Geometry: ONS LSOA 2021 BGC (V5) and LAD May 2026 BGC, EPSG:27700, GeoParquet. Seven published LSOA polygons
  that self-intersect are repaired with `make_valid` (polygonal parts kept); total area {lsoa['area_km2'].sum():,.0f} km².
  `lad21_ew` holds the latest (May 2026) LAD vintage, which is what `lad_cd` means in every output.

## Coverage per year

### Electricity (LSOA panel)

{md_table(fmt_cov(cov_e))}

### Gas (LSOA panel, weather-corrected)

{md_table(fmt_cov(cov_g))}

Column meanings: `missing` = no value published (electricity: mostly new-build 2021 LSOAs in the first years;
gas: LSOAs with no or too few gas meters, i.e. largely off-grid areas); `partial_source` = 2010-14 LSOA21 built
from several LSOA11s of which some were unpublished; `median_exact` = rows that keep a published median.

## Outlier rule

Applied identically to every year and both fuels (`02_energy_panels.py`). Nothing is dropped; `flag_outlier` is
the union of:

1. `small_n`: fewer than {energy.MIN_METERS} meters.
2. `low_mean`: mean below {energy.LOW_MEAN_FRAC} x that year's E&W median LSOA mean.
3. `extreme_mean`: within-year robust z-score of log(mean) beyond ±{energy.ROBUST_Z:g} (median / 1.4826·MAD).
4. `spike`: a one-year spike, log(mean) differs from both neighbouring years by more than log(1.75) in the same direction.

Informational only: `meter_jump` (|Δlog meters| > log 1.5 year on year; usually new-build growth or
re-allocation), `meters_per_hh2021` (meters ÷ Census 2021 households).

Electricity: {int(elec['flag_outlier'].sum()):,} flagged LSOA-years ({elec['flag_outlier'].mean()*100:.2f}%), in
{int(lsoa_any_outlier):,} LSOAs; {int((lsoa_out_years >= 3).sum()):,} LSOAs are flagged in 3+ years (persistent
high-consumption places such as all-electric rural areas rather than errors). Gas: {int(gas['flag_outlier'].sum()):,}
flagged LSOA-years ({gas['flag_outlier'].mean()*100:.2f}%), mostly `small_n` (few gas meters).

The 2018 minimum-mean outlier from the superseded 2010-21 release ({notes['old2018_min_code']},
{notes['old2018_min_mean']:,.0f} kWh/meter on {notes['old2018_min_meters']:,.0f} meters) is gone in the newest release:
its 2018 minimum LSOA mean is {notes['new2018_min_mean']:,.0f} kWh/meter.

## National sanity table (E&W)

{md_table(sanity, "{:,.2f}")}

The LSOA-sum and DESNZ LA-table means agree to within a few kWh. The 2022 fall (energy price crisis) and the
2020 rise (COVID-19 home working) are visible in both fuels. Gas is weather-corrected by DESNZ; electricity is not.

## Unallocated meters per year

{md_table(ua_tbl, "{:,.3f}")}

DESNZ cannot allocate some meters to an LSOA (they are in the LA tables but not the LSOA tables) and warns that
LSOA figures should not be summed to LA/regional level. The share is small (at most {ua_max_e:.2f}% of
electricity and {ua_max_g:.2f}% of gas meters in any year), so
`elec_lad_panel` sums LSOAs and records the gap: per LAD, the LSOA sums capture this share of LA-table meters:

{md_table(lad_cov, f4)}

For 2010-14 the LA tables use pre-2015 LA codes; only unchanged codes match, and the E&W LA-table total there
excludes the GB-wide unallocated row (so the "shortfall" can be slightly negative).

## 2011 → 2021 LSOA conversion (2010-2014 only)

The newest release publishes 2015-2024 on 2021 LSOAs and only 2010-2014 on 2011 LSOAs, so only those five years
are converted. ONS LSOA11→LSOA21 exact-fit lookup (change indicator per link):

{md_table(chg)}

Method (`01_geography.py`, `02_energy_panels.py`): each OA11→OA21 link gets the Census 2021 household count of its
OA21 (shared equally between parents when OA21s were merged), summed to LSOA11×LSOA21 pairs. Meters and total kWh
of each LSOA11 are apportioned by those household shares and summed per LSOA21; mean = total ÷ meters. U and M
LSOAs use share 1 (exact sums). Only S and X depend on the household weights ({int((w['weight_method']=='oa_households').sum()):,}
pairs; no zero-household fallback was needed: {int((w['weight_method']=='equal_fallback').sum())}). The OA-derived pair list matched
the LSOA lookup exactly (0 pairs missing either way). LSOA21 by class in converted years: {conv_rows.to_dict()}.
Converted meters are non-integer.

Medians are kept only where exact (unchanged 1:1 LSOAs): {int(elec.loc[elec.year.le(2014),'median_exact'].sum()):,} of
{int(elec.year.le(2014).sum()):,} electricity rows in 2010-14.

Conservation check (E&W totals before vs after conversion):

{md_table(cons, "{:,.6g}")}

### Validation of the conversion

The superseded 2010-21 release published 2015-2021 on 2011 LSOAs (old processing system). Converting it with the
same weights and comparing with the newest release's 2021-LSOA figures separates pure revision (class U) from
revision + conversion error (S, M, X). Class U therefore measures the old-vs-new DESNZ processing difference
(median {vu['median_abs_pct_diff'].median():.1f}%, p95 {vu['p95_abs_pct_diff'].median():.1f}%), and class S shows that household-weighted splitting
cannot recover within-parent differences (median {vs['median_abs_pct_diff'].median():.1f}%, correlation {vs['corr_mean'].min():.2f}-{vs['corr_mean'].max():.2f}):

{md_table(val_tbl, "{:,.3f}")}

## Release-to-release revisions (electricity)

Same LSOA codes and year, mean kWh/meter. (a) = superseded release, (b) = newest 2010-2024 release.

{md_table(rev_tbl, "{:,.3f}")}

2015-2021 were re-processed with DESNZ's new system in the 2010-22 release and revised again since; 2010-14
are identical in all three releases, i.e. they were never re-processed with the new system (see caveat 3).
{notes['lsoa21_codes_not_in_2011']:,} LSOA codes in the 2022 table do not exist in the 2011 geography and
{notes['lsoa11_codes_without_2022_row']:,} 2011 LSOA codes have no 2022 row, which is why joining 2015+ data to 2011
boundaries would produce spurious zeros and gaps.

## Covariates

Columns with missing values in `covariates_lsoa21` ({cov.shape[1]} columns):

{md_table(miss.rename('missing').reset_index().rename(columns={'index': 'column'}))}

- IMD (England) and WIMD (Wales) are missing by design for the other country; use
  `dep2025_pct_within_country` / `dep2019_pct_within_country` (rank percentile, 0 = most deprived, within country)
  when pooling, and do not compare English and Welsh scores directly.
- IMD 2019 / WIMD 2019: household-weighted averages of LSOA11 scores; ranks/deciles re-computed on the converted
  scores (`*_reranked`). Exact for unchanged LSOAs (verified: max score difference 0 for U).
- EPC (`epc25_*`, FYE2025, 2021 MSOAs) are MSOA values copied to every LSOA in the MSOA. `epc21_*`/`epc22_*` are
  the earlier ONS editions (2011 MSOAs) moved to 2021 MSOAs by household-weighted averaging. Missing EPC
  fuel percentages are ONS-suppressed small counts, not zeros.
- Population: mid-2011 to mid-2024 on 2021 LSOAs ({len(pop):,} rows). Mid-2010 is not published on 2021 LSOAs, so
  per-capita measures for 2010 are unavailable.

## Caveats the analysis must state

1. **Meters, not households or people.** DESNZ means are kWh per *meter*. Multi-rate (Economy 7) meters, landlord
   supplies, second homes, empty dwellings and communal heating distort small-area means.
2. **Do not sum LSOAs to LA/region as official totals.** Disclosure control and unallocated meters (at most
   {ua_max_e:.2f}% of electricity / {ua_max_g:.2f}% of gas meters per year) mean LSOA sums differ slightly from DESNZ LA tables; `elec_lad_panel` is computed as total kWh ÷ total
   meters from LSOA sums and records the gap.
3. **Two geographies and a processing break at 2014→2015.** 2010-14 are on 2011 LSOAs (converted here) and come
   from DESNZ's old processing system, never revised; 2015-24 are on 2021 LSOAs from the new system (introduced
   with the 2022 data, 2015-21 revised). For the same LSOAs and years the two systems differ by a median
   {vu['median_abs_pct_diff'].median():.1f}% (p95 {vu['p95_abs_pct_diff'].median():.1f}%). Treat 2014→2015 as a discontinuity in
   transition/mobility analyses and run key results on 2015-2024 as a robustness check.
4. **Split LSOAs share their parent's mean in 2010-14.** {int(conv_rows.get('S', 0)):,} LSOA21s created by splits (plus
   {int(conv_rows.get('X', 0))} complex) get identical 2010-14 means to their siblings, which mechanically inflates spatial
   autocorrelation and understates within-parent inequality in those years. Validation against 2015-21 shows split
   children really differ from the apportioned value by a median {vs['median_abs_pct_diff'].median():.1f}% (p95 ~{vs['p95_abs_pct_diff'].median():.0f}%).
   Medians are missing for all changed LSOAs in 2010-14. Consider excluding or down-weighting `lsoa21_change in
   (S, X)` for 2010-14, or starting LSOA-level dynamics in 2015.
5. **Electricity years run Feb–Jan** (non-half-hourly meters); gas years run Oct–Sep and are weather-corrected,
   electricity is not.
6. **Revisions.** DESNZ revises 2015 onwards every December; the panel is pinned to the 18 Dec 2025 release (2026
   release due December 2026).
7. **New-build LSOAs.** Some 2021 LSOAs were largely unbuilt in 2015-17 (missing values, `small_n`, `meter_jump`);
   their early means describe a handful of meters.
8. **Gas coverage.** About 1,000 LSOAs per year have no published gas value (off-grid or suppressed); gas means
   describe gas-connected households only and are not comparable to electricity across off-grid areas.
9. **Outliers are flagged, not removed.** Persistent `extreme_mean` areas are mostly real (electric heating, large
   dwellings); decide per analysis whether to exclude `flag_outlier` rows, and report sensitivity.
10. **Covariate timing.** Census, EPC and deprivation are cross-sections (2021, FYE2025, 2019/2025); they cannot
    explain year-to-year change on their own, and the 2019 deprivation indices are approximations for changed LSOAs.
11. **EPC is MSOA-level and selective.** Only dwellings with an EPC are covered (sales, lettings, new builds), and
    values are MSOA medians/percentages applied to all constituent LSOAs.
12. **England vs Wales deprivation indices are not comparable**; use within-country percentiles or country fixed
    effects.
13. **COVID-19 (2020-21) and the 2022 energy price shock** are large common shocks; year effects or detrending are
    needed before reading convergence/divergence.
"""
    (DATA / "DATA_AUDIT.md").write_text(md)
    log(f"wrote {DATA / 'DATA_AUDIT.md'}")
    log(rev_tbl.round(3).to_string(index=False))
    log(val_tbl.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
