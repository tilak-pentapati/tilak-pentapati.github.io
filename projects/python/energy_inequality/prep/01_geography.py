"""Geography: LSOA21 and LAD geometry for England & Wales, and LSOA11 -> LSOA21 weights.

Outputs
- processed/lsoa21_ew.parquet   GeoParquet, 35,672 LSOA21 with MSOA21, LAD22 and LAD26 codes
- processed/lad21_ew.parquet    GeoParquet, E&W local authority districts (May 2026 vintage, BGC)
- processed/lsoa11_lsoa21_weights.parquet  apportionment weights (one row per LSOA11 x LSOA21 pair)

Weights method (used to move 2010-14 DESNZ data from 2011 to 2021 LSOAs)
Each OA11 -> OA21 pair from the ONS exact-fit OA lookup is given the Census 2021
household count of its OA21 (TS041), divided equally among the OA21's OA11
parents when an OA21 was formed from several OA11s (merges / complex changes).
Pair households are summed to LSOA11 x LSOA21 pairs, and each LSOA11's values are
shared out in proportion to those household counts. For unchanged (U) and merged
(M) LSOA21s the share is exactly 1, so totals and meters are simply summed; only
split (S) and complex (X) cases depend on the household weights.
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # keep __pycache__ out of the repo
sys.path.insert(0, str(Path(__file__).resolve().parent))
import geopandas as gpd  # noqa: E402
import pandas as pd  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

from _common import AUDIT, CRS, GEO, PROCESSED, is_ew, log, read_census, read_lookup  # noqa: E402


def lsoa_geometry() -> gpd.GeoDataFrame:
    g = gpd.read_file(GEO / "LSOA_2021_EW_BGC_V5.geojson", engine="pyogrio")
    g = g.to_crs(CRS)[["LSOA21CD", "LSOA21NM", "BNG_E", "BNG_N", "LAT", "LONG", "geometry"]]
    g = g[is_ew(g["LSOA21CD"])]
    bad = ~g.geometry.is_valid
    if bad.any():
        # A handful of published BGC polygons self-intersect; keep only the polygonal parts after repair.
        fixed = g.loc[bad, "geometry"].make_valid()
        fixed = fixed.apply(lambda geom: unary_union([p for p in getattr(geom, "geoms", [geom])
                                                      if p.geom_type in ("Polygon", "MultiPolygon")]))
        g.loc[bad, "geometry"] = fixed
        log(f"repaired {int(bad.sum())} invalid LSOA21 geometries: {', '.join(g.loc[bad, 'LSOA21CD'])}")
        assert g.geometry.is_valid.all()

    oa = read_lookup("OA21_LSOA21_MSOA21_LAD22_EW_LU_v3.csv")
    msoa = oa.drop_duplicates("LSOA21CD")[["LSOA21CD", "MSOA21CD", "MSOA21NM", "LAD22CD", "LAD22NM"]]
    assert oa.groupby("LSOA21CD")["MSOA21CD"].nunique().max() == 1, "LSOA21 not nested in MSOA21"

    lad26 = read_lookup("LSOA21_WD26_LAD26_EW_LU.csv").drop_duplicates("LSOA21CD")[["LSOA21CD", "LAD26CD", "LAD26NM"]]
    g = g.merge(msoa, on="LSOA21CD", how="left", validate="1:1").merge(lad26, on="LSOA21CD", how="left", validate="1:1")
    if g[["MSOA21CD", "LAD22CD", "LAD26CD"]].isna().any().any():
        raise ValueError("LSOA21 without MSOA/LAD lookup")

    hh = read_census("ts041")["Number of households"].rename("households_2021")
    g = g.merge(hh, left_on="LSOA21CD", right_index=True, how="left")
    g["country"] = g["LSOA21CD"].str[0].map({"E": "England", "W": "Wales"})
    g["area_km2"] = g.geometry.area / 1e6
    g = g.rename(columns={
        "LSOA21CD": "lsoa21cd", "LSOA21NM": "lsoa21nm", "MSOA21CD": "msoa21cd", "MSOA21NM": "msoa21nm",
        "LAD22CD": "lad22cd", "LAD22NM": "lad22nm", "LAD26CD": "lad_cd", "LAD26NM": "lad_nm",
        "BNG_E": "bng_e", "BNG_N": "bng_n", "LAT": "lat", "LONG": "lon",
    })
    g["lad_vintage"] = "LAD May 2026"
    cols = ["lsoa21cd", "lsoa21nm", "msoa21cd", "msoa21nm", "lad_cd", "lad_nm", "lad_vintage", "lad22cd", "lad22nm",
            "country", "households_2021", "area_km2", "bng_e", "bng_n", "lat", "lon", "geometry"]
    return g[cols].sort_values("lsoa21cd").reset_index(drop=True)


def lad_geometry(lsoa: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    lad = gpd.read_file(GEO / "LAD_MAY_2026_UK_BGC.geojson", engine="pyogrio").to_crs(CRS)
    lad = lad[is_ew(lad["LAD26CD"])][["LAD26CD", "LAD26NM", "geometry"]]
    lad = lad.rename(columns={"LAD26CD": "lad_cd", "LAD26NM": "lad_nm"})
    have, need = set(lad["lad_cd"]), set(lsoa["lad_cd"])
    if have != need:
        raise ValueError(f"LAD boundary/lookup mismatch: {sorted(have ^ need)[:10]}")
    agg = lsoa.groupby("lad_cd").agg(n_lsoa21=("lsoa21cd", "size"), households_2021=("households_2021", "sum"))
    lad = lad.merge(agg, on="lad_cd", how="left")
    lad["country"] = lad["lad_cd"].str[0].map({"E": "England", "W": "Wales"})
    lad["lad_vintage"] = "LAD May 2026"
    lad["area_km2"] = lad.geometry.area / 1e6
    return lad.sort_values("lad_cd").reset_index(drop=True)


def lsoa11_to_lsoa21_weights() -> pd.DataFrame:
    oa = read_lookup("OA11_OA21_LAD22_EW_LU_v3.csv")[["OA11CD", "OA21CD", "CHNGIND"]]
    oa11 = read_lookup("OA11_LSOA11_MSOA11_LAD11_EW_LU.csv")[["OA11CD", "LSOA11CD"]]
    oa21 = read_lookup("OA21_LSOA21_MSOA21_LAD22_EW_LU_v3.csv")[["OA21CD", "LSOA21CD"]]
    hh_oa = read_census("ts041", "oa")["Number of households"].rename("hh_oa21")

    p = oa.merge(oa11, on="OA11CD", how="left").merge(oa21, on="OA21CD", how="left")
    if p[["LSOA11CD", "LSOA21CD"]].isna().any().any():
        raise ValueError("OA lookup rows without LSOA parents")
    p = p.merge(hh_oa, left_on="OA21CD", right_index=True, how="left")
    p["n_parents"] = p.groupby("OA21CD")["OA11CD"].transform("nunique")
    p["hh_pair"] = p["hh_oa21"].fillna(0) / p["n_parents"]
    w = p.groupby(["LSOA11CD", "LSOA21CD"], as_index=False)["hh_pair"].sum()

    lk = read_lookup("LSOA11_LSOA21_LAD22_EW_LU_v3.csv")[["LSOA11CD", "LSOA21CD", "CHGIND"]]
    w = lk.merge(w, on=["LSOA11CD", "LSOA21CD"], how="outer", indicator=True)
    only_oa = (w["_merge"] == "right_only").sum()
    only_lk = (w["_merge"] == "left_only").sum()
    log(f"weights: {len(lk):,} LSOA lookup pairs; pairs only in OA-derived set {only_oa}, only in LSOA lookup {only_lk}")
    # The LSOA exact-fit lookup is authoritative for which pairs exist.
    w = w[w["_merge"] != "right_only"].drop(columns="_merge")
    w["hh_pair"] = w["hh_pair"].fillna(0.0)
    tot = w.groupby("LSOA11CD")["hh_pair"].transform("sum")
    n = w.groupby("LSOA11CD")["LSOA21CD"].transform("size")
    w["share"] = (w["hh_pair"] / tot).where(tot > 0, 1.0 / n)
    w["weight_method"] = "oa_households"
    w.loc[n == 1, "weight_method"] = "whole"
    w.loc[(n > 1) & (tot <= 0), "weight_method"] = "equal_fallback"
    w.loc[n == 1, "share"] = 1.0

    # LSOA21-level change class: X if any complex link, else M if several sources, else S/U.
    src = w.groupby("LSOA21CD").agg(n_src=("LSOA11CD", "nunique"), any_x=("CHGIND", lambda s: (s == "X").any()),
                                   any_s=("CHGIND", lambda s: (s == "S").any()))
    cls = pd.Series("U", index=src.index)
    cls[src["any_s"]] = "S"
    cls[src["n_src"] > 1] = "M"
    cls[src["any_x"]] = "X"
    w = w.merge(cls.rename("lsoa21_change"), left_on="LSOA21CD", right_index=True)
    w = w.rename(columns={"LSOA11CD": "lsoa11cd", "LSOA21CD": "lsoa21cd", "CHGIND": "chgind"})
    return w[["lsoa11cd", "lsoa21cd", "chgind", "lsoa21_change", "hh_pair", "share", "weight_method"]]


def main() -> None:
    lsoa = lsoa_geometry()
    lsoa.to_parquet(PROCESSED / "lsoa21_ew.parquet", index=False)
    log(f"lsoa21_ew: {len(lsoa):,} rows, {lsoa['lad_cd'].nunique()} LADs, {lsoa['msoa21cd'].nunique():,} MSOAs")

    lad = lad_geometry(lsoa)
    lad.to_parquet(PROCESSED / "lad21_ew.parquet", index=False)
    log(f"lad21_ew: {len(lad)} LADs")

    w = lsoa11_to_lsoa21_weights()
    w.to_parquet(PROCESSED / "lsoa11_lsoa21_weights.parquet", index=False)
    chk = w.groupby("lsoa11cd")["share"].sum()
    assert (chk.sub(1).abs() < 1e-9).all(), "shares do not sum to 1 per LSOA11"
    missing21 = set(lsoa["lsoa21cd"]) - set(w["lsoa21cd"])
    assert not missing21, f"LSOA21 with no 2011 source: {sorted(missing21)[:5]}"

    summ = pd.DataFrame({
        "pairs_by_chgind": w["chgind"].value_counts(),
        "lsoa21_by_class": w.drop_duplicates("lsoa21cd")["lsoa21_change"].value_counts(),
        "lsoa11_by_chgind": w.drop_duplicates("lsoa11cd")["chgind"].value_counts(),
    }).fillna(0).astype(int)
    summ.index.name = "change"
    summ.to_csv(AUDIT / "lsoa11_lsoa21_change_summary.csv")
    log(summ.to_string())
    log(w["weight_method"].value_counts().to_string())


if __name__ == "__main__":
    main()
