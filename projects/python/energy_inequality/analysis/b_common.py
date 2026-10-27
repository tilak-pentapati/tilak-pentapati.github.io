"""Shared paths, data assembly, spatial weights and figure style for the
convergence and explanatory-modelling half of the energy-inequality analysis.

Run every script with the project env in isolated mode, e.g.

    /Users/jay/anaconda3/envs/envs563/bin/python -I projects/python/energy_inequality/analysis/b_run_models.py

`-I` stops Python adding the script folder to sys.path, so each script adds it
explicitly before importing this module.
"""
from __future__ import annotations

import json
import pickle
import time
from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
DATA = REPO / "data" / "data_python" / "energy_inequality"
PROCESSED = DATA / "processed"
OUT = DATA / "results" / "models"
FIG = OUT / "figures"
TAB = OUT / "tables"
CACHE = OUT / "cache"
for _d in (OUT, FIG, TAB, CACHE):
    _d.mkdir(parents=True, exist_ok=True)

CRS = "EPSG:27700"
YEARS = list(range(2010, 2025))
SEED = 20261026


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------
# Caching helpers: heavy fits are pickled so the page reads results only.
# --------------------------------------------------------------------------
def cached(name: str, fn, force: bool = False):
    path = CACHE / f"{name}.pkl"
    if path.exists() and not force:
        with open(path, "rb") as f:
            return pickle.load(f)
    obj = fn()
    with open(path, "wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    return obj


def save_table(df: pd.DataFrame, name: str, index: bool = False) -> Path:
    path = TAB / f"{name}.csv"
    df.to_csv(path, index=index, float_format="%.6g")
    return path


def save_json(obj, name: str) -> Path:
    path = OUT / f"{name}.json"
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return path


# --------------------------------------------------------------------------
# Figure style: large type so figures stay legible when scaled to a 375 px
# phone column; colour-blind-safe ramps only (cividis / viridis for ordered
# classes and magnitudes, a blue-grey-orange diverging map for signed values).
# --------------------------------------------------------------------------
INK = "#1f1f1f"
INK2 = "#52514e"
GRID = "#d9d8d4"
NS_GREY = "#e4e3df"  # non-significant local coefficients
# Okabe-Ito subset for <= 3 categorical series (validated CVD-safe set).
CAT3 = ["#0072B2", "#E69F00", "#009E73"]
DIVERGING = mpl.colors.LinearSegmentedColormap.from_list(
    "blue_grey_orange", ["#08519c", "#6baed6", "#f0efec", "#fdae6b", "#a63603"]
)


def set_style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "font.size": 13,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "lines.linewidth": 2.0,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
    })


def save_fig(fig, name: str) -> list[Path]:
    paths = []
    for ext in ("png", "svg"):
        p = FIG / f"{name}.{ext}"
        fig.savefig(p, bbox_inches="tight", facecolor="white")
        paths.append(p)
    plt.close(fig)
    return paths


def club_colours(n: int) -> list:
    """Ordered colours for clubs (club 1 = highest steady state = darkest)."""
    if n == 1:
        return [mpl.colormaps["viridis"](0.25)]
    return [mpl.colors.to_hex(mpl.colormaps["viridis"](v)) for v in np.linspace(0.0, 0.95, n)]


# --------------------------------------------------------------------------
# Data assembly
# --------------------------------------------------------------------------
def load_elec_lsoa() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "elec_lsoa_panel.parquet")


def load_cov_lsoa() -> pd.DataFrame:
    return pd.read_parquet(PROCESSED / "covariates_lsoa21.parquet")


def load_lsoa_geo() -> gpd.GeoDataFrame:
    return gpd.read_parquet(PROCESSED / "lsoa21_ew.parquet")


def load_lad_geo() -> gpd.GeoDataFrame:
    return gpd.read_parquet(PROCESSED / "lad21_ew.parquet")


def panel_wide(level: str, value: str = "mean_kwh") -> pd.DataFrame:
    """Unit x year table of kWh per meter at 'lsoa', 'msoa' or 'lad' level.

    Higher levels are always total kWh / total meters from LSOA sums (never an
    average of LSOA means)."""
    e = load_elec_lsoa()
    if level == "lsoa":
        return e.pivot(index="lsoa21cd", columns="year", values=value)
    if level == "msoa":
        cov = load_cov_lsoa()[["lsoa21cd", "msoa21cd"]]
        e = e.merge(cov, on="lsoa21cd")
        key = "msoa21cd"
    elif level == "lad":
        key = "lad_cd"
    else:
        raise ValueError(level)
    g = e.groupby([key, "year"])[["total_kwh", "meters"]].sum(min_count=1)
    g["mean_kwh"] = g.total_kwh / g.meters
    return g["mean_kwh"].unstack("year")


def window_mean(level: str, years) -> pd.Series:
    """Meter-weighted mean kWh per meter over a window of years."""
    e = load_elec_lsoa()
    e = e[e.year.isin(list(years))]
    if level == "lsoa":
        key = "lsoa21cd"
    elif level == "msoa":
        e = e.merge(load_cov_lsoa()[["lsoa21cd", "msoa21cd"]], on="lsoa21cd")
        key = "msoa21cd"
    else:
        key = "lad_cd"
    g = e.groupby(key)[["total_kwh", "meters"]].sum(min_count=1)
    n = e.groupby(key)["mean_kwh"].count()
    out = g.total_kwh / g.meters
    out[n < len(list(years))] = np.nan  # require every year in the window
    return out


# Covariates. Shares of households unless noted. Reference categories:
# gas central heating (heating), semi/terraced (dwelling type), owner-occupied (tenure).
HH_SHARES = [
    "share_ch_electric", "share_ch_oil", "share_ch_lpg", "share_ch_solid_wood",
    "share_ch_renewable", "share_ch_district", "share_ch_mains_gas", "share_no_central_heating",
    "share_detached", "share_semi", "share_terraced", "share_flat",
    "share_owned_outright", "share_owned_mortgage", "share_social_rent", "share_private_rent",
    "share_1person_hh", "share_overcrowded_bedrooms", "share_underoccupied_bedrooms",
    "share_no_car", "mean_hh_size", "mean_rooms", "mean_bedrooms",
]
PERSON_SHARES = [
    "share_age_0_14", "share_age_65plus", "share_age_75plus", "share_retired", "share_employed",
    "share_nssec_higher_mgr_prof", "share_nssec_routine", "share_eth_white_british",
]


def covariates(level: str) -> pd.DataFrame:
    """Analysis covariates at 'lsoa' or 'msoa' (or 'lad') level.

    LSOA values are aggregated with household weights (person weights for
    person-based shares). EPC measures are native MSOA values."""
    c = load_cov_lsoa().copy()
    c["share_ch_offgas"] = c[["share_ch_oil", "share_ch_lpg", "share_ch_solid_wood"]].sum(axis=1)
    c["rural"] = (c.urban_rural == "Rural").astype(float)
    c["wales"] = (c.country == "Wales").astype(float)
    c["dep19_pct"] = c["dep2019_pct_within_country"]  # 0 = most deprived, within country
    c["dep25_pct"] = c["dep2025_pct_within_country"]
    c["epc_ee"] = c["epc25_median_ee_existing_dwellings"]
    c["epc_c_plus"] = c["epc25_pct_band_c_plus_all_dwellings"]
    hh = HH_SHARES + ["share_ch_offgas", "rural", "wales", "dep19_pct", "dep25_pct"]
    if level == "lsoa":
        c["log_density"] = np.log(c.pop_density_2021_per_km2)
        return c.set_index("lsoa21cd")
    key = {"msoa": "msoa21cd", "lad": "lad_cd"}[level]
    w_h = c.households_2021
    w_p = c.usual_residents_2021
    agg = {}
    grp = c.groupby(key)
    for col in hh:
        agg[col] = (c[col] * w_h).groupby(c[key]).sum() / w_h.groupby(c[key]).sum()
    for col in PERSON_SHARES:
        agg[col] = (c[col] * w_p).groupby(c[key]).sum() / w_p.groupby(c[key]).sum()
    out = pd.DataFrame(agg)
    out["households_2021"] = grp.households_2021.sum()
    out["usual_residents_2021"] = grp.usual_residents_2021.sum()
    out["area_km2"] = grp.area_km2.sum()
    out["log_density"] = np.log(out.usual_residents_2021 / out.area_km2)
    out["country"] = grp.country.first()
    out["lad_cd"] = grp.lad_cd.first()
    if level == "msoa":
        for col in ["epc_ee", "epc_c_plus"]:
            out[col] = grp[col].first()
    else:  # LAD: household-weighted mean of MSOA EPC values
        for col in ["epc_ee", "epc_c_plus"]:
            out[col] = (c[col] * w_h).groupby(c[key]).sum() / w_h.groupby(c[key]).sum()
    # Rural = majority of households in rural LSOAs.
    out["rural"] = (out["rural"] > 0.5).astype(float)
    out["wales"] = (out["wales"] > 0.5).astype(float)
    return out


def unit_coords(level: str) -> pd.DataFrame:
    """Household-weighted centroids (BNG metres) of LSOA centroids."""
    g = load_lsoa_geo()[["lsoa21cd", "msoa21cd", "lad_cd", "bng_e", "bng_n", "households_2021"]]
    if level == "lsoa":
        return g.set_index("lsoa21cd")[["bng_e", "bng_n"]]
    key = {"msoa": "msoa21cd", "lad": "lad_cd"}[level]
    w = g.households_2021
    return pd.DataFrame({
        "bng_e": (g.bng_e * w).groupby(g[key]).sum() / w.groupby(g[key]).sum(),
        "bng_n": (g.bng_n * w).groupby(g[key]).sum() / w.groupby(g[key]).sum(),
    })


def geometry(level: str) -> gpd.GeoDataFrame:
    """Polygons indexed by unit code (MSOA polygons are dissolved once and cached)."""
    if level == "lsoa":
        return load_lsoa_geo().set_index("lsoa21cd")[["geometry"]]
    if level == "lad":
        return load_lad_geo().set_index("lad_cd")[["lad_nm", "geometry"]]
    path = CACHE / "msoa21_geom.parquet"
    if not path.exists():
        g = load_lsoa_geo()[["msoa21cd", "geometry"]]
        m = g.dissolve("msoa21cd")
        m["geometry"] = m.geometry.simplify(20)
        m.to_parquet(path)
    return gpd.read_parquet(path)


def display_geometry(level: str) -> gpd.GeoDataFrame:
    """Simplified polygons for maps (keeps SVGs small)."""
    path = CACHE / f"{level}_display_geom.parquet"
    if not path.exists():
        g = geometry(level).copy()
        tol = {"lsoa": 40, "msoa": 80, "lad": 150}[level]
        g["geometry"] = g.geometry.simplify(tol)
        g.to_parquet(path)
    return gpd.read_parquet(path)


def queen_weights(level: str):
    """Queen contiguity weights, islands attached to their nearest neighbour.

    MSOA and LAD contiguity is derived from LSOA contiguity (two higher units
    are neighbours if any of their LSOAs touch), which is exact and avoids
    dissolving polygons."""
    import libpysal

    path = CACHE / f"w_queen_{level}.pkl"
    if path.exists():
        with open(path, "rb") as f:
            return pickle.load(f)
    if level == "lsoa":
        g = load_lsoa_geo()[["lsoa21cd", "geometry"]].reset_index(drop=True)
        w = libpysal.weights.Queen.from_dataframe(g, ids=g.lsoa21cd.tolist(), use_index=False)
        if w.islands:
            w = libpysal.weights.util.attach_islands(w, libpysal.weights.KNN.from_dataframe(g, k=1, ids=g.lsoa21cd.tolist()))
    else:
        wl = queen_weights("lsoa")
        lk = load_lsoa_geo().set_index("lsoa21cd")
        key = {"msoa": "msoa21cd", "lad": "lad_cd"}[level]
        parent = lk[key].to_dict()
        neigh: dict = {u: set() for u in sorted(set(parent.values()))}
        for i, js in wl.neighbors.items():
            pi = parent[i]
            for j in js:
                pj = parent[j]
                if pi != pj:
                    neigh[pi].add(pj)
        ids = sorted(neigh)
        w = libpysal.weights.W({k: sorted(v) for k, v in neigh.items()}, id_order=ids, silence_warnings=True)
        if w.islands:
            c = unit_coords(level).loc[ids]
            knn = libpysal.weights.KNN(c.values, k=1, ids=ids)
            w = libpysal.weights.util.attach_islands(w, knn)
    w.transform = "r"
    with open(path, "wb") as f:
        pickle.dump(w, f)
    return w


def subset_weights(w, ids, coords: pd.DataFrame):
    """Row-standardised weights restricted to `ids` (order preserved);
    units left without neighbours get their nearest remaining neighbour."""
    import libpysal
    from scipy.spatial import cKDTree

    ids = list(ids)
    keep = set(ids)
    neigh = {i: [j for j in w.neighbors[i] if j in keep] for i in ids}
    lonely = [i for i, js in neigh.items() if not js]
    if lonely:
        xy = coords.loc[ids, ["bng_e", "bng_n"]].values
        tree = cKDTree(xy)
        _, nn = tree.query(coords.loc[lonely, ["bng_e", "bng_n"]].values, k=2)
        for i, (a, b) in zip(lonely, nn):
            j = ids[b] if ids[a] == i else ids[a]
            neigh[i] = [j]
            neigh[j] = sorted(set(neigh[j]) | {i})
    ws = libpysal.weights.W(neigh, id_order=ids, silence_warnings=True)
    ws.transform = "r"
    return ws
