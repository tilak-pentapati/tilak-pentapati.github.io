"""Shared paths, data loaders, spatial weights and figure style for the
descriptive and dynamics analysis (scripts prefixed ``a_``).

Run every script with the project env in isolated mode, e.g.

    /Users/jay/anaconda3/envs/envs563/bin/python -I \
        projects/python/energy_inequality/analysis/run_dynamics.py

`-I` stops Python adding the script folder to sys.path, so each script adds it
explicitly before importing this module.
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
DATA = REPO / "data" / "data_python" / "energy_inequality"
PROCESSED = DATA / "processed"
RAW = DATA / "raw"
RESULTS = DATA / "results" / "dynamics"
CACHE = RESULTS / "cache"
FIG = RESULTS / "figures"
TAB = RESULTS / "tables"
for _d in (RESULTS, CACHE, FIG, TAB):
    _d.mkdir(parents=True, exist_ok=True)

YEARS = list(range(2010, 2025))
CORE_YEARS = list(range(2015, 2025))   # consistent 2021-LSOA geography, new DESNZ system
SUBPERIODS = {"2010-14": (2010, 2014), "2015-19": (2015, 2019), "2020-24": (2020, 2024)}
SEED = 20260517
PERMS = 999
K_CLASSES = 5
CLASS_LABELS = ["Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"]

# --------------------------------------------------------------------------- palette
# Validated reference palette (dataviz skill): categorical slots 1-3 are safe for
# all-pairs use; blue sequential ramp; blue <-> red diverging with grey midpoint.
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]          # 5-class blue ramp
SEQ7 = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
RED = "#e34948"
BLUE = "#2a78d6"
LISA_COLORS = {  # HH, LH, LL, HL, ns (giddy/esda quadrant codes 1..4)
    "HH": "#c23b3a", "HL": "#f2a3a2", "LL": "#1c5cab", "LH": "#9ec5f4", "ns": "#e1e0d9",
}
SHOCK_YEARS = {2020: "COVID-19", 2022: "Price crisis"}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def set_style() -> None:
    """Matplotlib style: system sans, thin marks, recessive axes; sized so a
    figure exported at 7.5in x 300dpi stays legible when shrunk to 375px."""
    import matplotlib as mpl

    mpl.rcParams.update({
        # a family list enables per-glyph fallback (arrows, Greek) to DejaVu Sans
        "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"],
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 12,
        "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "figure.dpi": 100,
    })


def savefig(fig, name: str) -> list[str]:
    """Write PNG (300 dpi) and SVG; return the relative file names."""
    out = []
    for ext, kw in (("png", {"dpi": 300}), ("svg", {})):
        p = FIG / f"{name}.{ext}"
        fig.savefig(p, bbox_inches="tight", **kw)
        out.append(str(p.relative_to(RESULTS)))
    import matplotlib.pyplot as plt
    plt.close(fig)
    return out


def write_table(df: pd.DataFrame, name: str, csv: bool = True) -> None:
    df.to_parquet(TAB / f"{name}.parquet", index=False)
    if csv:
        df.to_csv(TAB / f"{name}.csv", index=False, float_format="%.6g")


def save_json(obj, name: str) -> None:
    (RESULTS / f"{name}.json").write_text(json.dumps(obj, indent=2, default=float))


# --------------------------------------------------------------------------- data
def load_lsoa_panel(fuel: str = "elec") -> pd.DataFrame:
    p = pd.read_parquet(PROCESSED / f"{fuel}_lsoa_panel.parquet")
    return p


def lad_region() -> pd.DataFrame:
    """LAD (May 2026) -> region (England's nine regions, Wales as one) from the
    DESNZ subnational table, which lists every LA under its region."""
    cache = CACHE / "lad_region.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    import geopandas as gpd

    t = pd.read_excel(RAW / "desnz" / "Subnational_electricity_consumption_statistics_2005-2024.xlsx",
                      "2024", header=4, usecols=[0, 1, 2])
    t.columns = ["code", "region", "la"]
    t = t[t["code"].astype(str).str.match(r"^[EW]0[6-9]")]
    t["region"] = t["region"].where(~t["code"].str.startswith("W"), "Wales")
    t["region"] = t["region"].replace({"Inner London": "London", "Outer London": "London"})
    lads = gpd.read_parquet(PROCESSED / "lad21_ew.parquet")[["lad_cd", "lad_nm", "country"]]
    out = lads.merge(t[["code", "region"]], left_on="lad_cd", right_on="code", how="left").drop(columns="code")
    miss = out["region"].isna()
    if miss.any():
        # LADs re-coded after the 2024 table: take the modal region of their
        # LSOAs' 2022 LADs.
        g = gpd.read_parquet(PROCESSED / "lsoa21_ew.parquet", columns=["lad_cd", "lad22cd"])
        g = g.merge(t[["code", "region"]], left_on="lad22cd", right_on="code", how="left")
        fill = g.groupby("lad_cd")["region"].agg(lambda s: s.mode().iat[0] if s.notna().any() else np.nan)
        out.loc[miss, "region"] = out.loc[miss, "lad_cd"].map(fill)
    assert out["region"].notna().all(), out[out["region"].isna()]
    out = pd.DataFrame(out)
    out.to_parquet(cache, index=False)
    return out


def wide(panel: pd.DataFrame, col: str = "mean_kwh", years=YEARS, ids=None) -> pd.DataFrame:
    w = panel.pivot(index="lsoa21cd", columns="year", values=col)[list(years)]
    if ids is not None:
        w = w.loc[ids]
    return w


def analysis_frame(fuel: str = "elec") -> dict:
    """The balanced analysis set and its metadata (cached).

    * ``ids``: LSOAs ordered as in the geometry file with a published value in
      every year 2010-2024 (the balanced panel used everywhere).
    * ``ever_flag``: LSOAs flagged as outliers in at least one year.
    * ``change``: ONS 2011->2021 change class (U/S/M/X).
    """
    import geopandas as gpd

    geo = gpd.read_parquet(PROCESSED / "lsoa21_ew.parquet")
    p = load_lsoa_panel(fuel)
    y = wide(p, "mean_kwh").reindex(geo["lsoa21cd"])
    complete = y.notna().all(axis=1)
    flag = p.groupby("lsoa21cd")["flag_outlier"].any().reindex(geo["lsoa21cd"])
    change = p[p["year"] == 2015].set_index("lsoa21cd")["lsoa21_change"].reindex(geo["lsoa21cd"])
    meters = wide(p, "meters").reindex(geo["lsoa21cd"])
    kwh = wide(p, "total_kwh").reindex(geo["lsoa21cd"])
    return {"geo": geo, "y": y, "complete": complete.values, "ever_flag": flag.fillna(False).values,
            "change": change.values, "meters": meters, "kwh": kwh}


# --------------------------------------------------------------------------- weights
def lsoa_weights(kind: str = "queen", ids_mask: np.ndarray | None = None, tag: str = "all"):
    """Row-standardised weights for the LSOAs in ``ids_mask`` (geometry order).

    Queen contiguity is built once on all 35,672 polygons and cached as a GAL
    adjacency list; islands (no contiguous neighbour) are joined to their
    nearest neighbour by centroid distance (libpysal ``attach_islands``) so
    every unit has a spatial lag. KNN uses k=8 on population-weighted
    centroids. Subsetting keeps neighbour relations among retained units, and
    any unit left without neighbours is again attached to its nearest
    retained unit.
    """
    import geopandas as gpd
    from libpysal import weights

    geo = gpd.read_parquet(PROCESSED / "lsoa21_ew.parquet", columns=["lsoa21cd", "bng_e", "bng_n", "geometry"])
    if ids_mask is None:
        ids_mask = np.ones(len(geo), bool)
    cache = CACHE / f"w_{kind}_{tag}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        w = _w_from_edges(z["i"], z["j"], int(z["n"]))
        w.transform = "r"
        return w
    sub = geo[ids_mask].reset_index(drop=True)
    xy = np.c_[sub["bng_e"].values, sub["bng_n"].values].astype(float)
    if kind == "queen":
        full_cache = CACHE / "w_queen_full_edges.npz"
        if full_cache.exists():
            z = np.load(full_cache)
            wf = _w_from_edges(z["i"], z["j"], int(z["n"]))
        else:
            t0 = time.time()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wf = weights.Queen.from_dataframe(geo, use_index=False, silence_warnings=True)
            i, j = _edges(wf)
            np.savez(full_cache, i=i, j=j, n=wf.n)
            log(f"queen contiguity built in {time.time() - t0:.0f}s: n={wf.n}, islands={len(wf.islands)}, "
                f"components={wf.n_components}")
        keep = np.flatnonzero(ids_mask)
        z = np.load(full_cache) if full_cache.exists() else None
        i, j = (z["i"], z["j"]) if z is not None else _edges(wf)
        if keep.size < len(geo):
            pos = np.full(len(geo), -1)
            pos[keep] = np.arange(keep.size)
            ok = (pos[i] >= 0) & (pos[j] >= 0)
            w = _w_from_edges(pos[i[ok]], pos[j[ok]], keep.size)
        else:
            w = wf
        if w.islands:
            knn1 = weights.KNN(xy, k=1)
            w = weights.attach_islands(w, knn1)
    elif kind.startswith("knn"):
        k = int(kind[3:])
        w = weights.KNN(xy, k=k)
    else:
        raise ValueError(kind)
    i, j = _edges(w)
    np.savez(cache, i=i, j=j, n=w.n)
    w = _w_from_edges(i, j, w.n)
    w.transform = "r"
    return w


def _edges(w):
    i, j = [], []
    for a, nb in w.neighbors.items():
        i.extend([a] * len(nb))
        j.extend(nb)
    return np.asarray(i, dtype=np.int64), np.asarray(j, dtype=np.int64)


def _w_from_edges(i, j, n):
    from libpysal import weights

    neighbors = {k: [] for k in range(n)}
    for a, b in zip(i.tolist(), j.tolist()):
        neighbors[a].append(b)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return weights.W(neighbors, silence_warnings=True)


def lad_weights():
    """Queen contiguity on the 318 LADs, islands attached to nearest LAD."""
    import geopandas as gpd
    from libpysal import weights

    lad = gpd.read_parquet(PROCESSED / "lad21_ew.parquet").sort_values("lad_cd").reset_index(drop=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = weights.Queen.from_dataframe(lad, use_index=False, silence_warnings=True)
        islands = [lad.loc[i, "lad_nm"] for i in w.islands]
        if w.islands:
            c = lad.geometry.representative_point()
            knn1 = weights.KNN(np.c_[c.x, c.y], k=1)
            w = weights.attach_islands(w, knn1)
    w.transform = "r"
    return lad, w, islands


def edge_arrays(w):
    """Directed neighbour pairs (i, j) and row-standardised weights."""
    i, j = _edges(w)
    return i, j


def relative(y: pd.DataFrame, meters: pd.DataFrame | None = None) -> pd.DataFrame:
    """Value / E&W mean per year. With ``meters`` the denominator is the
    meter-weighted E&W mean (total kWh / total meters) of the units supplied."""
    if meters is None:
        return y / y.mean(axis=0)
    return y / ((y * meters).sum(axis=0) / meters.sum(axis=0))


# --------------------------------------------------------------------------- core set
def core(fuel: str = "elec", exclude_flagged: bool = False, exclude_sx: bool = False) -> dict:
    """Balanced panel arrays for one analysis variant.

    Returns geometry-ordered ``ids`` (mask over the 35,672 LSOAs), the level
    (``y``, kWh/meter), relative consumption (``rel``: value / meter-weighted
    E&W mean of the retained LSOAs in that year), meters, LAD codes and the
    matching row-standardised queen weights.
    """
    A = analysis_frame(fuel)
    mask = A["complete"].copy()
    tag = f"{fuel}"
    if exclude_flagged:
        mask &= ~A["ever_flag"]
        tag += "_noflag"
    if exclude_sx:
        mask &= ~np.isin(A["change"], ["S", "X"])
        tag += "_noSX"
    y = A["y"][mask]
    meters = A["meters"][mask]
    rel = relative(y, meters)
    geo = A["geo"][mask].reset_index(drop=True)
    return {"mask": mask, "tag": tag, "y": y, "rel": rel, "meters": meters, "geo": geo,
            "lad": geo["lad_cd"].values, "change": A["change"][mask],
            "flag": A["ever_flag"][mask]}


# --------------------------------------------------------------------------- maps
def map_layers():
    """Simplified LSOA polygons (20 m tolerance, for drawing only) and region
    outlines, cached."""
    import geopandas as gpd

    c1, c2 = CACHE / "lsoa_simpl.parquet", CACHE / "region_outline.parquet"
    if c1.exists() and c2.exists():
        return gpd.read_parquet(c1), gpd.read_parquet(c2)
    g = gpd.read_parquet(PROCESSED / "lsoa21_ew.parquet", columns=["lsoa21cd", "geometry"])
    g["geometry"] = g.geometry.simplify(20, preserve_topology=True)
    g.to_parquet(c1)
    lad = gpd.read_parquet(PROCESSED / "lad21_ew.parquet").merge(lad_region()[["lad_cd", "region"]], on="lad_cd")
    reg = lad.dissolve("region").reset_index()[["region", "geometry"]]
    reg["geometry"] = reg.geometry.simplify(200)
    reg.to_parquet(c2)
    return g, reg


LONDON_BOX = (503000, 155000, 562000, 201000)


def draw_map(ax, gdf, colors, regions, title=None, inset=True):
    """Fill LSOAs with per-row colours; polygons rasterised so SVGs stay small
    while text remains vector. Adds a London inset."""
    gdf.plot(ax=ax, color=colors, linewidth=0, rasterized=True)
    regions.boundary.plot(ax=ax, color="white", linewidth=0.5)
    ax.set_axis_off()
    ax.set_aspect("equal")
    if title:
        ax.set_title(title)
    if inset:
        x0, y0, x1, y1 = LONDON_BOX
        ia = ax.inset_axes([0.62, 0.70, 0.38, 0.30])
        sub = gdf.cx[x0:x1, y0:y1]
        sub.plot(ax=ia, color=np.asarray(colors, dtype=object)[gdf.index.get_indexer(sub.index)].tolist(),
                 linewidth=0, rasterized=True)
        regions.boundary.plot(ax=ia, color="white", linewidth=0.6)
        ia.set_xlim(x0, x1)
        ia.set_ylim(y0, y1)
        ia.set_xticks([])
        ia.set_yticks([])
        ia.set_aspect("equal")
        for s in ia.spines.values():
            s.set_visible(True)
            s.set_color(MUTED)
            s.set_linewidth(0.6)
        ia.set_title("London", fontsize=10, loc="left", fontweight="normal", color=INK2, pad=2)
        ax.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle(
            (x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor=MUTED, linewidth=0.6))
    return ax
