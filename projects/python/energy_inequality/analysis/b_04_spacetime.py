"""Section 4: spatio-temporal models on the LSOA panel, 2010-2024.

GTWR is not used. With Census-2021-only covariates, temporal variation in
local coefficients is not identified from the covariates (it would only
smooth the outcome), and the genuinely time-varying covariates available
at LSOA level (meter counts; ONS population, age structure, people per
meter, 2011-2024) carry little within-area variation, so a local panel
estimator at 35k units x 15 years would be both very costly and noisy.
Two identified alternatives replace it:

A. Two-way fixed-effects panel (LSOA and year effects) with time-varying
   covariates, and an event-study extension in which year dummies are
   interacted with standardised 2021 baselines. The interactions are
   identified (time-invariant x year), and trace how the gap associated
   with each baseline evolved through the COVID-19 and price-crisis years.
   Reference year 2019; standard errors clustered by local authority.

B. A Fixed-Rank-Kriging-style trend surface: relative log kWh per meter
   (year mean removed) is projected onto a tensor product of three
   resolutions of bisquare spatial basis functions and piecewise-linear
   temporal basis functions, fitted by ridge regression with the penalty
   chosen by local-authority-blocked cross-validation. The fitted surface
   is the large-scale (>~25 km) spatio-temporal trend; residuals are
   fine-scale variation. This is a descriptive decomposition, not a causal
   model, and a Python re-implementation of the FRK idea rather than the R
   package.

LSOAs with split/irregular boundary changes (S/X; 2010-14 values
apportioned) and LSOAs with fewer than 50 meters in any year are excluded.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

import b_common as bc
import b_models as bmod

EVENT_Z = ["share_ch_electric", "share_ch_renewable", "share_ch_offgas", "dep19_pct", "share_detached",
           "mean_hh_size", "share_age_0_14", "share_private_rent"]
REF_YEAR = 2019


# --------------------------------------------------------------------------
# Panel assembly
# --------------------------------------------------------------------------
def panel() -> pd.DataFrame:
    def _build():
        e = bc.load_elec_lsoa()
        cov = bc.load_cov_lsoa().set_index("lsoa21cd")
        bad = set(e.loc[e.flag_small_n | e.flag_missing, "lsoa21cd"])
        bad |= set(cov.index[cov.lsoa21_change.isin(["S", "X"])])
        e = e[~e.lsoa21cd.isin(bad)].copy()
        e["log_kwh"] = np.log(e.mean_kwh)
        e["log_meters"] = np.log(e.meters)
        pop = pd.read_parquet(bc.PROCESSED / "population_lsoa_panel.parquet")
        e = e.merge(pop, on=["lsoa21cd", "year"], how="left")
        e["log_pop_per_meter"] = np.log(e.pop_total / e.meters)
        e["share_0_15"] = e.pop_0_15 / e.pop_total
        e["share_65plus_t"] = e.pop_65plus / e.pop_total
        return e[["lsoa21cd", "lad_cd", "year", "log_kwh", "log_meters", "log_pop_per_meter", "share_0_15",
                  "share_65plus_t", "flag_outlier"]].reset_index(drop=True)
    return bc.cached("st_panel", _build)


def _demean_twoway(df, cols, unit="lsoa21cd", time="year"):
    """Two-way within transformation for a balanced panel."""
    out = df[cols] - df.groupby(unit)[cols].transform("mean") - df.groupby(time)[cols].transform("mean") \
        + df[cols].mean()
    return out


def _cluster_ols(y, X, groups):
    """OLS with cluster-robust (CR1) covariance."""
    XtX_inv = np.linalg.inv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    g = pd.factorize(groups)[0]
    G = g.max() + 1
    S = np.zeros((G, X.shape[1]))
    np.add.at(S, g, X * e[:, None])
    n, k = X.shape
    adj = G / (G - 1) * (n - 1) / (n - k)
    V = adj * XtX_inv @ (S.T @ S) @ XtX_inv
    return b, np.sqrt(np.diag(V)), e


# --------------------------------------------------------------------------
# A. Two-way FE and event study
# --------------------------------------------------------------------------
def fe_models():
    p = panel()
    rows = []
    # A1: time-varying covariates, 2011-2024 (population starts in 2011).
    d = p[p.year >= 2011].dropna(subset=["log_pop_per_meter", "share_0_15", "share_65plus_t"]).copy()
    cnt = d.groupby("lsoa21cd").year.transform("count")
    d = d[cnt == cnt.max()]
    tv = ["log_meters", "log_pop_per_meter", "share_0_15", "share_65plus_t"]
    Z = _demean_twoway(d, ["log_kwh"] + tv)
    b, se, _ = _cluster_ols(Z["log_kwh"].values, Z[tv].values, d.lad_cd.values)
    ss_tot = (Z["log_kwh"] ** 2).sum()
    resid = Z["log_kwh"].values - Z[tv].values @ b
    within_r2 = 1 - (resid ** 2).sum() / ss_tot
    for t, bb, s in zip(tv, b, se):
        rows.append({"model": "twoway_fe_timevarying", "term": t, "coef": bb, "se": s, "n_obs": len(d),
                     "n_lsoa": d.lsoa21cd.nunique(), "within_r2": within_r2})
    # Within-LSOA variation available in each covariate (SD of two-way demeaned values)
    var_rows = [{"term": t, "sd_within": Z[t].std(), "sd_total": d[t].std(),
                 "share_within": Z[t].var() / d[t].var()} for t in tv + ["log_kwh"]]
    bc.save_table(pd.DataFrame(var_rows), "st_within_variation")
    coef = pd.DataFrame(rows)
    bc.save_table(coef, "st_fe_timevarying")
    return coef


def event_study():
    p = panel()
    cov = bc.covariates("lsoa")
    zs = (cov[EVENT_Z] - cov[EVENT_Z].mean()) / cov[EVENT_Z].std()
    d = p.join(zs, on="lsoa21cd").dropna(subset=EVENT_Z).copy()
    years = [y for y in sorted(d.year.unique()) if y != REF_YEAR]
    inter = {}
    for z in EVENT_Z:
        for y in years:
            inter[f"{z}_x_{y}"] = d[z].values * (d.year.values == y)
    I = pd.DataFrame(inter, index=d.index)
    cols = list(I.columns) + ["log_meters"]
    D = pd.concat([d[["lsoa21cd", "year", "log_kwh", "log_meters"]], I], axis=1)
    Z = _demean_twoway(D, ["log_kwh"] + cols)
    b, se, _ = _cluster_ols(Z["log_kwh"].values, Z[cols].values, d.lad_cd.values)
    res = pd.DataFrame({"name": cols, "coef": b, "se": se})
    out = []
    for z in EVENT_Z:
        out.append({"baseline": z, "year": REF_YEAR, "coef": 0.0, "se": 0.0})
        for y in years:
            r = res[res.name == f"{z}_x_{y}"].iloc[0]
            out.append({"baseline": z, "year": y, "coef": r.coef, "se": r.se})
    t = pd.DataFrame(out).sort_values(["baseline", "year"])
    t["n_lsoa"] = d.lsoa21cd.nunique()
    t["log_meters_coef"] = res[res.name == "log_meters"].coef.iloc[0]
    bc.save_table(t, "st_event_study")
    return t


def event_figure(t):
    sel = ["share_ch_electric", "share_ch_renewable", "dep19_pct", "share_detached"]
    fig, axes = plt.subplots(2, 2, figsize=(5.8, 5.6), sharex=True, sharey=True)
    for ax, z in zip(axes.ravel(), sel):
        d = t[t.baseline == z]
        ax.fill_between(d.year, 100 * (d.coef - 1.96 * d.se), 100 * (d.coef + 1.96 * d.se), color=bc.CAT3[0],
                        alpha=0.2, lw=0)
        ax.plot(d.year, 100 * d.coef, color=bc.CAT3[0], lw=2, marker="o", ms=3.5)
        ax.axhline(0, color=bc.INK2, lw=0.8)
        ax.axvline(2014.5, color=bc.INK2, lw=0.6, ls=":")
        ax.set_title(bmod.nice(z).replace(" (incl. heat pumps)", "\n(incl. heat pumps)"), fontsize=12, loc="left")
        ax.set_xticks([2010, 2015, 2019, 2024])
    for ax in axes[:, 0]:
        ax.set_ylabel("% vs 2019, per SD")
    fig.suptitle("How baseline gaps in kWh per meter evolved\n(LSOA two-way fixed effects, ref. 2019)",
                 x=0.0, ha="left", fontsize=13)
    fig.text(0.0, -0.04, "Year x standardised 2021 baseline interactions; 95% CIs clustered by local authority. "
                         "Dotted line: 2014/15 processing break.", fontsize=10, color=bc.INK2, wrap=True)
    fig.tight_layout()
    return bc.save_fig(fig, "st_event_study")


# --------------------------------------------------------------------------
# B. Basis-function trend surface
# --------------------------------------------------------------------------
SPACINGS = [100_000, 50_000, 25_000]
KNOTS = [2010, 2012, 2014, 2015, 2017, 2019, 2020, 2021, 2022, 2024]


def spatial_basis(xy: np.ndarray):
    """Bisquare basis functions on three regular grids (radius 1.5 x spacing),
    keeping centres within one radius of at least one LSOA."""
    from scipy.spatial import cKDTree

    tree = cKDTree(xy)
    blocks, meta = [], []
    for res, sp in enumerate(SPACINGS, start=1):
        r = 1.5 * sp
        xs = np.arange(xy[:, 0].min() - sp / 2, xy[:, 0].max() + sp, sp)
        ys = np.arange(xy[:, 1].min() - sp / 2, xy[:, 1].max() + sp, sp)
        C = np.array([(x, y) for x in xs for y in ys])
        near = np.array([len(tree.query_ball_point(c, r * 0.5)) > 0 for c in C])
        C = C[near]
        ct = cKDTree(C)
        pairs = ct.sparse_distance_matrix(tree, r, output_type="coo_matrix")
        d = pairs.data / r
        vals = (1 - d ** 2) ** 2
        B = sparse.csr_matrix((vals, (pairs.col, pairs.row)), shape=(len(xy), len(C)))
        blocks.append(B)
        meta += [{"res": res, "spacing_km": sp / 1000, "x": c[0], "y": c[1]} for c in C]
    return sparse.hstack(blocks).tocsr(), pd.DataFrame(meta)


def temporal_basis(years: np.ndarray):
    k = np.array(KNOTS, dtype=float)
    T = np.zeros((len(years), len(k)))
    for j in range(len(k)):
        left = k[j - 1] if j > 0 else None
        right = k[j + 1] if j < len(k) - 1 else None
        y = years.astype(float)
        v = np.zeros_like(y)
        if left is not None:
            m = (y >= left) & (y <= k[j])
            v[m] = (y[m] - left) / (k[j] - left)
        if right is not None:
            m = (y >= k[j]) & (y <= right)
            v[m] = np.maximum(v[m], (right - y[m]) / (right - k[j]))
        v[y == k[j]] = 1.0
        T[:, j] = v
    return T


def trend_surface():
    def _fit():
        p = panel()
        lsoa_ids = np.array(sorted(p.lsoa21cd.unique()))
        coords = bc.unit_coords("lsoa").loc[lsoa_ids][["bng_e", "bng_n"]].values
        Bs, meta = spatial_basis(coords)
        idx = pd.Series(np.arange(len(lsoa_ids)), index=lsoa_ids)
        p = p.copy()
        p["rel"] = p.log_kwh - p.groupby("year").log_kwh.transform("mean")
        Tt = temporal_basis(np.arange(2010, 2025))
        ii = idx.loc[p.lsoa21cd].values
        tt = (p.year - 2010).values
        # Tensor design: row i,t = Bs[i] (x) Tt[t]
        Bi = Bs[ii]
        Tm = sparse.csr_matrix(Tt[tt])
        X = sparse.hstack([Bi.multiply(Tm[:, [j]]) for j in range(Tt.shape[1])]).tocsr()
        y = p.rel.values
        # LAD-blocked 5-fold CV over the ridge penalty
        lads = p.lad_cd.values
        ulad = np.unique(lads)
        rng = np.random.default_rng(bc.SEED)
        fold_of = dict(zip(ulad, rng.integers(0, 5, len(ulad))))
        fold = np.array([fold_of[a] for a in lads])
        XtX = (X.T @ X).toarray()
        Xty = X.T @ y
        XtX_f = [(X[fold == f].T @ X[fold == f]).toarray() for f in range(5)]
        Xty_f = [X[fold == f].T @ y[fold == f] for f in range(5)]
        lambdas = [0.1, 1, 10, 100, 1000]
        cv = []
        for lam in lambdas:
            sse = 0.0
            for f in range(5):
                A = XtX - XtX_f[f] + lam * np.eye(XtX.shape[0])
                beta = np.linalg.solve(A, Xty - Xty_f[f])
                e = y[fold == f] - X[fold == f] @ beta
                sse += (e ** 2).sum()
            cv.append({"lambda": lam, "cv_mse": sse / len(y)})
            bc.log(f"trend surface: lambda={lam} cv_mse={sse / len(y):.6f}")
        cv = pd.DataFrame(cv)
        lam = float(cv.loc[cv.cv_mse.idxmin(), "lambda"])
        beta = np.linalg.solve(XtX + lam * np.eye(XtX.shape[0]), Xty)
        fitted = X @ beta
        p["trend"] = fitted
        p["fine"] = p.rel - p.trend
        return {"p": p[["lsoa21cd", "lad_cd", "year", "rel", "trend", "fine"]], "beta": beta,
                "meta": meta, "cv": cv, "lambda": lam, "n_basis_space": Bs.shape[1],
                "n_basis_time": Tt.shape[1]}
    return bc.cached("st_trend_surface", _fit)


def surface_outputs(ts):
    p = ts["p"]
    bc.save_table(ts["cv"], "st_trend_surface_cv")
    # Variance decomposition by year: large-scale trend vs fine-scale, and the
    # persistent (LSOA mean) part of the fine-scale component.
    p = p.copy()
    p["fine_persistent"] = p.groupby("lsoa21cd").fine.transform("mean")
    rows = []
    for y, g in p.groupby("year"):
        v = g.rel.var()
        rows.append({"year": y, "var_rel": v, "share_trend": g.trend.var() / v, "share_fine": g.fine.var() / v,
                     "cov_term": 2 * np.cov(g.trend, g.fine)[0, 1] / v,
                     "share_fine_persistent": g.fine_persistent.var() / v, "sd_trend": g.trend.std(),
                     "sd_rel": g.rel.std()})
    vd = pd.DataFrame(rows)
    bc.save_table(vd, "st_trend_surface_variance")
    meta = {"lambda": ts["lambda"], "n_basis_space": ts["n_basis_space"], "n_basis_time": ts["n_basis_time"],
            "n_coef": len(ts["beta"]), "n_obs": len(p), "n_lsoa": p.lsoa21cd.nunique(),
            "basis_per_resolution": ts["meta"].groupby("res").size().to_dict(),
            "spacings_km": [s / 1000 for s in SPACINGS], "knots": KNOTS}
    bc.save_json(meta, "st_trend_surface_meta")
    # MSOA-averaged surface for maps
    lk = bc.load_cov_lsoa().set_index("lsoa21cd").msoa21cd
    p["msoa21cd"] = p.lsoa21cd.map(lk)
    m = p.groupby(["msoa21cd", "year"]).trend.mean().unstack()
    m.columns = [str(c) for c in m.columns]
    m.reset_index().to_parquet(bc.TAB / "st_trend_surface_msoa.parquet")
    return vd, m


def surface_maps(m):
    import matplotlib as mpl

    g = bc.display_geometry("msoa")
    yrs = ["2010", "2019", "2022", "2024"]
    lim = np.nanquantile(np.abs(m[yrs].values), 0.99)
    norm = mpl.colors.TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    fig, axes = plt.subplots(2, 2, figsize=(5.8, 6.8))
    for ax, y in zip(axes.ravel(), yrs):
        gg = g.join(m[y].rename("v"), how="inner")
        gg.plot(ax=ax, column="v", cmap=bc.DIVERGING, norm=norm, linewidth=0)
        ax.set_axis_off()
        ax.set_title(y, fontsize=13, loc="left")
    cax = fig.add_axes([0.1, 0.06, 0.8, 0.02])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=bc.DIVERGING, norm=norm), cax=cax, orientation="horizontal")
    cb.set_label("Large-scale trend, log kWh vs E&W year mean", fontsize=11)
    cb.set_ticks([t for t in (-0.2, -0.1, 0.0, 0.1, 0.2) if abs(t) <= lim])
    cb.outline.set_visible(False)
    fig.suptitle("Spatio-temporal trend surface of relative consumption", x=0.02, ha="left", fontsize=13)
    fig.subplots_adjust(left=0.0, right=1.0, top=0.89, bottom=0.12, wspace=0.02, hspace=0.10)
    bc.save_fig(fig, "st_trend_surface_years")
    # Change surfaces
    for a, b, tag in [("2019", "2022", "shock"), ("2010", "2024", "long")]:
        d = (m[b] - m[a]).rename("v")
        gg = g.join(d, how="inner")
        import b_maps as bm
        bm.choropleth(gg, "v", f"Change in large-scale relative consumption\n{a} to {b}",
                      f"st_trend_surface_change_{tag}", cmap=bc.DIVERGING, norm=bm.symmetric_norm(d.values, 0.99),
                      cbar_label="Δ log kWh vs E&W mean", boundaries=bc.display_geometry("lad"),
                      note="Basis-function trend surface (3 resolutions: 100, 50, 25 km), MSOA averages.")


def main(force=False):
    bc.set_style()
    fe = fe_models()
    bc.log(f"\n{fe.to_string()}")
    ev = bc.cached("st_event_study_res", event_study, force)
    bc.save_table(ev, "st_event_study")
    event_figure(ev)
    ts = trend_surface()
    vd, m = surface_outputs(ts)
    bc.log(f"\n{vd.to_string()}")
    surface_maps(m)


if __name__ == "__main__":
    main(force="--force" in sys.argv)
