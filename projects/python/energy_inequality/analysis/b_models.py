"""Model library for the explanatory sections (cross-section and change):
design assembly, VIF screening, the global OLS / spatial-econometric suite
(spreg), GWR/MGWR (mgwr) with corrected local inference, and figures.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import b_common as bc

NICE = {
    "share_ch_electric": "Electric heating",
    "share_ch_offgas": "Oil/LPG/solid-fuel heating",
    "share_ch_renewable": "Renewable heating (incl. heat pumps)",
    "share_ch_district": "District/communal heating",
    "share_detached": "Detached houses",
    "share_flat": "Flats",
    "mean_rooms": "Mean rooms",
    "mean_hh_size": "Household size",
    "share_1person_hh": "One-person households",
    "share_social_rent": "Social rent",
    "share_private_rent": "Private rent",
    "share_age_0_14": "Aged 0-14",
    "share_age_65plus": "Aged 65+",
    "share_retired": "Retired",
    "dep19_pct": "Less deprived (IMD 2019 pct)",
    "dep25_pct": "Less deprived (IMD 2025 pct)",
    "epc_ee": "EPC efficiency (MSOA median)",
    "epc_c_plus": "EPC band C+ share (MSOA)",
    "log_density": "Log population density",
    "rural": "Rural",
    "wales": "Wales",
    "log_level_base": "Initial log kWh (2010-12)",
    "log_level_2019": "Initial log kWh (2019)",
    "dlog_meters": "Meter growth (log change)",
    "const": "Intercept",
}

CANDIDATES = ["share_ch_electric", "share_ch_offgas", "share_ch_renewable", "share_ch_district",
              "share_detached", "share_flat", "mean_rooms", "mean_hh_size", "share_1person_hh",
              "share_social_rent", "share_private_rent", "share_age_0_14", "share_age_65plus",
              "dep19_pct", "epc_ee", "log_density", "rural", "wales"]


# Covariates for the MSOA models (global and local alike, so AICc compares like
# with like). Dummies (rural, Wales) are excluded because they are constant
# inside most local kernels and make local designs singular; district heating
# is near-zero in most kernels; child share is dropped to keep MGWR tractable
# (it overlaps with household size at MSOA level).
XS_LOCAL = ["share_ch_electric", "share_ch_offgas", "share_ch_renewable", "share_detached", "share_flat",
            "mean_rooms", "mean_hh_size", "share_social_rent", "share_private_rent", "share_age_65plus",
            "dep19_pct", "epc_ee", "log_density"]


def nice(v: str) -> str:
    return NICE.get(v, v)


# --------------------------------------------------------------------------
# Design and screening
# --------------------------------------------------------------------------
def vif_table(X: pd.DataFrame) -> pd.Series:
    """Variance inflation factors (with intercept)."""
    Z = np.column_stack([np.ones(len(X)), X.values])
    out = {}
    for k, c in enumerate(X.columns, start=1):
        y = Z[:, k]
        others = np.delete(Z, k, axis=1)
        beta, *_ = np.linalg.lstsq(others, y, rcond=None)
        r2 = 1 - ((y - others @ beta) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        out[c] = 1.0 / max(1 - r2, 1e-12)
    return pd.Series(out)


def vif_screen(X: pd.DataFrame, threshold: float = 5.0, protect=()) -> tuple[list, pd.DataFrame]:
    """Drop the highest-VIF covariate until every VIF < threshold. Returns the
    kept list and a log of each step."""
    keep = list(X.columns)
    steps = []
    while True:
        v = vif_table(X[keep])
        cand = v.drop([p for p in protect if p in v.index])
        steps.append({"step": len(steps), "max_vif": v.max(), "max_var": v.idxmax(),
                      **{f"vif_{k}": val for k, val in v.items()}})
        if cand.max() < threshold:
            break
        keep.remove(cand.idxmax())
    return keep, pd.DataFrame(steps)


def standardise(X: pd.DataFrame) -> pd.DataFrame:
    return (X - X.mean()) / X.std()


# --------------------------------------------------------------------------
# Global models
# --------------------------------------------------------------------------
def moran_resid(e, w):
    import esda

    m = esda.Moran(np.asarray(e).ravel(), w, permutations=0)
    return m.I, m.z_norm, m.p_norm


def fit_global_suite(y: pd.Series, X: pd.DataFrame, w, name: str, ml: bool = True) -> dict:
    """OLS (White SEs + spatial diagnostics), SLX, SAR, SEM, SDM and SDEM.

    ML estimators use the sparse LU log-determinant, feasible at n = 35k.
    Returns dict with 'coefs' (long table), 'fit' (model comparison) and
    'diag' (LM tests) DataFrames, plus OLS residuals."""
    import spreg

    yv = y.values.reshape(-1, 1)
    Xv = X.values
    names = list(X.columns)
    coefs, fit = [], []

    def add_coefs(model, label, extra_names=None):
        nm = model.name_x if hasattr(model, "name_x") else ["CONSTANT"] + names
        if hasattr(model, "name_z"):
            nm = model.name_z
        se = np.sqrt(np.diag(model.vm))
        zt = model.z_stat if hasattr(model, "z_stat") else model.t_stat
        for k, n in enumerate(nm):
            coefs.append({"model": label, "term": n, "coef": float(model.betas[k, 0]), "se": float(se[k]),
                          "stat": float(zt[k][0]), "p": float(zt[k][1])})

    def add_fit(model, label, resid, k):
        ll = float(model.logll) if hasattr(model, "logll") else np.nan
        n = len(yv)
        aic = -2 * ll + 2 * k if np.isfinite(ll) else np.nan
        I, zI, pI = moran_resid(resid, w)
        pr2 = float(getattr(model, "pr2", getattr(model, "r2", np.nan)))
        fit.append({"model": label, "n": n, "k": k, "logL": ll, "AIC": aic,
                    "AICc": aic + 2 * k * (k + 1) / (n - k - 1) if np.isfinite(ll) else np.nan,
                    "pseudo_R2": pr2, "resid_moran_I": I, "resid_moran_z": zI, "resid_moran_p": pI})

    bc.log(f"{name}: OLS")
    ols = spreg.OLS(yv, Xv, w=w, robust="white", spat_diag=True, moran=True, name_x=names, name_y=y.name)
    add_coefs(ols, "OLS")
    add_fit(ols, "OLS", ols.u, Xv.shape[1] + 2)
    diag = {"model": "OLS", "lm_lag": ols.lm_lag[0], "lm_lag_p": ols.lm_lag[1],
            "lm_error": ols.lm_error[0], "lm_error_p": ols.lm_error[1],
            "rlm_lag": ols.rlm_lag[0], "rlm_lag_p": ols.rlm_lag[1],
            "rlm_error": ols.rlm_error[0], "rlm_error_p": ols.rlm_error[1],
            "lm_sarma": ols.lm_sarma[0], "lm_sarma_p": ols.lm_sarma[1],
            "moran_res_I": ols.moran_res[0], "moran_res_z": ols.moran_res[1],
            "jarque_bera": ols.jarque_bera["jb"], "breusch_pagan": ols.breusch_pagan["bp"],
            "r2": ols.r2, "adj_r2": ols.ar2}
    bc.log(f"{name}: SLX")
    slx = spreg.OLS(yv, Xv, w=w, slx_lags=1, robust="white", name_x=names, name_y=y.name)
    add_coefs(slx, "SLX")
    add_fit(slx, "SLX", slx.u, 2 * Xv.shape[1] + 2)
    out = {"ols_resid": pd.Series(ols.u.ravel(), index=y.index), "ols": ols}
    impacts = []
    if ml:
        # ML with dense log-determinant and information matrix: spreg's
        # variance step inverts (I - rho W) densely, so ML is used only for
        # n <= ~8,000 (MSOA, LAD); it gives log-likelihoods for AIC(c).
        specs = [("SAR", spreg.ML_Lag, {}), ("SEM", spreg.ML_Error, {}),
                 ("SDM", spreg.ML_Lag, {"slx_lags": 1}), ("SDEM", spreg.ML_Error, {"slx_lags": 1})]
        for label, cls, kw in specs:
            bc.log(f"{name}: {label} (ML)")
            m = cls(yv, Xv, w=w, method="full", name_x=names, name_y=y.name, **kw)
            add_coefs(m, label)
            resid = m.e_filtered if label in ("SEM", "SDEM") else m.u
            add_fit(m, label, resid, len(m.betas) + 1)
            if label in ("SAR", "SDM"):
                impacts.append(_impacts(m, names, w, label))
    else:
        # GMM / IV estimators (Kelejian-Prucha), feasible at n = 35k:
        # spatial 2SLS for lag models (White SEs), heteroskedasticity-robust
        # GMM for error models. No likelihood, so no AIC.
        specs = [("SAR", spreg.GM_Lag, {"robust": "white"}), ("SEM", spreg.GM_Error_Het, {}),
                 ("SDM", spreg.GM_Lag, {"slx_lags": 1, "robust": "white"}),
                 ("SDEM", spreg.GM_Error_Het, {"slx_lags": 1})]
        for label, cls, kw in specs:
            bc.log(f"{name}: {label} (GMM)")
            m = cls(yv, Xv, w=w, name_x=names, name_y=y.name, **kw)
            add_coefs(m, label)
            resid = m.e_filtered if label in ("SEM", "SDEM") else m.u
            add_fit(m, label, resid, len(m.betas))
            if label in ("SAR", "SDM"):
                m.rho = float(m.betas[-1, 0])
                impacts.append(_impacts(m, names, w, label))
    out["impacts"] = pd.concat(impacts)
    out["coefs"] = pd.DataFrame(coefs)
    out["fit"] = pd.DataFrame(fit)
    out["diag"] = pd.DataFrame([diag])
    return out


def _impacts(m, names, w, label):
    """Direct / indirect / total impacts for SAR and SDM (LeSage-Pace
    average impacts with the 'simple' row-standardised-W formulae)."""
    rho = float(m.rho)
    b = m.betas.ravel()
    nm = list(m.name_z if hasattr(m, "name_z") else m.name_x)
    rows = []
    # Under row-standardised W, average total multiplier = 1/(1-rho).
    # Average direct effect needs trace((I-rho W)^-1)/n, approximated by the
    # power series 1 + sum_k rho^k tr(W^k)/n, k = 2..10 (tr(W) = 0).
    W = w.sparse
    trs = []
    Wk = W.copy()
    for k in range(2, 11):
        Wk = Wk @ W
        trs.append(Wk.diagonal().mean())
    mdir = 1 + sum(rho ** k * t for k, t in zip(range(2, 11), trs))
    mtot = 1 / (1 - rho)
    for n in names:
        bj = b[nm.index(n)]
        tj = b[nm.index("W_" + n)] if (label == "SDM" and "W_" + n in nm) else 0.0
        total = (bj + tj) * mtot
        # Direct effect: tr((I-rho W)^-1 (I bj + W tj))/n
        direct = bj * mdir + (tj * (mdir - 1) / rho if rho != 0 else 0.0)
        rows.append({"model": label, "term": n, "direct": direct, "indirect": total - direct, "total": total,
                     "rho": rho})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# GWR / MGWR
# --------------------------------------------------------------------------
def fit_mgwr(coords: np.ndarray, y: np.ndarray, X: np.ndarray, n_chunks: int = 8, multi: bool = True,
             init_bw=None, sel_cache: str | None = None):
    """Adaptive bisquare (nearest-neighbour) GWR or MGWR with AICc-selected
    bandwidth(s); y and X must already be standardised. Returns a slim dict
    (pickle-friendly) with params, t-values, corrected-t mask, bandwidths,
    ENP, AICc, R2 and residuals."""
    from mgwr.gwr import GWR, MGWR
    from mgwr.sel_bw import Sel_BW

    y = y.reshape(-1, 1)
    if not multi:
        sel = Sel_BW(coords, y, X)
        bw = sel.search(criterion="AICc")
        res = GWR(coords, y, X, bw).fit()
        bws = np.full(X.shape[1] + 1, bw)
        local_cn = None
        local_vif = None
    else:
        # Vectorised MGWR (b_mgwr_fast): mgwr's own backfitting and golden-
        # section search with all-locations kernels, and an exact vectorised
        # version of mgwr's inference recursion; validated against mgwr.
        import b_mgwr_fast as mf

        def _search():
            s_, _, _ = mf.search_multi(coords, y.ravel(), X, init_bw, multi_bw_min=30, tol_multi=1e-5,
                                       max_iter_multi=200, verbose=True)
            return s_
        # Backfitting is cached on its own so that inference can be re-run.
        sel = bc.cached(sel_cache, _search) if sel_cache else _search()
        res, sel = mf.fit(coords, y.ravel(), X, init_bw, verbose=True, sel=sel)
        bws = np.asarray(sel.bw[0])
        local_cn = res.local_cn_fast  # = MGWRResults.local_collinearity()[0], vectorised
        local_vif = None
    filt = res.filter_tvals()  # zeroes t where |t| < corrected critical value
    crit = res.critical_tval() if hasattr(res, "critical_tval") else None
    resid = np.asarray(res.resid_response).ravel()
    out = {
        "params": np.asarray(res.params), "tvalues": np.asarray(res.tvalues),
        "sig": np.asarray(filt) != 0, "bws": bws, "crit_t": None if crit is None else np.asarray(crit),
        "aicc": float(res.aicc), "aic": float(res.aic), "R2": float(res.R2), "adj_R2": float(res.adj_R2),
        "enp": float(res.ENP), "enp_j": np.asarray(getattr(res, "ENP_j", [np.nan])),
        "resid": resid, "local_cn": None if local_cn is None else np.asarray(local_cn).ravel(),
        "local_vif": None if local_vif is None else np.asarray(local_vif), "sigma2": float(res.sigma2),
        "llf": float(res.llf),
        "n_backfit_iter": None if not multi else int(np.asarray(sel.bw[1]).shape[0]),
        "soc_final": None if not multi else float(np.asarray(sel.bw[2])[-1]),
        "n_inference_sweeps": getattr(res, "n_inference_sweeps", None),
    }
    return out


def local_cn_summary(fit: dict) -> pd.DataFrame:
    """Distribution of the MGWR local condition number (values above 30 are
    the usual warning threshold)."""
    cn = np.asarray(fit.get("local_cn"))
    q = np.percentile(cn, [5, 25, 50, 75, 95, 99])
    return pd.DataFrame([{"p05": q[0], "p25": q[1], "median": q[2], "p75": q[3], "p95": q[4], "p99": q[5],
                          "max": cn.max(), "share_above_30": float((cn > 30).mean()),
                          "n_backfit_iter": fit.get("n_backfit_iter"), "soc_final": fit.get("soc_final"),
                          "n_inference_sweeps": fit.get("n_inference_sweeps"), "enp": fit.get("enp"),
                          "aicc": fit.get("aicc"), "R2": fit.get("R2")}])


def local_summary(fit: dict, names: list, w=None) -> pd.DataFrame:
    """Per-covariate summary of local coefficients."""
    rows = []
    p = fit["params"]
    for j, n in enumerate(["const"] + names):
        v = p[:, j]
        s = fit["sig"][:, j]
        rows.append({"term": n, "bandwidth_nn": fit["bws"][j], "mean": v.mean(), "sd": v.std(),
                     "min": v.min(), "p25": np.quantile(v, 0.25), "median": np.median(v),
                     "p75": np.quantile(v, 0.75), "max": v.max(),
                     "share_sig": s.mean(), "share_sig_pos": (s & (v > 0)).mean(),
                     "share_sig_neg": (s & (v < 0)).mean(),
                     "crit_t": None if fit["crit_t"] is None else float(np.ravel(fit["crit_t"])[min(j, np.size(fit["crit_t"]) - 1)]),
                     "enp_j": float(fit["enp_j"][j]) if len(fit["enp_j"]) > j else np.nan})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def coef_forest(coefs: pd.DataFrame, models: list, title: str, name: str, note: str):
    """Point estimates with 95% intervals for up to three models."""
    import matplotlib.pyplot as plt

    terms = [t for t in coefs[coefs.model == models[0]].term if t not in ("CONSTANT", "const")
             and not t.startswith("W_") and t != "lambda" and t != "W_dep_var"]
    d = coefs[coefs.model == models[0]].set_index("term").loc[terms].sort_values("coef")
    order = list(d.index)
    fig, ax = plt.subplots(figsize=(5.6, 0.36 * len(order) + 1.6))
    off = np.linspace(-0.22, 0.22, len(models)) if len(models) > 1 else [0.0]
    markers = ["o", "s", "D"]
    for m, o, col, mk in zip(models, off, bc.CAT3, markers):
        dm = coefs[coefs.model == m].set_index("term").reindex(order)
        yy = np.arange(len(order)) + o
        ax.hlines(yy, dm.coef - 1.96 * dm.se, dm.coef + 1.96 * dm.se, color=col, lw=1.8)
        ax.plot(dm.coef, yy, mk, color=col, ms=6.5, mec="white", mew=1.0, label=m, linestyle="none")
    ax.axvline(0, color=bc.INK2, lw=0.9)
    ax.set_yticks(np.arange(len(order)), [nice(t) for t in order])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Coefficient (log kWh per SD of covariate)")
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_title(title, loc="left", fontsize=13)
    if len(models) > 1:
        ax.legend(loc="lower right", frameon=False)
    fig.text(0.0, -0.02, note, fontsize=10, color=bc.INK2, wrap=True, va="top")
    return bc.save_fig(fig, name)


def mgwr_maps(level: str, ids, fit: dict, names: list, prefix: str, title_suffix: str, terms=None):
    """One map per local coefficient, non-significant (corrected t) in grey."""
    import b_maps as bm

    g = bc.display_geometry(level)
    paths = []
    allnames = ["const"] + names
    for j, n in enumerate(allnames):
        if terms is not None and n not in terms:
            continue
        s = pd.Series(fit["params"][:, j], index=ids, name="b")
        ns = pd.Series(~fit["sig"][:, j], index=ids)
        gg = g.join(s, how="inner")
        norm = bm.symmetric_norm(s.values)
        bw = int(round(fit["bws"][j]))
        p = bm.choropleth(gg, "b", f"{nice(n)}\n{title_suffix}", f"{prefix}_{n}", cmap=bc.DIVERGING, norm=norm,
                          ns_mask=ns, cbar_label="Local coefficient (SD units)",
                          boundaries=bc.display_geometry("lad") if level != "lad" else None,
                          note=f"MGWR, adaptive bisquare bandwidth = {bw:,} nearest {level.upper()}s "
                               f"({100 * fit['sig'][:, j].mean():.0f}% significant).")
        paths.append(p)
    return paths
