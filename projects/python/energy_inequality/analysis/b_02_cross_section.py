"""Section 2: cross-sectional explanatory models of log kWh per meter.

Outcomes: log meter-weighted mean of 2022-2024 (recent, post-shock) and of
2019 (pre-COVID, pre-crisis). Covariates: Census 2021 shares, RUC 2021,
IMD/WIMD 2019 within-country percentile, MSOA EPC efficiency.

* LSOA (n ~ 35.6k): VIF screening, OLS with White SEs and LM diagnostics,
  SLX, SAR, SEM, SDM, SDEM (ML, sparse LU), impacts.
* MSOA (n = 7,264): the same suite plus GWR (both outcomes) and MGWR
  (2022-24). MSOA is the MGWR
  design because (i) the EPC covariate is native to MSOAs, so LSOA-level
  local models would replicate one EPC value across ~5 neighbours, (ii) an
  MGWR at n = 7,264 is computationally feasible on the full national
  surface (no sampling), and (iii) MSOA kWh per meter rests on ~3,600
  meters, so the outcome is less noisy.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import b_common as bc
import b_models as bmod

OUTCOMES = {"y2224": range(2022, 2025), "y2019": range(2019, 2020)}
OUT_LABEL = {"y2224": "2022–2024", "y2019": "2019"}
MAP_TERMS = ["const", "share_ch_electric", "share_ch_offgas", "share_detached", "mean_hh_size",
             "dep19_pct", "epc_ee", "share_private_rent", "share_flat", "share_age_65plus", "share_ch_renewable"]


def small_n_lsoas(years) -> set:
    e = bc.load_elec_lsoa()
    e = e[e.year.isin(list(years))]
    return set(e.loc[e.flag_small_n, "lsoa21cd"])


def outlier_lsoas(years) -> set:
    e = bc.load_elec_lsoa()
    e = e[e.year.isin(list(years))]
    return set(e.loc[e.flag_outlier, "lsoa21cd"])


def design(level: str, outcome: str, xs=None):
    yrs = OUTCOMES[outcome]
    y = np.log(bc.window_mean(level, yrs)).rename(f"log_kwh_{outcome}")
    cov = bc.covariates(level)
    xs = xs or bmod.CANDIDATES
    d = cov[xs].join(y, how="inner").dropna()
    if level == "lsoa":
        d = d.drop(index=list(small_n_lsoas(yrs) & set(d.index)))
    return d


def screen() -> list:
    """VIF screening on the LSOA 2022-24 design (threshold 5)."""
    d = design("lsoa", "y2224")
    keep, steps = bmod.vif_screen(d[bmod.CANDIDATES], 5.0)
    bc.save_table(steps, "xs_vif_steps")
    final = bmod.vif_table(d[keep]).rename("vif").reset_index().rename(columns={"index": "term"})
    msoa = design("msoa", "y2224")
    final["vif_msoa"] = final.term.map(bmod.vif_table(msoa[keep]))
    bc.save_table(final, "xs_vif_final")
    bc.log(f"VIF kept {len(keep)}: {keep}\n{final}")
    return keep


def global_models(level, outcome, xs, standardise_y=False):
    def _fit():
        d = design(level, outcome, xs)
        X = bmod.standardise(d[xs])
        y = d[f"log_kwh_{outcome}"]
        if standardise_y:
            y = (y - y.mean()) / y.std()
        w = bc.subset_weights(bc.queen_weights(level), d.index, bc.unit_coords(level))
        res = bmod.fit_global_suite(y, X, w, f"{level} {outcome}", ml=(level != "lsoa"))
        res.pop("ols")
        res["ids"] = list(d.index)
        res["y_mean"], res["y_sd"] = float(d[f"log_kwh_{outcome}"].mean()), float(d[f"log_kwh_{outcome}"].std())
        return res
    tag = f"xs_{level}_{outcome}" + ("_stdy" if standardise_y else "")
    return bc.cached(tag, _fit)


def outlier_sensitivity(xs):
    """OLS with and without LSOAs carrying any outlier flag in the window."""
    import statsmodels.api as sm

    rows = []
    for outcome in OUTCOMES:
        d = design("lsoa", outcome, xs)
        X = sm.add_constant(bmod.standardise(d[xs]))
        y = d[f"log_kwh_{outcome}"]
        drop = list(outlier_lsoas(OUTCOMES[outcome]) & set(d.index))
        for label, mask in [("all", np.ones(len(d), bool)), ("ex_outliers", ~d.index.isin(drop))]:
            r = sm.OLS(y[mask], X[mask]).fit(cov_type="HC1")
            for t in xs:
                rows.append({"outcome": outcome, "sample": label, "n": int(mask.sum()), "term": t,
                             "coef": r.params[t], "se": r.bse[t], "r2": r.rsquared})
    t = pd.DataFrame(rows)
    bc.save_table(t, "xs_lsoa_outlier_sensitivity")
    return t


def mgwr_models(outcome, xs, force=False, multi=True):
    """GWR and (optionally) MGWR at MSOA on standardised y and X; MGWR
    starts its backfitting from the GWR bandwidth."""
    def _fit(multi, init=None):
        d = design("msoa", outcome, xs)
        coords = bc.unit_coords("msoa").loc[d.index].values
        X = bmod.standardise(d[xs]).values
        y = d[f"log_kwh_{outcome}"].values
        y = (y - y.mean()) / y.std()
        bc.log(f"{'MGWR' if multi else 'GWR'} msoa {outcome}: n={len(y)}, k={X.shape[1]}")
        f = bmod.fit_mgwr(coords, y, X, multi=multi, init_bw=init,
                          sel_cache=f"mgwr_sel_msoa_{outcome}" if multi else None)
        f["ids"] = list(d.index)
        return f
    gwr = bc.cached(f"gwr_msoa_{outcome}", lambda: _fit(False), force)
    mg = bc.cached(f"mgwr_msoa_{outcome}", lambda: _fit(True, float(gwr["bws"][0])), force) if multi else None
    return gwr, mg


def comparison_table(outcome, glob, gwr, mg, tag):
    """AICc and residual Moran's I for global vs local models (same y scale)."""
    w = bc.subset_weights(bc.queen_weights("msoa"), gwr["ids"], bc.unit_coords("msoa"))
    rows = glob["fit"].to_dict("records")
    n = len(gwr["ids"])
    for label, f in [("GWR", gwr), ("MGWR", mg)]:
        if f is None:
            continue
        K = f["enp"] + 1
        I, z, p = bmod.moran_resid(f["resid"], w)
        rows.append({"model": label, "n": n, "k": K, "logL": f["llf"], "AIC": -2 * f["llf"] + 2 * K,
                     "AICc": -2 * f["llf"] + 2 * K + 2 * K * (K + 1) / (n - K - 1),
                     "pseudo_R2": f["R2"], "resid_moran_I": I, "resid_moran_z": z, "resid_moran_p": p,
                     "mgwr_package_aicc": f["aicc"]})
    t = pd.DataFrame(rows)
    t["delta_AICc"] = t.AICc - t.AICc.min()
    bc.save_table(t, f"{tag}_model_comparison")
    return t


def bandwidth_figure(mg, xs, name, title):
    import matplotlib.pyplot as plt

    terms = ["const"] + xs
    bw = pd.Series(mg["bws"], index=terms).sort_values()
    fig, ax = plt.subplots(figsize=(5.6, 0.36 * len(bw) + 1.4))
    ax.barh(np.arange(len(bw)), bw.values, color=bc.CAT3[0], height=0.6)
    ax.set_yticks(np.arange(len(bw)), [bmod.nice(t) for t in bw.index])
    ax.set_xscale("log")
    ax.axvline(len(mg["ids"]), color=bc.INK2, lw=0.9, ls="--")
    ax.text(len(mg["ids"]), len(bw) - 0.4, " global", fontsize=10, color=bc.INK2, va="top", ha="right")
    for k, v in enumerate(bw.values):
        ax.text(v * 1.05, k, f"{int(v):,}", va="center", fontsize=10.5, color=bc.INK)
    ax.set_xlabel("Bandwidth (nearest MSOAs, log scale)")
    ax.grid(axis="y", visible=False)
    ax.set_title(title, loc="left", fontsize=13)
    fig.text(0.0, -0.02, f"Adaptive bisquare kernels; n = {len(mg['ids']):,} MSOAs. Smaller = more local.",
             fontsize=10, color=bc.INK2, va="top")
    return bc.save_fig(fig, name)


def validate_mgwr_fast():
    """Check the vectorised MGWR against the mgwr package on 500 MSOAs."""
    import b_mgwr_fast as mf

    xs = ["share_ch_electric", "share_detached", "mean_hh_size", "dep19_pct", "epc_ee"]
    d = design("msoa", "y2224", xs)
    c = bc.unit_coords("msoa").loc[d.index].values
    r = mf.validate(c, d[f"log_kwh_y2224"].values, d[xs].values, n_sub=500)
    r["covariates"] = ", ".join(xs)
    t = pd.DataFrame([{k: (str(v) if isinstance(v, list) else v) for k, v in r.items()}])
    bc.save_table(t, "mgwr_fast_validation")
    return t


def main(force=False):
    bc.set_style()
    bc.log(f"\n{bc.cached('mgwr_fast_validation', validate_mgwr_fast, force).T}")
    xs = bc.cached("xs_vif_keep", screen, force)
    for outcome in OUTCOMES:
        res = global_models("lsoa", outcome, xs)
        for k in ("coefs", "fit", "diag", "impacts"):
            bc.save_table(res[k], f"xs_lsoa_{outcome}_{k}")
        bc.log(f"\nLSOA {outcome}\n{res['fit'].to_string()}\n{res['diag'].T.to_string()}")
        bmod.coef_forest(res["coefs"], ["OLS", "SLX", "SEM"],
                         f"Correlates of log kWh per meter, LSOAs, {OUT_LABEL[outcome]}",
                         f"xs_lsoa_{outcome}_forest",
                         f"n = {res['fit'].n.iloc[0]:,} LSOAs; covariates standardised; 95% intervals "
                         f"(OLS/SLX White; SEM heteroskedasticity-robust GMM). SLX: own-LSOA terms only.")
    outlier_sensitivity(xs)
    xl = [x for x in bmod.XS_LOCAL if x in xs]
    msoa = design("msoa", "y2224", xl)
    v = bmod.vif_table(msoa[xl]).rename("vif_msoa").reset_index().rename(columns={"index": "term"})
    bc.save_table(v, "xs_vif_msoa_local_set")
    for outcome in OUTCOMES:
        xs = xl  # MSOA models use the local-model covariate set
        glob = global_models("msoa", outcome, xs, standardise_y=True)
        for k in ("coefs", "fit", "diag", "impacts"):
            bc.save_table(glob[k], f"xs_msoa_{outcome}_{k}")
        # MGWR (about an hour per fit at n = 7,264, k = 18) is fitted for the
        # main 2022-24 outcome; 2019 gets GWR and the global suite only.
        gwr, mg = mgwr_models(outcome, xs, force, multi=(outcome == "y2224"))
        bc.save_table(bmod.local_summary(gwr, xs), f"xs_gwr_msoa_{outcome}_local_summary")
        cmp_ = comparison_table(outcome, glob, gwr, mg, f"xs_msoa_{outcome}")
        bc.log(f"\nMSOA {outcome} comparison\n{cmp_.to_string()}")
        if mg is None:
            continue
        summ = bmod.local_summary(mg, xs)
        bc.save_table(summ, f"xs_mgwr_msoa_{outcome}_local_summary")
        bc.save_table(bmod.local_cn_summary(mg), f"xs_mgwr_msoa_{outcome}_diagnostics")
        bc.log(f"\n{summ.to_string()}")
        local = pd.DataFrame(mg["params"], index=mg["ids"], columns=["const"] + xs)
        local.index.name = "msoa21cd"
        local.join(pd.DataFrame(mg["sig"], index=mg["ids"], columns=[f"sig_{c}" for c in ["const"] + xs])) \
            .reset_index().to_parquet(bc.TAB / f"xs_mgwr_msoa_{outcome}_local.parquet")
        bandwidth_figure(mg, xs, f"xs_mgwr_msoa_{outcome}_bandwidths",
                         f"MGWR bandwidths, log kWh per meter, {OUT_LABEL[outcome]}")
        if outcome == "y2224":
            bmod.mgwr_maps("msoa", mg["ids"], mg, xs, f"xs_mgwr_msoa_{outcome}_map",
                           f"MGWR local coefficient, {OUT_LABEL[outcome]}",
                           terms=[t for t in MAP_TERMS if t == "const" or t in xs])


if __name__ == "__main__":
    main(force="--force" in sys.argv)
