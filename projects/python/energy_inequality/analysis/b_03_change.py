"""Section 3: models of change in log kWh per meter on 2021 baselines.

Outcomes
* long:  log(mean 2022-24) - log(mean 2010-12)   (12-year change)
* shock: log(2022) - log(2019)                   (COVID + price-crisis shock;
         2019 is the last pre-COVID year, 2022 the first crisis meter-year)

Each outcome is regressed on the VIF-screened 2021 covariates, the log
change in meter counts over the same window (the one genuinely time-varying
covariate: new-build growth), and the initial log level (conditional
beta-convergence). Census 2021 values are treated as baselines: they are
not time-varying and do not identify *within*-area change mechanisms.

* LSOA: OLS / SLX / SAR / SEM / SDM / SDEM. LSOAs created by splits or
  irregular changes (S/X) are excluded from the long change because their
  2010-12 values are apportioned from the parent LSOA.
* MSOA: the same suite plus GWR (both changes) and MGWR (long change),
  on standardised variables.
* LAD: OLS and SEM adding domestic solar PV installations per household
  (Feed-in Tariff register, cumulative to March 2019; DESNZ sub-regional
  FiT statistics), which exist only at local-authority level.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import b_common as bc
import b_models as bmod

CHANGES = {"long": ((2010, 2012), (2022, 2024)), "shock": ((2019, 2019), (2022, 2022))}
CH_LABEL = {"long": "2010–12 to 2022–24", "shock": "2019 to 2022"}
FIT_XLS = bc.OUT / "external" / "December_2019_Sub-regional_Feed-in_Tariffs_confirmed_CFR.xls"
MAP_TERMS = ["const", "log_level_base", "share_ch_electric", "share_ch_renewable", "share_ch_offgas",
             "dep19_pct", "dlog_meters", "share_detached", "mean_hh_size", "share_private_rent", "epc_ee"]


def change_design(level, change, xs, initial=True):
    (a0, a1), (b0, b1) = CHANGES[change]
    base = bc.window_mean(level, range(a0, a1 + 1))
    end = bc.window_mean(level, range(b0, b1 + 1))
    e = bc.load_elec_lsoa()
    if level == "msoa":
        e = e.merge(bc.load_cov_lsoa()[["lsoa21cd", "msoa21cd"]], on="lsoa21cd")
    key = {"lsoa": "lsoa21cd", "msoa": "msoa21cd", "lad": "lad_cd"}[level]
    m0 = e[e.year.between(a0, a1)].groupby(key).meters.sum() / (a1 - a0 + 1)
    m1 = e[e.year.between(b0, b1)].groupby(key).meters.sum() / (b1 - b0 + 1)
    d = pd.DataFrame({"dlog_kwh": np.log(end) - np.log(base), "log_level_base": np.log(base),
                      "dlog_meters": np.log(m1) - np.log(m0)})
    cov = bc.covariates(level)
    d = cov[xs].join(d, how="inner").replace([np.inf, -np.inf], np.nan).dropna()
    if level == "lsoa":
        win = list(range(a0, a1 + 1)) + list(range(b0, b1 + 1))
        ee = e[e.year.isin(win)]
        bad = set(ee.loc[ee.flag_small_n, "lsoa21cd"])
        if change == "long":
            ch = bc.load_cov_lsoa().set_index("lsoa21cd").lsoa21_change
            bad |= set(ch[ch.isin(["S", "X"])].index)
        d = d.drop(index=list(bad & set(d.index)))
    xcols = xs + ["dlog_meters"] + (["log_level_base"] if initial else [])
    return d, xcols


def global_change(level, change, xs, standardise_y=False):
    def _fit():
        d, xcols = change_design(level, change, xs)
        X = bmod.standardise(d[xcols])
        y = d["dlog_kwh"].rename(f"dlog_kwh_{change}")
        if standardise_y:
            y = (y - y.mean()) / y.std()
        w = bc.subset_weights(bc.queen_weights(level), d.index, bc.unit_coords(level))
        res = bmod.fit_global_suite(y, X, w, f"{level} {change}", ml=(level != "lsoa"))
        res.pop("ols")
        res["ids"] = list(d.index)
        res["y_mean"], res["y_sd"] = float(d.dlog_kwh.mean()), float(d.dlog_kwh.std())
        return res
    return bc.cached(f"ch_{level}_{change}" + ("_stdy" if standardise_y else ""), _fit)


def beta_convergence(xs):
    """Unconditional and conditional beta-convergence (OLS, HC1) at all levels."""
    import statsmodels.api as sm

    rows = []
    for level in ("lad", "msoa", "lsoa"):
        for change in CHANGES:
            d, xcols = change_design(level, change, xs if level != "lad" else ["share_ch_electric"])
            for spec, cols in [("unconditional", ["log_level_base"]), ("conditional", xcols)]:
                X = sm.add_constant(d[cols])
                r = sm.OLS(d.dlog_kwh, X).fit(cov_type="HC1")
                b = r.params["log_level_base"]
                rows.append({"level": level, "change": change, "spec": spec, "n": len(d), "beta": b,
                             "se": r.bse["log_level_base"], "t": r.tvalues["log_level_base"], "r2": r.rsquared,
                             "mean_change": d.dlog_kwh.mean()})
    t = pd.DataFrame(rows)
    bc.save_table(t, "ch_beta_convergence")
    return t


def mgwr_change(change, xs, force=False, multi=True):
    def _fit(multi, init=None):
        d, xcols = change_design("msoa", change, xs)
        coords = bc.unit_coords("msoa").loc[d.index].values
        X = bmod.standardise(d[xcols]).values
        y = d.dlog_kwh.values
        y = (y - y.mean()) / y.std()
        bc.log(f"{'MGWR' if multi else 'GWR'} msoa {change}: n={len(y)}, k={X.shape[1]}")
        f = bmod.fit_mgwr(coords, y, X, multi=multi, init_bw=init,
                          sel_cache=f"mgwr_sel_msoa_ch_{change}" if multi else None)
        f["ids"] = list(d.index)
        f["xcols"] = xcols
        return f
    gwr = bc.cached(f"gwr_msoa_ch_{change}", lambda: _fit(False), force)
    mg = bc.cached(f"mgwr_msoa_ch_{change}", lambda: _fit(True, float(gwr["bws"][0])), force) if multi else None
    if gwr.get("xcols") is None:
        gwr["xcols"] = change_design("msoa", change, xs)[1]
    return gwr, mg


# --------------------------------------------------------------------------
# Solar PV (LAD level)
# --------------------------------------------------------------------------
def solar_lad() -> pd.Series:
    """Domestic PV installations per 100 households (2021 households) by LAD
    (May 2026), from the final DESNZ sub-regional FiT table (cumulative to
    end March 2019). 2019 LA codes are mapped to 2026 LADs directly where
    the code still exists, otherwise via the 2011 OA lookups."""
    raw = pd.read_excel(FIT_XLS, sheet_name="Latest Quarter - LA", header=None, engine="calamine")
    t = raw.iloc[7:, [1, 2, 4]].copy()
    t.columns = ["la_code", "la_name", "pv_domestic"]
    t = t[t.la_code.astype(str).str.match(r"^[EW]\d{8}$")]
    t["pv_domestic"] = pd.to_numeric(t.pv_domestic, errors="coerce")
    geo = bc.DATA / "raw" / "geography"
    lad26 = set(bc.load_lad_geo().lad_cd)
    oa11 = pd.read_csv(geo / "OA11_LSOA11_MSOA11_LAD11_EW_LU.csv", encoding="utf-8-sig", dtype=str)
    oa = pd.read_csv(geo / "OA11_OA21_LAD22_EW_LU_v3.csv", encoding="utf-8-sig", dtype=str)
    c11 = [c for c in oa11.columns if c.upper().startswith("OA11")][0]
    l11 = [c for c in oa11.columns if c.upper().startswith("LAD11CD")][0]
    o11 = [c for c in oa.columns if c.upper().startswith("OA11")][0]
    l22 = [c for c in oa.columns if c.upper().startswith("LAD22CD")][0]
    m = oa11[[c11, l11]].merge(oa[[o11, l22]], left_on=c11, right_on=o11)
    l22_26 = bc.load_lsoa_geo().groupby("lad22cd").lad_cd.agg(lambda s: s.mode().iloc[0])
    m["lad26"] = m[l22].map(l22_26)
    old_to_new = m.groupby(l11).lad26.agg(lambda s: s.mode().iloc[0])
    manual = {"E07000246": "E06000066"}  # Somerset West and Taunton -> Somerset (2023)

    def map_code(c):
        if c in lad26:
            return c
        if c in manual:
            return manual[c]
        return old_to_new.get(c, np.nan)

    t["lad_cd"] = t.la_code.map(map_code)
    miss = t[t.lad_cd.isna()]
    if len(miss):
        bc.log(f"FiT codes not mapped: {miss.la_code.tolist()}")
    pv = t.groupby("lad_cd").pv_domestic.sum()
    hh = bc.load_lsoa_geo().groupby("lad_cd").households_2021.sum()
    out = (100 * pv / hh).rename("pv_per_100hh").dropna()
    pd.DataFrame({"pv_domestic": pv, "households_2021": hh, "pv_per_100hh": out}).reset_index() \
        .rename(columns={"index": "lad_cd"}).pipe(bc.save_table, "ch_lad_solar_fit")
    return out


def lad_change_models():
    """LAD-level change models with solar PV (OLS HC1 and ML SEM)."""
    import spreg

    xs = ["share_ch_electric", "share_ch_offgas", "share_ch_renewable", "share_detached", "mean_hh_size",
          "share_private_rent", "dep19_pct", "epc_ee"]
    pv = solar_lad()
    rows, fits = [], []
    for change in CHANGES:
        d, xcols = change_design("lad", change, xs)
        d = d.join(pv, how="inner").dropna()
        xcols = xcols + ["pv_per_100hh"]
        X = bmod.standardise(d[xcols])
        y = d.dlog_kwh.values.reshape(-1, 1)
        w = bc.subset_weights(bc.queen_weights("lad"), d.index, bc.unit_coords("lad"))
        ols = spreg.OLS(y, X.values, w=w, robust="white", spat_diag=True, name_x=xcols)
        sem = spreg.ML_Error(y, X.values, w=w, method="full", name_x=xcols)
        for lab, mdl in [("OLS", ols), ("SEM", sem)]:
            se = np.sqrt(np.diag(mdl.vm))
            zt = mdl.z_stat if hasattr(mdl, "z_stat") else mdl.t_stat
            for k, n in enumerate(mdl.name_x):
                rows.append({"change": change, "model": lab, "term": n, "coef": mdl.betas[k, 0], "se": se[k],
                             "p": zt[k][1]})
        I, z, p = bmod.moran_resid(ols.u, w)
        fits.append({"change": change, "n": len(d), "ols_r2": ols.r2, "ols_resid_moran_I": I, "ols_resid_moran_p": p,
                     "lm_error": ols.lm_error[0], "lm_error_p": ols.lm_error[1], "sem_lambda": float(sem.lam),
                     "sem_logL": float(sem.logll), "ols_logL": float(ols.logll),
                     "corr_pv_dep": d.pv_per_100hh.corr(d.dep19_pct),
                     "corr_pv_detached": d.pv_per_100hh.corr(d.share_detached)})
    c = pd.DataFrame(rows)
    bc.save_table(c, "ch_lad_solar_coefs")
    bc.save_table(pd.DataFrame(fits), "ch_lad_solar_fit_stats")
    return c


def main(force=False):
    bc.set_style()
    xs = bc.cached("xs_vif_keep", lambda: (_ for _ in ()).throw(RuntimeError("run b_02 first")))
    beta = beta_convergence(xs)
    bc.log(f"\n{beta.to_string()}")
    for change in CHANGES:
        res = global_change("lsoa", change, xs)
        for k in ("coefs", "fit", "diag", "impacts"):
            bc.save_table(res[k], f"ch_lsoa_{change}_{k}")
        bc.log(f"\nLSOA change {change}\n{res['fit'].to_string()}\n{res['diag'].T.to_string()}")
        bmod.coef_forest(res["coefs"], ["OLS", "SLX", "SEM"],
                         f"Correlates of change in log kWh per meter\nLSOAs, {CH_LABEL[change]}",
                         f"ch_lsoa_{change}_forest",
                         f"n = {res['fit'].n.iloc[0]:,} LSOAs; covariates standardised; 95% intervals. "
                         f"Negative = larger fall. Census 2021 covariates are baselines.")
    lad = lad_change_models()
    bc.log(f"\n{lad.to_string()}")
    xs = [x for x in bmod.XS_LOCAL if x in xs]  # MSOA models use the local-model covariate set
    for change in CHANGES:
        glob = global_change("msoa", change, xs, standardise_y=True)
        for k in ("coefs", "fit", "diag", "impacts"):
            bc.save_table(glob[k], f"ch_msoa_{change}_{k}")
        # MGWR (1-2 h per fit here) for the long change only; the shock
        # change gets the global suite and GWR.
        gwr, mg = mgwr_change(change, xs, force, multi=(change == "long"))
        xcols = gwr["xcols"]
        bc.save_table(bmod.local_summary(gwr, xcols), f"ch_gwr_msoa_{change}_local_summary")
        import b_02_cross_section as b2
        cmp_ = b2.comparison_table(change, glob, gwr, mg, f"ch_msoa_{change}")
        bc.log(f"\nMSOA change {change}\n{cmp_.to_string()}")
        if mg is None:
            continue
        summ = bmod.local_summary(mg, xcols)
        bc.save_table(summ, f"ch_mgwr_msoa_{change}_local_summary")
        bc.save_table(bmod.local_cn_summary(mg), f"ch_mgwr_msoa_{change}_diagnostics")
        bc.log(f"\n{summ.to_string()}")
        local = pd.DataFrame(mg["params"], index=mg["ids"], columns=["const"] + xcols)
        local.index.name = "msoa21cd"
        local.join(pd.DataFrame(mg["sig"], index=mg["ids"], columns=[f"sig_{c}" for c in ["const"] + xcols])) \
            .reset_index().to_parquet(bc.TAB / f"ch_mgwr_msoa_{change}_local.parquet")
        b2.bandwidth_figure(mg, xcols, f"ch_mgwr_msoa_{change}_bandwidths",
                            f"MGWR bandwidths, change {CH_LABEL[change]}")
        bmod.mgwr_maps("msoa", mg["ids"], mg, xcols, f"ch_mgwr_msoa_{change}_map",
                       f"MGWR, change {CH_LABEL[change]}", terms=[t for t in MAP_TERMS if t == "const" or t in xcols])


if __name__ == "__main__":
    main(force="--force" in sys.argv)
