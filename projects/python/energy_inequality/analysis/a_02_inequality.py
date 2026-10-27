"""Section 2: inequality of LSOA mean kWh/meter.

* Theil T (GE(1)) and MLD (GE(0)), meter-weighted, with an exact nested
  decomposition region -> LAD -> LSOA; unweighted (place-based) Theil checked
  against ``inequality.theil.TheilD``; permutation test of the between-LAD
  component (random relabelling of LSOAs to LADs of the same sizes).
* Gini (place-based and meter-weighted) and the spatial Gini of Rey & Smith
  (2013) with permutation inference (binary queen contiguity).
* Within-LAD dispersion vs LAD mean: Spearman correlation with bootstrap CI,
  OLS of log dispersion on log mean with controls and HC3 errors, Moran's I of
  the residuals and a spatial-error (GMM, heteroskedasticity-robust) check.

Outputs: tables/ineq_*.{parquet,csv}, figures/ineq_*.{png,svg}.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import esda  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from inequality.gini import Gini_Spatial  # noqa: E402
from inequality.theil import TheilD  # noqa: E402
from libpysal import weights  # noqa: E402
from scipy import stats  # noqa: E402

from a_common import *  # noqa: E402,F401,F403


# --------------------------------------------------------------------------- Theil
def theil_t(e, p):
    """Meter-weighted Theil T: sum s_i ln(s_i / pi_i)."""
    s, pi = e / e.sum(), p / p.sum()
    return float(np.sum(s * np.log(s / pi)))


def mld(e, p):
    """Meter-weighted mean log deviation (Theil L): sum pi_i ln(pi_i / s_i)."""
    s, pi = e / e.sum(), p / p.sum()
    return float(np.sum(pi * np.log(pi / s)))


def theil_decomp(e, p, groups):
    """Exact T = between + within for a partition. Returns (total, between, within)."""
    g, inv = np.unique(groups, return_inverse=True)
    E, P = np.bincount(inv, e), np.bincount(inv, p)
    S, PI = E / E.sum(), P / P.sum()
    between = float(np.sum(S * np.log(S / PI)))
    s_in, pi_in = e / E[inv], p / P[inv]
    t_g = np.bincount(inv, s_in * np.log(s_in / pi_in))
    within = float(np.sum(S * t_g))
    return theil_t(e, p), between, within


def mld_decomp(e, p, groups):
    g, inv = np.unique(groups, return_inverse=True)
    E, P = np.bincount(inv, e), np.bincount(inv, p)
    S, PI = E / E.sum(), P / P.sum()
    between = float(np.sum(PI * np.log(PI / S)))
    s_in, pi_in = e / E[inv], p / P[inv]
    l_g = np.bincount(inv, pi_in * np.log(pi_in / s_in))
    return mld(e, p), between, float(np.sum(PI * l_g))


def theil_tables(C, reg_of_lad):
    rows, rng = [], np.random.default_rng(SEED)
    lad = pd.factorize(C["lad"])[0]                       # integer codes: fast np.unique
    region = pd.factorize(pd.Series(C["lad"]).map(reg_of_lad))[0]
    for t in YEARS:
        p = C["meters"][t].values
        e = C["y"][t].values * p
        T, b_reg, _ = theil_decomp(e, p, region)
        _, b_lad, w_lad = theil_decomp(e, p, lad)
        # nested: between LADs within regions = between-LAD total - between-region
        L, lb_reg, _ = mld_decomp(e, p, region)
        _, lb_lad, lw_lad = mld_decomp(e, p, lad)
        # unweighted (place-based) Theil, checked against PySAL
        y = C["y"][t].values
        td = TheilD(y, lad)
        T_unw = float(td.T[0] if np.ndim(td.T) else td.T)
        b_unw = float(np.ravel(td.bg)[0])
        ones = np.ones_like(y)
        assert abs(theil_decomp(y, ones, lad)[0] - T_unw) < 1e-9 and abs(theil_decomp(y, ones, lad)[1] - b_unw) < 1e-6
        # permutation test of between-LAD share: shuffle LAD labels
        perm = np.array([theil_decomp(e, p, rng.permutation(lad))[1] for _ in range(PERMS)])
        rows.append({
            "year": t, "theil_T": T, "between_region": b_reg, "between_lad_within_region": b_lad - b_reg,
            "within_lad": w_lad, "share_between_region": b_reg / T,
            "share_between_lad_within_region": (b_lad - b_reg) / T, "share_within_lad": w_lad / T,
            "share_between_lad": b_lad / T, "check_sum": b_lad + w_lad - T,
            "mld": L, "mld_share_between_region": lb_reg / L, "mld_share_between_lad": lb_lad / L,
            "mld_share_within_lad": lw_lad / L,
            "theil_T_unweighted": T_unw, "share_between_lad_unweighted": b_unw / T_unw,
            "between_lad_perm_mean": perm.mean(), "between_lad_perm_p": (np.sum(perm >= b_lad) + 1) / (PERMS + 1),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- Gini
def gini(x, w=None):
    x = np.asarray(x, float)
    if w is None:
        xs = np.sort(x)
        n = xs.size
        return float((2 * np.arange(1, n + 1) - n - 1) @ xs / (n * xs.sum()))
    o = np.argsort(x)
    x, w = x[o], np.asarray(w, float)[o]
    cw = np.cumsum(w)
    cxw = np.cumsum(x * w)
    # area under the weighted Lorenz curve
    W, X = cw[-1], cxw[-1]
    lor = np.r_[0, cxw / X]
    pop = np.r_[0, cw / W]
    return float(1 - np.sum((pop[1:] - pop[:-1]) * (lor[1:] + lor[:-1])))


def spatial_gini(x, i, j, perms, rng):
    """Rey & Smith (2013) spatial Gini on binary symmetric contiguity given as
    directed edge arrays. Vectorised equivalent of inequality.gini.Gini_Spatial."""
    x = np.asarray(x, float)
    n = x.size
    g = gini(x)
    d = g * 2 * n * n * x.mean()          # sum_i sum_j |x_i - x_j|
    wg = np.abs(x[i] - x[j]).sum()          # neighbour pairs (both directions)
    wcg = d - wg
    n_pairs = n * (n - 1) / 2
    nn = i.size / 2
    scale = nn / (n_pairs - nn)
    pol = wcg / wg * scale
    sims = np.empty(perms)
    for k in range(perms):
        xp = x[rng.permutation(n)]
        sims[k] = d - np.abs(xp[i] - xp[j]).sum()
    larger = (sims >= wcg).sum()
    larger = min(larger, perms - larger)
    pol_sim = sims / (d - sims) * scale
    return {"gini": g, "wg_neighbour": wg, "wcg_nonneighbour": wcg, "wcg_share": wcg / d,
            "neighbour_share": wg / d, "expected_neighbour_share_random": nn * 2 / (n * (n - 1)),
            "e_wcg": sims.mean(), "z_wcg": (wcg - sims.mean()) / sims.std(), "p_sim": (larger + 1) / (perms + 1),
            "polarization": pol, "polarization_p_sim": ((pol_sim >= pol).sum() + 1) / (perms + 1),
            "mean_abs_diff_neighbours": wg / i.size,
            "mean_abs_diff_all": d / (n * (n - 1))}


def gini_tables(C, wb):
    i, j = edge_arrays(wb)
    rng = np.random.default_rng(SEED)
    rows = []
    for t in YEARS:
        x = C["y"][t].values
        r = spatial_gini(x, i, j, PERMS, rng)
        r.update({"year": t, "gini_meter_weighted": gini(x, C["meters"][t].values)})
        rows.append(r)
    # validate against PySAL on one year (deterministic parts only)
    t = 2024
    ref = Gini_Spatial(C["y"][t].values, _binary_from_edges(i, j, C["y"].shape[0]), permutations=0)
    chk = rows[-1]
    assert abs(ref.g - chk["gini"]) < 1e-10 and abs(ref.wg - chk["wg_neighbour"]) / ref.wg < 1e-9, "spatial Gini mismatch"
    return pd.DataFrame(rows)


def _binary_from_edges(i, j, n):
    nb = {k: [] for k in range(n)}
    for a, b in zip(i.tolist(), j.tolist()):
        nb[a].append(b)
    return weights.W(nb, silence_warnings=True)


# --------------------------------------------------------------------------- dispersion vs mean
def lad_dispersion(C, Cx, reg_of_lad):
    lp = pd.read_parquet(PROCESSED / "elec_lad_panel.parquet")
    rows = []
    for name, CC in (("all", C), ("excl_flagged", Cx)):
        for t in YEARS:
            d = pd.DataFrame({"lad_cd": CC["lad"], "y": CC["y"][t].values, "m": CC["meters"][t].values})
            d["e"] = d.y * d.m
            g = d.groupby("lad_cd")
            out = pd.DataFrame({
                "n_lsoa": g.size(),
                "mean_lsoa_unweighted": g.y.mean(),
                "median_lsoa": g.y.median(),
                "mean_kwh": g.e.sum() / g.m.sum(),
                "cv": g.y.std(ddof=1) / g.y.mean(),
                "cv_meter_weighted": g.apply(lambda s: np.sqrt(np.average((s.y - np.average(s.y, weights=s.m)) ** 2,
                                                                     weights=s.m)) / np.average(s.y, weights=s.m)),
                "qcd": g.y.quantile(0.75).sub(g.y.quantile(0.25)) / g.y.quantile(0.75).add(g.y.quantile(0.25)),
                "theil_within": g.apply(lambda s: theil_t(s.e.values, s.m.values)),
            }).reset_index()
            out["year"], out["sample"] = t, name
            rows.append(out)
    D = pd.concat(rows, ignore_index=True)
    D["region"] = D.lad_cd.map(reg_of_lad)
    # cross-check LAD mean against the processed LAD panel (all LSOAs)
    chk = D[D["sample"] == "all"].merge(lp[["lad_cd", "year", "mean_kwh"]], on=["lad_cd", "year"], suffixes=("", "_panel"))
    D.attrs["max_rel_diff_vs_lad_panel"] = float((chk.mean_kwh / chk.mean_kwh_panel - 1).abs().max())
    return D


def dispersion_tests(D, lad_geo, w_lad):
    rng = np.random.default_rng(SEED)
    rows = []
    order = lad_geo["lad_cd"].values
    for (sample, t), d in D.groupby(["sample", "year"]):
        d = d.set_index("lad_cd").loc[order].reset_index()
        # dispersion needs several LSOAs: drop LADs with < 3 (Isles of Scilly has 1)
        keep = np.flatnonzero(d.n_lsoa.values >= 3)
        d = d.iloc[keep].reset_index(drop=True)
        wk = weights.w_subset(w_lad, list(keep), silence_warnings=True)
        wk = weights.W({keep.tolist().index(a): [keep.tolist().index(b) for b in nb]
                        for a, nb in wk.neighbors.items()}, silence_warnings=True)
        wk.transform = "r"
        for level, disp in (("mean_kwh", "cv"), ("mean_kwh", "cv_meter_weighted"), ("mean_kwh", "qcd"),
                            ("median_lsoa", "cv"), ("median_lsoa", "qcd")):
            # the median as the level guards against a mechanical link in which
            # a heavy upper tail raises both the LAD mean and its CV
            x, y = np.log(d[level].values), np.log(d[disp].values)
            rho = stats.spearmanr(x, y).statistic
            bs = []
            for _ in range(2000):
                k = rng.integers(0, len(d), len(d))
                bs.append(stats.spearmanr(x[k], y[k]).statistic)
            dd = d.assign(lcv=y, lmean=x, ln_n=np.log(d.n_lsoa))
            ols = smf.ols("lcv ~ lmean + ln_n + C(region)", dd).fit(cov_type="HC3")
            res = ols.resid.values
            np.random.seed(SEED)
            mi = esda.Moran(res, wk, permutations=PERMS)
            row = {"sample": sample, "year": t, "level": level, "dispersion": disp, "n_lad": len(d),
                   "spearman_rho": rho,
                   "rho_ci_lo": np.percentile(bs, 2.5), "rho_ci_hi": np.percentile(bs, 97.5),
                   "elasticity_ols": ols.params["lmean"], "elasticity_se_hc3": ols.bse["lmean"],
                   "elasticity_p": ols.pvalues["lmean"], "r2": ols.rsquared,
                   "resid_moran_I": mi.I, "resid_moran_p": mi.p_sim}
            if disp == "cv" and level == "mean_kwh":
                from spreg import GM_Error_Het
                X = pd.get_dummies(dd[["lmean", "ln_n", "region"]], columns=["region"], drop_first=True,
                                   dtype=float)
                m = GM_Error_Het(dd[["lcv"]].values, X.values, w=wk, name_x=list(X.columns))
                k = 1 + list(X.columns).index("lmean")
                row.update({"elasticity_sem": float(m.betas[k][0]), "elasticity_sem_se": float(m.std_err[k]),
                            "elasticity_sem_p": float(m.z_stat[k][1]), "sem_lambda": float(m.betas[-1][0])})
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- shocks
SHOCK_WINDOWS = [(2019, 2020), (2020, 2021), (2021, 2022), (2021, 2023), (2019, 2024)]


def shock_response(C):
    """Change in kWh/meter over shock windows by start-year consumption
    quintile and by deprivation quintile (2019 indices, within-country rank
    percentile; Q1 = most deprived). Aggregate change is meter-weighted
    (sum kWh / sum meters in each group); the mean LSOA log change is given too."""
    cov = pd.read_parquet(PROCESSED / "covariates_lsoa21.parquet",
                          columns=["lsoa21cd", "dep2019_pct_within_country"]).set_index("lsoa21cd")
    dep = cov.loc[C["geo"]["lsoa21cd"], "dep2019_pct_within_country"].values
    depq = np.minimum((dep * 5).astype(int), 4) + 1
    rows = []
    for a, b in SHOCK_WINDOWS:
        ya, yb = C["y"][a].values, C["y"][b].values
        ma, mb = C["meters"][a].values, C["meters"][b].values
        nat = (yb * mb).sum() / mb.sum() / ((ya * ma).sum() / ma.sum()) - 1
        startq = pd.qcut(C["rel"][a].values, 5, labels=False) + 1
        for grouping, grp in (("start-year consumption quintile", startq), ("deprivation quintile (2019)", depq)):
            for g in range(1, 6):
                k = grp == g
                agg = (yb[k] * mb[k]).sum() / mb[k].sum() / ((ya[k] * ma[k]).sum() / ma[k].sum()) - 1
                rows.append({"from_year": a, "to_year": b, "grouping": grouping, "group": g, "n": int(k.sum()),
                             "pct_change_aggregate": 100 * agg, "mean_log_change": float(np.mean(np.log(yb[k] / ya[k]))),
                             "start_mean_kwh": float((ya[k] * ma[k]).sum() / ma[k].sum()),
                             "pct_change_national": 100 * nat})
    return pd.DataFrame(rows)


def shock_figure(S):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), sharey=True)
    for ax, grouping, ttl, xl in ((axes[0], "start-year consumption quintile", "By 2021 consumption quintile",
                                   ["Q1\nlowest", "Q2", "Q3", "Q4", "Q5\nhighest"]),
                                  (axes[1], "deprivation quintile (2019)", "By deprivation quintile",
                                   ["Q1\nmost", "Q2", "Q3", "Q4", "Q5\nleast"])):
        d = S[(S.grouping == grouping) & (S.from_year == 2021) & (S.to_year == 2022)]
        ax.bar(d.group, d.pct_change_aggregate, color=CAT[0], width=0.65)
        ax.axhline(d.pct_change_national.iloc[0], color=INK, lw=1, ls="--", label="E&W change")
        ax.set_xticks(range(1, 6), xl)
        ax.set_title(ttl, fontsize=11)
        for x, v in zip(d.group, d.pct_change_aggregate):
            ax.text(x, v + 0.25, f"{v:.1f}", ha="center", va="bottom", fontsize=9, color="white",
                    fontweight="bold")
    axes[0].set_ylabel("Change in kWh per meter, 2021→22 (%)")
    axes[1].legend(loc="upper right", bbox_to_anchor=(1.0, 1.22), fontsize=9)
    fig.suptitle("Price-crisis fall in electricity use", x=0.02, ha="left", fontweight="bold")
    fig.tight_layout()
    return savefig(fig, "ineq_shock_2022_by_group")


# --------------------------------------------------------------------------- figures
def figures(T, G, Dt, D):
    files = []
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 6.4), sharex=True)
    ax = axes[0]
    ax.plot(G.year, G.gini, color=CAT[0], marker="o", ms=3.5, label="Gini (LSOAs)")
    ax.plot(T.year, T.theil_T, color=CAT[1], marker="s", ms=3.5, ls="--", label="Theil T (meter-weighted)")
    ax.set_ylabel("Index")
    ax.set_title("Inequality between LSOAs")
    ax.legend(fontsize=10)
    ax = axes[1]
    ax.stackplot(T.year, T.share_between_region * 100, T.share_between_lad_within_region * 100,
                 T.share_within_lad * 100, colors=[CAT[0], CAT[2], "#d9d8d1"],
                 labels=["Between regions", "Between LADs within regions", "Within LADs"],
                 edgecolor="white", linewidth=0.8)
    ax.set_ylabel("Share of Theil T (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Theil decomposition")
    ax.legend(loc="lower left", fontsize=9, framealpha=0.9, frameon=True)
    for a in axes:
        a.axvline(2014.5, color=MUTED, lw=0.8)
        a.set_xticks(range(2010, 2025, 2))
    fig.tight_layout()
    files += savefig(fig, "ineq_theil_gini")

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    d = D[(D["sample"] == "all") & (D.year == 2024)]
    ax.scatter(d.mean_kwh, d.cv, s=14, color=CAT[0], alpha=0.75, edgecolor="white", linewidth=0.4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    b = np.polyfit(np.log(d.mean_kwh), np.log(d.cv), 1)
    xx = np.linspace(d.mean_kwh.min(), d.mean_kwh.max(), 50)
    ax.plot(xx, np.exp(np.polyval(b, np.log(xx))), color=CAT[1], lw=1.5)
    ax.set_xlabel("LAD mean (kWh per meter, log scale)")
    ax.set_ylabel("Within-LAD CV of LSOA means (log)")
    ax.set_title("Within-LAD dispersion vs LAD mean, 2024")
    for lab in d.nlargest(3, "cv").itertuples():
        ax.annotate(lab.lad_cd, (lab.mean_kwh, lab.cv), fontsize=8, color=INK2, xytext=(3, 2),
                    textcoords="offset points")
    ax.grid(True, which="major", axis="both")
    fig.tight_layout()
    files += savefig(fig, "ineq_cv_vs_mean_2024")

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for (samp, disp), c, ls, lab in ((("all", "cv"), CAT[0], "-", "CV, all LSOAs"),
                                     (("excl_flagged", "cv"), CAT[1], "--", "CV, excl. flagged"),
                                     (("all", "qcd"), CAT[2], ":", "Quartile coef., all")):
        s = Dt[(Dt["sample"] == samp) & (Dt.dispersion == disp) & (Dt.level == "mean_kwh")]
        ax.plot(s.year, s.spearman_rho, color=c, ls=ls, marker="o", ms=3, label=lab)
        if disp == "cv" and samp == "all":
            ax.fill_between(s.year, s.rho_ci_lo, s.rho_ci_hi, color=c, alpha=0.15, lw=0)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_ylabel("Spearman ρ (dispersion, mean)")
    ax.set_title("Within-LAD dispersion vs LAD mean")
    ax.legend(fontsize=9)
    ax.set_xticks(range(2010, 2025, 2))
    fig.tight_layout()
    files += savefig(fig, "ineq_cv_mean_rho_series")
    return files


def run():
    set_style()
    log("Inequality: loading")
    reg = lad_region().set_index("lad_cd")["region"]
    C = core("elec")
    Cx = core("elec", exclude_flagged=True)
    wq = lsoa_weights("queen", C["mask"], C["tag"])
    T = theil_tables(C, reg)
    write_table(T, "ineq_theil")
    log("Inequality: Gini / spatial Gini")
    G = gini_tables(C, wq)
    write_table(G, "ineq_gini_spatial")
    log("Inequality: within-LAD dispersion")
    D = lad_dispersion(C, Cx, reg)
    write_table(D, "ineq_lad_dispersion")
    lad_geo, w_lad, islands = lad_weights()
    Dt = dispersion_tests(D, lad_geo, w_lad)
    write_table(Dt, "ineq_dispersion_vs_mean")
    Sh = shock_response(C)
    write_table(Sh, "ineq_shock_response")
    files = figures(T, G, Dt, D) + shock_figure(Sh)
    save_json({"lad_islands_attached": islands, "lad_mean_max_rel_diff_vs_panel": D.attrs["max_rel_diff_vs_lad_panel"],
               "figures": files}, "ineq_meta")
    log("Inequality done")


if __name__ == "__main__":
    run()
