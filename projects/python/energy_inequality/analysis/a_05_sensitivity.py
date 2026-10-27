"""Section 5: sensitivity analyses.

* Excluding the 216 LSOAs of the balanced panel flagged as outliers in any
  year (``flag_outlier``): Moran's I, Gini, Theil shares, classic and spatial
  Markov summaries.
* Gas (weather-corrected) as a comparison fuel on LSOAs with a complete gas
  series, alongside electricity on exactly the same LSOAs: Moran's I, Gini,
  Theil between-LAD share, Markov persistence and Van Kerm decomposition.
* Scale: Theil decomposition with MSOAs as the grouping (within-MSOA share).
* Permutation-based (spatially randomised) p-value for the spatial Markov x2
  statistic, which does not assume independent transitions.

Outputs: tables/sens_*.{parquet,csv}, figures/sens_*.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import esda  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from giddy.markov import Spatial_Markov  # noqa: E402

from a_common import *  # noqa: E402,F401,F403
from a_02_inequality import gini, theil_decomp  # noqa: E402
from a_03_markov import K, classify, counts_window, row_p, shorrocks  # noqa: E402
from a_04_mobility import vankerm  # noqa: E402


def gas_core():
    """Electricity and gas on LSOAs with complete series for both fuels."""
    E = analysis_frame("elec")
    G = analysis_frame("gas")
    mask = E["complete"] & G["complete"]
    out = {}
    for fuel, A in (("elec", E), ("gas", G)):
        y, m = A["y"][mask], A["meters"][mask]
        out[fuel] = {"y": y, "meters": m, "rel": relative(y, m)}
    geo = E["geo"][mask].reset_index(drop=True)
    return out, mask, geo


def variant_stats(name, Y, M, R, lad, w, years=YEARS):
    rows = []
    for t in years:
        v = Y[t].values
        np.random.seed(SEED + t)
        mi = esda.Moran(v, w, permutations=PERMS)
        T, b, wi = theil_decomp(v * M[t].values, M[t].values, lad)
        rows.append({"variant": name, "year": t, "n": v.size, "moran_I": mi.I, "moran_p": mi.p_sim,
                     "gini": gini(v), "theil_T": T, "share_between_lad": b / T, "cv": v.std(ddof=1) / v.mean()})
    return rows


def markov_stats(name, R, cuts, w=None):
    cls = classify(R[YEARS].values, cuts)
    yi = {t: k for k, t in enumerate(YEARS)}
    out = []
    for lab, (a, b) in {"2015-24": (2015, 2024), "2010-14": (2010, 2014), "2015-19": (2015, 2019),
                        "2020-24": (2020, 2024)}.items():
        T = counts_window(cls, yi[a], yi[b])
        P = row_p(T)
        out.append({"variant": name, "window": lab, "shorrocks": shorrocks(P),
                    **{f"diag_{c + 1}": P[c, c] for c in range(K)}})
    if w is not None:
        sm = Spatial_Markov(R[CORE_YEARS].values, w, k=K, m=K, fixed=True, cutoffs=cuts)
        out[0].update({"spatial_Q": sm.Q, "spatial_Q_p": sm.Q_p_value, "spatial_LR": sm.LR,
                       "p_up_Q1_lowlag": sm.P[0][0, 1:].sum(), "p_up_Q1_highlag": sm.P[K - 1][0, 1:].sum()})
    return out


def run():
    set_style()
    rows, mrows, files = [], [], []
    # ---------------- outliers
    log("Sensitivity: excluding flagged LSOAs")
    C = core("elec")
    Cx = core("elec", exclude_flagged=True)
    w = lsoa_weights("queen", C["mask"], C["tag"])
    wx = lsoa_weights("queen", Cx["mask"], Cx["tag"])
    cuts = np.quantile(C["rel"][YEARS].values.ravel(), [0.2, 0.4, 0.6, 0.8])
    rows += variant_stats("elec: all balanced", C["y"], C["meters"], C["rel"], C["lad"], w)
    rows += variant_stats("elec: excl. flagged", Cx["y"], Cx["meters"], Cx["rel"], Cx["lad"], wx)
    mrows += markov_stats("elec: all balanced", C["rel"], cuts, w)
    mrows += markov_stats("elec: excl. flagged", Cx["rel"], cuts, wx)

    # ---------------- gas comparison
    log("Sensitivity: gas comparison")
    F, gmask, ggeo = gas_core()
    wg = lsoa_weights("queen", gmask, "gas_and_elec")
    lad = ggeo["lad_cd"].values
    for fuel in ("elec", "gas"):
        d = F[fuel]
        name = f"{fuel}: gas-complete LSOAs"
        rows += variant_stats(name, d["y"], d["meters"], d["rel"], lad, wg)
        c = np.quantile(d["rel"][YEARS].values.ravel(), [0.2, 0.4, 0.6, 0.8])
        mrows += markov_stats(name, d["rel"], c, wg)
    V = []
    for fuel in ("elec", "gas"):
        Y = F[fuel]["y"]
        for a, b in [(t, t + 1) for t in YEARS[:-1]] + [(2015, 2024), (2019, 2022)]:
            r = vankerm(Y[a].values, Y[b].values, "log")
            V.append({"fuel": fuel, "from_year": a, "to_year": b, **{k: r[k] for k in ("total", "growth",
                                                                                    "dispersion", "exchange")}})
    V = pd.DataFrame(V)
    write_table(V, "sens_gas_vankerm")

    # ---------------- MSOA scale for Theil
    log("Sensitivity: MSOA grouping")
    msoa = C["geo"]["msoa21cd"].values
    th = []
    for t in YEARS:
        p = C["meters"][t].values
        e = C["y"][t].values * p
        T, b_m, w_m = theil_decomp(e, p, pd.factorize(msoa)[0])
        _, b_l, _ = theil_decomp(e, p, pd.factorize(C["lad"])[0])
        th.append({"year": t, "theil_T": T, "share_between_msoa": b_m / T, "share_within_msoa": w_m / T,
                   "share_between_lad": b_l / T, "share_between_msoa_within_lad": (b_m - b_l) / T})
    write_table(pd.DataFrame(th), "sens_theil_msoa")

    # ---------------- spatially randomised x2 for spatial Markov
    log("Sensitivity: spatial Markov permutation x2 (99 permutations)")
    np.random.seed(SEED)
    sm = Spatial_Markov(C["rel"][CORE_YEARS].values, w, k=K, m=K, fixed=True, cutoffs=cuts, permutations=99)
    perm = {"x2": sm.x2, "x2_dof": sm.x2_dof, "x2_p_analytic": sm.x2_pvalue, "x2_p_perm": sm.x2_rpvalue,
            "x2_perm_max": float(np.max(sm.x2_realizations)), "x2_perm_mean": float(np.mean(sm.x2_realizations))}

    S = pd.DataFrame(rows)
    write_table(S, "sens_cross_section")
    write_table(pd.DataFrame(mrows), "sens_markov")

    fig, axes = plt.subplots(2, 1, figsize=(5.2, 6.2), sharex=True)
    spec = [("elec: all balanced", CAT[0], "-", "Electricity, all"),
            ("elec: excl. flagged", CAT[0], ":", "Electricity, excl. flagged"),
            ("gas: gas-complete LSOAs", CAT[1], "--", "Gas (weather-corrected)")]
    for name, c, ls, lab in spec:
        d = S[S.variant == name]
        axes[0].plot(d.year, d.moran_I, color=c, ls=ls, marker="o", ms=3, label=lab)
        axes[1].plot(d.year, d.gini, color=c, ls=ls, marker="o", ms=3, label=lab)
    axes[0].set_ylabel("Moran's I (queen)")
    axes[0].set_title("Spatial autocorrelation")
    axes[1].set_ylabel("Gini (LSOAs)")
    axes[1].set_title("Inequality")
    axes[0].legend(fontsize=9)
    axes[1].set_xticks(range(2010, 2025, 2))
    for a in axes:
        a.axvline(2014.5, color=MUTED, lw=0.8)
    fig.tight_layout()
    files += savefig(fig, "sens_elec_gas_outliers")
    save_json({"spatial_markov_permutation": perm, "n_excl_flagged": int(Cx["mask"].sum()),
               "n_gas_complete": int(gmask.sum()), "figures": files}, "sens_meta")
    log("Sensitivity done")


if __name__ == "__main__":
    run()
