"""Section 1: exploratory spatial data analysis of LSOA mean kWh/meter.

Outputs (results/dynamics/):
  tables/esda_summary_stats.*        distribution summary per year
  tables/esda_kde_bandwidth.*        KDE bandwidths and modes
  tables/esda_fj_breaks.*            pooled sampled Fisher-Jenks breaks
  tables/esda_moran_global.*         global Moran's I (queen/KNN8, level/log, with/without S-X LSOAs)
  tables/esda_lisa_counts.*          LISA cluster counts per year, raw and FDR
  tables/esda_lisa_knn_agreement.*   queen vs KNN8 LISA agreement
  cache/lisa_queen.parquet           LSOA x year LISA (Is, quadrant, p_sim, cluster labels)
  figures/esda_*                     KDEs, maps, Moran series, LISA maps
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import esda  # noqa: E402
import mapclassify  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import stats  # noqa: E402

from a_common import *  # noqa: E402,F401,F403

LISA_PERMS = 9999
QUAD = {1: "HH", 2: "LH", 3: "LL", 4: "HL"}
MAP_YEARS = [2010, 2015, 2020, 2022, 2024]


def summary_stats(C):
    y, m = C["y"], C["meters"]
    rows = []
    for t in YEARS:
        v = y[t].values
        q = np.percentile(v, [5, 10, 25, 50, 75, 90, 95])
        rows.append({
            "year": t, "n": v.size, "mean": v.mean(), "mean_meter_weighted": (v * m[t]).sum() / m[t].sum(),
            "median": q[3], "sd": v.std(ddof=1), "cv": v.std(ddof=1) / v.mean(),
            "skewness": stats.skew(v), "kurtosis_excess": stats.kurtosis(v),
            "skewness_log": stats.skew(np.log(v)), "min": v.min(), "p5": q[0], "p10": q[1], "p25": q[2],
            "p75": q[4], "p90": q[5], "p95": q[6], "max": v.max(),
            "iqr_over_median": (q[4] - q[2]) / q[3], "p90_p10": q[5] / q[1],
            "bowley_skewness": (q[4] + q[2] - 2 * q[3]) / (q[4] - q[2]),
        })
    return pd.DataFrame(rows)


def kde_figure(C):
    rel = C["rel"]
    grid = np.linspace(0.2, 3.0, 800)
    rows, curves = [], {}
    for t in YEARS:
        v = rel[t].values
        k = stats.gaussian_kde(np.log(v), bw_method="silverman")
        # density of rel via change of variable from log scale
        d = k(np.log(grid)) / grid
        curves[t] = d
        peaks = np.flatnonzero((d[1:-1] > d[:-2]) & (d[1:-1] > d[2:])) + 1
        peaks = peaks[d[peaks] > 0.02 * d.max()]
        rows.append({"year": t, "kernel": "gaussian on log(relative)", "bandwidth_rule": "Silverman",
                     "bandwidth_log_scale": float(k.factor * np.std(np.log(v), ddof=1)),
                     "n_modes": int(peaks.size), "mode_relative": float(grid[peaks[np.argmax(d[peaks])]])})
    pd.DataFrame({"relative": grid, **{str(t): curves[t] for t in YEARS}}).to_parquet(TAB / "esda_kde_curves.parquet")

    fig, axes = plt.subplots(2, 1, figsize=(5.2, 6.4), sharex=True)
    sel = [2010, 2015, 2019]
    for ax, years, ttl in ((axes[0], [2010, 2015, 2019], "Before the shocks"),
                           (axes[1], [2019, 2020, 2022, 2024], "Shock years")):
        cols = {2010: SEQ[1], 2015: SEQ[2], 2019: SEQ[4], 2020: CAT[1], 2022: CAT[7], 2024: CAT[2]}
        for t in years:
            ax.plot(grid, curves[t], color=cols[t], lw=1.8, label=str(t))
        ax.set_title(ttl)
        ax.set_ylabel("Density")
        ax.legend(loc="upper right")
        ax.axvline(1, color=MUTED, lw=0.8, ls=":")
        ax.set_xlim(0.3, 2.5)
    axes[1].set_xlabel("Relative consumption (LSOA ÷ E&W mean)")
    fig.tight_layout()
    files = savefig(fig, "esda_kde_relative")
    del sel
    return pd.DataFrame(rows), files


def fj_maps(C):
    """Sampled Fisher-Jenks breaks on the pooled 2015-2024 relative values
    (fixed across years so maps are comparable), applied to selected years."""
    g, reg = map_layers()
    g = g[C["mask"]].reset_index(drop=True)
    pooled = C["rel"][CORE_YEARS].values.ravel()
    np.random.seed(SEED)
    fj = mapclassify.FisherJenksSampled(pooled, k=7, pct=0.05, truncate=True)
    bins = fj.bins.copy()
    bins[-1] = np.inf
    breaks = pd.DataFrame({"class": range(1, 8), "upper_bound_relative": fj.bins,
                           "count_pooled_2015_24": fj.counts})
    files = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 6.2))
    for ax, t in zip(axes, [2015, 2024]):
        cls = np.searchsorted(bins, C["rel"][t].values, side="left")
        draw_map(ax, g, [SEQ7[c] for c in cls], reg, title=str(t))
    labels = _bin_labels(fj.bins)
    fig.legend(handles=[Patch(color=SEQ7[i], label=labels[i]) for i in range(7)], loc="lower center",
               ncol=4, title="Relative consumption (Fisher–Jenks, pooled 2015–24)", fontsize=10,
               title_fontsize=10)
    fig.subplots_adjust(bottom=0.16, wspace=0.02)
    files += savefig(fig, "esda_map_fj_2015_2024")
    for t in MAP_YEARS:
        fig, ax = plt.subplots(figsize=(5.2, 6.2))
        cls = np.searchsorted(bins, C["rel"][t].values, side="left")
        draw_map(ax, g, [SEQ7[c] for c in cls], reg, title=f"Domestic electricity per meter, {t}")
        ax.legend(handles=[Patch(color=SEQ7[i], label=labels[i]) for i in range(7)], loc="upper left",
                  title="× E&W mean", fontsize=9, title_fontsize=9, bbox_to_anchor=(0.0, 0.93))
        files += savefig(fig, f"esda_map_fj_{t}")
    return breaks, files


def _bin_labels(b):
    lo = [0] + list(b[:-1])
    lab = [f"{a:.2f}–{c:.2f}" for a, c in zip(lo, b)]
    lab[0] = f"< {b[0]:.2f}"
    lab[-1] = f"> {b[-2]:.2f}"
    return lab


def moran_table(C, w_q, w_k, w_q_nosx, C_nosx):
    rows = []
    for t in YEARS:
        for wname, w, CC in (("queen", w_q, C), ("knn8", w_k, C), ("queen_excl_SX", w_q_nosx, C_nosx)):
            for var in ("level", "log"):
                v = CC["y"][t].values
                v = np.log(v) if var == "log" else v
                if wname == "knn8" and var == "log":
                    continue
                np.random.seed(SEED + t)
                m = esda.Moran(v, w, permutations=PERMS)
                rows.append({"year": t, "weights": wname, "variable": var, "n": v.size, "I": m.I, "EI": m.EI,
                             "z_sim": m.z_sim, "p_sim": m.p_sim, "z_norm": m.z_norm, "p_norm": m.p_norm})
    return pd.DataFrame(rows)


def lisa_all(C, w, name):
    cache = CACHE / f"lisa_{name}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    out = []
    ids = C["geo"]["lsoa21cd"].values
    for t in YEARS:
        v = C["rel"][t].values
        lm = esda.Moran_Local(v, w, permutations=LISA_PERMS, n_jobs=-1, seed=SEED + t)
        thr = esda.fdr(lm.p_sim, 0.05)
        q = pd.Series(lm.q).map(QUAD).values
        sig_raw = lm.p_sim < 0.05
        sig_fdr = lm.p_sim <= thr
        out.append(pd.DataFrame({"lsoa21cd": ids, "year": t, "Is": lm.Is, "quad": q, "p_sim": lm.p_sim,
                                 "fdr_threshold": thr,
                                 "cluster_raw": np.where(sig_raw, q, "ns"),
                                 "cluster_fdr": np.where(sig_fdr, q, "ns")}))
        log(f"  LISA {name} {t}: FDR threshold {thr:.5f}, {sig_fdr.sum()} significant")
    df = pd.concat(out, ignore_index=True)
    df.to_parquet(cache, index=False)
    return df


def lisa_counts(L):
    a = L.groupby(["year", "cluster_raw"]).size().unstack(fill_value=0).add_prefix("raw_")
    b = L.groupby(["year", "cluster_fdr"]).size().unstack(fill_value=0).add_prefix("fdr_")
    thr = L.groupby("year")["fdr_threshold"].first()
    return pd.concat([a, b], axis=1).assign(fdr_threshold=thr).reset_index()


def lisa_maps(C, L):
    g, reg = map_layers()
    g = g[C["mask"]].reset_index(drop=True)
    files = []
    order = ["HH", "LL", "HL", "LH", "ns"]
    names = {"HH": "High–High", "LL": "Low–Low", "HL": "High–Low", "LH": "Low–High", "ns": "Not significant"}
    for t in (2015, 2024, 2010, 2022):
        lab = L.loc[L["year"] == t, "cluster_fdr"].values
        fig, ax = plt.subplots(figsize=(5.2, 6.2))
        draw_map(ax, g, [LISA_COLORS[c] for c in lab], reg, title=f"LISA clusters (FDR 5%), {t}")
        ax.legend(handles=[Patch(color=LISA_COLORS[c], label=names[c]) for c in order], loc="upper left",
                  fontsize=9, bbox_to_anchor=(0.0, 0.93))
        files += savefig(fig, f"esda_lisa_map_{t}")
    return files


def moran_figure(M):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    styles = {("queen", "level"): (CAT[0], "-", "Queen, kWh"), ("knn8", "level"): (CAT[1], "--", "KNN-8, kWh"),
              ("queen", "log"): (CAT[2], ":", "Queen, log kWh"),
              ("queen_excl_SX", "level"): (MUTED, "-.", "Queen, excl. split LSOAs")}
    for (wn, var), (c, ls, lab) in styles.items():
        s = M[(M.weights == wn) & (M.variable == var)]
        ax.plot(s.year, s.I, color=c, ls=ls, label=lab, marker="o", ms=3.5)
    ax.axvline(2014.5, color=MUTED, lw=0.8)
    ax.text(2014.6, ax.get_ylim()[0], " geography/system break", color=INK2, fontsize=9, va="bottom")
    ax.set_ylabel("Global Moran's I")
    ax.set_title("Spatial autocorrelation of LSOA consumption")
    ax.legend(fontsize=9, loc="center left", bbox_to_anchor=(0, 0.42))
    ax.set_xticks(range(2010, 2025, 2))
    fig.tight_layout()
    return savefig(fig, "esda_moran_series")


def knn_agreement(C, Lq, w_k):
    rows = []
    for t in (2015, 2024):
        v = C["rel"][t].values
        lm = esda.Moran_Local(v, w_k, permutations=LISA_PERMS, n_jobs=-1, seed=SEED + t)
        thr = esda.fdr(lm.p_sim, 0.05)
        k = np.where(lm.p_sim <= thr, pd.Series(lm.q).map(QUAD).values, "ns")
        q = Lq.loc[Lq.year == t, "cluster_fdr"].values
        po = (k == q).mean()
        cats = np.unique(np.r_[k, q])
        pe = sum((k == c).mean() * (q == c).mean() for c in cats)
        rows.append({"year": t, "share_same_label": po, "cohen_kappa": (po - pe) / (1 - pe),
                     "n_sig_queen": int((q != "ns").sum()), "n_sig_knn8": int((k != "ns").sum()),
                     "HH_queen": int((q == "HH").sum()), "HH_knn8": int((k == "HH").sum()),
                     "LL_queen": int((q == "LL").sum()), "LL_knn8": int((k == "LL").sum())})
    return pd.DataFrame(rows)


def run():
    set_style()
    log("ESDA: loading data and weights")
    C = core("elec")
    w_q = lsoa_weights("queen", C["mask"], C["tag"])
    w_k = lsoa_weights("knn8", C["mask"], C["tag"])
    Cn = core("elec", exclude_sx=True)
    w_qn = lsoa_weights("queen", Cn["mask"], Cn["tag"])
    meta = {"n_lsoa_total": 35672, "n_balanced": int(C["mask"].sum()),
            "n_dropped_incomplete": int(35672 - C["mask"].sum()), "n_excl_SX": int(Cn["mask"].sum()),
            "queen_mean_neighbours": w_q.mean_neighbors, "queen_min_neighbours": int(w_q.min_neighbors),
            "queen_components": int(w_q.n_components), "lisa_permutations": LISA_PERMS,
            "global_permutations": PERMS}
    files = []
    S = summary_stats(C)
    write_table(S, "esda_summary_stats")
    K, f = kde_figure(C)
    write_table(K, "esda_kde_bandwidth")
    files += f
    B, f = fj_maps(C)
    write_table(B, "esda_fj_breaks")
    files += f
    log("ESDA: global Moran")
    M = moran_table(C, w_q, w_k, w_qn, Cn)
    write_table(M, "esda_moran_global")
    files += moran_figure(M)
    log("ESDA: LISA (queen, 9,999 permutations)")
    L = lisa_all(C, w_q, "queen")
    write_table(lisa_counts(L), "esda_lisa_counts")
    files += lisa_maps(C, L)
    write_table(knn_agreement(C, L, w_k), "esda_lisa_knn_agreement")
    save_json({**meta, "figures": files}, "esda_meta")
    log("ESDA done")


if __name__ == "__main__":
    run()
