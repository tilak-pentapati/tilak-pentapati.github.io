"""Section 4: place mobility and rank concordance.

* Van Kerm (2004) decomposition of mobility into growth, dispersion and
  exchange components, as applied to places by Rey and co-authors. Mobility is
  the mean absolute change across places, measured in kWh/meter and in log
  points. Each component is a Shapley value over the three factors, using
  counterfactual distributions that switch on: growth (the mean of the end
  year), dispersion (the end-year shape, i.e. quantile function divided by its
  mean) and exchange (end-year ranks). Exact: components sum to the total.
* Kendall's tau (tau-b) between years; spatial tau (Rey 2004 form: concordance
  among contiguous pairs) with permutation inference; local neighbour tau
  (LIMA) with conditional permutation inference for LSOAs (vectorised) and
  LADs (giddy, validated).

Outputs: tables/mobility_*, tables/rank_*, cache/lima_lsoa_*.parquet, figures/mobility_*, figures/rank_*.
"""
from __future__ import annotations

import sys
from itertools import permutations as iperms
from math import factorial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import esda  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from giddy.rank import SpatialTau, Tau_Local_Neighbor  # noqa: E402
from libpysal import weights  # noqa: E402
from scipy import stats  # noqa: E402

from a_common import *  # noqa: E402,F401,F403

FACTORS = ("growth", "dispersion", "exchange")
WINDOWS_MULTI = [(2010, 2015), (2015, 2019), (2019, 2024), (2015, 2024), (2010, 2024), (2019, 2022)]


# --------------------------------------------------------------------------- Van Kerm
def _ranks(v):
    r = np.empty(v.size, dtype=np.int64)
    r[np.argsort(v, kind="stable")] = np.arange(v.size)
    return r


def counterfactual(x, y, on: frozenset):
    mx, my = x.mean(), y.mean()
    shape = np.sort(y) / my if "dispersion" in on else np.sort(x) / mx
    mean = my if "growth" in on else mx
    rank = _ranks(y) if "exchange" in on else _ranks(x)
    return mean * shape[rank]


def vankerm(x, y, index="abs"):
    """Shapley decomposition of M(x, y) = mean |f(y) - f(x)|."""
    f = (lambda v: v) if index == "abs" else np.log

    def M(z):
        return float(np.mean(np.abs(f(z) - f(x))))

    vals = {}
    for k in range(4):
        for S in __import__("itertools").combinations(FACTORS, k):
            vals[frozenset(S)] = M(counterfactual(x, y, frozenset(S)))
    out = {}
    for fac in FACTORS:
        tot = 0.0
        for order in iperms(FACTORS):
            pre = frozenset(order[:order.index(fac)])
            tot += vals[pre | {fac}] - vals[pre]
        out[fac] = tot / factorial(3)
    out["total"] = vals[frozenset(FACTORS)]
    assert abs(out["total"] - M(y)) < 1e-9 and abs(sum(out[f] for f in FACTORS) - out["total"]) < 1e-9
    # sequential (growth -> dispersion -> exchange) for reference
    out["seq_exchange"] = vals[frozenset(FACTORS)] - vals[frozenset({"growth", "dispersion"})]
    return out


def mobility_table(panels):
    rows = []
    for scope, Y in panels.items():
        pairs = [(t, t + 1) for t in YEARS[:-1]] + WINDOWS_MULTI
        for a, b in pairs:
            x, y = Y[a].values, Y[b].values
            for idx in ("abs", "log"):
                d = vankerm(x, y, idx)
                rows.append({"scope": scope, "from_year": a, "to_year": b, "span": b - a, "index": idx,
                             "total": d["total"], **{f: d[f] for f in FACTORS},
                             **{f"share_{f}": d[f] / d["total"] for f in FACTORS},
                             "exchange_sequential": d["seq_exchange"],
                             "mean_from": x.mean(), "mean_to": y.mean()})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- rank concordance
def spatial_tau(x, y, i, j, perms, rng):
    """Spatial tau over contiguous pairs (i<j), as giddy.rank.SpatialTau."""
    keep = i < j
    a, b = i[keep], j[keep]

    def calc(xx, yy):
        dxdy = (xx[a] - xx[b]) * (yy[a] - yy[b])
        n1 = np.sum((dxdy != 0) | (xx[a] != xx[b]))
        n2 = np.sum((dxdy != 0) | (yy[a] != yy[b]))
        s = np.sum(np.sign(dxdy))
        return s / np.sqrt(n1 * n2), int((dxdy > 0).sum()), int((dxdy < 0).sum())

    tau_s, c, d = calc(x, y)
    sims = np.empty(perms)
    for k in range(perms):
        p = rng.permutation(x.size)
        sims[k] = calc(x[p], y[p])[0]
    larger = (sims >= tau_s).sum()
    larger = min(larger, perms - larger)
    out = {"tau_spatial": tau_s, "concordant_pairs": c, "discordant_pairs": d, "n_pairs": int(a.size),
           "tau_spatial_perm_mean": sims.mean(), "tau_spatial_p": (larger + 1) / (perms + 1)}
    out.update(matched_null(x, y, a, b, tau_s, rng))
    return out


def matched_null(x, y, a, b, tau_s, rng, reps=199, nbins=20):
    """Value-matched null for spatial tau. Neighbours have similar start values,
    and pairs that start close together swap ranks more easily whatever their
    location. Random (mostly non-neighbour) pairs are drawn and reweighted to
    the neighbour pairs' distribution of |x_i - x_j| (20 quantile bins); each
    replicate has as many pairs as there are neighbour pairs. If tau_spatial is
    below this distribution, neighbours re-rank against each other more than
    equally similar places elsewhere do."""
    n = x.size
    dn = np.abs(x[a] - x[b])
    edges = np.quantile(dn, np.linspace(0, 1, nbins + 1))
    edges[0], edges[-1] = 0, np.inf
    share = np.bincount(np.clip(np.searchsorted(edges, dn, side="right") - 1, 0, nbins - 1), minlength=nbins) / dn.size
    m = a.size
    need = np.maximum(np.round(share * m).astype(int), 1)
    taus = np.empty(reps)
    for r in range(reps):
        # oversample random pairs, then keep the first need[k] pairs in each bin
        i = rng.integers(0, n, 4 * m)
        j = rng.integers(0, n, 4 * m)
        ok = i != j
        i, j = i[ok], j[ok]
        bi = np.clip(np.searchsorted(edges, np.abs(x[i] - x[j]), side="right") - 1, 0, nbins - 1)
        o = np.argsort(bi, kind="stable")
        bs = bi[o]
        start = np.searchsorted(bs, np.arange(nbins))
        rank_in_bin = np.arange(bs.size) - start[bs]
        keep = o[rank_in_bin < need[bs]]
        sgn = np.sign((x[i[keep]] - x[j[keep]]) * (y[i[keep]] - y[j[keep]]))
        taus[r] = sgn.sum() / (sgn != 0).sum()
    return {"tau_matched_null_mean": taus.mean(), "tau_matched_null_lo": np.percentile(taus, 2.5),
            "tau_matched_null_hi": np.percentile(taus, 97.5),
            "tau_matched_p": ((np.sum(taus <= tau_s) if tau_s < taus.mean() else np.sum(taus >= tau_s)) + 1)
            / (reps + 1) * 2}


def local_tau(x, y, i, j, n, perms, rng, chunk=8):
    """Local neighbour tau (LIMA, Tau_Local_Neighbor): mean concordance sign of
    each unit with its contiguous neighbours; conditional permutation inference
    draws the same number of other units at random (with replacement, which is
    immaterial at n=35k)."""
    sgn = np.sign((x[i] - x[j]) * (y[i] - y[j]))
    card = np.bincount(i, minlength=n)
    tau = np.bincount(i, sgn, minlength=n) / card
    kmax = card.max()
    mask = np.arange(kmax)[None, :] < card[:, None]
    larger = np.zeros(n)
    smaller = np.zeros(n)
    for s in range(0, perms, chunk):
        m = min(chunk, perms - s)
        r = rng.integers(0, n - 1, size=(m, n, kmax))
        r += r >= np.arange(n)[None, :, None]           # skip self
        sg = np.sign((x[None, :, None] - x[r]) * (y[None, :, None] - y[r]))
        sim = (sg * mask[None]).sum(axis=2) / card[None]
        larger += (sim >= tau[None]).sum(axis=0)
        smaller += (sim <= tau[None]).sum(axis=0)
    p = (np.minimum(larger, smaller) + 1) / (perms + 1)
    return tau, p


def rank_tables(C, w, Dl, wl):
    rng = np.random.default_rng(SEED)
    i, j = edge_arrays(w)
    li, lj = edge_arrays(wl)
    rows = []
    for scope, Y, ii, jj in (("LSOA", C["rel"], i, j), ("LAD", Dl["rel"], li, lj)):
        pairs = [(t, t + 1) for t in YEARS[:-1]] + WINDOWS_MULTI
        for a, b in pairs:
            x, y = Y[a].values, Y[b].values
            kt = stats.kendalltau(x, y)
            sp = spatial_tau(x, y, ii, jj, PERMS, rng)
            rows.append({"scope": scope, "from_year": a, "to_year": b, "span": b - a, "kendall_tau_b": kt.statistic,
                         "kendall_p": kt.pvalue, "spearman": stats.spearmanr(x, y).statistic, **sp})
    T = pd.DataFrame(rows)
    # validate the vectorised spatial tau against giddy on the LAD panel
    wb = weights.W(dict(wl.neighbors), silence_warnings=True)
    ref = SpatialTau(Dl["rel"][2015].values, Dl["rel"][2024].values, wb, permutations=0)
    mine = T[(T.scope == "LAD") & (T.from_year == 2015) & (T.to_year == 2024)].iloc[0]
    assert abs(ref.tau_spatial - mine.tau_spatial) < 1e-12, (ref.tau_spatial, mine.tau_spatial)
    return T


def lima(C, w, Dl, wl):
    rng = np.random.default_rng(SEED + 1)
    i, j = edge_arrays(w)
    out, summ = [], []
    for a, b in ((2015, 2024), (2015, 2019), (2019, 2024)):
        x, y = C["rel"][a].values, C["rel"][b].values
        tau, p = local_tau(x, y, i, j, x.size, PERMS, rng)
        thr = esda.fdr(p, 0.05)
        out.append(pd.DataFrame({"lsoa21cd": C["geo"]["lsoa21cd"].values, "window": f"{a}-{b}", "tau_local": tau,
                                 "p_sim": p, "sig_fdr": p <= thr}))
        summ.append({"scope": "LSOA", "window": f"{a}-{b}", "mean_local_tau": tau.mean(),
                     "share_negative": float((tau < 0).mean()), "n_sig_raw": int((p < 0.05).sum()),
                     "n_sig_fdr": int((p <= thr).sum()), "fdr_threshold": thr,
                     "n_sig_fdr_negative": int(((p <= thr) & (tau < 0)).sum())})
    L = pd.concat(out, ignore_index=True)
    L.to_parquet(CACHE / "lima_lsoa.parquet", index=False)
    # LADs with giddy (and a check that the vectorised statistic matches it)
    wb = weights.W(dict(wl.neighbors), silence_warnings=True)
    np.random.seed(SEED)
    x, y = Dl["rel"][2015].values, Dl["rel"][2024].values
    tl = Tau_Local_Neighbor(x, y, wb, permutations=PERMS)
    li, lj = edge_arrays(wl)
    mine, _ = local_tau(x, y, li, lj, x.size, 9, rng)
    assert np.allclose(mine, tl.tau_ln), "local tau mismatch"
    thr = esda.fdr(tl.tau_ln_pvalues, 0.05)
    lad_out = pd.DataFrame({"lad_cd": Dl["geo"]["lad_cd"].values, "lad_nm": Dl["geo"]["lad_nm"].values,
                            "tau_local": tl.tau_ln, "p_sim": tl.tau_ln_pvalues,
                            "sig_raw": tl.tau_ln_pvalues < 0.05, "sig_fdr": tl.tau_ln_pvalues <= thr})
    summ.append({"scope": "LAD", "window": "2015-2024", "mean_local_tau": tl.tau_ln.mean(),
                 "share_negative": float((tl.tau_ln < 0).mean()), "n_sig_raw": int((tl.tau_ln_pvalues < 0.05).sum()),
                 "n_sig_fdr": int((tl.tau_ln_pvalues <= thr).sum()), "fdr_threshold": thr,
                 "n_sig_fdr_negative": int(((tl.tau_ln_pvalues <= thr) & (tl.tau_ln < 0)).sum())})
    return L, lad_out, pd.DataFrame(summ)


# --------------------------------------------------------------------------- figures
def figures(M, T, C, L):
    files = []
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 6.6), sharex=True)
    for ax, scope in zip(axes, ("LSOA", "LAD")):
        d = M[(M.scope == scope) & (M.span == 1) & (M["index"] == "log")].sort_values("to_year")
        bottom = np.zeros(len(d))
        for f, c in zip(FACTORS, (CAT[1], CAT[0], CAT[2])):
            v = d[f].values * 100
            ax.bar(d.to_year, v, bottom=bottom, color=c, width=0.75, label=f.capitalize(), edgecolor="white",
                   linewidth=0.8)
            bottom += v
        ax.set_ylabel("Mean |Δ log| (×100)")
        ax.set_title(f"{scope}s: annual mobility decomposition")
        ax.axhline(0, color=AXIS, lw=0.8)
    axes[0].legend(fontsize=9, ncol=3, loc="upper left")
    axes[1].set_xticks(range(2011, 2025, 2))
    axes[1].set_xlabel("Year (end of one-year window)")
    fig.tight_layout()
    files += savefig(fig, "mobility_vankerm_annual")

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for scope, c, mk in (("LSOA", CAT[0], "o"), ("LAD", CAT[1], "s")):
        d = T[(T.scope == scope) & (T.span == 1)].sort_values("to_year")
        ax.plot(d.to_year, d.kendall_tau_b, color=c, marker=mk, ms=3.5, label=f"{scope}: Kendall τ")
        ax.plot(d.to_year, d.tau_spatial, color=c, marker=mk, ms=3.5, ls="--", label=f"{scope}: spatial τ")
    ax.set_ylabel("Rank concordance, year t−1 to t")
    ax.set_title("Rank stability")
    ax.set_xticks(range(2011, 2025, 2))
    ax.legend(fontsize=9, ncol=2, loc="lower left")
    fig.tight_layout()
    files += savefig(fig, "rank_tau_annual")

    g, reg = map_layers()
    g = g[C["mask"]].reset_index(drop=True)
    d = L[L.window == "2015-2024"]
    # No LSOA survives FDR correction (see rank_lima_summary), so the map shows
    # the local statistic itself rather than significance.
    t = d.tau_local.values
    cls = np.select([t < 0, t < 0.5, t < 1.0], [0, 1, 2], 3)
    pal = [RED, "#e1e0d9", "#86b6ef", "#1c5cab"]
    labs = ["τ < 0: re-ranked against neighbours", "0 ≤ τ < 0.5", "0.5 ≤ τ < 1", "τ = 1: order kept with all"]
    fig, ax = plt.subplots(figsize=(5.2, 6.2))
    draw_map(ax, g, [pal[c] for c in cls], reg, title="Local rank concordance, 2015→2024")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=pal[k], label=labs[k]) for k in range(4)],
              loc="upper left", fontsize=8.5, bbox_to_anchor=(0.0, 0.93))
    files += savefig(fig, "rank_lima_map_2015_2024")
    return files


def run():
    set_style()
    log("Mobility: data")
    C = core("elec")
    w = lsoa_weights("queen", C["mask"], C["tag"])
    lp = pd.read_parquet(PROCESSED / "elec_lad_panel.parquet")
    lad, wl, _ = lad_weights()
    yl = lp.pivot(index="lad_cd", columns="year", values="mean_kwh").loc[lad["lad_cd"], YEARS]
    ml = lp.pivot(index="lad_cd", columns="year", values="meters").loc[lad["lad_cd"], YEARS]
    Dl = {"y": yl, "rel": relative(yl, ml), "geo": lad}
    Cx = core("elec", exclude_flagged=True)
    Cn = core("elec", exclude_sx=True)
    M = mobility_table({"LSOA": C["y"], "LAD": yl, "LSOA excl. flagged": Cx["y"], "LSOA excl. S/X": Cn["y"]})
    write_table(M, "mobility_vankerm")
    log("Mobility: rank concordance")
    T = rank_tables(C, w, Dl, wl)
    write_table(T, "rank_concordance")
    log("Mobility: LIMA")
    L, lad_lima, S = lima(C, w, Dl, wl)
    write_table(lad_lima, "rank_lima_lad_2015_2024")
    write_table(S, "rank_lima_summary")
    files = figures(M, T, C, L)
    save_json({"figures": files}, "mobility_meta")
    log("Mobility done")


if __name__ == "__main__":
    run()
