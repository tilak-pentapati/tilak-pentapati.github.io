"""Section 3: distribution dynamics.

Classic and spatial Markov chains on quintiles of relative consumption
(LSOA or LAD mean / meter-weighted E&W mean in each year), LISA Markov chains
and directional LISA rose diagrams, with:

* fixed class bounds: pooled quintiles of relative consumption over 2010-2024
  (the same cut-offs for every window, so matrices are comparable);
* spatial homogeneity tests (Q and LR, Bickenbach & Bode 2003) via giddy;
* time homogeneity across 2010-14, 2015-19, 2020-24 and per annual transition,
  with the 2014->15 geography/processing break examined separately and with
  split (S/X) LSOAs excluded;
* a test of the first-order Markov property (transitions t->t+1 conditioned
  on the class at t-1);
* ergodic distributions, mean first passage times and sojourn times.

Outputs: tables/markov_*, tables/rose_*, cache/markov_*.npz, figures/markov_*, figures/rose_*.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
from giddy import ergodic  # noqa: E402
from giddy.directional import Rose  # noqa: E402
from giddy.markov import LISA_Markov, Markov, Spatial_Markov, homogeneity  # noqa: E402
from libpysal import weights  # noqa: E402

from a_common import *  # noqa: E402,F401,F403

K = K_CLASSES
QNAMES = {1: "HH", 2: "LH", 3: "LL", 4: "HL"}


# --------------------------------------------------------------------------- helpers
def classify(R: np.ndarray, cuts: np.ndarray) -> np.ndarray:
    return np.searchsorted(cuts, R, side="left")  # 0..K-1; cuts holds K-1 inner bounds


def counts_between(c0, c1, k=K):
    return np.bincount(c0 * k + c1, minlength=k * k).reshape(k, k).astype(float)


def counts_window(cls: np.ndarray, j0: int, j1: int, k=K):
    """Pooled transition counts over consecutive pairs of columns j0..j1."""
    T = np.zeros((k, k))
    for j in range(j0, j1):
        T += counts_between(cls[:, j], cls[:, j + 1], k)
    return T


def row_p(T):
    return T / T.sum(axis=1, keepdims=True)


def shorrocks(P):
    return (P.shape[0] - np.trace(P)) / (P.shape[0] - 1)


def chain_summary(T, label):
    P = row_p(T)
    s = ergodic.steady_state(P)
    F = ergodic.mfpt(P)
    soj = 1 / (1 - np.diag(P))
    return {"label": label, "n_transitions": int(T.sum()), "P": P, "T": T, "steady_state": np.ravel(s),
            "mfpt": np.asarray(F), "sojourn": soj, "shorrocks": shorrocks(P),
            "mean_diag": float(np.mean(np.diag(P))), "share_staying": float(np.trace(T) / T.sum())}


def hom(mats, names):
    h = homogeneity(list(mats), regime_names=list(names))
    P_h1 = h.p_h1
    P_h0 = h.p_h0
    maxdiff = float(np.max(np.abs(P_h1 - P_h0[None])))
    meandiag = [float(np.mean(np.diag(p))) for p in P_h1]
    return {"Q": h.Q, "Q_p": h.Q_p_value, "LR": h.LR, "LR_p": h.LR_p_value, "dof": h.dof,
            "max_abs_diff_p": maxdiff, "regime_mean_diag": meandiag}


def mat_rows(P, label, kind="P", names=None):
    names = names or CLASS_LABELS
    out = []
    for a in range(P.shape[0]):
        for b in range(P.shape[1]):
            out.append({"chain": label, "matrix": kind, "from": names[a], "to": names[b], "value": float(P[a, b])})
    return out


# --------------------------------------------------------------------------- LSOA chains
def lsoa_chains(C, w, cuts, label_prefix="lsoa"):
    R = C["rel"][YEARS].values
    cls = classify(R, cuts)
    yi = {t: k for k, t in enumerate(YEARS)}
    chains, rows_mat, rows_sum = {}, [], []

    windows = {"2015-24 (core)": (2015, 2024), "2010-24 (full)": (2010, 2024),
               "2010-14": (2010, 2014), "2015-19": (2015, 2019), "2020-24": (2020, 2024)}
    for name, (a, b) in windows.items():
        T = counts_window(cls, yi[a], yi[b])
        ch = chain_summary(T, f"{label_prefix} {name}")
        chains[name] = ch
    # giddy cross-check on the core window
    m = Markov(cls[:, yi[2015]:yi[2024] + 1], classes=np.arange(K), summary=False)
    assert np.allclose(m.transitions, chains["2015-24 (core)"]["T"]), "Markov counts mismatch"

    # annual transitions
    annual = []
    for t in YEARS[:-1]:
        T = counts_between(cls[:, yi[t]], cls[:, yi[t] + 1])
        P = row_p(T)
        annual.append({"from_year": t, "to_year": t + 1, "shorrocks": shorrocks(P),
                       "mean_diag": float(np.mean(np.diag(P))), "share_staying": float(np.trace(T) / T.sum()),
                       "T": T})
    # time homogeneity: sub-periods
    sub = hom([chains[s]["T"] for s in ("2010-14", "2015-19", "2020-24")], ["2010-14", "2015-19", "2020-24"])
    sub_core = hom([chains[s]["T"] for s in ("2015-19", "2020-24")], ["2015-19", "2020-24"])
    # each annual transition vs the rest (2015-24 core): Q contribution
    ann_core = [a for a in annual if a["from_year"] >= 2015]
    h_ann = homogeneity([a["T"] for a in ann_core], regime_names=[str(a["from_year"]) for a in ann_core])
    for a, q in zip(ann_core, h_ann.Q_table):
        a["Q_contribution_core"] = float(q.sum())
    # break: is 2014->15 different from neighbouring transitions?
    tr = {a["from_year"]: a["T"] for a in annual}
    brk = hom([tr[2013], tr[2014], tr[2015]], ["2013-14", "2014-15", "2015-16"])
    adj = hom([tr[2013], tr[2015]], ["2013-14", "2015-16"])
    # first-order Markov property (core): condition t->t+1 on class at t-1
    so = np.zeros((K, K, K))
    for t in range(2016, 2024):
        prev, now, nxt = cls[:, yi[t - 1]], cls[:, yi[t]], cls[:, yi[t] + 1]
        for c in range(K):
            s = prev == c
            so[c] += counts_between(now[s], nxt[s])
    mk = hom(so, [f"prev {c + 1}" for c in range(K)])
    # how much does a first-order chain mis-predict 2-step persistence?
    P1 = chains["2015-24 (core)"]["P"]
    T2 = np.zeros((K, K))
    for t in range(2015, 2023):
        T2 += counts_between(cls[:, yi[t]], cls[:, yi[t] + 2])
    two_step = {"observed_2step_diag": np.diag(row_p(T2)).tolist(), "implied_P2_diag": np.diag(P1 @ P1).tolist()}
    T5 = counts_between(cls[:, yi[2015]], cls[:, yi[2020]]) + counts_between(cls[:, yi[2019]], cls[:, yi[2024]])
    five_step = {"observed_5step_diag": np.diag(row_p(T5)).tolist(),
                 "implied_P5_diag": np.diag(np.linalg.matrix_power(P1, 5)).tolist()}

    for name, ch in chains.items():
        rows_mat += mat_rows(ch["P"], ch["label"], "P")
        rows_mat += mat_rows(ch["mfpt"], ch["label"], "mfpt")
        rows_sum.append({"chain": ch["label"], "n_transitions": ch["n_transitions"], "shorrocks": ch["shorrocks"],
                         "mean_diag": ch["mean_diag"], "share_staying": ch["share_staying"],
                         **{f"steady_{c + 1}": ch["steady_state"][c] for c in range(K)},
                         **{f"sojourn_{c + 1}": ch["sojourn"][c] for c in range(K)},
                         **{f"diag_{c + 1}": ch["P"][c, c] for c in range(K)}})
    tests = [
        {"test": "time homogeneity: 2010-14 vs 2015-19 vs 2020-24", **sub},
        {"test": "time homogeneity: 2015-19 vs 2020-24", **sub_core},
        {"test": "annual transitions 2015-16 ... 2023-24", "Q": h_ann.Q, "Q_p": h_ann.Q_p_value, "LR": h_ann.LR,
         "LR_p": h_ann.LR_p_value, "dof": h_ann.dof,
         "max_abs_diff_p": float(np.max(np.abs(h_ann.p_h1 - h_ann.p_h0[None])))},
        {"test": "break: 2013-14 vs 2014-15 vs 2015-16", **brk},
        {"test": "control: 2013-14 vs 2015-16 (no break)", **adj},
        {"test": "first-order Markov property (regimes = class at t-1), 2015-24", **mk},
    ]
    return {"cls": cls, "chains": chains, "annual": annual, "tests": tests, "rows_mat": rows_mat,
            "rows_sum": rows_sum, "two_step": two_step, "five_step": five_step}


def spatial_markov(C, w, cuts, years, label):
    R = C["rel"][years].values
    L = np.column_stack([weights.lag_spatial(w, R[:, j]) for j in range(R.shape[1])])
    lag_cuts = np.quantile(L.ravel(), np.linspace(0, 1, K + 1)[1:-1])
    sm = Spatial_Markov(R, w, k=K, m=K, fixed=True, cutoffs=cuts, lag_cutoffs=lag_cuts)
    rows = []
    for m_ in range(K):
        rows += mat_rows(sm.P[m_], f"{label} | lag {CLASS_LABELS[m_]}", "P")
    S = np.vstack([np.ravel(s) for s in sm.S])
    test = {"chain": label, "Q": sm.Q, "Q_p": sm.Q_p_value, "LR": sm.LR, "LR_p": sm.LR_p_value,
            "dof": sm.dof_hom, "x2": sm.x2, "x2_p": sm.x2_pvalue, "x2_dof": sm.x2_dof,
            "max_abs_diff_p": float(np.max(np.abs(sm.P - sm.p[None]))),
            "n_transitions": int(sm.transitions.sum())}
    # upward-mobility summary: P(move up) from Q1-Q2 by lag class; P(move down) from Q4-Q5 by lag class
    up = {f"p_up_from_Q{c + 1}_lag{m_ + 1}": float(sm.P[m_][c, c + 1:].sum()) for m_ in range(K) for c in range(K - 1)}
    down = {f"p_down_from_Q{c + 1}_lag{m_ + 1}": float(sm.P[m_][c, :c].sum()) for m_ in range(K) for c in range(1, K)}
    per_lag = [{"chain": label, "lag_class": m_ + 1, "chi2": float(sm.chi2[m_][0]), "chi2_p": float(sm.chi2[m_][1]),
                "dof": int(sm.chi2[m_][2]), "n": int(sm.T[m_].sum()),
                **{f"steady_{c + 1}": S[m_, c] for c in range(K)}} for m_ in range(K)]
    return sm, rows, test, {**up, **down, "chain": label}, per_lag, lag_cuts


def lad_rel():
    lp = pd.read_parquet(PROCESSED / "elec_lad_panel.parquet")
    lad, w, islands = lad_weights()
    y = lp.pivot(index="lad_cd", columns="year", values="mean_kwh").loc[lad["lad_cd"], YEARS]
    m = lp.pivot(index="lad_cd", columns="year", values="meters").loc[lad["lad_cd"], YEARS]
    return {"y": y, "rel": relative(y, m), "meters": m, "geo": lad}, w, islands


# --------------------------------------------------------------------------- LISA Markov
def lisa_markov(C, w, years):
    Lc = pd.read_parquet(CACHE / "lisa_queen.parquet")
    ids = C["geo"]["lsoa21cd"].values
    q = Lc.pivot(index="lsoa21cd", columns="year", values="quad").loc[ids, years].values
    f = Lc.pivot(index="lsoa21cd", columns="year", values="cluster_fdr").loc[ids, years].values
    names4 = ["HH", "LH", "LL", "HL"]
    idx4 = {n: i for i, n in enumerate(names4)}
    T4 = np.zeros((4, 4))
    names5 = ["HH", "LH", "LL", "HL", "ns"]
    idx5 = {n: i for i, n in enumerate(names5)}
    T5 = np.zeros((5, 5))
    for j in range(len(years) - 1):
        a = np.vectorize(idx4.get)(q[:, j])
        b = np.vectorize(idx4.get)(q[:, j + 1])
        T4 += counts_between(a, b, 4)
        a = np.vectorize(idx5.get)(f[:, j])
        b = np.vectorize(idx5.get)(f[:, j + 1])
        T5 += counts_between(a, b, 5)
    lm = LISA_Markov(C["rel"][years].values, w)
    assert np.allclose(lm.transitions, T4), "LISA Markov counts mismatch"
    P4, P5 = row_p(T4), row_p(T5)
    rows = mat_rows(P4, "LISA Markov (quadrants)", "P", names4) + mat_rows(P5, "LISA Markov (FDR clusters + ns)", "P",
                                                                          names5)
    # move-type shares (Rey 2001 typology): type I own moves only, II lag only, III both, 0 none
    tt = {"none": np.trace(T4)}
    cross = T4.copy()
    np.fill_diagonal(cross, 0)
    # own changes quadrant along H/L axis but lag does not: HH<->LH, LL<->HL
    own = T4[0, 1] + T4[1, 0] + T4[2, 3] + T4[3, 2]
    lagm = T4[0, 3] + T4[3, 0] + T4[1, 2] + T4[2, 1]
    both = T4[0, 2] + T4[2, 0] + T4[1, 3] + T4[3, 1]
    tot = T4.sum()
    test = {"chi2_independence_own_vs_lag": float(lm.chi_2[0]), "chi2_p": float(lm.chi_2[1]),
            "chi2_dof": int(lm.chi_2[2]), "share_no_move": tt["none"] / tot, "share_own_only": own / tot,
            "share_lag_only": lagm / tot, "share_joint_move": both / tot,
            "steady_state_quadrants": dict(zip(names4, np.ravel(ergodic.steady_state(P4)).tolist())),
            "steady_state_fdr": dict(zip(names5, np.ravel(ergodic.steady_state(P5)).tolist())),
            "sojourn_fdr": dict(zip(names5, (1 / (1 - np.diag(P5))).tolist()))}
    return rows, test, P4, P5


# --------------------------------------------------------------------------- rose
ROSE_WINDOWS = [(2010, 2014), (2015, 2019), (2019, 2020), (2021, 2022), (2015, 2024)]


def roses(C, w):
    rows, objs = [], {}
    for a, b in ROSE_WINDOWS:
        Y = C["rel"][[a, b]].values
        r = Rose(Y, w, k=8)
        np.random.seed(SEED + a)
        r.permute(permutations=PERMS)
        for s in range(8):
            rows.append({"window": f"{a}-{b}", "sector": s + 1, "deg_from": s * 45, "deg_to": (s + 1) * 45,
                         "count": int(r.counts[s]), "expected_perm": float(r.expected_perm[s]),
                         "p_two_sided": float(r.p[s])})
        objs[(a, b)] = r
    return pd.DataFrame(rows), objs


# --------------------------------------------------------------------------- figures
def heat(ax, P, title, names=None, fmt="{:.2f}"):
    names = names or ["Q1", "Q2", "Q3", "Q4", "Q5"]
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("b", ["#ffffff"] + SEQ)
    ax.imshow(P, cmap=cmap, vmin=0, vmax=1)
    for a in range(P.shape[0]):
        for b in range(P.shape[1]):
            ax.text(b, a, fmt.format(P[a, b]), ha="center", va="center", fontsize=10,
                    color="white" if P[a, b] > 0.55 else INK)
    ax.set_xticks(range(len(names)), names)
    ax.set_yticks(range(len(names)), names)
    ax.set_xlabel("Class at t+1")
    ax.set_ylabel("Class at t")
    ax.set_title(title, fontsize=11)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def shorrocks_figure(a, b, d):
    """Annual Shorrocks index for LSOAs, LSOAs excl. split and LADs (rows of ``markov_annual``)."""
    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    ax.plot(a.to_year, a.shorrocks, color=CAT[0], marker="o", ms=3.5, label="LSOAs")
    ax.plot(b.to_year, b.shorrocks, color=CAT[0], ls=":", lw=1.5, label="LSOAs excl. split")
    ax.plot(d.to_year, d.shorrocks, color=CAT[1], marker="s", ms=3.5, ls="--", label="LADs")
    ax.set_ylabel("Shorrocks mobility index")
    ax.set_title("Year-to-year quintile mobility", pad=10)
    ax.set_xticks(range(2011, 2025, 2))
    # headroom above the data holds the event labels, inside the axes and clear of the title
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.18 * (hi - lo))
    for t, lab, ha in ((2015, "2015 break", "left"), (2020, "COVID-19", "right"), (2022, "price crisis", "left")):
        ax.axvline(t, color=AXIS, lw=0.8, zorder=0)
        ax.text(t + (0.15 if ha == "left" else -0.15), 0.97, lab, transform=ax.get_xaxis_transform(),
                fontsize=8.5, color=INK2, ha=ha, va="top")
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False,
              handlelength=2.2, columnspacing=1.2)
    fig.tight_layout()
    return savefig(fig, "markov_shorrocks_annual")


def figures(LS, LD, SM, annual_core_noSX, rose_objs, P5):
    files = []
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    heat(ax, LS["chains"]["2015-24 (core)"]["P"], "LSOA quintile transitions, 2015–24")
    fig.tight_layout()
    files += savefig(fig, "markov_classic_lsoa_core")

    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    heat(ax, LD["chains"]["2015-24 (core)"]["P"], "LAD quintile transitions, 2015–24")
    fig.tight_layout()
    files += savefig(fig, "markov_classic_lad_core")

    sm = SM
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
    heat(axes[0], sm.P[0], "Neighbours in lowest quintile")
    heat(axes[1], sm.P[K - 1], "Neighbours in highest quintile")
    fig.suptitle("Spatial Markov, LSOAs 2015–24", x=0.02, ha="left", fontweight="bold")
    fig.tight_layout()
    files += savefig(fig, "markov_spatial_lsoa_extremes")

    a = pd.DataFrame([{k: v for k, v in x.items() if k != "T"} for x in LS["annual"]])
    b = pd.DataFrame([{k: v for k, v in x.items() if k != "T"} for x in annual_core_noSX])
    d = pd.DataFrame([{k: v for k, v in x.items() if k != "T"} for x in LD["annual"]])
    files += shorrocks_figure(a, b, d)

    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    heat(ax, P5, "LISA Markov (FDR clusters), 2015–24", names=["HH", "LH", "LL", "HL", "ns"])
    fig.tight_layout()
    files += savefig(fig, "markov_lisa_fdr")

    sel = [(2015, 2019), (2019, 2020), (2021, 2022), (2015, 2024)]
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 7.6), subplot_kw={"projection": "polar"})
    for ax, key in zip(axes.ravel(), sel):
        r = rose_objs[key]
        th = np.deg2rad(np.arange(8) * 45 + 22.5)
        sig = r.p < 0.05
        cols = [CAT[1] if (s and c > e) else (CAT[0] if s else "#c3c2b7")
                for s, c, e in zip(sig, r.counts, r.expected_perm)]
        ax.bar(th, r.counts, width=np.deg2rad(43), color=cols, edgecolor="white", linewidth=1)
        ax.plot(np.r_[th, th[0]], np.r_[r.expected_perm, r.expected_perm[0]], color=INK, lw=1, ls="--")
        ax.set_title(f"{key[0]}–{key[1]}", fontsize=11, pad=14)
        ax.set_xticks(np.deg2rad([0, 90, 180, 270]), ["own up", "neighbours\nup", "own\ndown", "neighbours\ndown"],
                      fontsize=9)
        top = max(r.counts.max(), r.expected_perm.max())
        ax.set_yticks([round(top / 2, -3)], [f"{round(top / 2, -3):,.0f}"])
        ax.set_rlabel_position(112)
        ax.tick_params(axis="y", labelsize=8, colors=MUTED)
        ax.grid(True, color=GRID)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig.legend(handles=[Patch(color=CAT[1], label="More than expected (p<0.05)"),
                        Patch(color=CAT[0], label="Fewer than expected (p<0.05)"),
                        Patch(color="#c3c2b7", label="Not significant"),
                        Line2D([], [], color=INK, ls="--", lw=1, label="Permutation mean")],
               loc="lower center", ncol=2, fontsize=9)
    fig.suptitle("Directional LISA: co-movement of LSOAs and neighbours", x=0.02, y=1.0, ha="left",
                 fontweight="bold")
    fig.subplots_adjust(bottom=0.13, top=0.86, hspace=0.42, wspace=0.45)
    files += savefig(fig, "rose_directional_lisa")
    return files


def run():
    set_style()
    log("Markov: data")
    C = core("elec")
    w = lsoa_weights("queen", C["mask"], C["tag"])
    pooled = C["rel"][YEARS].values.ravel()
    cuts = np.quantile(pooled, [0.2, 0.4, 0.6, 0.8])
    LS = lsoa_chains(C, w, cuts)
    Cn = core("elec", exclude_sx=True)
    LSn = lsoa_chains(Cn, None, cuts, "lsoa excl. S/X")
    Dd, wl, islands = lad_rel()
    cuts_lad = np.quantile(Dd["rel"].values.ravel(), [0.2, 0.4, 0.6, 0.8])
    LD = lsoa_chains(Dd, None, cuts_lad, "lad")  # same machinery on the LAD panel
    log("Markov: spatial Markov")
    SM, rows_sm, test_sm, upd, per_lag, lag_cuts = spatial_markov(C, w, cuts, CORE_YEARS, "LSOA 2015-24")
    sm_full = spatial_markov(C, w, cuts, YEARS, "LSOA 2010-24")
    sm_lad = spatial_markov(Dd, wl, cuts_lad, CORE_YEARS, "LAD 2015-24")
    sm_sub = {k: spatial_markov(C, w, cuts, list(range(a, b + 1)), f"LSOA {k}") for k, (a, b) in SUBPERIODS.items()}
    log("Markov: LISA Markov and roses")
    rows_lm, test_lm, P4, P5 = lisa_markov(C, w, CORE_YEARS)
    R, rose_objs = roses(C, w)

    # ---- tables
    write_table(pd.DataFrame([{"scope": "LSOA", "cuts_relative": list(cuts)},
                              {"scope": "LAD", "cuts_relative": list(cuts_lad)},
                              {"scope": "LSOA spatial lag (2015-24)", "cuts_relative": list(lag_cuts)}])
                .assign(cuts_relative=lambda d: d.cuts_relative.map(lambda v: ", ".join(f"{x:.4f}" for x in v))),
                "markov_class_bounds")
    write_table(pd.DataFrame(LS["rows_mat"] + LSn["rows_mat"] + LD["rows_mat"] + rows_sm + sm_lad[1] + sm_full[1]
                             + sum((v[1] for v in sm_sub.values()), []) + rows_lm), "markov_matrices")
    write_table(pd.DataFrame(LS["rows_sum"] + LSn["rows_sum"] + LD["rows_sum"]), "markov_summary")
    ann = []
    for scope, X in (("LSOA", LS), ("LSOA excl. S/X", LSn), ("LAD", LD)):
        for a in X["annual"]:
            ann.append({"scope": scope, **{k: v for k, v in a.items() if k != "T"},
                        **{f"diag_{c + 1}": row_p(a["T"])[c, c] for c in range(K)}})
    write_table(pd.DataFrame(ann), "markov_annual")
    tests = []
    for scope, X in (("LSOA", LS), ("LSOA excl. S/X", LSn), ("LAD", LD)):
        for t in X["tests"]:
            tests.append({"scope": scope, **{k: (str(v) if isinstance(v, list) else v) for k, v in t.items()}})
    write_table(pd.DataFrame(tests), "markov_time_homogeneity_tests")
    sp = [test_sm, sm_full[2], sm_lad[2]] + [v[2] for v in sm_sub.values()]
    write_table(pd.DataFrame(sp), "markov_spatial_homogeneity_tests")
    write_table(pd.DataFrame(per_lag + sm_lad[4] + sm_full[4]), "markov_spatial_by_lag")
    write_table(pd.DataFrame([upd, sm_full[3], sm_lad[3]] + [v[3] for v in sm_sub.values()]), "markov_spatial_updown")
    write_table(R, "rose_counts")
    save_json({"lisa_markov": test_lm, "lsoa_markov_order": {"two_step": LS["two_step"], "five_step": LS["five_step"]},
               "lad_islands_attached": islands, "lad_markov_order": {"two_step": LD["two_step"]}},
              "markov_meta")
    np.savez(CACHE / "markov_lsoa_classes.npz", cls=LS["cls"], ids=C["geo"]["lsoa21cd"].values, cuts=cuts)
    files = figures(LS, LD, SM, LSn["annual"], rose_objs, P5)
    meta = __import__("json").loads((RESULTS / "markov_meta.json").read_text())
    meta["figures"] = files
    save_json(meta, "markov_meta")
    log("Markov done")


if __name__ == "__main__":
    run()
