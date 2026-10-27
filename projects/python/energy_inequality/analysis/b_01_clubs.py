"""Section 1: club convergence in domestic electricity per meter.

Phillips-Sul log t test, clustering and merging on HP-filtered log kWh per
meter, at LAD (main), MSOA and LSOA (robustness), for 2010-2024 and the
pre-shock window 2010-2019 (plus 2015-2024, the post-processing-break
window, at LAD and MSOA). Then: transition-path plots, club maps, join-count
tests of spatial clustering of membership, and ordered/multinomial logit
profiles of club membership on 2021 covariates.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import b_common as bc
import b_maps as bm
import b_phillips_sul as ps

WINDOWS = {"full": (2010, 2024), "pre": (2010, 2019), "post15": (2015, 2024), "levels": (2010, 2024)}
RUNS = [("lad", "full"), ("lad", "pre"), ("lad", "post15"),
        ("msoa", "full"), ("msoa", "pre"), ("msoa", "post15"),
        ("lsoa", "full"), ("lsoa", "pre"), ("lad", "levels"), ("msoa", "levels")]
WINDOW_LABEL = {"full": "2010–2024", "pre": "2010–2019 (pre-shock)", "post15": "2015–2024",
                "levels": "2010–2024, levels"}
PROFILE_X = ["share_ch_electric", "share_ch_offgas", "share_ch_renewable", "share_detached", "share_flat",
             "mean_hh_size", "share_social_rent", "share_private_rent", "share_age_65plus",
             "dep19_pct", "epc_ee", "rural", "wales"]


def prepared_panel(level: str, window: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Raw kWh/meter and HP trend of log kWh/meter for complete units.

    Window 'levels' filters kWh/meter itself instead of its log: h_it is then
    a ratio to the cross-sectional mean, so common proportional shocks
    (COVID-19, the 2022 price fall) cancel exactly."""
    y0, y1 = WINDOWS[window]
    raw = bc.panel_wide(level).loc[:, y0:y1]
    raw = raw.dropna()
    vals = raw.values if window == "levels" else np.log(raw.values)
    trend = pd.DataFrame(ps.hp_trend(vals), index=raw.index, columns=raw.columns)
    return raw, trend


def fit_all(force=False) -> dict:
    out = {}
    for level, window in RUNS:
        def _fit(level=level, window=window):
            raw, X = prepared_panel(level, window)
            res = ps.run(X.values)
            return {"ids": X.index.tolist(), "years": list(X.columns), "res": res,
                    "n_dropped": int(bc.panel_wide(level).shape[0] - X.shape[0])}
        bc.log(f"clubs: {level} {window}")
        out[(level, window)] = bc.cached(f"clubs_{level}_{window}", _fit, force)
    return out


def summary_tables(fits):
    rows, merges, memb = [], [], {}
    for (level, window), f in fits.items():
        raw, X = prepared_panel(level, window)
        res = f["res"]
        lab = pd.Series(res.labels, index=f["ids"])
        memb[(level, window)] = lab
        g = res.global_test
        rows.append({"level": level, "window": window, "club": "all", "n": g["n"], "beta": g["beta"],
                     "t": g["t"], "converges": g["t"] > ps.CRIT,
                     "kwh_first": raw.iloc[:, 0].mean(), "kwh_last": raw.iloc[:, -1].mean(),
                     "n_initial_clubs": len(set(res.initial_labels) - {0}),
                     "n_final_clubs": len(set(res.labels) - {0}), "n_divergent": int((res.labels == 0).sum()),
                     "n_dropped_incomplete": f["n_dropped"]})
        for c in res.clubs:
            m = (lab == c["club"]).values
            rows.append({"level": level, "window": window, "club": c["club"], "n": c["n"], "beta": c["beta"],
                         "t": c["t"], "converges": c["t"] > ps.CRIT,
                         "kwh_first": raw.iloc[m, 0].mean(), "kwh_last": raw.iloc[m, -1].mean()})
        for d in res.merge_log:
            merges.append({"level": level, "window": window, **d})
    tab = pd.DataFrame(rows)
    bc.save_table(tab, "clubs_summary")
    bc.save_table(pd.DataFrame(merges), "clubs_merge_log")
    m = pd.DataFrame({f"{lv}_{w}": s for (lv, w), s in memb.items() if lv == "lad"})
    m.index.name = "lad_cd"
    bc.save_table(m.reset_index(), "clubs_membership_lad")
    m = pd.DataFrame({f"{lv}_{w}": s for (lv, w), s in memb.items() if lv == "msoa"})
    m.index.name = "msoa21cd"
    bc.save_table(m.reset_index(), "clubs_membership_msoa")
    m = pd.DataFrame({f"{lv}_{w}": s for (lv, w), s in memb.items() if lv == "lsoa"})
    m.index.name = "lsoa21cd"
    m.reset_index().to_parquet(bc.TAB / "clubs_membership_lsoa.parquet")
    return tab, memb


def monte_carlo_table():
    def _mc():
        out = []
        for T in (10, 15):
            for N in (300, 2000):
                d = ps.monte_carlo(reps=200, N=N, T=T, seed=T * N)
                d["T"], d["N"] = T, N
                out.append(d)
        return pd.concat(out)
    mc = bc.cached("ps_monte_carlo", _mc)
    s = mc.groupby(["N", "T"]).agg(
        size_reject_convergent=("reject_under_convergence", "mean"),
        power_reject_divergent=("reject_under_divergence", "mean"),
        power_reject_3club=("reject_global_3club", "mean"),
        share_3_clubs_found=("n_clubs_found", lambda v: (v == 3).mean()),
        mean_adjusted_rand=("ari_3club", "mean")).reset_index()
    bc.save_table(s, "clubs_monte_carlo")
    return s


def transition_figure(level, window, fits):
    raw, X = prepared_panel(level, window)
    lab = fits[(level, window)]["res"].labels
    h = ps.transition_paths(X.values)
    clubs = sorted(set(lab) - {0})
    cols = bc.club_colours(len(clubs))
    years = np.array(X.columns)
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    if level == "lad":
        for c, col in zip(clubs, cols):
            for row in h[lab == c]:
                ax.plot(years, row, color=col, alpha=0.15, lw=0.6)
    meds = []
    for c, col in zip(clubs, cols):
        hc = h[lab == c]
        med = np.median(hc, axis=0)
        meds.append(med[-1])
        if level != "lad" and len(clubs) <= 12:
            ax.fill_between(years, np.quantile(hc, 0.25, 0), np.quantile(hc, 0.75, 0), color=col, alpha=0.18, lw=0)
        ax.plot(years, med, color=col, lw=2.4 if len(clubs) <= 12 else 1.6)
    # Direct labels (every club when <= 12, otherwise every third plus the last).
    show = list(range(len(clubs))) if len(clubs) <= 12 else sorted(set(range(0, len(clubs), 3)) | {len(clubs) - 1})
    ylim = ax.get_ylim()
    gap = (ylim[1] - ylim[0]) * 0.045
    ypos = bm.repel([meds[k] for k in show], gap)
    for k, yp in zip(show, ypos):
        n = int((lab == clubs[k]).sum())
        ax.annotate(f"{clubs[k]} (n={n:,})", (years[-1], meds[k]), xytext=(years[-1] + 0.4, yp),
                    textcoords="data", va="center", fontsize=11, color=bc.INK,
                    arrowprops=dict(arrowstyle="-", color=bc.GRID, lw=0.6))
    if (lab == 0).any():
        ax.plot(years, np.median(h[lab == 0], 0), color=bc.INK, lw=1.4, ls="--")
    ax.axhline(1, color=bc.INK2, lw=0.8)
    ax.set_xlim(years[0], years[-1])
    ax.set_xticks([y for y in years if ((y - years[0]) % 4 == 0 and years[-1] - y >= 3) or y == years[-1]])
    ax.set_ylabel("Relative transition path $h_{it}$")
    lvl = {"lad": "local authorities", "msoa": "MSOAs", "lsoa": "LSOAs"}[level]
    ax.set_title(f"Club transition paths, {lvl}\n{WINDOW_LABEL[window]}", loc="left", fontsize=14)
    ax.text(1.01, 1.0, "Club", transform=ax.transAxes, fontsize=11, color=bc.INK2, va="bottom")
    var = "kWh per meter" if window == "levels" else "log kWh per meter"
    band = "" if (level == "lad" or len(clubs) > 12) else " Band = interquartile range."
    div = f" Dashed = divergent units (n={(lab == 0).sum()})." if (lab == 0).any() else ""
    fig.text(0.0, -0.06, f"HP-filtered {var} (λ = 400) relative to the cross-sectional mean. "
                         f"Bold line = club median.{band}{div}", fontsize=10, color=bc.INK2, wrap=True)
    fig.subplots_adjust(right=0.74)
    return bc.save_fig(fig, f"clubs_paths_{level}_{window}")


def club_map(level, window, fits):
    lab = pd.Series(fits[(level, window)]["res"].labels, index=fits[(level, window)]["ids"], name="club")
    g = bc.display_geometry(level).join(lab, how="left")
    clubs = sorted(set(lab) - {0})
    lvl = {"lad": "local authorities", "msoa": "MSOAs", "lsoa": "LSOAs"}[level]
    bounds = bc.display_geometry("lad") if level != "lad" else None
    title = f"Convergence clubs, {lvl}\n{WINDOW_LABEL[window]}"
    n_div, n_inc = int((lab == 0).sum()), int(g.club.isna().sum())
    extra = (f" Divergent: {n_div}." if n_div else "") + (f" Incomplete series (blank): {n_inc}." if n_inc else "")
    lines = " White lines: local authorities." if bounds is not None else ""
    if len(clubs) <= 12:
        cols = dict(zip(clubs, bc.club_colours(len(clubs))))
        labels = {c: f"Club {c} ({(lab == c).sum():,})" for c in clubs}
        if n_div:
            cols[0] = "#ffffff"
            labels[0] = f"Divergent ({n_div})"
        g["club"] = g["club"].fillna(-1).astype(int)
        return bm.choropleth(g, "club", title, f"clubs_map_{level}_{window}", colours=cols, labels=labels,
                             boundaries=bounds, note=f"Club 1 = highest steady-state kWh per meter.{lines}"
                             + (f" Incomplete series (blank): {n_inc}." if n_inc else ""))
    import matplotlib as mpl
    g.loc[g.club == 0, "club"] = np.nan
    g = g.dropna(subset=["club"])
    cmap = mpl.colors.ListedColormap(bc.club_colours(len(clubs)))
    norm = mpl.colors.BoundaryNorm(np.arange(0.5, len(clubs) + 1.5), len(clubs))
    return bm.choropleth(g, "club", title, f"clubs_map_{level}_{window}", cmap=cmap, norm=norm,
                         cbar_label=f"Club (1 = highest of {len(clubs)})", boundaries=bounds,
                         cbar_ticks=sorted({1, *range(5, len(clubs) + 1, 5), len(clubs)}),
                         note=f"Divergent and incomplete units left blank.{extra}{lines}")


def join_counts(memb, perms=999):
    """Spatial clustering of club membership.

    For each club, the BB join count (links joining two members) against
    random permutation; and, over all clubs, the share of neighbour links
    joining units in the same club (a multi-category join-count statistic).
    The same label permutations serve every club, so the BB counts are
    computed for all clubs at once from the link list. Moran's I of the
    ordinal club number is added for comparison."""
    import esda

    rows = []
    rng = np.random.default_rng(bc.SEED)
    for (level, window), lab in memb.items():
        w = bc.queen_weights(level)
        ids = [i for i in w.id_order if i in lab.index]
        ws = bc.subset_weights(w, ids, bc.unit_coords(level))
        ws.transform = "b"
        y = lab.loc[ids].values.astype(int)
        K = y.max()
        coo = ws.sparse.tocoo()
        r_, c_ = coo.row, coo.col  # each undirected link appears twice

        def stats(lbl):
            same = lbl[r_] == lbl[c_]
            bb = np.bincount(lbl[r_][same], minlength=K + 1) / 2.0
            return same.mean(), bb

        same_obs, bb_obs = stats(y)
        sims = [stats(rng.permutation(y)) for _ in range(perms)]
        same_sim = np.array([s[0] for s in sims])
        bb_sim = np.vstack([s[1] for s in sims])
        rows.append({"level": level, "window": window, "club": "all", "statistic": "same-club link share",
                     "observed": same_obs, "expected_random": same_sim.mean(),
                     "ratio": same_obs / same_sim.mean(), "z_sim": (same_obs - same_sim.mean()) / same_sim.std(),
                     "p_sim": (1 + (same_sim >= same_obs).sum()) / (1 + perms)})
        for c in range(1, K + 1):
            sim = bb_sim[:, c]
            rows.append({"level": level, "window": window, "club": c, "statistic": "BB joins",
                         "observed": bb_obs[c], "expected_random": sim.mean(),
                         "ratio": bb_obs[c] / sim.mean() if sim.mean() > 0 else np.nan,
                         "z_sim": (bb_obs[c] - sim.mean()) / sim.std() if sim.std() > 0 else np.nan,
                         "p_sim": (1 + (sim >= bb_obs[c]).sum()) / (1 + perms)})
        keep = y > 0
        wr = bc.subset_weights(w, list(np.array(ids)[keep]), bc.unit_coords(level))
        mi = esda.Moran(y[keep].astype(float), wr, permutations=perms)
        rows.append({"level": level, "window": window, "club": "ordinal", "statistic": "Moran's I of club number",
                     "observed": mi.I, "expected_random": mi.EI, "ratio": np.nan, "z_sim": mi.z_sim,
                     "p_sim": mi.p_sim})
    tab = pd.DataFrame(rows)
    bc.save_table(tab, "clubs_join_counts")
    return tab


def profile_models(memb):
    """Ordered logit (main) and multinomial logit (robustness) of club
    membership on standardised 2021 covariates. MSOA: cluster-robust SEs by
    LAD; LAD: a reduced covariate set because n = 318."""
    import statsmodels.api as sm
    from statsmodels.miscmodels.ordinal_model import OrderedModel

    rows, ame_rows, fit_rows = [], [], []
    specs = [("msoa", "full", PROFILE_X), ("msoa", "pre", PROFILE_X),
             ("lad", "full", ["share_ch_electric", "share_ch_offgas", "share_detached", "mean_hh_size",
                              "dep19_pct", "epc_ee"])]
    for level, window, xs in specs:
        lab = memb[(level, window)]
        cov = bc.covariates(level)
        d = cov.join(lab.rename("club"), how="inner")
        d = d[d.club > 0].dropna(subset=xs)
        K = d.club.max()
        y = (K + 1 - d.club).astype(int)  # higher = higher-consumption club
        Z = (d[xs] - d[xs].mean()) / d[xs].std()
        cat = pd.Series(pd.Categorical(y, ordered=True), index=d.index)
        om = OrderedModel(cat, Z, distr="logit")
        kw = {}
        if level == "msoa":
            kw = {"cov_type": "cluster", "cov_kwds": {"groups": pd.factorize(d.lad_cd)[0]}}
        r = om.fit(method="bfgs", disp=False, maxiter=2000, **kw)
        for v in xs:
            rows.append({"level": level, "window": window, "model": "ordered_logit", "term": v,
                         "coef": r.params[v], "se": r.bse[v], "z": r.tvalues[v], "p": r.pvalues[v],
                         "odds_ratio_per_sd": np.exp(r.params[v])})
        # pseudo R2 against thresholds-only model
        p0 = y.value_counts(normalize=True)
        ll0 = (y.map(np.log(p0))).sum()
        fit_rows.append({"level": level, "window": window, "model": "ordered_logit", "n": len(d),
                         "n_clubs": int(K), "llf": r.llf, "ll_null": ll0, "mcfadden_r2": 1 - r.llf / ll0,
                         "aic": r.aic})
        # Multinomial logit robustness (no ordering assumption). With ~20 clubs
        # a 20-outcome MNL is not estimable, so clubs are grouped into three
        # tiers by their position in the club order: 'high' = clubs that start
        # within the top 20% of units, 'low' = clubs that end in the bottom 20%.
        if level == "msoa":
            sizes = d.club.value_counts().sort_index()
            cum_end = sizes.cumsum() / sizes.sum()
            cum_start = cum_end - sizes / sizes.sum()
            tier_of = {c: ("high" if cum_start[c] < 0.2 and cum_end[c] <= 0.4 else
                           "low" if cum_end[c] > 0.8 and cum_start[c] >= 0.6 else "middle") for c in sizes.index}
            tier = d.club.map(tier_of)
            codes = tier.map({"middle": 0, "high": 1, "low": 2}).values
            mn = sm.MNLogit(codes, sm.add_constant(Z)).fit(method="bfgs", disp=False, maxiter=2000)
            me = mn.get_margeff(at="overall")
            sf = me.summary_frame().reset_index()
            sf.columns = ["outcome", "term", "dy_dx", "se", "z", "p", "ci_lo", "ci_hi"]
            sf["tier"] = sf.outcome.str.extract(r"=(\d)").iloc[:, 0].astype(int).map({0: "middle", 1: "high", 2: "low"})
            sf["n_tier"] = sf.tier.map(tier.value_counts())
            sf["level"], sf["window"] = level, window
            ame_rows.append(sf)
            fit_rows.append({"level": level, "window": window, "model": "multinomial_logit", "n": len(d),
                             "n_clubs": int(K), "llf": mn.llf, "ll_null": mn.llnull,
                             "mcfadden_r2": mn.prsquared, "aic": mn.aic})
    coef = pd.DataFrame(rows)
    bc.save_table(coef, "clubs_profile_ordered_logit")
    bc.save_table(pd.DataFrame(fit_rows), "clubs_profile_fit")
    if ame_rows:
        bc.save_table(pd.concat(ame_rows), "clubs_profile_mnlogit_ame")
    return coef


def profile_figure(coef):
    d = coef[(coef.level == "msoa") & (coef.window == "full")].copy()
    nice = {"share_ch_electric": "Electric heating", "share_ch_offgas": "Oil/LPG/solid heating",
            "share_ch_renewable": "Renewable heating", "share_detached": "Detached", "share_flat": "Flats",
            "mean_hh_size": "Household size", "share_social_rent": "Social rent",
            "share_private_rent": "Private rent", "share_age_65plus": "Aged 65+",
            "dep19_pct": "Less deprived (IMD pct)", "epc_ee": "EPC efficiency (median)",
            "rural": "Rural", "wales": "Wales"}
    d["label"] = d.term.map(nice)
    d = d.sort_values("coef")
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    yy = np.arange(len(d))
    lo, hi = d.coef - 1.96 * d.se, d.coef + 1.96 * d.se
    ax.hlines(yy, lo, hi, color=bc.CAT3[0], lw=2)
    ax.plot(d.coef, yy, "o", color=bc.CAT3[0], ms=8, mec="white", mew=1.5)
    ax.axvline(0, color=bc.INK2, lw=0.9)
    ax.set_yticks(yy, d.label)
    ax.set_xlabel("Log-odds per SD (higher club = higher consumption)")
    ax.set_title("Who is in the high-consumption clubs?\nOrdered logit, MSOAs, 2010–2024 clubs", loc="left",
                 fontsize=13)
    ax.grid(axis="y", visible=False)
    fig.text(0.0, -0.02, "95% intervals; standard errors clustered by local authority. 2021 covariates.",
             fontsize=9.5, color=bc.INK2)
    return bc.save_fig(fig, "clubs_profile_ordered_logit")


def lad_cross_level(memb):
    """Agreement between scales and windows."""
    from scipy.stats import spearmanr
    from sklearn.metrics import adjusted_rand_score

    rows = []
    lad = memb[("lad", "full")]
    for a, b in [(("lad", "full"), ("lad", "pre")), (("lad", "full"), ("lad", "post15")),
                 (("lad", "full"), ("lad", "levels")), (("msoa", "full"), ("msoa", "levels")),
                 (("msoa", "full"), ("msoa", "pre")), (("lsoa", "full"), ("lsoa", "pre"))]:
        s = pd.concat([memb[a].rename("a"), memb[b].rename("b")], axis=1).dropna()
        rows.append({"comparison": f"{a[0]} {a[1]} vs {b[1]}", "n": len(s),
                     "adjusted_rand": adjusted_rand_score(s.a, s.b),
                     "spearman": spearmanr(s.a, s.b).statistic})
    # MSOA club vs the club of its LAD (rank agreement only; club counts differ)
    cov = bc.covariates("msoa")[["lad_cd"]].join(memb[("msoa", "full")].rename("msoa_club"))
    cov["lad_club"] = cov.lad_cd.map(lad)
    cov = cov.dropna()
    cov = cov[(cov.msoa_club > 0) & (cov.lad_club > 0)]
    rows.append({"comparison": "msoa full vs parent lad full", "n": len(cov), "adjusted_rand": np.nan,
                 "spearman": spearmanr(cov.msoa_club, cov.lad_club).statistic})
    t = pd.DataFrame(rows)
    bc.save_table(t, "clubs_agreement")
    ct = pd.crosstab(memb[("lad", "full")].rename("full"), memb[("lad", "pre")].rename("pre"))
    bc.save_table(ct, "clubs_crosstab_lad_full_pre", index=True)
    return t


def dispersion_table():
    """Sigma-convergence diagnostics: cross-sectional SD of log kWh/meter and
    CV of kWh/meter per year and level (complete units)."""
    rows = []
    for level in ("lad", "msoa", "lsoa"):
        w = bc.panel_wide(level).dropna()
        for y in w.columns:
            rows.append({"level": level, "year": y, "n": len(w), "sd_log": np.log(w[y]).std(),
                         "cv": w[y].std() / w[y].mean(), "mean_kwh": w[y].mean()})
    t = pd.DataFrame(rows)
    bc.save_table(t, "clubs_sigma_dispersion")
    return t


def main(force=False):
    bc.set_style()
    dispersion_table()
    mc = monte_carlo_table()
    bc.log(f"Monte Carlo:\n{mc}")
    fits = fit_all(force)
    tab, memb = summary_tables(fits)
    bc.log(f"\n{tab.to_string()}")
    for level, window in RUNS:
        transition_figure(level, window, fits)
        if level != "lsoa" or window == "full":
            club_map(level, window, fits)
    jc = join_counts(memb)
    bc.log(f"\n{jc[jc.club == 'all'].to_string()}")
    coef = profile_models(memb)
    profile_figure(coef)
    agree = lad_cross_level(memb)
    bc.log(f"\n{agree.to_string()}")


if __name__ == "__main__":
    main(force="--force" in sys.argv)
