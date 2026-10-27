"""Phillips & Sul (2007, 2009) log t convergence test, club clustering and
club merging, vectorised so that it runs on every LSOA (N ~ 35,000).

Notation follows Phillips & Sul (2007, Econometrica 75:1771-1855) and
Phillips & Sul (2009, J. Applied Econometrics 24:1153-1185), and the
algorithm mirrors the R package ConvergenceClubs (Sichera & Pizzuto 2019):

* X_it   positive panel (here the HP-filtered trend of log kWh per meter).
* h_it   relative transition path X_it / (N^-1 sum_i X_it).
* H_t    cross-sectional variance of h_it, N^-1 sum_i (h_it - 1)^2.
* log t  regression  log(H_1/H_t) - 2 log(log t) = a + b log t + u_t,
         t = [rT], ..., T, after discarding the first r fraction (r = 1/3).
         Convergence (b >= 0) is rejected when t_b < -1.65 (one-sided 5%),
         using a Newey-West HAC standard error. b/2 estimates the speed alpha.

Because H_t for a group only needs n, sum_i X_it and sum_i X_it^2, the
statistics for thousands of candidate groups (core-group prefixes, sieve
candidates) are computed at once from cumulative sums.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve

CRIT = -1.65


# --------------------------------------------------------------------------
# Filtering
# --------------------------------------------------------------------------
def hp_trend(y: np.ndarray, lamb: float = 400.0) -> np.ndarray:
    """Hodrick-Prescott trend of each row of y (N x T). lambda = 400 is the
    value Phillips & Sul (2007) use for annual data."""
    T = y.shape[1]
    d = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(T - 2, T))
    A = (sparse.eye(T) + lamb * (d.T @ d)).tocsc()
    Ainv = np.linalg.inv(A.toarray())
    return y @ Ainv.T


# --------------------------------------------------------------------------
# log t regression, vectorised over groups
# --------------------------------------------------------------------------
def _trim_index(T: int, r: float) -> np.ndarray:
    """Time index t (1-based) kept in the regression: t = [rT], ..., T."""
    start = max(int(np.floor(r * T)), 2)  # log(log 1) undefined, so t >= 2
    return np.arange(start, T + 1)


def _nw_lag(n: int) -> int:
    return int(np.floor(4 * (n / 100.0) ** (2.0 / 9.0)))


def logt_from_H(H: np.ndarray, r: float = 1 / 3) -> tuple[np.ndarray, np.ndarray]:
    """log t regression for each row of H (m x T, cross-sectional variance
    paths). Returns (b_hat, t_stat) arrays of length m."""
    H = np.atleast_2d(H)
    m, T = H.shape
    t = _trim_index(T, r)
    with np.errstate(divide="ignore", invalid="ignore"):
        Y = np.log(H[:, [0]] / H[:, t - 1]) - 2.0 * np.log(np.log(t))[None, :]
    x = np.log(t)
    z = x - x.mean()
    sxx = (z ** 2).sum()
    b = (Y - Y.mean(axis=1, keepdims=True)) @ z / sxx
    a = Y.mean(axis=1) - b * x.mean()
    e = Y - a[:, None] - b[:, None] * x[None, :]
    ez = e * z[None, :]
    S = (ez ** 2).sum(axis=1)
    L = _nw_lag(len(t))
    for lag in range(1, L + 1):
        wl = 1.0 - lag / (L + 1.0)
        S = S + 2.0 * wl * (ez[:, lag:] * ez[:, :-lag]).sum(axis=1)
    se = np.sqrt(np.maximum(S, 0.0)) / sxx
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = b / se
    # A group whose H_t is exactly zero (identical paths) has converged.
    tstat = np.where(np.isfinite(tstat), tstat, np.inf)
    return b, tstat


def H_from_sums(n, s, q):
    """H_t from group size n, sum and sum of squares of X over the group."""
    n = np.asarray(n, dtype=float)
    if n.ndim == 1:
        n = n[:, None]
    return n * q / s ** 2 - 1.0


def logt(X: np.ndarray, r: float = 1 / 3) -> dict:
    """Log t test on a single group (rows of X are units)."""
    X = np.asarray(X, dtype=float)
    H = H_from_sums(np.array([X.shape[0]]), X.sum(0)[None], (X ** 2).sum(0)[None])
    b, t = logt_from_H(H, r)
    return {"n": X.shape[0], "beta": float(b[0]), "alpha": float(b[0] / 2), "t": float(t[0])}


# --------------------------------------------------------------------------
# Clustering (Phillips & Sul 2009, steps 1-4) and merging
# --------------------------------------------------------------------------
@dataclass
class ClubResult:
    labels: np.ndarray            # club number (1 = highest), 0 = divergent
    clubs: list = field(default_factory=list)   # dicts with n, beta, t per club
    initial_labels: np.ndarray | None = None
    initial_clubs: list = field(default_factory=list)
    merge_log: list = field(default_factory=list)
    global_test: dict | None = None


def _core_group(Xs: np.ndarray, r: float):
    """Indices (into the sorted block Xs) of the core group, or None."""
    N = Xs.shape[0]
    cs = np.cumsum(Xs, axis=0)
    cq = np.cumsum(Xs ** 2, axis=0)
    # Pairwise tests (i, i+1) to find the first converging pair.
    s2 = Xs[:-1] + Xs[1:]
    q2 = Xs[:-1] ** 2 + Xs[1:] ** 2
    _, tpair = logt_from_H(H_from_sums(np.full(N - 1, 2.0), s2, q2), r)
    ok = np.where(tpair > CRIT)[0]
    if len(ok) == 0:
        return None
    i0 = ok[0]
    # Prefix groups i0..i0+k-1, k = 2..N-i0, all at once.
    base_s = cs[i0 - 1] if i0 > 0 else 0.0
    base_q = cq[i0 - 1] if i0 > 0 else 0.0
    s = cs[i0 + 1:] - base_s
    q = cq[i0 + 1:] - base_q
    k = np.arange(2, N - i0 + 1)
    _, tk = logt_from_H(H_from_sums(k.astype(float), s, q), r)
    fail = np.where(tk <= CRIT)[0]
    run = tk[: fail[0]] if len(fail) else tk
    kstar = int(k[int(np.argmax(np.where(np.isinf(run), 1e12, run)))])
    return np.arange(i0, i0 + kstar)


def _sieve(Xs: np.ndarray, core: np.ndarray, r: float, cstar: float = 0.0):
    """Add non-core units whose log t with the core exceeds c*, then shrink
    the candidate set (raising c*) until the whole club passes the test."""
    N = Xs.shape[0]
    rest = np.setdiff1d(np.arange(N), core)
    sc, qc, nc = Xs[core].sum(0), (Xs[core] ** 2).sum(0), len(core)
    if len(rest) == 0:
        return core
    s = sc[None, :] + Xs[rest]
    q = qc[None, :] + Xs[rest] ** 2
    _, tj = logt_from_H(H_from_sums(np.full(len(rest), nc + 1.0), s, q), r)
    cand = rest[tj > cstar]
    if len(cand) == 0:
        return core
    order = np.argsort(-tj[tj > cstar], kind="stable")
    cand = cand[order]
    # Club = core + first m candidates; raising c* = shortening the prefix.
    s = sc[None, :] + np.cumsum(Xs[cand], axis=0)
    q = qc[None, :] + np.cumsum(Xs[cand] ** 2, axis=0)
    m = np.arange(1, len(cand) + 1)
    _, tm = logt_from_H(H_from_sums(nc + m.astype(float), s, q), r)
    good = np.where(tm > CRIT)[0]
    if len(good) == 0:
        return core
    return np.concatenate([core, cand[: good[-1] + 1]])


def find_clubs(X: np.ndarray, r: float = 1 / 3, sort_cols: int = 1, cstar: float = 0.0):
    """Phillips & Sul (2009) clustering. X: N x T positive panel.

    Units are ordered by the mean of their last `sort_cols` observations.
    Returns integer labels (1..K, 0 = divergent) in descending club order."""
    X = np.asarray(X, dtype=float)
    N = X.shape[0]
    labels = np.zeros(N, dtype=int)
    remaining = np.arange(N)
    club = 0
    while len(remaining) >= 2:
        Xr = X[remaining]
        # Step 4 shortcut: if the remaining units converge, they form a club.
        if logt(Xr, r)["t"] > CRIT:
            club += 1
            labels[remaining] = club
            break
        order = np.argsort(-Xr[:, -sort_cols:].mean(1), kind="stable")
        Xs = Xr[order]
        core = _core_group(Xs, r)
        if core is None:
            break  # no converging pair left: the rest diverge
        members = _sieve(Xs, core, r, cstar)
        club += 1
        idx = remaining[order[members]]
        labels[idx] = club
        remaining = np.setdiff1d(remaining, idx)
    return labels


def club_stats(X: np.ndarray, labels: np.ndarray, r: float = 1 / 3) -> list:
    out = []
    for c in sorted(set(labels) - {0}):
        d = logt(X[labels == c], r)
        d["club"] = int(c)
        out.append(d)
    if (labels == 0).sum() >= 2:
        d = logt(X[labels == 0], r)
        d["club"] = 0
        out.append(d)
    return out


def merge_clubs(X: np.ndarray, labels: np.ndarray, r: float = 1 / 3):
    """Phillips & Sul (2009) merging: test adjacent clubs jointly and merge
    whenever the joint log t statistic exceeds -1.65, repeating until no
    adjacent pair can be merged. Returns (new_labels, log)."""
    lab = labels.copy()
    logs = []
    changed = True
    while changed:
        changed = False
        ids = sorted(set(lab) - {0})
        for a, b in zip(ids[:-1], ids[1:]):
            d = logt(X[(lab == a) | (lab == b)], r)
            logs.append({"club_a": a, "club_b": b, "t": d["t"], "beta": d["beta"], "merged": d["t"] > CRIT})
            if d["t"] > CRIT:
                lab[lab == b] = a
                changed = True
                break
        # renumber consecutively
        ids = sorted(set(lab) - {0})
        remap = {c: i + 1 for i, c in enumerate(ids)}
        lab = np.array([remap.get(c, 0) for c in lab])
    return lab, logs


def run(X: np.ndarray, r: float = 1 / 3, merge: bool = True) -> ClubResult:
    lab0 = find_clubs(X, r)
    res = ClubResult(labels=lab0, initial_labels=lab0.copy())
    res.global_test = logt(X, r)
    res.initial_clubs = club_stats(X, lab0, r)
    if merge:
        lab, logs = merge_clubs(X, lab0, r)
        res.labels, res.merge_log = lab, logs
    res.clubs = club_stats(X, res.labels, r)
    return res


def transition_paths(X: np.ndarray) -> np.ndarray:
    """Relative transition paths h_it = X_it / cross-sectional mean."""
    return X / X.mean(axis=0, keepdims=True)


# --------------------------------------------------------------------------
# Monte Carlo check of the implementation (size and club recovery)
# --------------------------------------------------------------------------
def simulate_panel(N, T, thetas, shares, sigma=0.3, alpha=0.5, rho=0.4, seed=0):
    """PS (2009, section 4) style DGP: X_it = theta_it * mu_t with
    theta_it = theta_c + sigma * xi_it / (log(t+1) * t^alpha), xi_it AR(1).
    Returns X (positive) and the true club labels."""
    rng = np.random.default_rng(seed)
    counts = np.round(np.array(shares) * N).astype(int)
    counts[-1] = N - counts[:-1].sum()
    lab = np.repeat(np.arange(1, len(thetas) + 1), counts)
    theta = np.repeat(thetas, counts)
    xi = np.zeros((N, T))
    eps = rng.normal(size=(N, T))
    xi[:, 0] = eps[:, 0]
    for t in range(1, T):
        xi[:, t] = rho * xi[:, t - 1] + eps[:, t]
    tt = np.arange(1, T + 1)
    th = theta[:, None] + sigma * xi / (np.log(tt + 1) * tt ** alpha)[None, :]
    mu = np.exp(0.02 * tt)
    return th * mu[None, :] + 5.0, lab  # shift keeps X positive


def monte_carlo(reps=200, N=300, T=15, seed=1) -> pd.DataFrame:
    """Rejection rate of the log t test under convergence / divergence and
    club recovery (adjusted Rand index) for a three-club DGP."""
    from sklearn.metrics import adjusted_rand_score

    rows = []
    rng = np.random.default_rng(seed)
    for rep in range(reps):
        s = int(rng.integers(1e9))
        Xc, _ = simulate_panel(N, T, [1.0], [1.0], seed=s)
        Xd, _ = simulate_panel(N, T, [1.0], [1.0], alpha=0.0, seed=s)
        Xd = Xd + np.random.default_rng(s + 1).normal(0, 0.3, size=(N, 1)) * np.linspace(0, 1, T)[None]
        X3, lab = simulate_panel(N, T, [1.6, 1.0, 0.4], [0.3, 0.4, 0.3], seed=s)
        lab_hat = find_clubs(X3)
        lab_m, _ = merge_clubs(X3, lab_hat)
        rows.append({
            "rep": rep,
            "reject_under_convergence": logt(Xc)["t"] < CRIT,
            "reject_under_divergence": logt(Xd)["t"] < CRIT,
            "reject_global_3club": logt(X3)["t"] < CRIT,
            "n_clubs_found": len(set(lab_m) - {0}),
            "ari_3club": adjusted_rand_score(lab, lab_m),
        })
    return pd.DataFrame(rows)
