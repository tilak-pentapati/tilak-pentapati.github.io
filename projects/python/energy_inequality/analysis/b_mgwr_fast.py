"""Vectorised MGWR for n of several thousand, numerically equivalent to the
`mgwr` package (Oshan et al. 2019) but feasible at n = 7,264.

Why: mgwr builds every kernel row with a Python-level call (~2-4 ms each).
Its inference step (Li & Fotheringham 2020) replays every backfitting
iteration for every covariate and every row, so a 14-covariate model with
~130 backfitting iterations at n = 7,264 would take more than 100 hours here.

What is reused from mgwr, unchanged:
* `mgwr.search.multi_bw` (backfitting loop, score-of-change criterion,
  `bws_same_times` rule) and `mgwr.search.golden_section` (bandwidth search,
  same bounds, tolerance and integer rounding as `Sel_BW`);
* `mgwr.gwr.MGWRResults` for every diagnostic (sigma2, AICc, ENP, standard
  errors, t-values, the da Silva-Fotheringham corrected critical t).

What is replaced: the kernel construction and the local fits (computed for
all locations at once from a pre-sorted distance matrix), and the inference
recursion of `MGWR._chunk_compute_R`, written as dense matrix products. The
recursion is the same; it is stopped early only when the remaining
bandwidth history is constant and five further sweeps change ENP_j and CCT
by less than 1e-10 (relative), i.e. the recursion has reached its fixed point.

`validate()` compares bandwidths, coefficients, ENP_j and standard errors
with the mgwr package on a subsample.
"""
from __future__ import annotations

import time

import numpy as np
from spglm.family import Gaussian

EPS = 1.0000001  # mgwr's adaptive-bandwidth inflation factor


class KernelCache:
    """Adaptive bisquare kernels (nearest-neighbour bandwidth), all rows at once."""

    def __init__(self, coords: np.ndarray, keep: int = 3):
        c = np.asarray(coords, dtype=float)
        sq = (c ** 2).sum(1)
        D2 = np.maximum(sq[:, None] + sq[None, :] - 2 * c @ c.T, 0.0)
        self.D = np.sqrt(D2)
        np.fill_diagonal(self.D, 0.0)
        self.Ds = np.sort(self.D, axis=1)
        self.n = len(c)
        self._cache: dict = {}
        self._order: list = []
        self.keep = keep

    def W(self, bw) -> np.ndarray:
        b = int(bw)
        if b in self._cache:
            return self._cache[b]
        t = self.Ds[:, b - 1] * EPS
        # (1 - (d/t)^2)^2 for d < t, else 0, computed in place.
        W = self.D * (1.0 / t)[:, None]
        np.square(W, out=W)
        np.subtract(1.0, W, out=W)
        np.maximum(W, 0.0, out=W)
        np.square(W, out=W)
        self._cache[b] = W
        self._order.append(b)
        if len(self._order) > self.keep:
            self._cache.pop(self._order.pop(0), None)
        return W


class _Lite:
    """The attributes of a lite GWR fit that multi_bw uses."""

    def __init__(self, params, predy, resid):
        self.params = params
        self.predy = predy
        self.resid_response = resid


def _uni_fit(K, y, x, bw):
    """Univariate local regression without intercept (all locations)."""
    W = K.W(bw)
    den = W @ (x * x)
    beta = (W @ (x * y)) / den
    predy = x * beta
    influ = x * x / den  # W_ii = 1
    return beta, predy, y - predy, influ


def _multi_fit(K, y, X, bw):
    """Multivariate local regression (all locations); returns params, predy,
    resid and the per-location inverse cross-product matrices."""
    W = K.W(bw)
    n, k = X.shape
    P = (X[:, :, None] * X[:, None, :]).reshape(n, k * k)
    XtWX = (W @ P).reshape(n, k, k)
    XtWy = W @ (X * y[:, None])
    M = np.linalg.inv(XtWX)
    beta = np.einsum("nab,nb->na", M, XtWy)
    predy = (X * beta).sum(1)
    return beta, predy, y - predy, M


def _aicc(y, predy, tr_s):
    n = len(y)
    llf = Gaussian().loglike(y.reshape(-1, 1), predy.reshape(-1, 1))
    return -2.0 * llf + 2.0 * n * (tr_s + 1.0) / (n - tr_s - 2.0)


class _FastSel:
    """Stand-in for the univariate Sel_BW created inside multi_bw."""

    def __init__(self, K, y, X):
        self.K, self.y, self.x = K, y.ravel(), X.ravel()
        self.sel_hist = []

    def search(self, bw_min=None, bw_max=None, tol=1e-6, max_iter=200):
        from mgwr.search import golden_section

        n = self.K.n
        a = 40 + 2 * 1 if bw_min is None else bw_min  # Sel_BW._init_section, one local variable
        c = n if (bw_max is None or bw_max is np.inf) else bw_max

        def f(bw):
            _, predy, _, influ = _uni_fit(self.K, self.y, self.x, bw)
            return _aicc(self.y, predy, influ.sum())

        out = golden_section(a, c, 0.38197, f, tol, max_iter, bw_max, True, False)
        self.sel_hist = out[2]
        return out[0]


class FastSelector:
    """Results of the multiscale bandwidth search, laid out like Sel_BW."""

    def __init__(self, bw_tuple):
        self.bw = bw_tuple
        self.params = bw_tuple[3]
        self.sel_hist = bw_tuple[-2]
        self.bw_init = bw_tuple[-1]


def search_multi(coords, y, X, init_bw, multi_bw_min=30, tol_multi=1e-5, max_iter_multi=200,
                 bws_same_times=5, tol=1e-6, max_iter=200, verbose=True, K=None):
    """MGWR backfitting with mgwr.search.multi_bw. X excludes the constant."""
    from mgwr.search import multi_bw

    K = K or KernelCache(coords)
    y = np.asarray(y, float).ravel()
    Xc = np.column_stack([np.ones(len(y)), np.asarray(X, float)])
    n, k = Xc.shape

    def gwr_func(yy, XX, bw):
        yy = yy.ravel()
        if XX.shape[1] == 1:
            beta, predy, resid, _ = _uni_fit(K, yy, XX.ravel(), bw)
            return _Lite(beta.reshape(-1, 1), predy.reshape(-1, 1), resid.reshape(-1, 1))
        beta, predy, resid, _ = _multi_fit(K, yy, XX, bw)
        return _Lite(beta, predy.reshape(-1, 1), resid.reshape(-1, 1))

    def bw_func(yy, XX):
        return _FastSel(K, yy, XX)

    def sel_func(sel, bw_min=None, bw_max=None):
        return sel.search(bw_min=bw_min, bw_max=bw_max, tol=tol, max_iter=max_iter)

    t0 = time.time()
    res = multi_bw(init_bw, y.reshape(-1, 1), Xc, n, k, Gaussian(), tol_multi, max_iter_multi, False, gwr_func,
                   bw_func, sel_func, [multi_bw_min] * k, [None] * k, bws_same_times, verbose=verbose)
    if verbose:
        print(f"backfitting: {len(res[1])} iterations, {time.time() - t0:.0f} s", flush=True)
    return FastSelector(res), K, Xc


def inference(sel: FastSelector, K: KernelCache, Xc: np.ndarray, verbose=True, rtol=1e-10):
    """ENP_j and CCT (n x k) by the recursion of MGWR._chunk_compute_R,
    computed for all columns at once with dense products."""
    n, k = Xc.shape
    t0 = time.time()
    W0 = K.W(sel.bw_init)
    P = (Xc[:, :, None] * Xc[:, None, :]).reshape(n, k * k)
    M = np.linalg.inv((W0 @ P).reshape(n, k, k))  # (X'W_iX)^-1
    pR = np.empty((k, n, n))
    for j in range(k):
        # pR_j[i, m] = X[i, j] * [(X'W_iX)^-1 X' W_i]_{j, m}
        pR[j] = Xc[:, j][:, None] * ((M[:, j, :] @ Xc.T) * W0)
    err = np.eye(n) - pR.sum(0)
    hist = np.asarray(sel.bw[1])
    prev = None
    used = 0
    for it in range(hist.shape[0]):
        for j in range(k):
            pRj_old = pR[j] + err
            W = K.W(hist[it, j])
            x = Xc[:, j]
            A = (x[:, None] * W * x[None, :]) / (W @ (x * x))[:, None]
            pR[j] = A @ pRj_old
            err = pRj_old - pR[j]
        used = it + 1
        last = it == hist.shape[0] - 1
        if not (last or (it + 1) % 5 == 0):
            continue  # ENP_j / CCT are only needed at the end and for the stopping check
        enp = np.array([np.trace(pR[j]) for j in range(k)])
        cct = np.stack([((pR[j] / Xc[:, j][:, None]) ** 2).sum(1) for j in range(k)], axis=1)
        if verbose:
            print(f"inference sweep {it + 1}/{hist.shape[0]}: ENP={enp.sum():.3f} ({time.time() - t0:.0f} s)",
                  flush=True)
        stable_after = bool(np.all(hist[it:] == hist[-1]))
        if prev is not None and stable_after:
            d = max(np.max(np.abs(enp - prev[0]) / np.abs(prev[0])), np.max(np.abs(cct - prev[1]) / prev[1]))
            if d < rtol:
                break
        prev = (enp, cct)
    return enp, cct, used


def local_cn(K: KernelCache, Xc: np.ndarray, bws) -> np.ndarray:
    """Local condition number, identical to MGWRResults.local_collinearity()[0]
    (columns x_j * w_ij with covariate-specific kernels, scaled to unit
    length), computed from the k x k Gram matrices instead of one full-matrix
    SVD per location."""
    n, k = Xc.shape
    Ws = [K.W(b) for b in bws]
    K._cache.clear()
    K._order.clear()
    G = np.empty((n, k, k))
    for a in range(k):
        for b in range(a, k):
            g = (Ws[a] * Ws[b]) @ (Xc[:, a] * Xc[:, b])
            G[:, a, b] = g
            G[:, b, a] = g
    d = np.sqrt(np.einsum("nii->ni", G))
    G = G / d[:, :, None] / d[:, None, :]
    ev = np.linalg.eigvalsh(G)
    return np.sqrt(ev[:, -1] / np.maximum(ev[:, 0], 1e-300))


def fit(coords, y, X, init_bw, name_x=None, verbose=True, sel=None):
    """Full MGWR: backfitting + inference, returning mgwr's MGWRResults."""
    from mgwr.gwr import MGWR, MGWRResults

    K = KernelCache(coords)
    y = np.asarray(y, float).ravel()
    Xc = np.column_stack([np.ones(len(y)), np.asarray(X, float)])
    if sel is None:
        sel, K, Xc = search_multi(coords, y, X, init_bw, verbose=verbose, K=K)
    enp, cct, used = inference(sel, K, Xc, verbose=verbose)
    model = MGWR(np.asarray(coords, float), y.reshape(-1, 1), np.asarray(X, float), sel, kernel="bisquare",
                 fixed=False, constant=True, name_x=name_x)
    params = sel.params
    predy = (Xc * params).sum(1).reshape(-1, 1)
    res = MGWRResults(model, params, predy, cct, enp, np.ones(len(y)), None, name_x)
    res.n_inference_sweeps = used
    res.local_cn_fast = local_cn(K, Xc, np.asarray(sel.bw[0]))
    return res, sel


def validate(coords, y, X, n_sub=600, seed=0, init_bw=None):
    """Fit the same subsample with mgwr and with this module; return the
    largest discrepancies."""
    from mgwr.gwr import MGWR
    from mgwr.sel_bw import Sel_BW

    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(len(y), n_sub, replace=False))
    c, yy, XX = coords[idx], y[idx].reshape(-1, 1), X[idx]
    yy = (yy - yy.mean()) / yy.std()
    XX = (XX - XX.mean(0)) / XX.std(0)
    gsel = Sel_BW(c, yy, XX)
    gbw = gsel.search() if init_bw is None else init_bw
    t = time.time()
    s = Sel_BW(c, yy, XX, multi=True)
    s.search(multi_bw_min=[30], init_multi=gbw)
    ref = MGWR(c, yy, XX, s).fit()
    t_ref = time.time() - t
    t = time.time()
    mine, msel = fit(c, yy.ravel(), XX, gbw, verbose=False)
    t_mine = time.time() - t
    return {
        "n": n_sub, "k": XX.shape[1] + 1, "bws_mgwr": list(map(float, s.bw[0])),
        "bws_fast": list(map(float, msel.bw[0])), "bws_equal": bool(np.all(s.bw[0] == msel.bw[0])),
        "iters_mgwr": int(len(s.bw[1])), "iters_fast": int(len(msel.bw[1])),
        "max_abs_diff_params": float(np.max(np.abs(ref.params - mine.params))),
        "max_rel_diff_enp_j": float(np.max(np.abs(ref.ENP_j - mine.ENP_j) / np.abs(ref.ENP_j))),
        "max_rel_diff_bse": float(np.max(np.abs(ref.bse - mine.bse) / ref.bse)),
        "aicc_mgwr": float(ref.aicc), "aicc_fast": float(mine.aicc),
        "seconds_mgwr": t_ref, "seconds_fast": t_mine,
    }
