"""Small statistics toolkit shared by the text, hedge and event scripts (the Python half of
include/rcmm/stats.hpp).  Everything here is deterministic given the seed; the C++ uses its own
xoshiro stream so bootstrap CIs agree in width, not in the third decimal."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np


# ---- descriptive (same conventions as the C++: sample variance, linear-interpolated quantiles) ----------
def mean(v: Sequence[float]) -> float:
    v = np.asarray(v, dtype=float)
    return float(v.mean()) if v.size else 0.0


def variance(v: Sequence[float]) -> float:
    v = np.asarray(v, dtype=float)
    return float(v.var(ddof=1)) if v.size >= 2 else 0.0


def stdev(v: Sequence[float]) -> float:
    return float(np.sqrt(variance(v)))


def quantile(v: Sequence[float], q: float) -> float:
    v = np.asarray(v, dtype=float)
    return float(np.quantile(v, q)) if v.size else float("nan")   # numpy's default = C++ linear interpolation


def median(v: Sequence[float]) -> float:
    return quantile(v, 0.5)


def kurtosis(v: Sequence[float]) -> float:
    """Excess kurtosis (population moments), 0 for n < 4."""
    v = np.asarray(v, dtype=float)
    if v.size < 4:
        return 0.0
    d = v - v.mean()
    s2 = (d * d).mean()
    return float((d ** 4).mean() / (s2 * s2) - 3.0) if s2 > 0 else 0.0


def autocorr(v: Sequence[float], lag: int) -> float:
    v = np.asarray(v, dtype=float)
    if v.size <= lag + 1:
        return 0.0
    d = v - v.mean()
    den = float((d * d).sum())
    return float((d[lag:] * d[:-lag]).sum() / den) if den > 0 and lag > 0 else (1.0 if lag == 0 else 0.0)


# ---- distances ------------------------------------------------------------------------------------------
def ks_distance(a: Sequence[float], b: Sequence[float]) -> float:
    """Two-sample Kolmogorov-Smirnov distance sup |F_a - F_b|."""
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    if a.size == 0 or b.size == 0:
        return 1.0
    grid = np.concatenate([a, b])
    fa = np.searchsorted(a, grid, side="right") / a.size
    fb = np.searchsorted(b, grid, side="right") / b.size
    return float(np.abs(fa - fb).max())


def wasserstein1(a: Sequence[float], b: Sequence[float], grid: int = 512) -> float:
    """1-Wasserstein distance between empirical distributions, both quantile functions on a common grid."""
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    if a.size == 0 or b.size == 0:
        return float("nan")
    q = (np.arange(grid) + 0.5) / grid
    qa = a[np.minimum(a.size - 1, (q * a.size).astype(int))]
    qb = b[np.minimum(b.size - 1, (q * b.size).astype(int))]
    return float(np.abs(qa - qb).mean())


# ---- resampling -----------------------------------------------------------------------------------------
@dataclass
class CI:
    point: float
    lo: float
    hi: float
    n: int

    def as_list(self) -> list[float]:
        return [self.lo, self.hi]

    def __str__(self) -> str:
        return f"{self.point:.4g} [{self.lo:.4g}, {self.hi:.4g}] (n={self.n})"


def bootstrap_ci(x: Sequence[float], stat: Callable[[np.ndarray], float] = np.mean, B: int = 2000,
                 alpha: float = 0.05, seed: int = 7) -> CI:
    """Percentile bootstrap CI of `stat` over i.i.d. resamples of x."""
    x = np.asarray(x, dtype=float)
    ci = CI(float(stat(x)) if x.size else float("nan"), float("nan"), float("nan"), int(x.size))
    if x.size < 2:
        return ci
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(B, x.size))
    s = np.array([stat(x[i]) for i in idx])
    ci.lo, ci.hi = float(np.quantile(s, alpha / 2)), float(np.quantile(s, 1 - alpha / 2))
    return ci


def bootstrap_mean_ci(x: Sequence[float], B: int = 2000, alpha: float = 0.05, seed: int = 7) -> CI:
    return bootstrap_ci(x, np.mean, B, alpha, seed)


def paired_bootstrap_ci(x: Sequence[float], y: Sequence[float], stat: Callable[[np.ndarray, np.ndarray], float],
                        B: int = 2000, alpha: float = 0.05, seed: int = 7) -> CI:
    """Bootstrap CI of a statistic of two paired samples (rows resampled together), e.g. a std ratio or a
    correlation - the CIs quoted in the README for the hedge ratios and the text nowcast."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = min(x.size, y.size)
    ci = CI(float(stat(x[:n], y[:n])) if n else float("nan"), float("nan"), float("nan"), n)
    if n < 3:
        return ci
    rng = np.random.default_rng(seed)
    s = np.array([stat(x[i], y[i]) for i in rng.integers(0, n, size=(B, n))])
    s = s[np.isfinite(s)]
    ci.lo, ci.hi = float(np.quantile(s, alpha / 2)), float(np.quantile(s, 1 - alpha / 2))
    return ci


def block_bootstrap_ci(x: Sequence[float], stat: Callable[[np.ndarray], float] = np.mean, block: int = 10,
                       B: int = 2000, alpha: float = 0.05, seed: int = 7) -> CI:
    """Moving-block bootstrap for serially dependent series (daily P&L, daily hedge errors)."""
    x = np.asarray(x, dtype=float)
    n = x.size
    ci = CI(float(stat(x)) if n else float("nan"), float("nan"), float("nan"), int(n))
    if n < 2 * block:
        return bootstrap_ci(x, stat, B, alpha, seed)
    rng = np.random.default_rng(seed)
    nb = int(np.ceil(n / block))
    starts = rng.integers(0, n - block + 1, size=(B, nb))
    offs = np.arange(block)
    s = np.array([stat(x[(st[:, None] + offs).ravel()[:n]]) for st in starts])
    ci.lo, ci.hi = float(np.quantile(s, alpha / 2)), float(np.quantile(s, 1 - alpha / 2))
    return ci


def permutation_p_value(x: Sequence[float], y: Sequence[float], stat: Callable[[np.ndarray, np.ndarray], float]
                        = lambda a, b: np.corrcoef(a, b)[0, 1], n_perm: int = 5000, seed: int = 7,
                        two_sided: bool = False) -> tuple[float, float]:
    """Permutation test of `stat(x, y)` against shuffles of y.  Returns (p, 95th percentile of the null):
    the placebo used for the text nowcast - how often a random pairing of statements and moves gives a
    correlation this large."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    rng = np.random.default_rng(seed)
    obs = float(stat(x, y))
    null = np.array([stat(x, rng.permutation(y)) for _ in range(n_perm)])
    p = float(np.mean(np.abs(null) >= abs(obs))) if two_sided else float(np.mean(null >= obs))
    return p, float(np.quantile(null, 0.95))


# ---- regression -----------------------------------------------------------------------------------------
@dataclass
class Ols:
    a: float
    b: float
    r2: float
    n: int


def ols(x: Sequence[float], y: Sequence[float]) -> Ols:
    """y = a + b x; {0, 0, 0} for n < 3, as the C++."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = min(x.size, y.size)
    if n < 3:
        return Ols(0.0, 0.0, 0.0, n)
    x, y = x[:n], y[:n]
    dx, dy = x - x.mean(), y - y.mean()
    sxx, sxy, syy = float(dx @ dx), float(dx @ dy), float(dy @ dy)
    b = sxy / sxx if sxx > 0 else 0.0
    r2 = sxy * sxy / (sxx * syy) if sxx > 0 and syy > 0 else 0.0
    return Ols(float(y.mean() - b * x.mean()), b, r2, n)


def ridge(X: np.ndarray, y: np.ndarray, lam: float = 1.0, intercept: bool = True) -> np.ndarray:
    """Ridge coefficients (intercept unpenalised).  Used for the learned hedge ratios and the text nowcast."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    if intercept:
        X = np.column_stack([np.ones(X.shape[0]), X])
    pen = lam * np.eye(X.shape[1])
    if intercept:
        pen[0, 0] = 0.0
    return np.linalg.solve(X.T @ X + pen, X.T @ y)


def walk_forward(X: np.ndarray, y: np.ndarray, window: int, fit: Callable[[np.ndarray, np.ndarray], np.ndarray]
                 = lambda a, b: ridge(a, b, 1.0), intercept: bool = True) -> np.ndarray:
    """One-step-ahead predictions: the model is fitted on the trailing `window` rows and applied to the
    next row, so no prediction ever sees its own target.  NaN for the first `window` rows."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    out = np.full(y.shape[0], np.nan)
    for t in range(window, y.shape[0]):
        w = fit(X[t - window:t], y[t - window:t])
        xt = np.concatenate([[1.0], X[t]]) if intercept else X[t]
        out[t] = float(xt @ w)
    return out


def logistic_fit(X: np.ndarray, y: np.ndarray, ridge_lam: float = 1e-3, iters: int = 50) -> np.ndarray:
    """Logistic regression by Newton-Raphson with a ridge penalty; X rows already include the intercept 1."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    w = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-X @ w))
        g = X.T @ (p - y) + ridge_lam * w
        H = (X * (p * (1 - p))[:, None]).T @ X + ridge_lam * np.eye(X.shape[1])
        step = np.linalg.solve(H, g)
        w -= step
        if np.abs(step).max() < 1e-10:
            break
    return w


def logistic_predict(w: np.ndarray, x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.asarray(x, float) @ w))
