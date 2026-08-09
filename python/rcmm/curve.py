"""Yield curves in numpy: a line-for-line port of include/rcmm/curve.hpp.

Hagan-West (2006) monotone-convex forward interpolation with the closed-form sector integrals, the global
bootstrap of par instruments (bills, semi-annual par bonds, OIS), key-rate DV01s, the meeting-date step
front end, and the kernel-ridge discount curve of Filipovic, Pelger & Ye.  Times are years (ACT/365.25),
rates decimals, DV01s per 1 bp on unit notional in price units - the same conventions as the C++, so the
two implementations can be compared to 1e-9 in the tests.

The C++ is the production path (it does the 2,927-day walk-forward in seconds); this module exists so
that every number the C++ produces can be re-derived independently, and so the notebooks can build a
curve without a compiler.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Sequence

import numpy as np

EPS = 1e-14


def _llround(x: float) -> int:
    """std::llround for non-negative x (half away from zero)."""
    return int(math.floor(x + 0.5))


# ---- Hagan & West monotone convex interpolation of discrete forwards --------------------------------------
class MonotoneConvex:
    """Knots t_0 = 0 < t_1 < ... < t_n with discrete forwards fd_i on (t_{i-1}, t_i]."""

    def __init__(self) -> None:
        self.t = np.zeros(0)
        self.fd = np.zeros(0)
        self.f = np.zeros(0)

    def build(self, knots: Sequence[float], disc_fwd: Sequence[float]) -> "MonotoneConvex":
        t = np.asarray(knots, dtype=float)
        fd = np.asarray(disc_fwd, dtype=float)
        n = fd.size
        f = np.zeros(n + 1)
        if n == 1:
            f[:] = fd[0]
        else:
            f[1:n] = ((t[1:n] - t[0:n - 1]) * fd[1:n] + (t[2:n + 1] - t[1:n]) * fd[0:n - 1]) / (t[2:n + 1] - t[0:n - 1])
            f[0] = fd[0] - 0.5 * (f[1] - fd[0])
            f[n] = fd[n - 1] - 0.5 * (f[n - 1] - fd[n - 1])
        self.t, self.fd, self.f = t, fd, f
        return self

    # g(x) on one interval, x in [0, 1]: the four sectors of Hagan & West (2006), degenerate cases guarded
    @staticmethod
    def g(g0: float, g1: float, x: float) -> float:
        if abs(g0) < EPS and abs(g1) < EPS:
            return 0.0
        if (g0 > 0 and -0.5 * g0 >= g1 >= -2 * g0) or (g0 < 0 and -0.5 * g0 <= g1 <= -2 * g0):
            return g0 * (1 - 4 * x + 3 * x * x) + g1 * (-2 * x + 3 * x * x)
        if (g0 < 0 and g1 > -2 * g0) or (g0 > 0 and g1 < -2 * g0):
            eta = (g1 + 2 * g0) / (g1 - g0)
            if x <= eta or 1 - eta < EPS:
                return g0
            return g0 + (g1 - g0) * ((x - eta) / (1 - eta)) ** 2
        if (g0 > 0 > g1 > -0.5 * g0) or (g0 < 0 < g1 < -0.5 * g0):
            eta = 3 * g1 / (g1 - g0)
            if x >= eta or eta < EPS:
                return g1
            return g1 + (g0 - g1) * ((eta - x) / eta) ** 2
        # same sign (or one of them zero)
        if abs(g0 + g1) < EPS:
            return 0.0
        eta = g1 / (g0 + g1)
        A = -g0 * g1 / (g0 + g1)
        if x < eta:
            return A if eta < EPS else A + (g0 - A) * ((eta - x) / eta) ** 2
        return A if 1 - eta < EPS else A + (g1 - A) * ((x - eta) / (1 - eta)) ** 2

    # int_0^x g, sector by sector (closed form; int_0^1 g = 0 so each interval reproduces its discrete forward)
    @staticmethod
    def G(g0: float, g1: float, x: float) -> float:
        cube = lambda v: v * v * v
        if abs(g0) < EPS and abs(g1) < EPS:
            return 0.0
        if (g0 > 0 and -0.5 * g0 >= g1 >= -2 * g0) or (g0 < 0 and -0.5 * g0 <= g1 <= -2 * g0):
            return g0 * (x - 2 * x * x + x * x * x) + g1 * (-x * x + x * x * x)
        if (g0 < 0 and g1 > -2 * g0) or (g0 > 0 and g1 < -2 * g0):
            eta = (g1 + 2 * g0) / (g1 - g0)
            if 1 - eta < EPS:
                return g0 * x
            return g0 * x + ((g1 - g0) * (1 - eta) / 3 * cube((x - eta) / (1 - eta)) if x > eta else 0.0)
        if (g0 > 0 > g1 > -0.5 * g0) or (g0 < 0 < g1 < -0.5 * g0):
            eta = 3 * g1 / (g1 - g0)
            if eta < EPS:
                return g1 * x
            return g1 * x + (g0 - g1) * eta / 3 * (1 - cube((eta - min(x, eta)) / eta))
        if abs(g0 + g1) < EPS:
            return 0.0
        eta = g1 / (g0 + g1)
        A = -g0 * g1 / (g0 + g1)
        acc = A * x
        if eta >= EPS:
            acc += (g0 - A) * eta / 3 * (1 - cube((eta - min(x, eta)) / eta))
        if 1 - eta >= EPS and x > eta:
            acc += (g1 - A) * (1 - eta) / 3 * cube((x - eta) / (1 - eta))
        return acc

    def forward(self, tt: float) -> float:
        n = self.fd.size
        if n == 0:
            return 0.0
        if tt <= 0:
            return float(self.f[0])
        if tt >= self.t[n]:
            return float(self.f[n])
        i = max(1, int(np.searchsorted(self.t, tt, side="left")))   # t[i-1] < tt <= t[i]
        x = (tt - self.t[i - 1]) / (self.t[i] - self.t[i - 1])
        return float(self.fd[i - 1] + self.g(self.f[i - 1] - self.fd[i - 1], self.f[i] - self.fd[i - 1], x))

    def integral(self, tt: float) -> float:
        """int_0^tt f, exact on each knot interval; flat extrapolation past the last knot."""
        n = self.fd.size
        if n == 0:
            return 0.0
        acc = 0.0
        t, fd, f = self.t, self.fd, self.f
        for i in range(1, n + 1):
            if not t[i - 1] < tt:
                break
            w = t[i] - t[i - 1]
            x = min(1.0, (tt - t[i - 1]) / w)
            acc += w * (fd[i - 1] * x + self.G(f[i - 1] - fd[i - 1], f[i] - fd[i - 1], x))
        if tt > t[n]:
            acc += f[n] * (tt - t[n])
        return float(acc)

    def discount(self, tt: float) -> float:
        return math.exp(-self.integral(tt))


# ---- instruments -----------------------------------------------------------------------------------------
class InstKind(Enum):
    BILL = 0        # bond-equivalent simple yield, prices to 1
    PAR_BOND = 1    # semi-annual coupon = par rate
    OIS = 2         # annual fixed


@dataclass(frozen=True)
class Instrument:
    kind: InstKind
    T: float
    rate: float
    freq: int = 2

    def bumped(self, bp: float = 1.0) -> "Instrument":
        return Instrument(self.kind, self.T, self.rate + bp * 1e-4, self.freq)


DiscountFn = Callable[[float], float]


def price_instrument(ins: Instrument, P: DiscountFn) -> float:
    """Price of a par instrument off a discount function; equals 1 when the curve reprices it."""
    if ins.kind is InstKind.BILL:
        return P(ins.T) * (1 + ins.rate * ins.T)
    if ins.kind is InstKind.PAR_BOND:
        n = _llround(ins.T * ins.freq)
        v = sum(ins.rate / ins.freq * P(k / ins.freq) for k in range(1, n + 1))
        return v + P(ins.T)
    n = _llround(ins.T)
    if n < 1:
        return P(ins.T) * (1 + ins.rate * ins.T)
    return sum(ins.rate * P(float(k)) for k in range(1, n + 1)) + P(ins.T)


# ---- bootstrap -------------------------------------------------------------------------------------------
class Bootstrap:
    """Discrete forwards on the instrument maturities so that every instrument prices at par under the
    monotone-convex interpolant.  The interpolant couples neighbouring intervals, so the sequential secant
    solves are swept to a global fixed point (tol on the max forward change per sweep; 1e-10 = 1e-6 bp)."""

    def __init__(self) -> None:
        self.mc = MonotoneConvex()
        self.inst: list[Instrument] = []
        self.iterations = 0
        self.max_error_bp = 0.0

    def fit(self, instruments: Sequence[Instrument], tol: float = 1e-10, max_iter: int = 60) -> "Bootstrap":
        inst = sorted(instruments, key=lambda i: i.T)
        n = len(inst)
        knots = np.zeros(n + 1)
        knots[1:] = [i.T for i in inst]
        fd = np.array([i.rate for i in inst], dtype=float)   # start from the par rates
        mc = self.mc.build(knots, fd)

        def err(i: int, x: float) -> float:
            fd[i] = x
            mc.build(knots, fd)
            return price_instrument(inst[i], mc.discount) - 1.0

        self.iterations = 0
        while self.iterations < max_iter:
            maxchg = 0.0
            for i in range(n):
                start = x0 = float(fd[i])
                e0 = err(i, x0)
                if abs(e0) <= 1e-14:
                    continue
                x1 = x0 + 1e-4
                e1 = err(i, x1)
                for _ in range(30):
                    if abs(e1) <= 1e-14 or abs(e1 - e0) < 1e-18:
                        break
                    x2 = x1 - e1 * (x1 - x0) / (e1 - e0)
                    x0, e0 = x1, e1
                    x1 = x2
                    e1 = err(i, x1)
                maxchg = max(maxchg, abs(x1 - start))
            self.iterations += 1
            if maxchg < tol:
                break
        self.inst = inst
        self.max_error_bp = max(abs(price_instrument(i, mc.discount) - 1.0) * 1e4 for i in inst)
        return self

    def P(self, tt: float) -> float:
        return self.mc.discount(tt)

    def zero(self, tt: float) -> float:
        return -math.log(self.P(tt)) / tt if tt > 0 else self.mc.forward(0.0)

    def forward(self, tt: float) -> float:
        return self.mc.forward(tt)

    def par_yield(self, T: float, freq: int = 2) -> float:
        n = _llround(T * freq)
        a = sum(self.P(k / freq) / freq for k in range(1, n + 1))
        return (1 - self.P(T)) / a


def key_rate_dv01(inst: Sequence[Instrument], T: float, coupon: float, freq: int = 2) -> np.ndarray:
    """Price change of a par bond (unit notional) per 1 bp bump of each input par rate, rebootstrapping
    the curve for each bump.  Negative for rising rates; the entries sum to the parallel DV01."""
    base = Bootstrap().fit(inst)
    bond = Instrument(InstKind.PAR_BOND, T, coupon, freq)
    p0 = price_instrument(bond, base.P)
    kr = np.zeros(len(inst))
    for k in range(len(inst)):
        up = list(inst)
        up[k] = up[k].bumped(1.0)
        kr[k] = p0 - price_instrument(bond, Bootstrap().fit(up).P)
    return kr


# ---- meeting-date step-function front end ----------------------------------------------------------------
@dataclass
class StepObs:
    """An average-forward observation over [a, b]: a bill (0, T, ln(1 + yT)/T) or a futures strip."""
    a: float
    b: float
    avg: float


class StepFrontEnd:
    """Instantaneous forward piecewise constant between policy meeting dates, least squares to a set of
    average-rate observations with a small ridge on the step changes (keeps unobserved steps continuous)
    and today's fixing anchoring the first step."""

    def __init__(self) -> None:
        self.meeting_t = np.zeros(0)
        self.level = np.zeros(1)

    def fit(self, meetings: Sequence[float], obs: Sequence[StepObs], sofr_today: float, lam: float = 1e-3) -> "StepFrontEnd":
        m = np.asarray(meetings, dtype=float)
        K = m.size + 1
        bounds = np.concatenate([[0.0], m, [1e9]])
        rows, y = [], []
        for o in obs:
            row = np.zeros(K)
            for k in range(K):
                lo, hi = max(bounds[k], o.a), min(bounds[k + 1], o.b)
                if hi > lo:
                    row[k] += (hi - lo) / (o.b - o.a)
            rows.append(row)
            y.append(o.avg)
        anchor = np.zeros(K)
        anchor[0] = 1
        rows.append(anchor)
        y.append(sofr_today)
        for k in range(1, K):
            row = np.zeros(K)
            row[k], row[k - 1] = lam, -lam
            rows.append(row)
            y.append(0.0)
        A, yv = np.array(rows), np.array(y)
        self.meeting_t = m
        self.level = np.linalg.solve(A.T @ A, A.T @ yv)   # normal equations, as the C++
        return self

    def forward(self, tt: float) -> float:
        k = int(np.searchsorted(self.meeting_t, tt, side="right"))
        return float(self.level[min(k, self.level.size - 1)])

    def integral(self, tt: float) -> float:
        acc, a = 0.0, 0.0
        for k, lev in enumerate(self.level):
            b = self.meeting_t[k] if k < self.meeting_t.size else 1e9
            if tt <= a:
                break
            acc += lev * (min(b, tt) - a)
            a = b
        return acc

    def P(self, tt: float) -> float:
        return math.exp(-self.integral(tt))


# ---- kernel ridge discount curve (Filipovic, Pelger & Ye) --------------------------------------------------
def sobolev_kernel(s: float, t: float) -> float:
    """Reproducing kernel of int d''(t)^2 dt with d(0) = 1: min^2 (3 max - min) / 6."""
    a, b = min(s, t), max(s, t)
    return a * a * (3 * b - a) / 6.0


class KernelRidge:
    """d(t) = 1 - beta t - sum_j w_j k(t, t_j) on the cash-flow grid; the linear term is the unpenalised
    null space of the curvature penalty.  Weights minimise sum_i (price_i - 1)^2 / T_i + lambda w'Kw;
    pricing is linear in d, so the fit is one linear solve."""

    def __init__(self) -> None:
        self.tj = np.zeros(0)
        self.w = np.zeros(0)
        self.lambda_ = 1e-6
        self.beta = 0.0

    @staticmethod
    def cashflow_grid(inst: Sequence[Instrument]) -> np.ndarray:
        times: list[float] = []
        for ins in inst:
            if ins.kind is InstKind.PAR_BOND:
                times += [k / ins.freq for k in range(1, _llround(ins.T * ins.freq) + 1)]
            elif ins.kind is InstKind.OIS:
                times += [float(k) for k in range(1, _llround(ins.T) + 1)] + [ins.T]
            else:
                times.append(ins.T)
        times.sort()
        out: list[float] = []
        for x in times:
            if not out or abs(x - out[-1]) >= 1e-9:
                out.append(x)
        return np.array(out)

    def fit(self, inst: Sequence[Instrument], lam: float = 1e-6) -> "KernelRidge":
        self.lambda_ = lam
        tj = self.cashflow_grid(inst)
        m, n = tj.size, len(inst)
        idx = lambda tt: int(np.searchsorted(tj, tt - 1e-9, side="left"))
        C = np.zeros((n, m))
        wt = np.ones(n)
        for i, ins in enumerate(inst):
            wt[i] = 1.0 / max(ins.T, 0.05)
            if ins.kind is InstKind.PAR_BOND:
                for k in range(1, _llround(ins.T * ins.freq) + 1):
                    C[i, idx(k / ins.freq)] += ins.rate / ins.freq
                C[i, idx(ins.T)] += 1
            elif ins.kind is InstKind.OIS:
                for k in range(1, _llround(ins.T) + 1):
                    C[i, idx(float(k))] += ins.rate
                C[i, idx(ins.T)] += 1
            else:
                C[i, idx(ins.T)] += 1 + ins.rate * ins.T
        b = C.sum(axis=1) - 1.0                          # residual when d == 1
        ka, kb = np.minimum.outer(tj, tj), np.maximum.outer(tj, tj)
        K = ka * ka * (3 * kb - ka) / 6.0
        G = np.column_stack([C @ K, C @ tj])             # [C K | C t], unknowns (w, beta)
        M = (G * wt[:, None]).T @ G
        M[:m, :m] += lam * K
        v = G.T @ (wt * b)
        try:
            u = np.linalg.solve(M, v)
        except np.linalg.LinAlgError:
            u = np.linalg.lstsq(M, v, rcond=None)[0]
        self.tj, self.w, self.beta = tj, u[:m], float(u[m])
        return self

    def P(self, tt: float) -> float:
        a = np.minimum(tt, self.tj)
        b = np.maximum(tt, self.tj)
        return float(1 - self.beta * tt - np.dot(self.w, a * a * (3 * b - a) / 6.0))

    def zero(self, tt: float) -> float:
        return -math.log(self.P(tt)) / tt

    def forward(self, tt: float, h: float = 1e-4) -> float:
        lo = max(tt - h, 0.0)
        return -(math.log(self.P(tt + h)) - math.log(self.P(lo))) / (tt + h - lo)

    def max_error_bp(self, inst: Sequence[Instrument]) -> float:
        return max(abs(price_instrument(ins, self.P) - 1.0) * 1e4 for ins in inst)


# ---- diagnostics ------------------------------------------------------------------------------------------
def forward_roughness(fwd: Callable[[float], float], T: float, steps: int = 600) -> float:
    """int_0^T f'(t)^2 dt by central differences on a uniform grid (same grid as the C++)."""
    h = T / steps
    acc = 0.0
    for k in range(1, steps):
        d = (fwd((k + 1) * h) - fwd((k - 1) * h)) / (2 * h)
        acc += d * d * h
    return acc


def leave_one_out(inst: Sequence[Instrument], lam: float = 1e-4) -> list[dict]:
    """Drop each instrument, rebuild both curves from the rest, reprice it: error in price bp."""
    out = []
    for k, held in enumerate(inst):
        sub = [x for j, x in enumerate(inst) if j != k]
        bb = Bootstrap().fit(sub)
        kk = KernelRidge().fit(sub, lam)
        out.append({"T": held.T, "loo_err_bp_boot": (price_instrument(held, bb.P) - 1) * 1e4,
                    "loo_err_bp_kr": (price_instrument(held, kk.P) - 1) * 1e4})
    return out


def price_error_yield_bp(held: Instrument, P: DiscountFn) -> float:
    """Price error of `held` off P in yield bp: (price - 1) / (dP/dc), with dP/dc = -dP/dy the price
    change per 1 bp of coupon on the same curve.  Positive when the curve prices the instrument rich."""
    p = price_instrument(held, P)
    dpdc = price_instrument(held.bumped(1.0), P) - p
    return (p - 1.0) / dpdc if dpdc != 0 else float("nan")
