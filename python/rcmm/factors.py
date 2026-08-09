"""Factor risk (port of include/rcmm/factors.hpp): PCA of daily par-yield changes into level / slope /
curvature, the ladder covariance used by the quoters, and the loadings interpolated to any maturity.
The C++ uses a Jacobi sweep; here numpy's eigh - same eigenpairs up to sign, which the sign conventions
below then fix the same way."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd

from .data import CMT_T, FACTOR_TENORS, daily_changes_bp


@dataclass
class FactorModel:
    K: int = 3
    tenor_idx: list[int] = field(default_factory=lambda: list(FACTOR_TENORS))
    T: np.ndarray = field(default_factory=lambda: CMT_T[FACTOR_TENORS].copy())
    cov: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))         # bp^2 / day
    evals: np.ndarray = field(default_factory=lambda: np.zeros(0))            # descending
    loadings: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))    # [factor, tenor], unit norm
    n_days: int = 0
    frm: str = ""
    to: str = ""

    # ---- fit ---------------------------------------------------------------------------------------------
    def fit(self, cmt: pd.DataFrame, d0: str, d1: str, k: int = 3) -> "FactorModel":
        X = daily_changes_bp(cmt, self.tenor_idx, d0, d1).values
        self.K, self.n_days, self.frm, self.to = k, X.shape[0], d0, d1
        self.cov = np.cov(X, rowvar=False, ddof=1) if X.shape[0] > 1 else np.zeros((X.shape[1],) * 2)
        return self._decompose()

    def from_cov(self, cov: np.ndarray, k: int = 3) -> "FactorModel":
        self.cov = np.asarray(cov, dtype=float)
        self.K = k
        return self._decompose()

    def _decompose(self) -> "FactorModel":
        evals, V = np.linalg.eigh(self.cov)
        order = np.argsort(evals)[::-1]
        self.evals, V = evals[order], V[:, order]
        n = self.T.size
        L = V[:, :self.K].T.copy()
        for f in range(self.K):   # level positive; slope rising with maturity; curvature positive in the belly
            if f == 0:
                s = L[f].sum()
            elif f == 1:
                s = L[f, n - 1] - L[f, 0]
            else:
                s = L[f, n // 2] - 0.5 * (L[f, 0] + L[f, n - 1])
            if s < 0:
                L[f] = -L[f]
        self.loadings = L
        return self

    # ---- queries -----------------------------------------------------------------------------------------
    def explained(self, k: int) -> float:
        tot = self.evals.sum()
        return float(self.evals[:k].sum() / tot) if tot > 0 else 0.0

    def loading_at(self, f: int, Tq: float) -> float:
        """Loading of factor f at maturity Tq, linear between the CMT tenors, flat outside."""
        return float(np.interp(Tq, self.T, self.loadings[f]))

    def ladder_cov(self, ladder_T: Sequence[float]) -> np.ndarray:
        """K-factor covariance of yield changes at the ladder maturities (bp^2 / day):
        Sigma = sum_k lambda_k L_k L_k^T with the loadings interpolated to the ladder."""
        L = np.array([[self.loading_at(f, T) for T in ladder_T] for f in range(self.K)])
        return (L.T * self.evals[:self.K]) @ L

    def factor_scores(self, dy_bp: np.ndarray) -> np.ndarray:
        """Project daily changes (n x tenors, bp) onto the K factors."""
        return np.asarray(dy_bp) @ self.loadings.T

    def simulate(self, n_days: int, seed: int = 0) -> np.ndarray:
        """Daily yield changes (n_days x tenors, bp) from the K-factor model with Gaussian factors."""
        rng = np.random.default_rng(seed)
        z = rng.standard_normal((n_days, self.K)) * np.sqrt(self.evals[:self.K])
        return z @ self.loadings

    def to_dict(self) -> dict:
        """Same layout as FactorModel::to_json in the C++ (results/factors.json entries)."""
        return {"from": self.frm, "to": self.to, "n_days": int(self.n_days), "T": self.T.tolist(),
                "evals_bp2_per_day": self.evals.tolist(), "loadings": self.loadings.tolist(),
                "explained": [self.explained(k) for k in range(1, self.K + 1)], "cov": self.cov.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> "FactorModel":
        m = cls(K=len(d["loadings"]))
        m.T = np.asarray(d["T"], dtype=float)
        m.cov = np.asarray(d["cov"], dtype=float)
        m.evals = np.asarray(d["evals_bp2_per_day"], dtype=float)
        m.loadings = np.asarray(d["loadings"], dtype=float)
        m.n_days, m.frm, m.to = int(d["n_days"]), d["from"], d["to"]
        return m


def inventory_risk(q_dv01: Sequence[float], sigma: np.ndarray) -> float:
    """q' Sigma q: daily P&L variance ($^2 / day) of a ladder inventory q in $/bp."""
    q = np.asarray(q_dv01, dtype=float)
    return float(q @ sigma @ q)
