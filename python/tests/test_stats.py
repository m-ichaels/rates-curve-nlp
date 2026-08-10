"""Statistics toolkit: C++ conventions, resampling CIs cover the truth, permutation nulls behave."""
import numpy as np
import pytest

from rcmm import stats as S


def test_descriptives_match_cpp_conventions():
    v = [1.0, 2.0, 4.0, 8.0]
    assert S.mean(v) == 3.75 and S.variance(v) == pytest.approx(np.var(v, ddof=1))
    assert S.quantile(v, 0.5) == 3.0 and S.median([5.0]) == 5.0
    assert S.quantile([], 0.5) != S.quantile([], 0.5)       # NaN
    assert S.variance([1.0]) == 0.0 and S.kurtosis([1, 2, 3]) == 0.0
    assert S.kurtosis(np.random.default_rng(0).standard_normal(200_000)) == pytest.approx(0.0, abs=0.05)


def test_autocorr_of_ar1():
    rng = np.random.default_rng(1)
    x = np.zeros(50_000)
    for t in range(1, x.size):
        x[t] = 0.5 * x[t - 1] + rng.standard_normal()
    assert S.autocorr(x, 1) == pytest.approx(0.5, abs=0.02)
    assert S.autocorr(x, 2) == pytest.approx(0.25, abs=0.02)
    assert S.autocorr([1, 2], 5) == 0.0


def test_ks_and_wasserstein():
    rng = np.random.default_rng(2)
    a, b = rng.standard_normal(5000), rng.standard_normal(5000) + 1.0
    assert S.ks_distance(a, a) == 0.0 and S.ks_distance(a, []) == 1.0
    assert 0.3 < S.ks_distance(a, b) < 0.45                  # sup |Phi(x) - Phi(x - 1)| = 0.383
    assert S.wasserstein1(a, b) == pytest.approx(1.0, abs=0.05)


def test_bootstrap_ci_covers_the_mean():
    rng = np.random.default_rng(3)
    hits = 0
    for k in range(60):
        x = rng.standard_normal(80) + 0.3
        ci = S.bootstrap_mean_ci(x, B=500, seed=k)
        hits += ci.lo <= 0.3 <= ci.hi
        assert ci.lo <= ci.point <= ci.hi and ci.n == 80
    assert hits >= 48                                         # ~95 % nominal; 48/60 is a 4-sigma floor


def test_bootstrap_ci_degenerate():
    ci = S.bootstrap_mean_ci([1.0])
    assert ci.point == 1.0 and np.isnan(ci.lo) and ci.n == 1


def test_paired_bootstrap_std_ratio():
    rng = np.random.default_rng(4)
    x, y = rng.standard_normal(500), 2 * rng.standard_normal(500)
    ci = S.paired_bootstrap_ci(x, y, lambda a, b: a.std() / b.std(), B=400)
    assert ci.lo < 0.5 < ci.hi and ci.hi - ci.lo < 0.2


def test_block_bootstrap_wider_for_dependent_series():
    rng = np.random.default_rng(5)
    x = np.zeros(2000)
    for t in range(1, x.size):
        x[t] = 0.8 * x[t - 1] + rng.standard_normal()
    iid, blk = S.bootstrap_mean_ci(x, B=400), S.block_bootstrap_ci(x, block=20, B=400)
    assert (blk.hi - blk.lo) > 1.5 * (iid.hi - iid.lo)
    assert S.block_bootstrap_ci([1.0, 2.0, 3.0], block=10).n == 3   # falls back to iid


def test_permutation_p_value():
    rng = np.random.default_rng(6)
    x = rng.standard_normal(60)
    p_null, q95 = S.permutation_p_value(x, rng.standard_normal(60), n_perm=1000)
    p_sig, _ = S.permutation_p_value(x, x + 0.5 * rng.standard_normal(60), n_perm=1000)
    assert p_sig < 0.01 < p_null and 0.1 < q95 < 0.35


def test_ols_and_ridge():
    rng = np.random.default_rng(7)
    x = rng.standard_normal(500)
    y = 1.0 + 2.0 * x + 0.1 * rng.standard_normal(500)
    o = S.ols(x, y)
    assert o.a == pytest.approx(1.0, abs=0.02) and o.b == pytest.approx(2.0, abs=0.02) and o.r2 > 0.99
    assert S.ols([1, 2], [1, 2]) == S.Ols(0.0, 0.0, 0.0, 2)
    w = S.ridge(x[:, None], y, lam=0.0)
    np.testing.assert_allclose(w, [o.a, o.b], atol=1e-10)
    assert abs(S.ridge(x[:, None], y, lam=1e4)[1]) < abs(o.b)   # shrinks the slope, not the intercept


def test_walk_forward_never_sees_its_target():
    rng = np.random.default_rng(8)
    X = rng.standard_normal((400, 2))
    y = X @ [1.0, -1.0]
    pred = S.walk_forward(X, y, window=50, fit=lambda a, b: S.ridge(a, b, 0.0))
    assert np.isnan(pred[:50]).all() and np.isfinite(pred[50:]).all()
    np.testing.assert_allclose(pred[50:], y[50:], atol=1e-6)
    y2 = y.copy()
    y2[300] = 1e6                                             # an outlier is only felt after its date
    pred2 = S.walk_forward(X, y2, window=50, fit=lambda a, b: S.ridge(a, b, 0.0))
    np.testing.assert_allclose(pred2[:301], pred[:301], atol=1e-6)
    assert np.abs(pred2[301:351] - pred[301:351]).max() > 1.0


def test_logistic_recovers_hit_curve():
    """The tier hit model f(delta) = 1 / (1 + exp(alpha + beta delta)) with alpha = -1, beta = 2."""
    rng = np.random.default_rng(9)
    delta = rng.uniform(0, 2, 20_000)
    p = 1 / (1 + np.exp(-1 + 2 * delta))
    hit = (rng.uniform(size=delta.size) < p).astype(float)
    w = S.logistic_fit(np.column_stack([np.ones_like(delta), delta]), hit, ridge_lam=1e-6)
    assert w[0] == pytest.approx(1.0, abs=0.1) and w[1] == pytest.approx(-2.0, abs=0.1)
    assert S.logistic_predict(w, [1.0, 0.5]) == pytest.approx(1 / (1 + np.exp(0)), abs=0.03)
