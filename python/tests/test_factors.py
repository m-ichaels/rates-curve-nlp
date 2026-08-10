"""PCA factor model: orthonormal loadings, sign conventions, ladder covariance, and the same
decomposition as the C++ Jacobi sweep on the checked-in covariances (results/factors.json)."""
import numpy as np
import pytest

from rcmm import data as D
from rcmm.factors import FactorModel, inventory_risk


@pytest.fixture(scope="module")
def model(cmt):
    return FactorModel().fit(cmt, "2024-01-01", "2026-12-31")


def test_daily_changes_match_manual_diff(cmt):
    d = D.daily_changes_bp(cmt, D.FACTOR_TENORS, "2024-01-01", "2026-12-31")
    cols = [D.CMT_COLUMNS[i] for i in D.FACTOR_TENORS]
    manual = cmt[cols].diff().dropna() * 1e4
    assert len(d) == len(manual) - 0 or len(d) == len(manual)     # the sample has no gaps
    np.testing.assert_allclose(d.values, manual.loc[d.index].values)


def test_loadings_orthonormal_and_signed(model):
    L = model.loadings
    np.testing.assert_allclose(L @ L.T, np.eye(model.K), atol=1e-10)
    assert L[0].sum() > 0                                    # level positive
    assert L[1, -1] - L[1, 0] > 0                            # slope rises with maturity
    n = L.shape[1]
    assert L[2, n // 2] - 0.5 * (L[2, 0] + L[2, -1]) > 0     # curvature positive in the belly


def test_explained_variance_is_monotone_and_large(model):
    ex = [model.explained(k) for k in (1, 2, 3)]
    assert ex == sorted(ex) and ex[2] > 0.9 and ex[0] > 0.6
    assert model.evals[0] >= model.evals[1] >= model.evals[2] > 0


def test_ladder_cov_is_the_k_factor_reconstruction(model):
    S = model.ladder_cov(model.T)
    recon = sum(model.evals[k] * np.outer(model.loadings[k], model.loadings[k]) for k in range(model.K))
    np.testing.assert_allclose(S, recon, atol=1e-10)
    full = model.cov
    assert np.trace(S) / np.trace(full) == pytest.approx(model.explained(3), rel=1e-10)


def test_loading_at_interpolates_and_clamps(model):
    assert model.loading_at(0, 0.001) == model.loadings[0, 0]
    assert model.loading_at(0, 50.0) == model.loadings[0, -1]
    a, b = model.loading_at(1, 5.0), model.loading_at(1, 7.0)
    assert model.loading_at(1, 6.0) == pytest.approx(0.5 * (a + b))


def test_inventory_risk_quadratic(model):
    S = model.ladder_cov([2, 5, 10, 30])
    q = np.array([100.0, -50.0, 20.0, 0.0])
    assert inventory_risk(q, S) == pytest.approx(q @ S @ q)
    assert inventory_risk(2 * q, S) == pytest.approx(4 * inventory_risk(q, S))


def test_simulated_changes_have_the_model_covariance(model):
    X = model.simulate(200_000, seed=1)
    np.testing.assert_allclose(np.cov(X, rowvar=False), model.ladder_cov(model.T), atol=0.15 * model.evals[0] ** 0.5)


def test_round_trip_dict(model):
    m2 = FactorModel.from_dict(model.to_dict())
    np.testing.assert_allclose(m2.loadings, model.loadings)
    assert m2.n_days == model.n_days and m2.frm == model.frm


def test_cpp_eigendecomposition_matches(cpp_factors):
    """Same eigenvalues, loadings and signs as the Jacobi sweep, window by window."""
    for w, d in cpp_factors.items():
        m = FactorModel().from_cov(np.array(d["cov"]), k=len(d["loadings"]))
        np.testing.assert_allclose(m.evals, d["evals_bp2_per_day"], rtol=1e-8, atol=1e-6)
        np.testing.assert_allclose(m.loadings, d["loadings"], atol=1e-8)
        np.testing.assert_allclose([m.explained(k) for k in (1, 2, 3)], d["explained"], atol=1e-10)
