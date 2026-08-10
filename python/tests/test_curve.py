"""Curve engine: closed-form identities, repricing, DV01 sums, and a cross-check of every number in the
checked-in C++ output (results/curve_2025-04-14.json) against the numpy port."""
import math

import numpy as np
import pytest

from rcmm import curve as C
from rcmm.cli import build_curve

BOOT_TOL = 1e-6      # bp; the two bootstraps agree to ~5e-8
KR_TOL = 1e-3        # bp; the kernel-ridge normal equations are ill-conditioned once the 30y is dropped


# ---- Hagan-West sectors -------------------------------------------------------------------------------
@pytest.mark.parametrize("g0,g1", [(0.01, -0.003), (0.01, -0.03), (-0.01, 0.002), (0.01, 0.02), (-0.01, -0.005),
                                   (0.0, 0.01), (0.01, 0.0), (0.0, 0.0), (0.01, -0.01)])
def test_sector_integral_matches_quadrature(g0, g1):
    """G(x) = int_0^x g for every sector, and int_0^1 g = 0 so each interval reproduces its discrete forward."""
    xs = np.linspace(0, 1, 2001)
    gs = np.array([C.MonotoneConvex.g(g0, g1, x) for x in xs])
    for x in (0.2, 0.5, 0.8, 1.0):
        k = int(round(x * 2000))
        quad = np.trapezoid(gs[:k + 1], xs[:k + 1])
        assert C.MonotoneConvex.G(g0, g1, x) == pytest.approx(quad, abs=2e-7)
    assert abs(C.MonotoneConvex.G(g0, g1, 1.0)) < 1e-12


def test_interpolant_reproduces_discrete_forwards():
    mc = C.MonotoneConvex().build([0, 0.5, 1, 2, 5, 10], [0.04, 0.042, 0.045, 0.047, 0.05])
    for i in range(1, 6):
        avg = (mc.integral(mc.t[i]) - mc.integral(mc.t[i - 1])) / (mc.t[i] - mc.t[i - 1])
        assert avg == pytest.approx(mc.fd[i - 1], abs=1e-12)
    assert mc.forward(0.0) == mc.f[0] and mc.forward(20.0) == mc.f[-1]


def test_single_interval_is_flat():
    mc = C.MonotoneConvex().build([0, 1], [0.03])
    assert mc.forward(0.3) == 0.03 and mc.discount(1.0) == pytest.approx(math.exp(-0.03))


# ---- bootstrap ------------------------------------------------------------------------------------------
def test_bootstrap_reprices_all_inputs(inst):
    b = C.Bootstrap().fit(inst)
    assert b.max_error_bp < 1e-6
    for i in inst:
        assert C.price_instrument(i, b.P) == pytest.approx(1.0, abs=1e-10)
    assert 1 <= b.iterations <= 10


def test_bootstrap_par_yields_are_the_inputs(inst):
    b = C.Bootstrap().fit(inst)
    for i in inst:
        if i.kind is C.InstKind.PAR_BOND:
            assert b.par_yield(i.T) == pytest.approx(i.rate, abs=1e-12)


def test_ois_leg_prices_to_par():
    inst = [C.Instrument(C.InstKind.OIS, T, r, 1) for T, r in ((1, 0.04), (2, 0.041), (5, 0.043), (10, 0.045))]
    b = C.Bootstrap().fit(inst)
    assert b.max_error_bp < 1e-6


def test_key_rate_dv01s_sum_to_parallel_dv01(inst):
    b = C.Bootstrap().fit(inst)
    for T in (2.0, 10.0):
        c = b.par_yield(T)
        kr = C.key_rate_dv01(inst, T, c)
        up = [i.bumped(1.0) for i in inst]
        p0 = C.price_instrument(C.Instrument(C.InstKind.PAR_BOND, T, c), b.P)
        p1 = C.price_instrument(C.Instrument(C.InstKind.PAR_BOND, T, c), C.Bootstrap().fit(up).P)
        assert kr.sum() == pytest.approx(p0 - p1, rel=1e-3)      # equal up to the cross-convexity of the bumps
        # an on-the-run par bond loads on its own tenor only
        own = [k for k, i in enumerate(inst) if i.T == T][0]
        assert kr[own] == pytest.approx(kr.sum(), rel=1e-6)
        assert np.abs(np.delete(kr, own)).max() < 1e-12


# ---- kernel ridge -------------------------------------------------------------------------------------
def test_kernel_is_the_sobolev_reproducing_kernel():
    """k(s, t) = int_0^min (s - u)(t - u) du: the reproducing kernel of int d''^2 with d(0) = 1."""
    for s, t in ((0.5, 2.0), (3.0, 3.0), (10.0, 1.0)):
        u = np.linspace(0, min(s, t), 20001)
        assert C.sobolev_kernel(s, t) == pytest.approx(np.trapezoid((s - u) * (t - u), u), rel=1e-8)


def test_kernel_ridge_reprices_within_rule(inst):
    kr = C.KernelRidge().fit(inst, 1e-4)
    assert kr.max_error_bp(inst) < 1.0
    assert kr.P(0.0) == pytest.approx(1.0)
    ps = [kr.P(t) for t in np.linspace(0.05, 30, 100)]
    assert all(p1 < p0 for p0, p1 in zip(ps, ps[1:]))   # discount falls in maturity on this day


def test_smoother_lambda_is_rougher_fit(inst):
    """More penalty: rougher repricing, smoother forwards - the trade-off the selection rule walks."""
    errs, rough = [], []
    for lam in (1e-6, 1e-4, 1e-2):
        kr = C.KernelRidge().fit(inst, lam)
        errs.append(kr.max_error_bp(inst))
        rough.append(C.forward_roughness(kr.forward, 30.0))
    assert errs == sorted(errs) and rough == sorted(rough, reverse=True)


# ---- cross-check against the C++ output --------------------------------------------------------------
@pytest.fixture(scope="module")
def py_curve(cmt, cpp_curve):
    return build_curve(cmt.loc[cpp_curve["date"]], cpp_curve["kernel_ridge_lambda"])


def test_cpp_bootstrap_matches(py_curve, cpp_curve):
    assert py_curve["bootstrap_iterations"] == cpp_curve["bootstrap_iterations"]
    assert py_curve["bootstrap_max_reprice_error_bp"] < 1e-6
    for a, b in zip(py_curve["instruments"], cpp_curve["instruments"]):
        assert a["zero_boot"] == pytest.approx(b["zero_boot"], abs=1e-10)
        assert a["fwd_boot"] == pytest.approx(b["fwd_boot"], abs=1e-10)
    for a, b in zip(py_curve["grid"], cpp_curve["grid"]):
        assert a["fwd_boot"] == pytest.approx(b["fwd_boot"], abs=1e-10)
    assert py_curve["roughness_boot"] == pytest.approx(cpp_curve["roughness_boot"], rel=1e-8)


def test_cpp_kernel_ridge_matches(py_curve, cpp_curve):
    assert py_curve["kernel_ridge_max_reprice_error_bp"] == pytest.approx(cpp_curve["kernel_ridge_max_reprice_error_bp"], abs=1e-6)
    for a, b in zip(py_curve["instruments"], cpp_curve["instruments"]):
        assert a["zero_kr"] == pytest.approx(b["zero_kr"], abs=1e-8)
        assert a["fwd_kr"] == pytest.approx(b["fwd_kr"], abs=1e-6)
    assert py_curve["roughness_kr"] == pytest.approx(cpp_curve["roughness_kr"], rel=1e-5)


def test_cpp_leave_one_out_matches(py_curve, cpp_curve):
    for a, b in zip(py_curve["leave_one_out"], cpp_curve["leave_one_out"]):
        assert a["T"] == pytest.approx(b["T"])
        assert a["loo_err_bp_boot"] == pytest.approx(b["loo_err_bp_boot"], abs=BOOT_TOL)
        assert a["loo_err_bp_kr"] == pytest.approx(b["loo_err_bp_kr"], abs=KR_TOL)


def test_cpp_key_rate_dv01s_match(py_curve, cpp_curve):
    for a, b in zip(py_curve["ladder_krdv01"], cpp_curve["ladder_krdv01"]):
        assert a["coupon"] == pytest.approx(b["coupon"], abs=1e-12)
        assert a["dv01_per_1"] == pytest.approx(b["dv01_per_1"], abs=1e-12)
        np.testing.assert_allclose(a["krdv01_per_1"], b["krdv01_per_1"], atol=1e-11)


def test_cpp_front_end_matches(py_curve, cpp_curve):
    np.testing.assert_allclose(py_curve["front_end"]["levels"], cpp_curve["front_end"]["levels"], atol=1e-10)


# ---- step front end -----------------------------------------------------------------------------------
def test_step_front_end_recovers_a_step():
    """Two bills spanning one meeting with a 25 bp hike: the fitted levels are the true steps."""
    meetings = [0.5]
    true = C.StepFrontEnd()
    true.meeting_t, true.level = np.array(meetings), np.array([0.04, 0.0425])
    obs = [C.StepObs(0.0, T, true.integral(T) / T) for T in (0.25, 0.75, 1.0)]
    fe = C.StepFrontEnd().fit(meetings, obs, 0.04, lam=1e-6)
    np.testing.assert_allclose(fe.level, true.level, atol=1e-6)
    assert fe.forward(0.49) == pytest.approx(0.04, abs=1e-6) and fe.forward(0.51) == pytest.approx(0.0425, abs=1e-6)


def test_price_error_yield_bp_sign(inst):
    """A curve that prices a bond rich (price > 1) reports a positive yield error and the right magnitude."""
    b = C.Bootstrap().fit(inst)
    bond = C.Instrument(C.InstKind.PAR_BOND, 10.0, b.par_yield(10.0) - 1e-4)   # coupon 1 bp below par
    e = C.price_error_yield_bp(bond, b.P)
    assert e == pytest.approx(-1.0, abs=1e-6)
