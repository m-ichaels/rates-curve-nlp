"""Results loaders read the checked-in results/ and the derived tables, and the headline numbers the
README quotes come out of them (so the summary and the report cannot drift from the JSON)."""
import numpy as np
import pytest

from rcmm import data as D
from rcmm import results as R


def test_curve_day_tables():
    c = R.curve_day("2025-04-14")
    assert c is not None and c.date == "2025-04-14" and c.bootstrap_iterations == 5
    assert list(c.instruments.columns) == ["T", "par", "zero_boot", "fwd_boot", "zero_kr", "fwd_kr"]
    assert len(c.grid) == 600 and c.grid["t"].iloc[-1] == pytest.approx(30.0)
    assert len(c.loo_interior) == 8 and c.loo_interior["T"].between(0.25, 10).all()
    assert c.krdv01.shape == (4, 2 + 12) and c.krdv01.index.tolist() == ["2y", "5y", "10y", "30y"]
    np.testing.assert_allclose(c.krdv01.iloc[:, 2:].sum(axis=1), c.krdv01["dv01_per_1"], atol=1e-12)
    assert len(c.front_end) == 9 and c.front_end["meeting"].isna().sum() == 1
    assert c.reprice_boot_bp < 1e-6 < c.reprice_kr_bp < 1.0


def test_curvehist_tables_and_selection_rule():
    m = R.curvehist_methods()
    assert m is not None and "bootstrap" in m.index and m.loc["bootstrap", "lambda"] == 0
    lam = R.selected_lambda()
    assert lam is not None and m.loc[f"kr {lam:.0e}", "max_reprice_bp"] <= 1.0
    bigger = [l for l in m["lambda"] if l > lam]
    assert all(m[m["lambda"] == l]["max_reprice_bp"].iloc[0] > 1.0 for l in bigger)
    by = R.curvehist_by_year()
    assert by is not None and by.shape[0] == 12 and "bootstrap" in by.columns
    days = R.curvehist_days()
    assert days is not None and len(days) == 2927 and days.index.is_monotonic_increasing


def test_hedge_test_tables():
    h = R.hedge_test()
    assert h is not None and ("3Yr", "learned") in h.index and h.attrs["window"] == 250
    for T in ("3Yr", "7Yr"):
        assert h.loc[(T, "learned"), "std_bp"] < h.loc[(T, "model_boot"), "std_bp"] < h.loc[(T, "unhedged"), "std_bp"]
    r = R.hedge_test_ratio_ci()
    assert r is not None and (r["lo"] < r["ratio"]).all() and (r["ratio"] < r["hi"]).all()
    assert R.hedge_test_by_year().shape[0] >= 10


def test_factor_tables(cpp_factors):
    L = R.factor_loadings(window="2015-01-01:2026-12-31")
    assert L is not None and L.shape == (12, 3) and list(L.columns) == ["level", "slope", "curvature"]
    assert (L["level"] > 0).all()


def test_frontier_events_markouts():
    f = R.frontier()
    assert f is not None and len(f) == 42 and {"quoter", "gamma", "pnl_per_day", "decomp_spread"} <= set(f.columns)
    assert f.attrs["days"] > 0 and len(f.attrs["ladder_dv01"]) == 4
    e = R.events()
    assert e is not None and len(e) == 8 and e.attrs["n_event_days"] > 0
    ex = R.events("events_ex2020.json")
    assert ex is not None and ex.attrs["exclude_2020_22"] is True
    mk = R.markouts()
    assert mk is not None and {"tier", "pre_event", "markout_30m_bp"} <= set(mk.columns) and len(mk) == 42 * 6


def test_nlp_and_text_signal_tables():
    n = R.nlp_signal()
    assert n is not None and ("roberta", "r2") in n.index
    r = n.loc[("roberta", "r2")]
    assert r["corr_ci_lo"] < r["corr"] < r["corr_ci_hi"] and r["n"] == 44
    assert n.attrs["train_end"] == "2023-12-31"
    t = R.text_signal_table()
    assert t is not None and len(t) >= 2 and "corr_test" in t.columns
    assert R.text_signal("dict")["scorer"] != R.text_signal("lorenzo")["scorer"]


def test_headline_numbers_match_the_readme():
    h = R.headline()
    assert -0.05 < h["curve"]["loo_rel"] < -0.02                          # -3.6 %
    assert -0.25 < h["curve"]["hedge_3y_rel"] < -0.15                     # -20 %
    assert -0.45 < h["hedge"]["3Yr"]["learned_vs_boot"] < -0.30           # -38 %
    assert -0.30 < h["hedge"]["7Yr"]["learned_vs_boot"] < -0.20           # -25 %
    assert 0.30 < h["text"]["corr"] < 0.38 and h["text"]["p"] < 0.05      # corr 0.34, p 0.022
    assert "events" in h and len(h["events"]) == 8


def test_derived_tables_load():
    ev = D.read_events()
    assert ev.index.is_monotonic_increasing and {"is_event", "type", "alpha_bp", "jump2y"} <= set(ev.columns)
    assert ev["is_event"].sum() > 100
    ts = D.read_text_scores()
    assert {"bank", "date", "fomc_rob_h", "dict_h"} <= set(ts.columns) and set(ts["bank"]) == {"FOMC", "ECB"}
    nf = D.read_nlp_features()
    assert {"rob_redline", "r2", "novelty"} <= set(nf.columns) and len(nf) > 100


def test_cmt_sample_loader(cmt):
    assert cmt.index[0].strftime("%Y-%m-%d") == "2024-01-02" and cmt.shape[1] == 13
    inst = D.cmt_instruments(cmt.loc["2025-04-14"])
    assert len(inst) == 12 and [i.T for i in inst] == sorted(i.T for i in inst)
    assert sum(i.kind.name == "BILL" for i in inst) == 5
    assert len(D.cmt_instruments(cmt.loc["2025-04-14"], drop_4m=False)) == 13
