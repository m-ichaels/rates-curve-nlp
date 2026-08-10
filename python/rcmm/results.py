"""Loaders for everything the pipeline writes to results/: each returns tidy pandas frames so the
notebook, the summary and the report all read the same tables.  Keys follow the JSON the C++ CLI and
the tools/*.py scripts emit; nothing is recomputed here."""
from __future__ import annotations

import glob
import json
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .data import ROOT


def _load(name: str, results_dir: str | None = None) -> dict | None:
    p = os.path.join(results_dir or os.path.join(ROOT, "results"), name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


# ---- rcmm curve ---------------------------------------------------------------------------------------
@dataclass
class CurveDay:
    date: str
    lambda_: float
    instruments: pd.DataFrame     # T, par, zero_boot, fwd_boot, zero_kr, fwd_kr
    grid: pd.DataFrame            # t, fwd_boot, fwd_kr, zero_boot, zero_kr
    leave_one_out: pd.DataFrame   # T, loo_err_bp_boot, loo_err_bp_kr
    krdv01: pd.DataFrame          # rows = ladder bond, columns = bumped input tenor (per 1 notional per bp)
    front_end: pd.DataFrame       # meeting, level
    reprice_boot_bp: float
    reprice_kr_bp: float
    roughness_boot: float
    roughness_kr: float
    bootstrap_iterations: int

    @property
    def loo_interior(self) -> pd.DataFrame:
        """Tenors with 0.2 < T <= 10 (3m-10y), the filter scripts/summarize.py applies: the endpoints are
        extrapolation and excluded from the headline."""
        return self.leave_one_out[(self.leave_one_out["T"] > 0.2) & (self.leave_one_out["T"] <= 10)]


def curve_day(date: str | None = None, results_dir: str | None = None) -> CurveDay | None:
    """results/curve_<date>.json (latest when date is None)."""
    d = results_dir or os.path.join(ROOT, "results")
    files = sorted(glob.glob(os.path.join(d, "curve_*.json")))
    if date:
        files = [f for f in files if f.endswith(f"curve_{date}.json")]
    if not files:
        return None
    with open(files[-1]) as f:
        j = json.load(f)
    tenors = [i["T"] for i in j["instruments"]]
    kr = pd.DataFrame([k["krdv01_per_1"] for k in j["ladder_krdv01"]], columns=[f"{t:g}y" for t in tenors],
                      index=[f"{k['T']:g}y" for k in j["ladder_krdv01"]])
    kr.insert(0, "dv01_per_1", [k["dv01_per_1"] for k in j["ladder_krdv01"]])
    kr.insert(0, "coupon", [k["coupon"] for k in j["ladder_krdv01"]])
    fe = j["front_end"]
    front = pd.DataFrame({"meeting": list(fe["meetings"]) + [np.nan], "level": fe["levels"]})
    return CurveDay(j["date"], j["kernel_ridge_lambda"], pd.DataFrame(j["instruments"]), pd.DataFrame(j["grid"]),
                    pd.DataFrame(j["leave_one_out"]), kr, front, j["bootstrap_max_reprice_error_bp"],
                    j["kernel_ridge_max_reprice_error_bp"], j["roughness_boot"], j["roughness_kr"], j["bootstrap_iterations"])


# ---- rcmm curvehist -----------------------------------------------------------------------------------
def curvehist_methods(results_dir: str | None = None) -> pd.DataFrame | None:
    """One row per method (bootstrap, kernel ridge x lambda) with the walk-forward scores."""
    j = _load("curvehist.json", results_dir)
    if not j:
        return None
    rows = []
    for m in j["methods"]:
        r = {k: v for k, v in m.items() if k not in ("hedge_test", "loo_yield_bp_by_year")}
        for T, h in m["hedge_test"].items():
            r[f"hedged_{T}_std_bp"] = h["hedged_pnl_bp_std"]
            r[f"unhedged_{T}_std_bp"] = h["unhedged_bp_std"]
        rows.append(r)
    df = pd.DataFrame(rows)
    df["label"] = [m["method"] if m["lambda"] == 0 else f"kr {m['lambda']:.0e}" for m in j["methods"]]
    return df.set_index("label")


def curvehist_by_year(results_dir: str | None = None) -> pd.DataFrame | None:
    """Mean LOO yield error by year, one column per method."""
    j = _load("curvehist.json", results_dir)
    if not j:
        return None
    cols = {(m["method"] if m["lambda"] == 0 else f"kr {m['lambda']:.0e}"): m["loo_yield_bp_by_year"] for m in j["methods"]}
    return pd.DataFrame(cols).rename_axis("year")


def curvehist_days(results_dir: str | None = None) -> pd.DataFrame | None:
    """Per-day record: date, y10, dv01 / LOO / reprice per method (the hedge weight vectors are dropped)."""
    j = _load("curvehist.json", results_dir)
    if not j:
        return None
    rows = [{k: v for k, v in d.items() if not isinstance(v, list)} for d in j["days"]]
    df = pd.DataFrame(rows)
    df["date"] = pd.to_datetime(df["date"])
    return df.set_index("date")


def selected_lambda(results_dir: str | None = None, max_reprice_bp: float = 1.0) -> float | None:
    """The rule fixed in advance: the smoothest kernel ridge (largest lambda) whose max repricing error
    over all days is within `max_reprice_bp`."""
    j = _load("curvehist.json", results_dir)
    if not j:
        return None
    ok = [m["lambda"] for m in j["methods"] if m["method"] == "kernel_ridge" and m["max_reprice_bp"] <= max_reprice_bp]
    return max(ok) if ok else None


# ---- tools/hedge_test.py ------------------------------------------------------------------------------
def hedge_test(results_dir: str | None = None) -> pd.DataFrame | None:
    """Rows = (target tenor, hedge method): std of next-day hedged P&L in bp, n, reduction vs unhedged."""
    j = _load("hedge_test.json", results_dir)
    if not j:
        return None
    rows = []
    for T, tgt in j["targets"].items():
        for method, r in tgt.items():
            if isinstance(r, dict) and "std_bp" in r:
                rows.append({"target": T, "method": method, **r})
    df = pd.DataFrame(rows).set_index(["target", "method"])
    df.attrs.update(window=j["window"], test_from=j["test_from"], ridge=j["ridge"])
    return df


def hedge_test_ratio_ci(results_dir: str | None = None) -> pd.DataFrame | None:
    """learned / model-bootstrap std ratio with its bootstrap CI, per target."""
    j = _load("hedge_test.json", results_dir)
    if not j:
        return None
    rows = [{"target": T, "ratio": r["point"], "lo": r["ci95"][0], "hi": r["ci95"][1]}
            for T, tgt in j["targets"].items() if (r := tgt.get("learned_over_model_boot_std_ratio"))]
    return pd.DataFrame(rows).set_index("target")


def hedge_test_by_year(results_dir: str | None = None) -> pd.DataFrame | None:
    j = _load("hedge_test.json", results_dir)
    if not j:
        return None
    frames = {T: pd.DataFrame(tgt["std_by_year"]) for T, tgt in j["targets"].items()}
    return pd.concat(frames, axis=1).rename_axis("year")


# ---- rcmm factors -------------------------------------------------------------------------------------
def factors(results_dir: str | None = None) -> dict[str, dict] | None:
    """window label -> FactorModel::to_json dict; use FactorModel.from_dict to rehydrate."""
    return _load("factors.json", results_dir)


def factor_loadings(results_dir: str | None = None, window: str | None = None) -> pd.DataFrame | None:
    j = factors(results_dir)
    if not j:
        return None
    w = j[window or next(iter(j))]
    return pd.DataFrame(np.array(w["loadings"]).T, index=[f"{t:g}y" for t in w["T"]],
                        columns=["level", "slope", "curvature"][:len(w["loadings"])])


# ---- rcmm mm / events ---------------------------------------------------------------------------------
def frontier(results_dir: str | None = None) -> pd.DataFrame | None:
    """Closed-form frontier on synthetic days: one row per (quoter, gamma) with P&L, risk, decomposition."""
    j = _load("mm_frontier.json", results_dir)
    if not j:
        return None
    rows = []
    for c in j["cells"]:
        r = {k: v for k, v in c.items() if not isinstance(v, (list, dict))}
        r.update({f"decomp_{k}": v for k, v in c["decomposition"].items()})
        for i, h in enumerate(c["hit_ratio_by_tier"]):
            r[f"hit_tier{i}"] = h
        rows.append(r)
    df = pd.DataFrame(rows)
    df.attrs.update(days=j["days"], seeds=j["seeds"], ladder_dv01=j["ladder_dv01"], Sigma_ladder=j["Sigma_ladder"])
    return df


def markouts(name: str = "mm_frontier.json", results_dir: str | None = None) -> pd.DataFrame | None:
    """Mark-outs by (quoter, gamma / variant, tier, pre-event flag)."""
    j = _load(name, results_dir)
    if not j:
        return None
    rows = []
    for c in j["cells"]:
        key = {k: c[k] for k in ("quoter", "gamma", "variant") if k in c}
        rows += [{**key, **m} for m in c["markouts"]]
    return pd.DataFrame(rows)


def events(name: str = "events.json", results_dir: str | None = None) -> pd.DataFrame | None:
    """Real-day replay: one row per (quoter, variant)."""
    j = _load(name, results_dir)
    if not j:
        return None
    rows = []
    for c in j["cells"]:
        r = {k: v for k, v in c.items() if not isinstance(v, (list, dict))}
        r.update({f"decomp_{k}": v for k, v in c["decomposition"].items()})
        rows.append(r)
    df = pd.DataFrame(rows).set_index(["quoter", "variant"])
    df.attrs.update({k: j[k] for k in ("n_days", "n_event_days", "n_measured_jumps", "from", "to", "exclude_2020_22",
                                       "gamma", "events_by_type")})
    return df


# ---- tools/nlp_features.py / text_features.py --------------------------------------------------------
def nlp_signal(results_dir: str | None = None) -> pd.DataFrame | None:
    """Out-of-sample nowcast scores: rows = (feature set, target window), columns corr / CI / p / R2 ..."""
    j = _load("nlp_signal.json", results_dir)
    if not j:
        return None
    rows = []
    for fs, targets in j["feature_sets"].items():
        for tgt, r in targets.items():
            rr = {k: v for k, v in r.items() if not isinstance(v, list)}
            for k in ("corr_ci", "sign_hit_ci"):
                if k in r:
                    rr[f"{k}_lo"], rr[f"{k}_hi"] = r[k]
            rows.append({"features": fs, "target": tgt, **rr})
    df = pd.DataFrame(rows).set_index(["features", "target"])
    df.attrs.update(train_end=j["train_end"], n_train=j["n_train"], n_test=j["n_test"],
                    dated_checkpoint=j["dated_checkpoint"], size_model=j["size_model"], models=j["models"])
    return df


def text_signal(scorer: str = "", results_dir: str | None = None) -> dict | None:
    """results/text_signal[_<scorer>].json: the single-feature hawkishness regression with its placebo."""
    return _load(f"text_signal{'_' + scorer if scorer else ''}.json", results_dir)


def text_signal_table(results_dir: str | None = None) -> pd.DataFrame | None:
    """All scorers side by side."""
    d = results_dir or os.path.join(ROOT, "results")
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "text_signal*.json"))):
        with open(f) as fh:
            j = json.load(fh)
        rows.append({k: v for k, v in j.items() if not isinstance(v, (dict, list))} | {"b1_lo": j["b1_ci"][0], "b1_hi": j["b1_ci"][1]})
    return pd.DataFrame(rows).set_index("scorer") if rows else None


def headline(results_dir: str | None = None) -> dict:
    """The numbers the README quotes, in one dict, so the summary and the report cannot drift apart."""
    out: dict = {}
    m = curvehist_methods(results_dir)
    lam = selected_lambda(results_dir)
    if m is not None and lam is not None:
        boot, kr = m.loc["bootstrap"], m.loc[f"kr {lam:.0e}"]
        out["curve"] = {"lambda": lam, "n_days": None,
                        "loo_rel": kr["loo_yield_bp_mean"] / boot["loo_yield_bp_mean"] - 1,
                        "loo_median_rel": kr["loo_yield_bp_median"] / boot["loo_yield_bp_median"] - 1,
                        "roughness_rel": kr["roughness_mean"] / boot["roughness_mean"] - 1,
                        "hedge_3y_rel": kr["hedged_3y_std_bp"] / boot["hedged_3y_std_bp"] - 1,
                        "hedge_7y_rel": kr["hedged_7y_std_bp"] / boot["hedged_7y_std_bp"] - 1}
    h = hedge_test(results_dir)
    if h is not None:
        out["hedge"] = {T: {"learned_vs_boot": h.loc[(T, "learned"), "std_bp"] / h.loc[(T, "model_boot"), "std_bp"] - 1,
                            "learned_vs_kr": h.loc[(T, "learned"), "std_bp"] / h.loc[(T, "model_kr"), "std_bp"] - 1}
                        for T in ("3Yr", "7Yr") if (T, "model_kr") in h.index}
    n = nlp_signal(results_dir)
    if n is not None and ("roberta", "r2") in n.index:
        r = n.loc[("roberta", "r2")]
        out["text"] = {"corr": r["corr"], "ci": [r["corr_ci_lo"], r["corr_ci_hi"]], "n": int(r["n"]),
                       "p": r.get("placebo_p_value"), "r2": r["r2_oos"]}
    e = events(results_dir=results_dir)
    if e is not None:
        out["events"] = e[["pnl_per_day", "pnl_per_var"]].to_dict("index")
    return out
