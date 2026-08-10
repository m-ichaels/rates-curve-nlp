"""python -m rcmm <curve|factors|compare|summary> [--key value ...]

  curve    --date D [--cmt file] [--lambda 1e-4] [--out file]   the same JSON as `rcmm curve`, from numpy
  factors  [--cmt file] [--windows a:b,c:d] [--out file]         the same JSON as `rcmm factors`
  compare  --cpp results/curve_2025-04-14.json [--cmt file]      rebuild that day in Python, report max |diff|
  summary  [--results dir]                                       the headline numbers from results/*.json
"""
from __future__ import annotations

import json
import sys

import numpy as np

from . import curve as C
from . import data as D
from .factors import FactorModel
from .results import headline


def _args(argv: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    i = 0
    while i < len(argv):
        if argv[i].startswith("--"):
            key = argv[i][2:]
            if i + 1 < len(argv) and not argv[i + 1].startswith("--"):
                out[key] = argv[i + 1]
                i += 2
            else:
                out[key] = "1"
                i += 1
        else:
            i += 1
    return out


def build_curve(row, lam: float = 1e-4, meetings: list[float] | None = None) -> dict:
    """Everything `rcmm curve` writes for one CMT row, computed in numpy."""
    inst = D.cmt_instruments(row)
    b = C.Bootstrap().fit(inst)
    kr = C.KernelRidge().fit(inst, lam)
    out = {"bootstrap_iterations": b.iterations, "bootstrap_max_reprice_error_bp": b.max_error_bp,
           "kernel_ridge_max_reprice_error_bp": kr.max_error_bp(inst), "kernel_ridge_lambda": lam}
    out["instruments"] = [{"T": i.T, "par": i.rate, "zero_boot": b.zero(i.T), "fwd_boot": b.forward(i.T),
                           "zero_kr": kr.zero(i.T), "fwd_kr": kr.forward(i.T)} for i in inst]
    grid = []
    t = 0.05
    while t <= 30.0001:
        grid.append({"t": t, "fwd_boot": b.forward(t), "fwd_kr": kr.forward(t), "zero_boot": b.zero(t), "zero_kr": kr.zero(t)})
        t += 0.05
    out["grid"] = grid
    out["roughness_boot"] = C.forward_roughness(b.forward, 30.0)
    out["roughness_kr"] = C.forward_roughness(kr.forward, 30.0)
    out["leave_one_out"] = C.leave_one_out(inst, lam)
    kr01 = []
    for T in (2.0, 5.0, 10.0, 30.0):
        c = b.par_yield(T)
        v = C.key_rate_dv01(inst, T, c)
        kr01.append({"T": T, "coupon": c, "krdv01_per_1": v.tolist(), "dv01_per_1": float(v.sum())})
    out["ladder_krdv01"] = kr01
    meetings = meetings or [0.12, 0.25, 0.37, 0.5, 0.62, 0.75, 0.87, 1.0]
    obs = [C.StepObs(0.0, i.T, np.log(1 + i.rate * i.T) / i.T) for i in inst if i.kind is C.InstKind.BILL]
    fe = C.StepFrontEnd().fit(meetings, obs, inst[0].rate)
    out["front_end"] = {"meetings": meetings, "levels": fe.level.tolist(), "source": "bills + fixing (proxy)"}
    return out


def _diff(a, b, path: str = "") -> list[tuple[str, float]]:
    """Flatten two JSON trees into (path, |a - b|) for every numeric leaf present in both."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a.keys() & b.keys():
            out += _diff(a[k], b[k], f"{path}.{k}" if path else k)
    elif isinstance(a, list) and isinstance(b, list):
        for i, (x, y) in enumerate(zip(a, b)):
            out += _diff(x, y, f"{path}[{i}]")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        out.append((path, abs(float(a) - float(b))))
    return out


def cmd_curve(a: dict[str, str]) -> int:
    cmt = D.read_cmt(a.get("cmt"))
    date = a.get("date", cmt.index[-1].strftime("%Y-%m-%d"))
    row = cmt.loc[date]
    out = {"date": date, **build_curve(row, float(a.get("lambda", 1e-4)))}
    with open(a.get("out", f"results/curve_{date}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"date {date} bootstrap max reprice {out['bootstrap_max_reprice_error_bp']:.3g} bp ({out['bootstrap_iterations']} it); "
          f"KR max {out['kernel_ridge_max_reprice_error_bp']:.3g} bp; roughness boot {out['roughness_boot']:.3g} kr {out['roughness_kr']:.3g}")
    for e in out["leave_one_out"]:
        print(f"  LOO T={e['T']:g} boot {e['loo_err_bp_boot']:+.3f} bp, kr {e['loo_err_bp_kr']:+.3f} bp")
    return 0


def cmd_factors(a: dict[str, str]) -> int:
    cmt = D.read_cmt(a.get("cmt"))
    windows = a.get("windows", "2015-01-01:2026-12-31,2023-01-01:2026-12-31,2022-01-01:2022-12-31,2020-01-01:2022-12-31").split(",")
    out = {}
    for w in windows:
        d0, d1 = w.split(":")
        m = FactorModel().fit(cmt, d0, d1)
        out[w] = m.to_dict()
        print(f"{w}: {m.n_days} days, explained {', '.join(f'{100 * m.explained(k):.1f}%' for k in (1, 2, 3))}")
    with open(a.get("out", "results/factors.json"), "w") as f:
        json.dump(out, f, indent=1)
    return 0


def cmd_compare(a: dict[str, str]) -> int:
    with open(a["cpp"]) as f:
        cpp = json.load(f)
    cmt = D.read_cmt(a.get("cmt"))
    py = {"date": cpp["date"], **build_curve(cmt.loc[cpp["date"]], cpp["kernel_ridge_lambda"])}
    diffs = _diff(cpp, py)
    worst = sorted(diffs, key=lambda d: -d[1])[:8]
    print(f"{len(diffs)} numeric leaves compared for {cpp['date']}; max |diff| {worst[0][1]:.3g} at {worst[0][0]}")
    for p, d in worst[1:]:
        print(f"  {d:.3g}  {p}")
    tol = float(a.get("tol", 1e-8))
    return 0 if worst[0][1] <= tol else 1


def cmd_summary(a: dict[str, str]) -> int:
    print(json.dumps(headline(a.get("results")), indent=1, default=float))
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 1
    cmd, a = argv[0], _args(argv[1:])
    return {"curve": cmd_curve, "factors": cmd_factors, "compare": cmd_compare, "summary": cmd_summary}[cmd](a)


if __name__ == "__main__":
    sys.exit(main())
