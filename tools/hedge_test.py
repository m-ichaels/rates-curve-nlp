"""Which curve should a desk hedge with?  Walk-forward hedge test on the daily Treasury CMT curve, 2015-2026.
A bond at tenor T (3y, 7y, 20y) is hedged with its two neighbouring benchmark tenors.  Hedge ratios (yield-bp of T per bp of
each neighbour) come from
  model    : key-rate DV01s of the T bond off the curve built WITHOUT the T point (bootstrap / kernel ridge, from rcmm curvehist)
  learned  : ridge regression of dy_T on the neighbours' dy over the trailing window (rolling, no look-ahead)
  pca      : neutralise the level and slope exposures using the trailing-window PCA loadings
  naive    : 50/50 on the two neighbours (duration-matched)
Scored on the realised next-day hedged P&L in yield bp: std, and the reduction versus unhedged, on the test years.
Writes results/hedge_test.json.
"""
import os, sys, json
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
WINDOW = int(args.get("window", 250)); TEST_FROM = args.get("test_from", "2016-01-04"); RIDGE = float(args.get("ridge", 1.0))
TEN = ["1Mo", "2Mo", "3Mo", "6Mo", "1Yr", "2Yr", "3Yr", "5Yr", "7Yr", "10Yr", "20Yr", "30Yr"]   # cmt_instruments order (4Mo dropped)
NB = {"3Yr": ("2Yr", "5Yr"), "7Yr": ("5Yr", "10Yr"), "20Yr": ("10Yr", "30Yr")}

def main():
    cmt = pd.read_csv(os.path.join(ROOT, "data/raw/treasury_cmt.csv")); cmt["date"] = pd.to_datetime(cmt["date"]); cmt = cmt.sort_values("date").set_index("date")[TEN[4:]].dropna()   # 1y-30y (the 2m tenor only starts in 2018)
    dy = cmt.diff().dropna() * 100.0   # bp
    ch = json.load(open(os.path.join(ROOT, "results/curvehist.json"))) if os.path.exists(os.path.join(ROOT, "results/curvehist.json")) else None
    # model hedges: yesterday's key-rate weights (ex-T curve) applied to today's input changes, both from the curvehist record
    model_pnl = {}
    if ch:
        days = [d for d in ch["days"] if len(d.get("rates", [])) in (len(TEN), len(TEN) - 1)]
        for d0, d1 in zip(days[:-1], days[1:]):
            if len(d0["rates"]) != len(d1["rates"]): continue
            ten = TEN if len(d0["rates"]) == len(TEN) else [t for t in TEN if t != "2Mo"]   # the 2m tenor only exists from Oct 2018
            r0, r1 = np.array(d0["rates"]), np.array(d1["rates"]); dr = (r1 - r0) * 1e4
            for T in ("3y", "7y"):
                iT = ten.index(T[:-1] + "Yr"); idx = [j for j in range(len(ten)) if j != iT]
                for m in ("boot", "kr_2"):
                    key = f"hw_{T}_{m}"
                    if key in d0: model_pnl.setdefault((T, m), {})[pd.Timestamp(d1["date"])] = dr[iT] - sum(w * dr[j] for w, j in zip(d0[key], idx))
    out = {"window": WINDOW, "test_from": TEST_FROM, "ridge": RIDGE, "targets": {}}
    for T, (a, b) in NB.items():
        y = dy[T].values; X = dy[[a, b]].values; dates = dy.index
        res = {k: [] for k in ["unhedged", "naive", "learned", "pca", "model_boot", "model_kr"]}; used = []
        for i in range(WINDOW, len(y)):
            d = dates[i]
            if d < pd.Timestamp(TEST_FROM): continue
            mk = [(T[:-2] + "y", "boot"), (T[:-2] + "y", "kr_2")]
            if ch and mk[0] in model_pnl and not all(d in model_pnl[k] for k in mk): continue   # keep every method on the same days
            Xw, yw = X[i - WINDOW:i], y[i - WINDOW:i]
            w_l = np.linalg.solve(Xw.T @ Xw + RIDGE * np.eye(2), Xw.T @ yw)                       # ridge, fitted on the trailing window only
            W = dy.iloc[i - WINDOW:i].values; W = W - W.mean(0); C = W.T @ W / len(W); ev, V = np.linalg.eigh(C); L = V[:, ::-1][:, :2]   # level, slope loadings (trailing window)
            cols = list(dy.columns); ia, ib, iT = cols.index(a), cols.index(b), cols.index(T); A = np.array([[L[ia, 0], L[ib, 0]], [L[ia, 1], L[ib, 1]]]); w_p = np.linalg.solve(A, np.array([L[iT, 0], L[iT, 1]]))
            res["unhedged"].append(y[i]); res["naive"].append(y[i] - 0.5 * (X[i, 0] + X[i, 1])); res["learned"].append(y[i] - X[i] @ w_l); res["pca"].append(y[i] - X[i] @ w_p)
            for m, k in (("model_boot", (T[:-2] + "y", "boot")), ("model_kr", (T[:-2] + "y", "kr_2"))):
                if k in model_pnl and d in model_pnl[k]: res[m].append(model_pnl[k][d])
            used.append(d)
        tr = {k: {"std_bp": float(np.std(v)), "n": len(v), "reduction_vs_unhedged": float(1 - np.std(v) / np.std(res["unhedged"]))} for k, v in res.items() if v}
        # bootstrap CI on the std ratio learned / model_boot
        if res["model_boot"] and len(res["model_boot"]) == len(res["learned"]):   # same days by construction (both need a full record)
            rng = np.random.default_rng(0); l = np.array(res["learned"]); mb = np.array(res["model_boot"]); ratios = [np.std(l[ix]) / np.std(mb[ix]) for ix in (rng.integers(0, len(l), len(l)) for _ in range(1000))]
            tr["learned_over_model_boot_std_ratio"] = {"point": float(np.std(l) / np.std(mb)), "ci95": [float(np.percentile(ratios, 2.5)), float(np.percentile(ratios, 97.5))]}
        # by year
        by = {}
        for k in ("unhedged", "learned", "model_boot", "pca"):
            if res[k] and len(res[k]) == len(used):
                s = pd.Series(res[k], index=used); by[k] = {str(yr): float(g.std()) for yr, g in s.groupby(s.index.year)}
        tr["std_by_year"] = by; out["targets"][T] = tr
        print(T, {k: round(v["std_bp"], 3) for k, v in tr.items() if isinstance(v, dict) and "std_bp" in v})
    json.dump(out, open(os.path.join(ROOT, "results/hedge_test.json"), "w"), indent=1)

if __name__ == "__main__": main()
