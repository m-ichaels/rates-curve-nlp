#!/usr/bin/env python3
"""Figures from results/*.json.   python scripts/plots.py [results] [results/figures]"""
import glob, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1] if len(sys.argv) > 1 else "results"; F = sys.argv[2] if len(sys.argv) > 2 else os.path.join(R, "figures"); os.makedirs(F, exist_ok=True)
def load(n): p = os.path.join(R, n); return json.load(open(p)) if os.path.exists(p) else None
def save(fig, n): fig.tight_layout(); fig.savefig(os.path.join(F, n), dpi=130); plt.close(fig); print("wrote", n)
COL = {"Symmetric": "C7", "Bergault_oracle": "C8", "Bergault": "C0", "BarzykinCiceri": "C2", "BarzykinCiceri_QAHR": "C3", "CarteaWang": "C1"}

def curve():
    fs = sorted(glob.glob(os.path.join(R, "curve_*.json")))
    if not fs: return
    d = json.load(open(fs[-1])); g = d["grid"]; t = [x["t"] for x in g]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    ax[0, 0].plot(t, [100 * x["fwd_boot"] for x in g], label="monotone-convex bootstrap"); ax[0, 0].plot(t, [100 * x["fwd_kr"] for x in g], label="kernel ridge"); ax[0, 0].plot([i["T"] for i in d["instruments"]], [100 * i["par"] for i in d["instruments"]], "k.", label="par inputs (CMT)")
    ax[0, 0].set_xlabel("maturity (y)"); ax[0, 0].set_ylabel("%"); ax[0, 0].set_title(f"instantaneous forward, {d['date']}"); ax[0, 0].legend(fontsize=8)
    ax[0, 1].plot(t, [100 * x["zero_boot"] for x in g], label="bootstrap"); ax[0, 1].plot(t, [100 * x["zero_kr"] for x in g], label="kernel ridge"); ax[0, 1].set_title("zero rate"); ax[0, 1].set_xlabel("maturity (y)"); ax[0, 1].set_ylabel("%"); ax[0, 1].legend(fontsize=8)
    loo = [x for x in d["leave_one_out"] if 0.2 < x["T"] <= 10]; xs = np.arange(len(loo)); w = 0.4
    ax[1, 0].bar(xs - w / 2, [x["loo_err_bp_boot"] for x in loo], w, label="bootstrap"); ax[1, 0].bar(xs + w / 2, [x["loo_err_bp_kr"] for x in loo], w, label="kernel ridge"); ax[1, 0].set_xticks(xs); ax[1, 0].set_xticklabels([f"{x['T']:g}y" for x in loo]); ax[1, 0].set_ylabel("price error (bp)"); ax[1, 0].set_title("leave-one-out pricing error (interior tenors)"); ax[1, 0].legend(fontsize=8); ax[1, 0].axhline(0, color="k", lw=0.5)
    fe = d["front_end"]; m = [0] + fe["meetings"]; lv = fe["levels"]
    for k in range(len(lv)): ax[1, 1].plot([m[k], m[k + 1] if k + 1 < len(m) else m[-1] + 0.15], [100 * lv[k]] * 2, "C3", lw=2)
    ax[1, 1].plot([x["t"] for x in g if x["t"] <= 1.2], [100 * x["fwd_boot"] for x in g if x["t"] <= 1.2], "C0", alpha=0.6, label="bootstrap forward")
    for mm in fe["meetings"]: ax[1, 1].axvline(mm, color="grey", lw=0.5, ls=":")
    ax[1, 1].set_title(f"meeting-date step front end ({fe['source']})"); ax[1, 1].set_xlabel("years"); ax[1, 1].set_ylabel("%"); ax[1, 1].legend(fontsize=8)
    save(fig, "curve.png")

def curvehist():
    d = load("curvehist.json")
    if not d: return
    days = d["days"]; x = np.arange(len(days)); lab = [r["date"] for r in days]; tick = list(range(0, len(days), max(1, len(days) // 8)))
    ms = d["methods"]; kr = [m for m in ms if m["method"] == "kernel_ridge"]; ok = [m for m in kr if m["max_reprice_bp"] <= 1.0]; best = max(ok, key=lambda m: m["lambda"]) if ok else kr[-1]; ki = kr.index(best); tag = f"kr_{ki}"   # selection rule: smoothest KR that reprices within 1 bp
    fig, ax = plt.subplots(2, 2, figsize=(12, 7.5))
    ax[0, 0].plot(x, [r["dv01_10y_boot"] for r in days], lw=0.7, label="bootstrap"); ax[0, 0].plot(x, [r["dv01_10y_" + tag] for r in days], lw=0.7, label=f"kernel ridge (lambda {best['lambda']:.0e})"); ax[0, 0].set_ylabel("10y par DV01 (USD per bp per USD M)"); ax[0, 0].legend(fontsize=8); ax[0, 0].set_title("daily builds: 10y DV01", fontsize=10); ax[0, 0].set_xticks(tick); ax[0, 0].set_xticklabels([lab[i][:4] for i in tick], fontsize=7)
    yrs = sorted(ms[0]["loo_yield_bp_by_year"]); w = 0.8 / (1 + len(kr)); xs = np.arange(len(yrs))
    for k, m in enumerate(ms): ax[0, 1].bar(xs + (k - len(ms) / 2) * w + w / 2, [m["loo_yield_bp_by_year"][y] for y in yrs], w, label="bootstrap" if m["method"] == "bootstrap" else f"KR {m['lambda']:.0e}", color="k" if m["method"] == "bootstrap" else f"C{k}")
    ax[0, 1].set_xticks(xs); ax[0, 1].set_xticklabels(yrs, fontsize=7); ax[0, 1].set_ylabel("mean |LOO error| (yield bp)"); ax[0, 1].set_title("leave-one-out pricing error by year, interior tenors 6m-10y", fontsize=10); ax[0, 1].legend(fontsize=6, ncol=2)
    names = ["bootstrap"] + [f"KR {m['lambda']:.0e}" for m in kr]; xs = np.arange(len(names))
    ax[1, 0].bar(xs - 0.2, [m["hedge_test"]["3y"]["hedged_pnl_bp_std"] for m in ms], 0.4, label="3y bond hedged with the curve ex-3y"); ax[1, 0].bar(xs + 0.2, [m["hedge_test"]["7y"]["hedged_pnl_bp_std"] for m in ms], 0.4, label="7y bond hedged with the curve ex-7y")
    ax[1, 0].set_xticks(xs); ax[1, 0].set_xticklabels(names, fontsize=7, rotation=20); ax[1, 0].set_ylabel("std of next-day hedged P&L (yield bp)"); ax[1, 0].set_title(f"hedge test with each curve's key-rate DV01s (unhedged {ms[0]['hedge_test']['7y']['unhedged_bp_std']:.1f} bp)", fontsize=10); ax[1, 0].legend(fontsize=7)
    ax2 = ax[1, 1]; ax2.plot([m["roughness_mean"] for m in ms], [m["loo_yield_bp_mean"] for m in ms], "k.")
    for m, nme in zip(ms, names): ax2.annotate(nme, (m["roughness_mean"], m["loo_yield_bp_mean"]), fontsize=7)
    ax2.set_xlabel("mean forward roughness  int f'(t)^2 dt"); ax2.set_ylabel("mean |LOO error| (yield bp)"); ax2.set_title("smoothness vs out-of-sample pricing error", fontsize=10)
    save(fig, "curvehist.png")

def hedge():
    d = load("hedge_test.json")
    if not d: return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    keys = ["naive", "model_boot", "model_kr", "pca", "learned"]; labels = {"naive": "50/50 neighbours", "model_boot": "bootstrap key-rate", "model_kr": "kernel-ridge key-rate", "pca": "PCA level+slope neutral", "learned": "learned (ridge, trailing 250d)"}
    Ts = list(d["targets"]); xs = np.arange(len(Ts)); w = 0.8 / len(keys)
    for k, key in enumerate(keys): ax[0].bar(xs + (k - len(keys) / 2) * w + w / 2, [d["targets"][T].get(key, {}).get("std_bp", np.nan) for T in Ts], w, label=labels[key])
    ax[0].set_xticks(xs); ax[0].set_xticklabels([f"{T[:-2]}y bond" for T in Ts]); ax[0].set_ylabel("std of next-day hedged P&L (yield bp)"); ax[0].set_title(f"hedging with two neighbouring benchmarks, walk-forward from {d['test_from'][:4]}", fontsize=10); ax[0].legend(fontsize=7)
    T = "7Yr"; by = d["targets"][T]["std_by_year"]
    for key in ("model_boot", "pca", "learned"):
        if key in by: ax[1].plot(list(by[key]), list(by[key].values()), marker=".", label=labels[key])
    ax[1].set_ylabel("hedged P&L std (yield bp)"); ax[1].set_title("7y bond: hedge error by year", fontsize=10); ax[1].legend(fontsize=7); ax[1].tick_params(axis="x", labelsize=7)
    save(fig, "hedge.png")

def factors():
    d = load("factors.json")
    if not d: return
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
    for w, v in d.items():
        for k in range(3): ax[k].plot(v["T"], v["loadings"][k], marker=".", label=f"{w} (sd {np.sqrt(v['evals_bp2_per_day'][k]):.1f} bp/d)")
    for k, n in enumerate(["level", "slope", "curvature"]): ax[k].set_title(n); ax[k].set_xscale("log"); ax[k].set_xlabel("maturity (y)"); ax[k].axhline(0, color="k", lw=0.5); ax[k].legend(fontsize=6)
    save(fig, "factors.png")

def frontier():
    d = load("mm_frontier.json")
    if not d: return
    cells = d["cells"]; qs = sorted(set(c["quoter"] for c in cells), key=lambda q: list(COL).index(q) if q in COL else 9)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for q in qs:
        cs = sorted([c for c in cells if c["quoter"] == q], key=lambda c: c["gamma"])
        ax[0].errorbar([c["pnl_5min_std"] / 1e3 for c in cs], [c["pnl_per_day"] / 1e3 for c in cs], yerr=[[(c["pnl_per_day"] - c["pnl_lo"]) / 1e3 for c in cs], [(c["pnl_hi"] - c["pnl_per_day"]) / 1e3 for c in cs]], marker="o", ms=4, color=COL.get(q), label=q, capsize=2)
        for c in cs: ax[0].annotate(f"{c['gamma']:.0e}", (c["pnl_5min_std"] / 1e3, c["pnl_per_day"] / 1e3), fontsize=6, alpha=0.7)
    ax[0].set_xlabel("5-minute P&L std ($k)"); ax[0].set_ylabel("P&L per day ($k)"); ax[0].set_title("frontier over risk aversion (synthetic days)", fontsize=10); ax[0].legend(fontsize=7)
    g0 = sorted(set(c["gamma"] for c in cells))[len(set(c["gamma"] for c in cells)) // 2]; cs = [c for c in cells if abs(c["gamma"] - g0) < 1e-12]; x = np.arange(len(cs)); w = 0.27
    ax[1].bar(x - w, [c["decomposition"]["spread"] / 1e3 for c in cs], w, label="spread"); ax[1].bar(x, [c["decomposition"]["inventory"] / 1e3 for c in cs], w, label="inventory"); ax[1].bar(x + w, [-c["decomposition"]["hedge_cost"] / 1e3 for c in cs], w, label="-hedge cost")
    ax[1].plot(x, [c["pnl_per_day"] / 1e3 for c in cs], "k_", ms=18, label="total"); ax[1].set_xticks(x); ax[1].set_xticklabels([c["quoter"] for c in cs], rotation=20, fontsize=7); ax[1].set_ylabel("$k per day"); ax[1].set_title(f"P&L decomposition, gamma = {g0:.0e}"); ax[1].legend(fontsize=7); ax[1].axhline(0, color="k", lw=0.5)
    for i, c in enumerate(cs):
        mk = c["markouts"]; tiers = sorted(set(m["tier"] for m in mk))
        for j, tr in enumerate(tiers):
            pre = [m for m in mk if m["tier"] == tr and m["pre_event"]]; nor = [m for m in mk if m["tier"] == tr and not m["pre_event"]]
            ax[2].bar(i + (j - 1) * 0.25 - 0.06, nor[0]["markout_30m_bp"] if nor else 0, 0.12, color=f"C{j}", alpha=0.5); ax[2].bar(i + (j - 1) * 0.25 + 0.06, pre[0]["markout_30m_bp"] if pre else 0, 0.12, color=f"C{j}", label=f"{tr} (pre-event, dark)" if i == 0 else None)
    ax[2].set_xticks(x); ax[2].set_xticklabels([c["quoter"] for c in cs], rotation=20, fontsize=7); ax[2].set_ylabel("30-min mark-out against dealer (bp)"); ax[2].set_title("30-min mark-outs: normal (light) vs pre-event (dark)", fontsize=10); ax[2].legend(fontsize=6); ax[2].axhline(0, color="k", lw=0.5)
    save(fig, "frontier.png")

def policy():
    d = load("mm_frontier.json")
    if not d: return
    cells = [c for c in d["cells"] if "quote_structure" in c and c["quote_structure"]]; gs = sorted(set(c["gamma"] for c in cells)); g0 = gs[len(gs) // 2]
    qs = ["Symmetric", "Bergault", "BarzykinCiceri_QAHR", "Bergault_oracle"]; tiers = ["informed", "real_money", "other"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for pnl, (title, pre) in enumerate([("normal RFQs", False), ("30 min before an event", True)]):
        x = np.arange(len(tiers)); w = 0.8 / (2 * len(qs))
        for k, q in enumerate(qs):
            cs = [c for c in cells if c["quoter"] == q and c["gamma"] == g0]
            if not cs: continue
            for r, ra in enumerate([False, True]):
                v = [np.mean([[e["mean_delta_bp"] for e in c["quote_structure"] if e["tier"] == t and e["pre_event"] == pre and e["risk_adding"] == ra][0] for c in cs]) for t in tiers]
                ax[pnl].bar(x + (2 * k + r - len(qs)) * w + w / 2, v, w, color=COL.get(q, "C6"), alpha=1.0 if ra else 0.45, label=f"{q} ({'risk-adding' if ra else 'risk-reducing'} side)" if pnl == 0 else None)
        ax[pnl].set_xticks(x); ax[pnl].set_xticklabels(tiers); ax[pnl].set_ylabel("mean half-spread quoted (bp)"); ax[pnl].set_title(f"{title}, gamma = {g0:.0e}", fontsize=10)
    ax[0].legend(fontsize=6, ncol=2)
    save(fig, "policy.png")

def text_signal():
    d = load("text_signal.json")
    if not d: return
    ev = d["events"]; tr = [e for e in ev if e["date"] <= d["train_end"]]; te = [e for e in ev if e["date"] > d["train_end"]]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for bank, mk in [("FOMC", "o"), ("ECB", "s")]:
        ax[0].scatter([e["dh"] for e in tr if e["bank"] == bank], [e["j10"] for e in tr if e["bank"] == bank], s=14, marker=mk, alpha=0.6, label=f"{bank} train")
        ax[0].scatter([e["dh"] for e in te if e["bank"] == bank], [e["j10"] for e in te if e["bank"] == bank], s=22, marker=mk, edgecolor="k", label=f"{bank} test")
    xs = np.linspace(min(e["dh"] for e in ev), max(e["dh"] for e in ev), 10); ax[0].plot(xs, d["b0"] + d["b1_bp_per_unit_dh"] * xs, "k-", lw=1, label=f"fit (train): b1 = {d['b1_bp_per_unit_dh']:.1f} bp")
    ax[0].set_xlabel(f"change in hawkishness ({d['scorer']})"); ax[0].set_ylabel("10y reaction (bp)"); ax[0].set_title(f"text signal: train corr {d['corr_train']:.2f} (n={d['n_train']}), test corr {d['corr_test'] if d['corr_test'] is None else round(d['corr_test'], 2)} (n={d['n_test']})"); ax[0].legend(fontsize=6); ax[0].axhline(0, color="k", lw=0.4); ax[0].axvline(0, color="k", lw=0.4)
    import datetime as dt; D = lambda s: dt.datetime.strptime(s, "%Y-%m-%d")
    fo = [e for e in ev if e["bank"] == "FOMC"]; ax[1].plot([D(e["date"]) for e in fo], [e["h"] for e in fo], lw=1, label="FOMC"); ec = [e for e in ev if e["bank"] == "ECB"]; ax[1].plot([D(e["date"]) for e in ec], [e["h"] for e in ec], lw=1, label="ECB")
    ax[1].axvline(D(d["placebo"]["checkpoint_date"]), color="C3", ls="--", label="FOMC-RoBERTa checkpoint"); ax[1].axvline(D(d["train_end"]), color="k", ls=":", label="train / test split"); ax[1].tick_params(axis="x", labelsize=7, rotation=30); ax[1].set_ylabel("hawkishness, mean P(hawk) - P(dove)"); ax[1].set_title("statement hawkishness"); ax[1].legend(fontsize=7)
    pl = d["placebo"]; ax[2].bar(["in model training\nperiod (<= Oct 2022)", "after checkpoint\n(> Sep 2023)"], [pl["corr_in_training_period"] or 0, pl["corr_after_checkpoint"] or 0], color=["C3", "C2"]); ax[2].set_ylabel("corr(dh, 10y reaction)"); ax[2].set_title(f"look-ahead placebo (FOMC, n = {pl['n_in']} / {pl['n_after']})"); ax[2].axhline(0, color="k", lw=0.5)
    save(fig, "text_signal.png")

def nlp():
    d = load("nlp_signal.json")
    if not d: return
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
    fs = list(d["feature_sets"]); T = ["r2", "r5", "r10", "r30"]; x = np.arange(len(T)); w = 0.8 / len(fs)
    for k, n in enumerate(fs):
        v = [d["feature_sets"][n][t]["corr"] for t in T]; lo = [d["feature_sets"][n][t]["corr"] - d["feature_sets"][n][t]["corr_ci"][0] for t in T]; hi = [d["feature_sets"][n][t]["corr_ci"][1] - d["feature_sets"][n][t]["corr"] for t in T]
        ax[0].bar(x + (k - len(fs) / 2) * w + w / 2, v, w, yerr=[lo, hi], capsize=1.5, label=n)
        ax[0].plot(x + (k - len(fs) / 2) * w + w / 2, [d["feature_sets"][n][t]["placebo_corr_p95"] for t in T], "k_", ms=6)
    ax[0].axhline(0, color="k", lw=0.5); ax[0].set_xticks(x); ax[0].set_xticklabels(["2y", "5y", "10y", "30y"]); ax[0].set_ylabel("out-of-sample corr (pred, reaction)"); ax[0].set_title(f"text -> event-window move, test n = {d['n_test']} (ticks: placebo 95th pct)", fontsize=9); ax[0].legend(fontsize=6)
    import pandas as pd; f = pd.read_csv(os.path.join(R, "..", "data", "derived", "nlp_features.csv")) if os.path.exists(os.path.join(R, "..", "data", "derived", "nlp_features.csv")) else None
    if f is not None:
        te = f[f["date"] > d["train_end"]]; tr = f[f["date"] <= d["train_end"]]
        ax[1].scatter(tr["alpha_r2"], tr["r2"], s=10, alpha=0.4, label="train (leave-one-out)"); ax[1].scatter(te["alpha_r2"], te["r2"], s=18, edgecolor="k", label="test (out of sample)")
        for bank, mk in [("FOMC", "o"), ("ECB", "s")]: pass
        ax[1].axhline(0, color="k", lw=0.4); ax[1].axvline(0, color="k", lw=0.4); ax[1].set_xlabel("nowcast of the 2y move (bp)"); ax[1].set_ylabel("2y reaction (bp)"); ax[1].set_title(f"2y nowcast ({d['alpha_source']} features)", fontsize=9); ax[1].legend(fontsize=7)
    if "dated_checkpoint" in d:
        c = d["dated_checkpoint"]; ax[2].bar(["redline\ninside", "redline\nafter", "doc dh\ninside", "doc dh\nafter"], [c["corr_redline_r2_inside"], c["corr_redline_r2_after"], c["corr_dh_r2_inside"], c["corr_dh_r2_after"]], color=["C3", "C2", "C3", "C2"])
        ax[2].axhline(0, color="k", lw=0.5); ax[2].set_ylabel("corr with 2y reaction"); ax[2].set_title(f"FOMC-RoBERTa: inside its training period (n={c['n_inside']}) vs after its checkpoint (n={c['n_after']})", fontsize=8)
    sm = d["size_model"]; names = list(sm); ax[3].bar(names, [sm[n]["corr"] for n in names], yerr=[[sm[n]["corr"] - sm[n]["corr_ci"][0] for n in names], [sm[n]["corr_ci"][1] - sm[n]["corr"] for n in names]], capsize=2, color=["C1", "C7", "C2"])
    ax[3].axhline(0, color="k", lw=0.5); ax[3].set_ylabel("OOS corr with |2y reaction|"); ax[3].set_title("size of the move: text (novelty, dispersion) vs realised vol", fontsize=9)
    save(fig, "nlp.png")

def events():
    d = load("events.json")
    if not d: return
    cells = d["cells"]; fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    short = {"text+calendar": "text+cal", "placebo_text+calendar": "placebo", "calendar_only": "cal only", "none": "none"}; labels = [c["quoter"].replace("BarzykinCiceri_QAHR", "BC_QAHR").replace("Bergault_oracle", "oracle") + ("" if c["quoter"] != "CarteaWang" else "\n" + short[c["variant"]]) for c in cells]; x = np.arange(len(cells))
    ax[0].bar(x, [c["pnl_per_day"] / 1e3 for c in cells], yerr=[[(c["pnl_per_day"] - c["pnl_lo"]) / 1e3 for c in cells], [(c["pnl_hi"] - c["pnl_per_day"]) / 1e3 for c in cells]], color=[COL.get(c["quoter"], "C6") for c in cells], alpha=[1.0 if c["variant"] == "signal" else 0.45 for c in cells][0] if False else None, capsize=2)
    for i, c in enumerate(cells):
        if c["variant"] != "text+calendar": ax[0].patches[i].set_alpha(0.45)
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=6, rotation=45); ax[0].set_ylabel("P&L per day ($k)"); ax[0].set_title(f"replay {d['from']}..{d['to']}: {d['n_days']} days, {d['n_event_days']} event days", fontsize=10)
    ax[1].bar(x, [c["pnl_per_var"] * 1e3 for c in cells], color=[COL.get(c["quoter"], "C6") for c in cells])
    for i, c in enumerate(cells):
        if c["variant"] != "text+calendar": ax[1].patches[i].set_alpha(0.45)
    ax[1].set_xticks(x); ax[1].set_xticklabels(labels, fontsize=6, rotation=45); ax[1].set_ylabel("P&L per unit 5-min variance (x1e-3)"); ax[1].set_title("P&L per unit variance (faded: permuted text / calendar only / neither)", fontsize=9)
    base = [c for c in cells if c["variant"] == "text+calendar"]
    for i, c in enumerate(base):
        mk = c["markouts"]; tiers = sorted(set(m["tier"] for m in mk))
        for j, tr in enumerate(tiers):
            pre = [m for m in mk if m["tier"] == tr and m["pre_event"]]; nor = [m for m in mk if m["tier"] == tr and not m["pre_event"]]
            ax[2].bar(i + (j - 1) * 0.25 - 0.06, nor[0]["markout_30m_bp"] if nor else 0, 0.12, color=f"C{j}", alpha=0.5); ax[2].bar(i + (j - 1) * 0.25 + 0.06, pre[0]["markout_30m_bp"] if pre else 0, 0.12, color=f"C{j}", label=f"{tr} (pre-event, dark)" if i == 0 else None)
    ax[2].set_xticks(np.arange(len(base))); ax[2].set_xticklabels([c["quoter"] for c in base], rotation=20, fontsize=7); ax[2].set_ylabel("30-min mark-out against dealer (bp)"); ax[2].set_title("30-min mark-outs: normal (light) vs pre-event (dark)", fontsize=10); ax[2].legend(fontsize=6); ax[2].axhline(0, color="k", lw=0.5)
    save(fig, "events.png")

if __name__ == "__main__":
    for f in (curve, curvehist, hedge, factors, frontier, policy, text_signal, nlp, events):
        try: f()
        except Exception as e: print(f.__name__, "failed:", repr(e))
