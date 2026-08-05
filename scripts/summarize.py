#!/usr/bin/env python3
"""results/summary.md from results/*.json.   python scripts/summarize.py [results]"""
import glob, json, os, sys
import numpy as np

R = sys.argv[1] if len(sys.argv) > 1 else "results"
def load(n): p = os.path.join(R, n); return json.load(open(p)) if os.path.exists(p) else None
out = []
def P(s=""): out.append(s)
def table(cols, rows): P("| " + " | ".join(cols) + " |"); P("|" + "---|" * len(cols)); [P("| " + " | ".join(str(c) for c in r) + " |") for r in rows]; P()

P("# Results summary\n")
fs = sorted(glob.glob(os.path.join(R, "curve_*.json")))
if fs:
    d = json.load(open(fs[-1])); P(f"## Curve construction ({d['date']})\n")
    P(f"Bootstrap (Hagan-West monotone convex, global iteration): max repricing error {d['bootstrap_max_reprice_error_bp']:.2e} bp in {d['bootstrap_iterations']} sweeps. Kernel ridge (lambda = {d['kernel_ridge_lambda']:.0e}): max repricing error {d['kernel_ridge_max_reprice_error_bp']:.3f} bp. Forward roughness int f'(t)^2 dt: bootstrap {d['roughness_boot']:.2e}, kernel ridge {d['roughness_kr']:.2e}.\n")
    loo = d["leave_one_out"]; inter = [x for x in loo if 0.2 < x["T"] <= 10]
    P(f"Leave-one-out pricing error, interior tenors 6m-10y (mean |err|): bootstrap {np.mean([abs(x['loo_err_bp_boot']) for x in inter]):.2f} bp, kernel ridge {np.mean([abs(x['loo_err_bp_kr']) for x in inter]):.2f} bp. Endpoints are extrapolation and excluded (20y: {[x for x in loo if x['T']==20][0]['loo_err_bp_boot']:.0f} / {[x for x in loo if x['T']==20][0]['loo_err_bp_kr']:.0f} bp; 30y: {[x for x in loo if x['T']==30][0]['loo_err_bp_boot']:.0f} / {[x for x in loo if x['T']==30][0]['loo_err_bp_kr']:.0f} bp).\n")
    table(["tenor", "par %", "LOO boot (bp)", "LOO KR (bp)"], [[f"{x['T']:g}y", f"{100 * i['par']:.3f}", f"{x['loo_err_bp_boot']:+.2f}", f"{x['loo_err_bp_kr']:+.2f}"] for x, i in zip(loo, d["instruments"])])
    P("Ladder key-rate DV01s (per 1 notional per bp; row = bond, columns = bumped input):\n")
    table(["bond", "coupon %", "DV01"] + [f"{i['T']:g}y" for i in d["instruments"]], [[f"{k['T']:g}y", f"{100 * k['coupon']:.3f}", f"{k['dv01_per_1']:.2e}"] + [f"{v:.1e}" for v in k["krdv01_per_1"]] for k in d["ladder_krdv01"]])
    fe = d["front_end"]; P(f"Meeting-date step front end ({fe['source']}): levels " + ", ".join(f"{100 * l:.2f}%" for l in fe["levels"]) + " on " + ", ".join(f"{m:.2f}y" for m in fe["meetings"]) + ".\n")
d = load("curvehist.json")
if d:
    P(f"## Daily curve builds, {d['days'][0]['date']} .. {d['days'][-1]['date']} ({d['n']} days)\n")
    ms = d["methods"]; boot = ms[0]; kr = [m for m in ms if m["method"] == "kernel_ridge"]; ok = [m for m in kr if m["max_reprice_bp"] <= 1.0]; sel = max(ok, key=lambda m: m["lambda"]) if ok else kr[-1]
    P("Selection rule fixed in advance: the smoothest kernel ridge (largest lambda) whose max repricing error over all days is within 1 bp; that is lambda = %.0e. Leave-one-out = each interior tenor (6m-10y) dropped and repriced off the curve built from the rest; hedge test = a 3y (7y) par bond priced off the curve built without that point and hedged with the remaining inputs using that curve's key-rate DV01s, realised next-day hedged P&L in yield bp.\n" % sel["lambda"])
    table(["method", "lambda", "max reprice bp", "LOO mean |err| yield bp", "LOO median", "10y DV01 day-to-day std", "forward roughness", "hedged 3y std bp", "hedged 7y std bp"],
          [[m["method"], f"{m['lambda']:.0e}" if m["lambda"] else "-", f"{m['max_reprice_bp']:.2g}", f"{m['loo_yield_bp_mean']:.3f}", f"{m['loo_yield_bp_median']:.3f}", f"{m['dv01_10y_daily_change_std']:.3f}", f"{m['roughness_mean']:.2e}", f"{m['hedge_test']['3y']['hedged_pnl_bp_std']:.3f}", f"{m['hedge_test']['7y']['hedged_pnl_bp_std']:.3f}"] for m in ms])
    P(f"Selected kernel ridge vs bootstrap: LOO error {100 * (sel['loo_yield_bp_mean'] / boot['loo_yield_bp_mean'] - 1):+.1f}% (median {100 * (sel['loo_yield_bp_median'] / boot['loo_yield_bp_median'] - 1):+.1f}%), roughness {100 * (sel['roughness_mean'] / boot['roughness_mean'] - 1):+.1f}%, hedged 3y / 7y std {100 * (sel['hedge_test']['3y']['hedged_pnl_bp_std'] / boot['hedge_test']['3y']['hedged_pnl_bp_std'] - 1):+.1f}% / {100 * (sel['hedge_test']['7y']['hedged_pnl_bp_std'] / boot['hedge_test']['7y']['hedged_pnl_bp_std'] - 1):+.1f}%, 10y DV01 day-to-day std {100 * (sel['dv01_10y_daily_change_std'] / boot['dv01_10y_daily_change_std'] - 1):+.1f}%. Unhedged 7y daily std {boot['hedge_test']['7y']['unhedged_bp_std']:.2f} bp.\n")
    yrs = sorted(boot["loo_yield_bp_by_year"]); table(["year"] + yrs, [[("bootstrap" if m["method"] == "bootstrap" else f"KR {m['lambda']:.0e}")] + [f"{m['loo_yield_bp_by_year'][y]:.2f}" for y in yrs] for m in ms])
d = load("hedge_test.json")
if d:
    P(f"## Learned vs model hedge ratios (walk-forward from {d['test_from']}, trailing window {d['window']} days)\n")
    labels = {"unhedged": "unhedged", "naive": "50/50 neighbours", "model_boot": "bootstrap key-rate DV01s", "model_kr": "kernel-ridge key-rate DV01s", "pca": "PCA level+slope neutral (trailing window)", "learned": "ridge regression on neighbours (trailing window)"}
    rows = []
    for T, v in d["targets"].items():
        for k in ["unhedged", "naive", "model_boot", "model_kr", "pca", "learned"]:
            if k in v: rows.append([f"{T[:-2]}y", labels[k], v[k]["n"], f"{v[k]['std_bp']:.3f}", f"{100 * v[k]['reduction_vs_unhedged']:.1f}%"])
    table(["bond", "hedge ratios", "days", "hedged P&L std (yield bp)", "reduction vs unhedged"], rows)
    for T, v in d["targets"].items():
        if "learned_over_model_boot_std_ratio" in v: r = v["learned_over_model_boot_std_ratio"]; P(f"{T[:-2]}y: learned / bootstrap-key-rate std ratio {r['point']:.3f} (bootstrap 95% CI [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]).")
    P()
d = load("factors.json")
if d:
    P("## Factor model (PCA of daily par-yield changes, 12 tenors)\n")
    table(["window", "days", "explained 1/2/3", "factor sd bp/day (level, slope, curvature)"], [[w, v["n_days"], " / ".join(f"{e:.3f}" for e in v["explained"]), ", ".join(f"{np.sqrt(e):.1f}" for e in v["evals_bp2_per_day"][:3])] for w, v in d.items()])
d = load("mm_frontier.json")
if d:
    P(f"## Closed-form quoters on synthetic days ({d['days']} days x {d['seeds']} seeds, common random numbers)\n")
    P("P&L per day with 95% bootstrap CI over seeds; 5-minute P&L std; decomposition (spread + inventory - hedge cost = total, residual is machine zero); hit ratios by tier (informed / real money / other).\n")
    table(["quoter", "gamma", "P&L/day $k", "CI", "5-min std $k", "P&L/var x1e-3", "spread", "inventory", "hedge", "max |resid|", "hit ratios"],
          [[c["quoter"], f"{c['gamma']:.0e}", f"{c['pnl_per_day'] / 1e3:.1f}", f"[{c['pnl_lo'] / 1e3:.1f}, {c['pnl_hi'] / 1e3:.1f}]", f"{c['pnl_5min_std'] / 1e3:.1f}", f"{1e3 * c['pnl_per_var']:.2f}", f"{c['decomposition']['spread'] / 1e3:.1f}", f"{c['decomposition']['inventory'] / 1e3:.1f}", f"{-c['decomposition']['hedge_cost'] / 1e3:.1f}", f"{c['decomposition']['max_abs_residual']:.1e}", "/".join(f"{h:.2f}" for h in c["hit_ratio_by_tier"])] for c in d["cells"]])
    gs = sorted(set(c["gamma"] for c in d["cells"])); g0 = gs[len(gs) // 2]; P(f"Mark-outs at 30 minutes (bp against the dealer), gamma = {g0:.0e}:\n")
    rows = []
    for c in [c for c in d["cells"] if abs(c["gamma"] - g0) < 1e-12]:
        for m in c["markouts"]: rows.append([c["quoter"], m["tier"], "pre-event" if m["pre_event"] else "normal", m["n"], f"{m['mean_delta_bp']:.2f}", f"{m['markout_5m_bp']:+.2f}", f"{m['markout_30m_bp']:+.2f}", f"{m['markout_60m_bp']:+.2f}"])
    table(["quoter", "tier", "window", "fills", "half-spread bp", "5m", "30m", "60m"], rows)
d = load("text_signal.json")
if d:
    P(f"## Text signal ({d['scorer']}), training window to {d['train_end']}\n")
    P(f"Regression of the 10y reaction on the change in statement hawkishness, pooled FOMC + ECB: b1 = {d['b1_bp_per_unit_dh']:.1f} bp per unit (bootstrap 95% CI [{d['b1_ci'][0]:.1f}, {d['b1_ci'][1]:.1f}]), n = {d['n_train']}, in-sample corr {d['corr_train']:.3f}, R2 {d['r2_train']:.3f}. Out of sample (n = {d['n_test']}): corr {d['corr_test']}, R2 {d['r2_test_oos']}, sign hit rate {d['sign_hit_test']}.\n")
    P("By bank (training window): " + "; ".join(f"{k}: n {v['n']}, corr(dh) {v['corr_dh']:.3f}, corr(level) {v['corr_level']:.3f}" for k, v in d["by_bank_train"].items()) + ".\n")
    pl = d["placebo"]; P(f"Look-ahead placebo (FOMC-RoBERTa checkpoint {pl['checkpoint_date']}): corr(dh, reaction) inside the model's training period {pl['corr_in_training_period']} (n = {pl['n_in']}) vs after the checkpoint {pl['corr_after_checkpoint']} (n = {pl['n_after']}).\n")
d = load("nlp_signal.json")
if d:
    P(f"## NLP: statement text -> event-window move (train to {d['train_end']}, n = {d['n_train']}; test n = {d['n_test']})\n")
    P("Ridge on hand-scored features (lexicon and FOMC-RoBERTa, whole document and redline = sentences added minus removed vs the previous statement, novelty, dispersion) or on the 384-d MiniLM embedding of the redline; lambda by leave-one-out on the training set; placebo = permutation test of the out-of-sample correlation (95th percentile of |corr| under permutation, and the two-sided p-value).\n")
    rows = []
    for n, r in d["feature_sets"].items():
        for t in ["r2", "r5", "r10", "r30"]: x = r[t]; rows.append([n, t[1:] + "y", x["lambda"], f"{x['corr_loo_train']:.2f}", f"{x['corr']:.2f} [{x['corr_ci'][0]:.2f}, {x['corr_ci'][1]:.2f}]", f"{x['r2_oos']:.3f}" if x["r2_oos"] is not None else "", f"{x['sign_hit']:.2f} [{x['sign_hit_ci'][0]:.2f}, {x['sign_hit_ci'][1]:.2f}]", f"{x['placebo_corr_p95']:.2f}", f"{x['placebo_p_value']:.3f}"])
    table(["features", "target", "lambda", "LOO corr (train)", "OOS corr [CI]", "OOS R2", "sign hit [CI]", "placebo corr p95", "placebo p"], rows)
    if "dated_checkpoint" in d: c = d["dated_checkpoint"]; P(f"Dated-checkpoint arm (FOMC-RoBERTa checkpoint {c['checkpoint']}, training data to {c['training_data_end']}): corr(redline score, 2y reaction) inside the training period {c['corr_redline_r2_inside']:.2f} (n = {c['n_inside']}) vs after the checkpoint {c['corr_redline_r2_after']:.2f} (n = {c['n_after']}); document-level dh {c['corr_dh_r2_inside']:.2f} vs {c['corr_dh_r2_after']:.2f}.\n")
    P("Size of the move (|2y reaction|), OOS corr: " + "; ".join(f"{k}: {v['corr']:.2f} [{v['corr_ci'][0]:.2f}, {v['corr_ci'][1]:.2f}]" for k, v in d["size_model"].items()) + f". Replay alphas come from the {d['alpha_source']} feature set.\n")
for n, title in [("events.json", "Real-day replay with event features"), ("events_ex2020.json", "Real-day replay, factor model fitted without 2020-22")]:
    d = load(n)
    if d:
        P(f"## {title}: {d['from']}..{d['to']}, {d['n_days']} days, {d['n_event_days']} event days ({d['n_measured_jumps']} with measured USMPD jumps), gamma {d['gamma']:.0e}\n")
        P("Events by type: " + ", ".join(f"{k} {int(v)}" for k, v in d["events_by_type"].items()) + ".\n")
        table(["quoter", "variant", "P&L/day $k", "CI", "5-min std $k", "P&L/var x1e-3", "spread", "inventory", "hedge", "hit ratios"],
              [[c["quoter"], c["variant"], f"{c['pnl_per_day'] / 1e3:.1f}", f"[{c['pnl_lo'] / 1e3:.1f}, {c['pnl_hi'] / 1e3:.1f}]", f"{c['pnl_5min_std'] / 1e3:.1f}", f"{1e3 * c['pnl_per_var']:.2f}", f"{c['decomposition']['spread'] / 1e3:.1f}", f"{c['decomposition']['inventory'] / 1e3:.1f}", f"{-c['decomposition']['hedge_cost'] / 1e3:.1f}", "/".join(f"{h:.2f}" for h in c["hit_ratio_by_tier"])] for c in d["cells"]])
        rows = []
        for c in [c for c in d["cells"] if c["variant"] == "text+calendar"]:
            for m in c["markouts"]: rows.append([c["quoter"], m["tier"], "pre-event" if m["pre_event"] else "normal", m["n"], f"{m['mean_delta_bp']:.2f}", f"{m['markout_30m_bp']:+.2f}"])
        table(["quoter", "tier", "window", "fills", "half-spread bp", "30m mark-out bp"], rows)
t = os.path.join(R, "tests.txt")
if os.path.exists(t): P("## Tests\n"); P(open(t).read().strip().splitlines()[-1] + "\n")
open(os.path.join(R, "summary.md"), "w").write("\n".join(out)); print("wrote", os.path.join(R, "summary.md"))
