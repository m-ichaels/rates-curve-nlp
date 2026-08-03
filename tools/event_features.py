"""Event calendar and features for the real-day replay.
Inputs: data/derived/text_scores.csv (tools/text_features.py), data/raw/USMPD.xlsx (SF Fed event-study reactions),
data/raw/treasury_cmt.csv, data/raw/treasury_auctions.csv, data/raw/cftc/*.txt, data/raw/nyfed_rates.csv.
Outputs:
  data/derived/events.csv   date,is_event,type,alpha_bp,info_mult,jump2y,jump5y,jump10y,jump30y
  results/text_signal.json  regression of the FOMC/ECB reaction on the change in hawkishness (training window), OOS test,
                            dated-checkpoint arm and the look-ahead placebo
The signal alpha (expected 10y yield change at the event, bp) is b * dh with b fitted on the training window only, dh the
change in the statement's hawkishness versus the previous statement, from text available at the release time.
"""
import os, sys, json, glob, datetime as dt
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRAIN_END = "2023-12-31"
args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
SCORER = args.get("scorer", "fomc_rob")          # dict | fomc_rob | lorenzo
TRAIN_END = args.get("train_end", TRAIN_END)

def nfp_dates(y0, y1):
    """BLS rule: released the third Friday after the reference week (the week containing the 12th)."""
    out = []
    for y in range(y0, y1 + 1):
        for m in range(1, 13):
            d12 = dt.date(y, m, 12); sat = d12 + dt.timedelta(days=(5 - d12.weekday()) % 7); out.append(sat + dt.timedelta(days=20))
    return set(out)

def zscore(x, w=52):
    return (x - x.rolling(w, min_periods=20).mean()) / x.rolling(w, min_periods=20).std()

def main():
    cmt = pd.read_csv(os.path.join(ROOT, "data/raw/treasury_cmt.csv")); cmt["date"] = pd.to_datetime(cmt["date"]); cmt = cmt.sort_values("date").set_index("date")
    days = cmt.index
    ts = pd.read_csv(os.path.join(ROOT, "data/derived/text_scores.csv")); ts["date"] = pd.to_datetime(ts["date"])
    hcol = {"dict": "dict_h", "fomc_rob": "fomc_rob_h", "lorenzo": "lorenzo_h"}[SCORER]; dcol = hcol.replace("_h", "_disp")
    if hcol not in ts.columns: hcol, dcol = "dict_h", "dict_disp"; print("scorer not available, using dict")
    ts = ts.sort_values(["bank", "date"]); ts["dh"] = ts.groupby("bank")[hcol].diff(); ts["disp"] = ts[dcol]
    # reactions: USMPD statements sheet (pp -> bp) for FOMC; ECB uses the daily CMT change (no intraday source)
    us = pd.read_excel(os.path.join(ROOT, "data/raw/USMPD.xlsx"), sheet_name="Statements"); us["Date"] = pd.to_datetime(us["Date"])
    us = us[["Date", "UST2Y", "UST5Y", "UST10Y", "UST30Y", "MP1", "Unscheduled"]].set_index("Date") * 1.0
    for c in ["UST2Y", "UST5Y", "UST10Y", "UST30Y", "MP1"]: us[c] = us[c] * 100.0
    dy = cmt[["2Yr", "5Yr", "10Yr", "30Yr"]].diff() * 100.0
    ev = []
    for _, r in ts.iterrows():
        d = r["date"]
        if d not in days:   # statement on a non-business day: attach to the next trading day
            nxt = days[days >= d]; d = nxt[0] if len(nxt) else None
        if d is None or pd.isna(r["dh"]): continue
        if r["bank"] == "FOMC" and r["date"] in us.index: j = us.loc[r["date"], ["UST2Y", "UST5Y", "UST10Y", "UST30Y"]].values.astype(float); src = "usmpd"
        else: j = dy.loc[d].values.astype(float) * 0.6 if d in dy.index else np.full(4, np.nan); src = "daily x 0.6"
        ev.append({"date": d, "bank": r["bank"], "dh": r["dh"], "disp": r["disp"], "h": r[hcol], "j2": j[0], "j5": j[1], "j10": j[2], "j30": j[3], "src": src})
    ev = pd.DataFrame(ev).dropna(subset=["j10"])
    # signal regression on the training window, pooled Fed + ECB, with bootstrap CIs (n is small: ~8 meetings per bank per year)
    tr = ev[ev["date"] <= TRAIN_END]; te = ev[ev["date"] > TRAIN_END]
    def fit(x, y):
        X = np.c_[np.ones(len(x)), x]; b = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ b; r2 = 1 - res.var() / y.var() if len(y) > 2 else 0; return b, r2
    b, r2 = fit(tr["dh"].values, tr["j10"].values); rng = np.random.default_rng(0)
    bs = [fit(tr["dh"].values[i], tr["j10"].values[i])[0][1] for i in (rng.integers(0, len(tr), len(tr)) for _ in range(2000))]
    corr_tr = np.corrcoef(tr["dh"], tr["j10"])[0, 1]; corr_te = np.corrcoef(te["dh"], te["j10"])[0, 1] if len(te) > 3 else np.nan
    pred_te = b[0] + b[1] * te["dh"].values; oos_r2 = 1 - ((te["j10"].values - pred_te) ** 2).mean() / te["j10"].var() if len(te) > 3 else np.nan
    hit_te = np.mean(np.sign(pred_te) == np.sign(te["j10"].values)) if len(te) else np.nan
    # dated-checkpoint arm and Gao-Jiang-Yan placebo: is the text-reaction correlation stronger inside the model's training
    # period (contaminated) than after its checkpoint date?
    ck = "2023-09-26"; fo = ev[ev["bank"] == "FOMC"]
    pre = fo[fo["date"] <= "2022-10-31"]; post = fo[fo["date"] > ck]
    placebo = {"scorer": SCORER, "checkpoint_date": ck, "corr_in_training_period": float(np.corrcoef(pre["dh"], pre["j10"])[0, 1]) if len(pre) > 3 else None, "n_in": int(len(pre)),
               "corr_after_checkpoint": float(np.corrcoef(post["dh"], post["j10"])[0, 1]) if len(post) > 3 else None, "n_after": int(len(post))}
    # per-bank in-sample correlations, and level (not change) as a control
    by_bank = {bk: {"n": int((tr["bank"] == bk).sum()), "corr_dh": float(np.corrcoef(tr[tr["bank"] == bk]["dh"], tr[tr["bank"] == bk]["j10"])[0, 1]), "corr_level": float(np.corrcoef(tr[tr["bank"] == bk]["h"], tr[tr["bank"] == bk]["j10"])[0, 1])} for bk in ["FOMC", "ECB"] if (tr["bank"] == bk).sum() > 3}
    sig = {"scorer": SCORER, "train_end": TRAIN_END, "n_train": int(len(tr)), "n_test": int(len(te)), "b0": float(b[0]), "b1_bp_per_unit_dh": float(b[1]), "b1_ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
           "r2_train": float(r2), "corr_train": float(corr_tr), "corr_test": float(corr_te) if corr_te == corr_te else None, "r2_test_oos": float(oos_r2) if oos_r2 == oos_r2 else None, "sign_hit_test": float(hit_te) if hit_te == hit_te else None,
           "by_bank_train": by_bank, "placebo": placebo, "events": ev.assign(date=ev["date"].dt.strftime("%Y-%m-%d")).to_dict("records")}
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True); json.dump(sig, open(os.path.join(ROOT, "results/text_signal.json"), "w"), indent=1, default=float)
    print(f"signal b1 = {b[1]:.2f} bp per unit dh (CI {sig['b1_ci'][0]:.2f}..{sig['b1_ci'][1]:.2f}), train corr {corr_tr:.3f} n={len(tr)}, test corr {corr_te:.3f} n={len(te)}, OOS R2 {oos_r2:.3f}, sign hit {hit_te:.2f}")
    print("placebo", placebo)
    # other features
    au = pd.read_csv(os.path.join(ROOT, "data/raw/treasury_auctions.csv")); au["auction_date"] = pd.to_datetime(au["auction_date"])
    big = au[(au["security_type"].isin(["Note", "Bond"])) & (au["security_term"].str.contains("^(?:9-Year|10-Year|29-Year|30-Year|19-Year|20-Year)", regex=True, na=False))]
    auction_days = set(big["auction_date"].dt.normalize())
    nfp = nfp_dates(2015, 2026)
    cf = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(ROOT, "data/raw/cftc/fut_fin_*.txt"))]); cf["Report_Date_as_YYYY-MM-DD"] = pd.to_datetime(cf["Report_Date_as_YYYY-MM-DD"])
    t10 = cf[cf["Market_and_Exchange_Names"].str.startswith("10-YEAR U.S. TREASURY NOTES")].sort_values("Report_Date_as_YYYY-MM-DD").set_index("Report_Date_as_YYYY-MM-DD")
    lev = (t10["Lev_Money_Positions_Long_All"] - t10["Lev_Money_Positions_Short_All"]) / t10["Open_Interest_All"]; lev_z = zscore(lev)
    lev_daily = lev_z.reindex(days, method="ffill").shift(3)   # report published Friday for Tuesday positions: use with a 3-business-day lag
    ny = pd.read_csv(os.path.join(ROOT, "data/raw/nyfed_rates.csv")); ny["Effective Date"] = pd.to_datetime(ny["Effective Date"])
    sofr = ny[ny["Rate Type"] == "SOFR"].set_index("Effective Date").sort_index(); effr = ny[ny["Rate Type"] == "EFFR"].set_index("Effective Date").sort_index()
    stress = ((sofr["Rate (%)"] - effr["Rate (%)"].reindex(sofr.index)) * 100).reindex(days, method="ffill"); tail = ((sofr["99th Percentile (%)"] - sofr["Rate (%)"]) * 100).reindex(days, method="ffill")
    stress_z = zscore(stress.fillna(0), 250); tail_z = zscore(tail.fillna(0), 250)
    # per-tenor nowcasts from tools/nlp_features.py when available (redline / RoBERTa ridge models, training-window fit, LOO inside it)
    nlp = None; nlp_path = os.path.join(ROOT, "data/derived/nlp_features.csv")
    if os.path.exists(nlp_path): nlp = pd.read_csv(nlp_path); nlp["date"] = pd.to_datetime(nlp["date"]); nlp = nlp.set_index(["bank", "date"]); print("using NLP alphas from", nlp_path)
    # assemble the daily calendar
    evd = ev.set_index("date")
    rows = []
    for d in days:
        typ = ""; alpha = 0.0; info = 1.0; jumps = ["", "", "", ""]; alphas = ["", "", "", ""]
        if d in evd.index:
            r = evd.loc[d]; r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
            typ = r["bank"]; alpha = float(b[0] + b[1] * r["dh"]); info = 3.0 * (1.0 + 0.5 * min(2.0, r["disp"] / max(1e-9, ts["disp"].median())))   # dispersion of the statement raises informed intensity
            if nlp is not None:
                key = None
                for cand in nlp.index:
                    if cand[0] == typ and abs((cand[1] - d).days) <= 3: key = cand; break
                if key is not None: alphas = [f"{nlp.loc[key, 'alpha_' + t]:.3f}" for t in ("r2", "r5", "r10", "r30")]; alpha = float(nlp.loc[key, "alpha_r10"])
            if r["src"] == "usmpd": jumps = [f"{r['j2']:.2f}", f"{r['j5']:.2f}", f"{r['j10']:.2f}", f"{r['j30']:.2f}"]
        elif d.date() in nfp: typ = "NFP"; info = 2.5
        elif d.normalize() in auction_days: typ = "AUCTION"; info = 1.5
        z = lev_daily.get(d, np.nan); s = stress_z.get(d, np.nan); t = tail_z.get(d, np.nan)
        if typ: info *= 1.0 + 0.25 * min(2.0, abs(z) if z == z else 0.0) + 0.25 * max(0.0, s if s == s else 0.0) + 0.25 * max(0.0, t if t == t else 0.0)   # positioning and funding stress scale the informed intensity
        rows.append({"date": d.strftime("%Y-%m-%d"), "is_event": int(bool(typ)), "type": typ, "alpha_bp": round(alpha, 3), "info_mult": round(info, 3), "jump2y": jumps[0], "jump5y": jumps[1], "jump10y": jumps[2], "jump30y": jumps[3], "alpha2y": alphas[0], "alpha5y": alphas[1], "alpha10y": alphas[2], "alpha30y": alphas[3],
                     "cftc_lev_z": round(z, 3) if z == z else "", "sofr_effr_z": round(s, 3) if s == s else "", "sofr_p99_z": round(t, 3) if t == t else ""})
    out = pd.DataFrame(rows); os.makedirs(os.path.join(ROOT, "data/derived"), exist_ok=True); out.to_csv(os.path.join(ROOT, "data/derived/events.csv"), index=False)
    print("events.csv:", len(out), "days;", out["type"].value_counts().to_dict())

if __name__ == "__main__": main()
