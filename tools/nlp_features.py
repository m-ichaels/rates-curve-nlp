"""NLP arm: what does the statement text predict about the event-window move, out of sample?
Features per statement (FOMC, ECB), all computable at release time from the statement and its predecessor:
  redline      : sentences added / removed versus the previous statement (difflib); hawkishness of added minus removed
                 sentences under FOMC-RoBERTa (tim9510019 re-upload, commit 2023-09-26) and the lexicon
  document     : hawkishness change of the whole document, level, sentence dispersion
  novelty      : fraction of sentences changed, number added, length
  embedding    : mean MiniLM embedding of the added sentences minus the removed ones (384-d), ridge-regressed
Targets: USMPD event-window reactions (FOMC: 2y/5y/10y/30y, MP1) and daily changes for ECB days (US curve).
Models trained on statements to TRAIN_END (pooled banks), tested after; bootstrap CIs; permutation placebo; the dated-checkpoint
split for FOMC-RoBERTa (inside its training period vs after its checkpoint).
Writes results/nlp_signal.json and data/derived/nlp_features.csv (per-event alpha vector for the replay).
"""
import os, sys, glob, json, difflib, hashlib
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("HF_HOME", os.path.join(ROOT, "storage", "hf"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from text_features import split_sentences, fomc_text, ecb_text, dict_score, HF
args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
TRAIN_END = args.get("train_end", "2023-12-31"); USE_MODELS = args.get("models", "1") == "1"
TARGETS = ["r2", "r5", "r10", "r30"]

def load_docs():
    docs = []
    for bank, reader, pat in [("FOMC", fomc_text, "data/raw/fomc/*.txt"), ("ECB", ecb_text, "data/raw/ecb/*.txt")]:
        for f in sorted(glob.glob(os.path.join(ROOT, pat))): docs.append((bank, pd.Timestamp(os.path.basename(f)[:8]), split_sentences(reader(f))))
    return docs

def reactions(bank, date, us, dy):
    if bank == "FOMC" and date in us.index: x = us.loc[date]; return {"r2": x["UST2Y"] * 100, "r5": x["UST5Y"] * 100, "r10": x["UST10Y"] * 100, "r30": x["UST30Y"] * 100, "mp1": x["MP1"] * 100, "src": "usmpd"}
    nxt = dy.index[dy.index >= date]
    if len(nxt): x = dy.loc[nxt[0]]; return {"r2": x["2Yr"], "r5": x["5Yr"], "r10": x["10Yr"], "r30": x["30Yr"], "mp1": np.nan, "src": "daily"}
    return None

def ridge_fit(X, y, lam):
    mu, sd = X.mean(0), X.std(0) + 1e-12; Z = (X - mu) / sd; ym = y.mean()
    w = np.linalg.solve(Z.T @ Z + lam * np.eye(Z.shape[1]), Z.T @ (y - ym)); return (mu, sd, ym, w)
def ridge_pred(m, X): mu, sd, ym, w = m; return ym + ((X - mu) / sd) @ w
def score(y, p):
    y, p = np.asarray(y), np.asarray(p); c = np.corrcoef(y, p)[0, 1] if len(y) > 2 else np.nan
    return {"n": int(len(y)), "corr": float(c), "r2_oos": float(1 - ((y - p) ** 2).mean() / y.var()) if len(y) > 2 else None, "sign_hit": float(np.mean(np.sign(p) == np.sign(y))), "mae_bp": float(np.abs(y - p).mean())}
def boot_ci(fn, n, B=2000, seed=0):
    rng = np.random.default_rng(seed); v = [fn(rng.integers(0, n, n)) for _ in range(B)]; return [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))]

def main():
    us = pd.read_excel(os.path.join(ROOT, "data/raw/USMPD.xlsx"), sheet_name="Statements"); us["Date"] = pd.to_datetime(us["Date"]); us = us.set_index("Date")
    cmt = pd.read_csv(os.path.join(ROOT, "data/raw/treasury_cmt.csv")); cmt["date"] = pd.to_datetime(cmt["date"]); cmt = cmt.set_index("date").sort_index(); dy = cmt[["2Yr", "5Yr", "10Yr", "30Yr"]].diff() * 100
    docs = load_docs(); models = {}; emb = None
    if USE_MODELS:
        try:
            models["rob"] = HF("tim9510019/FOMC-RoBERTa")
            from sentence_transformers import SentenceTransformer; from huggingface_hub import model_info
            ename = "sentence-transformers/all-MiniLM-L6-v2"; emb = SentenceTransformer(ename); emb_sha = model_info(ename).sha
        except Exception as e: print("models unavailable:", e); models = {}; emb = None
    cache = os.path.join(ROOT, "data/derived/nlp_cache.json"); rows = []; prev = {}
    if os.path.exists(cache) and args.get("recompute", "0") != "1":
        rows = json.load(open(cache)); print("using cached document features", cache); docs = []
        for r in rows: r["date"] = pd.Timestamp(r["date"])
    for bank, date, sents in docs:
        p = prev.get(bank); prev[bank] = sents
        if p is None: continue
        ps = set(p); ss = set(sents); added = [x for x in sents if x not in ps]; removed = [x for x in p if x not in ss]
        sm = difflib.SequenceMatcher(a=p, b=sents, autojunk=False); same = sum(b.size for b in sm.get_matching_blocks())
        r = reactions(bank, date, us, dy)
        if r is None: continue
        row = {"bank": bank, "date": date, "n_sent": len(sents), "n_added": len(added), "n_removed": len(removed), "novelty": 1 - same / max(len(sents), 1), **r}
        dall = np.array([dict_score(x) for x in sents]); dprev = np.array([dict_score(x) for x in p])
        row["dict_h"] = dall.mean(); row["dict_dh"] = dall.mean() - dprev.mean(); row["dict_disp"] = dall.std()
        row["dict_redline"] = (np.mean([dict_score(x) for x in added]) if added else 0.0) - (np.mean([dict_score(x) for x in removed]) if removed else 0.0)
        if models:
            ra = np.array(models["rob"].scores(sents)); rp = np.array(models["rob"].scores(p))
            row["rob_h"] = ra.mean(); row["rob_dh"] = ra.mean() - rp.mean(); row["rob_disp"] = ra.std()
            row["rob_redline"] = (np.mean(models["rob"].scores(added)) if added else 0.0) - (np.mean(models["rob"].scores(removed)) if removed else 0.0)
        if emb is not None:
            ea = emb.encode(added).mean(0) if added else np.zeros(384); er = emb.encode(removed).mean(0) if removed else np.zeros(384)
            row["emb"] = (ea - er).tolist()
        rows.append(row); print(bank, date.date(), f"added {len(added)} removed {len(removed)} novelty {row['novelty']:.2f}", flush=True)
    if docs: json.dump([{**r, "date": r["date"].strftime("%Y-%m-%d")} for r in rows], open(cache, "w"))
    df = pd.DataFrame(rows).sort_values("date").reset_index(drop=True)
    tr = df[df["date"] <= TRAIN_END]; te = df[df["date"] > TRAIN_END]
    out = {"train_end": TRAIN_END, "n_train": int(len(tr)), "n_test": int(len(te)), "models": {}, "feature_sets": {}}
    if models: out["models"]["fomc_roberta"] = {"name": models["rob"].name, "commit": models["rob"].sha}
    if emb is not None: out["models"]["embedding"] = {"name": "sentence-transformers/all-MiniLM-L6-v2", "commit": emb_sha}
    # feature sets: hand-scored (small, interpretable) and embedding (384-d, ridge)
    hand = ["dict_dh", "dict_redline", "dict_h", "novelty", "dict_disp"] + (["rob_dh", "rob_redline", "rob_h", "rob_disp"] if models else [])
    sets = {"lexicon": ["dict_dh", "dict_redline", "dict_h", "novelty", "dict_disp"], "roberta": ["rob_dh", "rob_redline", "rob_h", "rob_disp", "novelty"] if models else None, "all_hand": hand}
    if emb is not None: sets["embedding"] = "emb"
    if models: sets["roberta_redline_only"] = ["rob_redline"]; sets["roberta_doc_only"] = ["rob_dh"]
    def X_of(d, fs): return np.array(d["emb"].tolist()) if fs == "emb" else d[fs].values.astype(float)
    preds = {}
    for name, fs in sets.items():
        if fs is None: continue
        res = {}
        for t in TARGETS:
            Xtr, ytr, Xte, yte = X_of(tr, fs), tr[t].values, X_of(te, fs), te[t].values
            # lambda by leave-one-out on the training set
            best = None
            for lam in ([1, 3, 10, 30, 100, 300, 1000] if fs == "emb" else [0.1, 1, 3, 10, 30, 100]):
                loo = np.array([ridge_pred(ridge_fit(np.delete(Xtr, i, 0), np.delete(ytr, i), lam), Xtr[i:i + 1])[0] for i in range(len(ytr))])
                c = np.corrcoef(loo, ytr)[0, 1]
                if best is None or c > best[0]: best = (c, lam, loo)
            m = ridge_fit(Xtr, ytr, best[1]); p = ridge_pred(m, Xte)
            s = score(yte, p); s["lambda"] = best[1]; s["corr_loo_train"] = float(best[0])
            s["corr_ci"] = boot_ci(lambda ix: np.corrcoef(yte[ix], p[ix])[0, 1], len(yte)); s["sign_hit_ci"] = boot_ci(lambda ix: np.mean(np.sign(p[ix]) == np.sign(yte[ix])), len(yte))
            # permutation placebo: refit with shuffled training targets, score on the test set (distribution of corr under no signal)
            rng = np.random.default_rng(1); perm = np.array([np.corrcoef(yte, rng.permutation(p))[0, 1] for _ in range(5000)])   # permutation test of the OOS correlation
            s["placebo_corr_p95"] = float(np.nanpercentile(np.abs(perm), 95)); s["placebo_p_value"] = float(np.mean(np.abs(perm) >= abs(s["corr"])))
            refit = [np.corrcoef(yte, ridge_pred(ridge_fit(Xtr, rng.permutation(ytr), best[1]), Xte))[0, 1] for _ in range(200)]   # refit on permuted training targets (a random-weight model)
            s["refit_placebo_corr_p95"] = float(np.nanpercentile(np.abs(refit), 95))
            res[t] = s; preds[(name, t)] = (m, p)
        out["feature_sets"][name] = res
        print(name, {t: (round(res[t]["corr"], 2), round(res[t]["sign_hit"], 2)) for t in TARGETS})
    # dated-checkpoint arm for the RoBERTa features: correlation inside the model's training period vs after its checkpoint (FOMC only)
    if models:
        fo = df[df.bank == "FOMC"]; ins = fo[fo["date"] <= "2022-10-31"]; post = fo[fo["date"] > "2023-09-26"]
        out["dated_checkpoint"] = {"checkpoint": "2023-09-26", "training_data_end": "2022-10", "corr_redline_r2_inside": float(np.corrcoef(ins["rob_redline"], ins["r2"])[0, 1]), "n_inside": int(len(ins)), "corr_redline_r2_after": float(np.corrcoef(post["rob_redline"], post["r2"])[0, 1]), "n_after": int(len(post)),
                                  "corr_dh_r2_inside": float(np.corrcoef(ins["rob_dh"], ins["r2"])[0, 1]), "corr_dh_r2_after": float(np.corrcoef(post["rob_dh"], post["r2"])[0, 1])}
    # size model: |r2| from novelty, n_added, dispersion, length (+ recent realised vol of the 2y as the non-text control)
    vol = (dy["2Yr"].rolling(20).std()).reindex(df["date"], method="ffill").values
    df["rvol20"] = vol; size_sets = {"text_only": ["novelty", "n_added", "n_sent", "dict_disp"] + (["rob_disp"] if models else []), "vol_only": ["rvol20"], "text_plus_vol": ["novelty", "n_added", "n_sent", "dict_disp", "rvol20"] + (["rob_disp"] if models else [])}
    out["size_model"] = {}
    tr = df[df["date"] <= TRAIN_END]; te = df[df["date"] > TRAIN_END]
    for name, fs in size_sets.items():
        Xtr, Xte = tr[fs].values.astype(float), te[fs].values.astype(float); ytr, yte = np.abs(tr["r2"].values), np.abs(te["r2"].values)
        m = ridge_fit(Xtr, ytr, 3.0); p = ridge_pred(m, Xte); out["size_model"][name] = {"corr": float(np.corrcoef(yte, p)[0, 1]), "corr_ci": boot_ci(lambda ix: np.corrcoef(yte[ix], p[ix])[0, 1], len(yte)), "n": int(len(yte)), "mae_bp": float(np.abs(yte - p).mean())}
        print("size", name, round(out["size_model"][name]["corr"], 2))
    # per-event alpha vector for the replay: the hand-scored set with the best training-window LOO correlation, prediction per tenor, fitted on the training window only (test events) and
    # leave-one-out on training events, so no event's own reaction enters its alpha
    best_set = max((k for k in out["feature_sets"] if k != "embedding"), key=lambda k: out["feature_sets"][k]["r2"]["corr_loo_train"])   # chosen on the training window (leave-one-out), never on the test set
    fs = sets[best_set]; alpha = {t: np.zeros(len(df)) for t in TARGETS}
    for t in TARGETS:
        lam = out["feature_sets"][best_set][t]["lambda"]; Xa = X_of(df, fs); ya = df[t].values; is_tr = (df["date"] <= TRAIN_END).values
        m = ridge_fit(Xa[is_tr], ya[is_tr], lam)
        for i in range(len(df)):
            if is_tr[i]: mi = ridge_fit(np.delete(Xa[is_tr], np.where(np.where(is_tr)[0] == i)[0], 0), np.delete(ya[is_tr], np.where(np.where(is_tr)[0] == i)[0]), lam); alpha[t][i] = ridge_pred(mi, Xa[i:i + 1])[0]
            else: alpha[t][i] = ridge_pred(m, Xa[i:i + 1])[0]
    feat = df[["bank", "date", "novelty", "n_added", "dict_dh", "dict_redline"] + (["rob_dh", "rob_redline", "rob_disp"] if models else []) + TARGETS + ["src"]].copy()
    for t in TARGETS: feat["alpha_" + t] = alpha[t]
    feat["date"] = feat["date"].dt.strftime("%Y-%m-%d"); feat.to_csv(os.path.join(ROOT, "data/derived/nlp_features.csv"), index=False)
    out["alpha_source"] = best_set; json.dump(out, open(os.path.join(ROOT, "results/nlp_signal.json"), "w"), indent=1, default=float)
    print("alpha from", best_set, "; wrote nlp_features.csv and nlp_signal.json")

if __name__ == "__main__": main()
