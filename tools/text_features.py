"""Hawkishness scores for FOMC and ECB statements, sentence by sentence, with three scorers:
  dict     : a fixed hawkish/dovish lexicon with negation (no model, no look-ahead)
  fomc_rob : FOMC-RoBERTa (Shah, Paturi, Chava 2023; open re-upload tim9510019/FOMC-RoBERTa, commit 2023-09-26,
             trained on statements up to Oct 2022) -> the dated-checkpoint arm: only statements after 2023-09-26 are
             out of sample for it
  lorenzo  : LorenzoAleCon29/roberta-base-{fomc,ecb}-hawkish-dovish (commit 2026-08-13; training range unknown, so
             every score before that date may be contaminated - used as a second opinion only)
Each document gets h = mean(P(hawk) - P(dove)) over sentences and the sentence-level dispersion.  Model commit hashes
are written to the output so the checkpoint used is on record.
Writes data/derived/text_scores.csv.
"""
import os, re, sys, glob, json, time, hashlib
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("HF_HOME", os.path.join(ROOT, "storage", "hf"))

HAWK = ["raise", "raising", "raised", "increase in the target", "tighten", "tightening", "tighter", "higher rates", "inflation pressures", "inflation remains elevated",
        "elevated", "upside risks", "overheating", "strong", "strengthen", "robust", "solid", "firming", "above target", "persistent", "restrictive", "reduce its holdings",
        "balance sheet reduction", "additional firming", "further increases", "hike", "vigilant", "price stability risks", "wage pressures", "too high"]
DOVE = ["lower", "lowering", "lowered", "cut", "cuts", "reduce the target", "ease", "easing", "accommodative", "accommodation", "downside risks", "weak", "weaken", "weakness",
        "slow", "slowing", "slowed", "soft", "softening", "subdued", "below target", "moderate", "moderated", "moderating", "patient", "uncertainty", "deteriorat", "asset purchases",
        "purchase", "support the economy", "stimulus", "disinflation", "eased", "recession", "unemployment rate has risen", "cooling"]
NEG = ["not", "no", "less", "without", "neither", "nor"]

def split_sentences(text):
    text = re.sub(r"\s+", " ", text)
    s = re.split(r"(?<=[.!?])\s+(?=[A-Z])", text)
    return [x.strip() for x in s if 25 <= len(x.strip()) <= 600]

def dict_score(sent):
    s = sent.lower(); h = sum(s.count(w) for w in HAWK); d = sum(s.count(w) for w in DOVE)
    neg = any(re.search(r"\b" + n + r"\b", s) for n in NEG)
    if neg: h, d = d * 0.5 + h * 0.5, h * 0.5 + d * 0.5   # a negated sentence carries mixed sign; halve the polarity
    return (h - d) / (h + d + 1.0)

def fomc_text(path):
    return open(path, encoding="utf-8", errors="ignore").read()

def ecb_text(path):
    s = open(path, encoding="utf-8", errors="ignore").read()
    starts = [s.find(m) for m in ["Good afternoon", "Ladies and gentlemen"]]; starts = [x for x in starts if x > 0]
    i = min(starts) if starts else 0
    ends = [s.find(m, i) for m in ["ready to take your questions", "at your disposal for questions", "We are now ready", "Question:"]]; ends = [x for x in ends if x > 0]
    j = min(ends) if ends else i + 8000
    return s[i:j]

class HF:
    def __init__(self, name):
        from transformers import AutoTokenizer, AutoModelForSequenceClassification
        from huggingface_hub import model_info
        import torch
        self.torch = torch; self.name = name; self.sha = model_info(name).sha
        self.tok = AutoTokenizer.from_pretrained(name, revision=self.sha); self.m = AutoModelForSequenceClassification.from_pretrained(name, revision=self.sha).eval()
        labels = {v.lower(): k for k, v in self.m.config.id2label.items()}
        # gtfintechlab convention for LABEL_i: 0 dovish, 1 hawkish, 2 neutral
        self.i_dove = labels.get("dovish", 0); self.i_hawk = labels.get("hawkish", 1)
    def scores(self, sents, bs=16):
        out = []
        for k in range(0, len(sents), bs):
            x = self.tok(sents[k:k + bs], return_tensors="pt", padding=True, truncation=True, max_length=256)
            with self.torch.no_grad(): p = self.torch.softmax(self.m(**x).logits, -1).numpy()
            out.extend((p[:, self.i_hawk] - p[:, self.i_dove]).tolist())
        return out

def main():
    docs = []
    for f in sorted(glob.glob(os.path.join(ROOT, "data/raw/fomc/*.txt"))): docs.append(("FOMC", os.path.basename(f)[:8], fomc_text(f)))
    for f in sorted(glob.glob(os.path.join(ROOT, "data/raw/ecb/*.txt"))): docs.append(("ECB", os.path.basename(f)[:8], ecb_text(f)))
    use_models = "--no-models" not in sys.argv
    models = {}
    if use_models:
        try:
            models["fomc_rob"] = HF("tim9510019/FOMC-RoBERTa")
            models["lorenzo_fomc"] = HF("LorenzoAleCon29/roberta-base-fomc-hawkish-dovish")
            models["lorenzo_ecb"] = HF("LorenzoAleCon29/roberta-base-ecb-hawkish-dovish")
        except Exception as e:
            print("models unavailable, dictionary only:", e); models = {}
    rows = []; t0 = time.time()
    for bank, d, text in docs:
        sents = split_sentences(text)
        if not sents: continue
        row = {"bank": bank, "date": f"{d[:4]}-{d[4:6]}-{d[6:]}", "n_sentences": len(sents), "text_sha1": hashlib.sha1(text.encode("utf-8", "ignore")).hexdigest()[:12]}
        ds = np.array([dict_score(s) for s in sents]); row["dict_h"] = ds.mean(); row["dict_disp"] = ds.std()
        if models:
            r = np.array(models["fomc_rob"].scores(sents)); row["fomc_rob_h"] = r.mean(); row["fomc_rob_disp"] = r.std()
            l = np.array(models["lorenzo_fomc" if bank == "FOMC" else "lorenzo_ecb"].scores(sents)); row["lorenzo_h"] = l.mean(); row["lorenzo_disp"] = l.std()
        rows.append(row); print(f"{bank} {row['date']} n={len(sents)} dict={row['dict_h']:+.3f}" + (f" rob={row['fomc_rob_h']:+.3f} lor={row['lorenzo_h']:+.3f}" if models else "") + f"  [{time.time() - t0:.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["bank", "date"])
    os.makedirs(os.path.join(ROOT, "data/derived"), exist_ok=True); df.to_csv(os.path.join(ROOT, "data/derived/text_scores.csv"), index=False)
    meta = {"models": {k: {"name": m.name, "commit": m.sha} for k, m in models.items()}, "n_docs": len(df), "lexicon_hawk": HAWK, "lexicon_dove": DOVE,
            "notes": "FOMC-RoBERTa training statements end Oct 2022, checkpoint 2023-09-26: scores before that are in-sample for the model; LorenzoAleCon29 checkpoints dated 2026-08-13 with unknown training range"}
    json.dump(meta, open(os.path.join(ROOT, "data/derived/text_models.json"), "w"), indent=1)
    print("wrote", len(df), "documents")

if __name__ == "__main__": main()
