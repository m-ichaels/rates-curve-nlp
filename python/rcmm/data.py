"""Loaders for the downloaded data (port of include/rcmm/data.hpp): Treasury CMT par yields as bootstrap
instruments, NY Fed SOFR / EFFR fixings, and the FOMC / ECB statement tables written by tools/download.py.
Yields come back as decimals; missing CMT tenors as NaN."""
from __future__ import annotations

import os
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .curve import InstKind, Instrument

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CMT_T = np.array([1 / 12, 2 / 12, 0.25, 4 / 12, 0.5, 1, 2, 3, 5, 7, 10, 20, 30], dtype=float)
CMT_NAMES = ["1M", "2M", "3M", "4M", "6M", "1Y", "2Y", "3Y", "5Y", "7Y", "10Y", "20Y", "30Y"]
CMT_COLUMNS = ["1Mo", "2Mo", "3Mo", "4Mo", "6Mo", "1Yr", "2Yr", "3Yr", "5Yr", "7Yr", "10Yr", "20Yr", "30Yr"]
FACTOR_TENORS = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12]     # all but the 4M column


def path(*parts: str) -> str:
    return os.path.join(ROOT, *parts)


def read_cmt(csv: str | None = None) -> pd.DataFrame:
    """Daily CMT par yields (decimals), indexed by date, columns CMT_COLUMNS.  Default: data/raw/treasury_cmt.csv."""
    df = pd.read_csv(csv or path("data", "raw", "treasury_cmt.csv"))
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values("date").set_index("date")
    for c in CMT_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    return df[CMT_COLUMNS].astype(float) / 100.0


def cmt_instruments(row: pd.Series | Sequence[float], drop_4m: bool = True) -> list[Instrument]:
    """Instruments from one CMT row: bills (<= 1y) as bond-equivalent simple yields, the rest as
    semi-annual par bonds.  The 4M tenor is dropped by default (it only exists from late 2022)."""
    y = np.asarray(row, dtype=float)
    out = []
    for k in range(13):
        if np.isnan(y[k]) or (drop_4m and k == 3):
            continue
        kind = InstKind.BILL if CMT_T[k] <= 1.0 + 1e-9 else InstKind.PAR_BOND
        out.append(Instrument(kind, float(CMT_T[k]), float(y[k]), 2))
    return out


def daily_changes_bp(cmt: pd.DataFrame, tenors: Iterable[int] = FACTOR_TENORS, d0: str | None = None,
                     d1: str | None = None) -> pd.DataFrame:
    """Day-to-day changes in bp on the chosen tenors, matching FactorModel::fit in the C++: a day is
    kept when it and the preceding row in the file both have every tenor, and its date lies in [d0, d1]
    (the first in-window day differences against the last row before the window)."""
    cols = [CMT_COLUMNS[i] for i in tenors]
    full = cmt[cols]
    ok = full.notna().all(axis=1)
    d = (full - full.shift(1)) * 1e4
    keep = ok & ok.shift(1, fill_value=False)
    if d0 is not None:
        keep &= d.index >= pd.Timestamp(d0)
    if d1 is not None:
        keep &= d.index <= pd.Timestamp(d1)
    return d[keep]


def read_nyfed(csv: str | None = None) -> pd.DataFrame:
    """NY Fed reference rates -> one row per date with sofr, effr, sofr_p99 (decimals) and sofr_vol ($bn)."""
    p = csv or path("data", "raw", "nyfed_rates.csv")
    if not os.path.exists(p):
        return pd.DataFrame(columns=["sofr", "effr", "sofr_p99", "sofr_vol"])
    raw = pd.read_csv(p)
    raw.columns = [c.strip() for c in raw.columns]
    date_col, type_col, rate_col = raw.columns[0], raw.columns[1], raw.columns[2]
    raw["date"] = pd.to_datetime(raw[date_col])
    out = pd.DataFrame(index=sorted(raw["date"].unique()))
    for name, key in (("sofr", "SOFR"), ("effr", "EFFR")):
        s = raw[raw[type_col] == key].set_index("date")
        out[name] = s[rate_col].astype(float) / 100.0
        if key == "SOFR":
            out["sofr_p99"] = s[raw.columns[6]].astype(float) / 100.0
            out["sofr_vol"] = s[raw.columns[7]].astype(float)
    return out


def read_statements(bank: str = "FOMC", raw_dir: str | None = None) -> pd.DataFrame:
    """Statement texts downloaded by tools/download.py (data/raw/<bank>/YYYYMMDD.txt): bank, date, text."""
    import glob
    d = raw_dir or path("data", "raw", bank.lower())
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "*.txt"))):
        with open(f, encoding="utf-8") as fh:
            rows.append({"bank": bank.upper(), "date": pd.Timestamp(os.path.basename(f)[:8]), "text": fh.read()})
    return pd.DataFrame(rows, columns=["bank", "date", "text"])


def read_events(csv: str | None = None) -> pd.DataFrame:
    """data/derived/events.csv (the replay's feature table): one row per business day."""
    df = pd.read_csv(csv or path("data", "derived", "events.csv"))
    df["date"] = pd.to_datetime(df["date"])
    return df.set_index("date")


def read_text_scores(csv: str | None = None) -> pd.DataFrame:
    """data/derived/text_scores.csv: per-statement hawkishness under each scorer."""
    df = pd.read_csv(csv or path("data", "derived", "text_scores.csv"))
    df["date"] = pd.to_datetime(df["date"])
    return df


def read_nlp_features(csv: str | None = None) -> pd.DataFrame:
    """data/derived/nlp_features.csv: redline / novelty / dispersion features and event-window moves."""
    df = pd.read_csv(csv or path("data", "derived", "nlp_features.csv"))
    df["date"] = pd.to_datetime(df["date"])
    return df
