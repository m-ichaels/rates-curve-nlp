#!/usr/bin/env python3
"""Free data for the curve, the events and the positioning features.  python tools/download.py data/raw

  Treasury CMT par yields (daily, per year)        home.treasury.gov
  SOFR / EFFR fixings and volumes                  markets.newyorkfed.org
  SRF (standing repo facility) usage               markets.newyorkfed.org (repo operations)
  FOMC statements (text)                           federalreserve.gov
  ECB monetary policy statements (text)            ecb.europa.eu (best effort)
  SF Fed US Monetary Policy Event-Study Database   frbsf.org (USMPD.xlsx)
  CFTC Traders in Financial Futures                cftc.gov (yearly zips)
  Treasury auction results                         api.fiscaldata.treasury.gov
CME SOFR futures settlements are NOT downloaded: cmegroup.com blocks non-browser clients; the
front-end builder accepts their `stlint` format if you have it (see curve.py).
"""
import datetime as dt
import io
import json
import os
import re
import sys
import time
import zipfile

import requests

S = requests.Session(); S.headers.update({"User-Agent": "Mozilla/5.0 (research; rates-curve-nlp)"})


def get(url, **kw):
    for k in range(3):
        try:
            r = S.get(url, timeout=60, **kw)
            if r.status_code == 200:
                return r
            print("  ", url[:90], r.status_code)
        except Exception as e:  # noqa: BLE001
            print("  ", url[:90], "ERR", str(e)[:60])
        time.sleep(2 * (k + 1))
    return None


def treasury_cmt(out, years):
    rows = []
    for y in years:
        r = get(f"https://home.treasury.gov/resource-center/data-chart-center/interest-rates/daily-treasury-rates.csv/{y}/all?type=daily_treasury_yield_curve&field_tdr_date_value={y}&page&_format=csv")
        if r is None:
            continue
        lines = r.text.strip().splitlines()
        hdr = [h.strip('"') for h in lines[0].split(",")]
        for line in lines[1:]:
            rows.append(dict(zip(hdr, [c.strip('"') for c in line.split(",")])))
        print("  CMT", y, len(lines) - 1, "days")
    cols = ["1 Mo", "2 Mo", "3 Mo", "4 Mo", "6 Mo", "1 Yr", "2 Yr", "3 Yr", "5 Yr", "7 Yr", "10 Yr", "20 Yr", "30 Yr"]
    with open(os.path.join(out, "treasury_cmt.csv"), "w") as f:
        f.write("date," + ",".join(c.replace(" ", "") for c in cols) + "\n")
        for r in sorted(rows, key=lambda r: dt.datetime.strptime(r["Date"], "%m/%d/%Y")):
            d = dt.datetime.strptime(r["Date"], "%m/%d/%Y").date().isoformat()
            f.write(d + "," + ",".join(r.get(c, "") or "" for c in cols) + "\n")


def nyfed_rates(out):
    r = get("https://markets.newyorkfed.org/api/rates/all/search.csv?startDate=2018-04-01&endDate=2030-01-01&type=rate")
    if r: open(os.path.join(out, "nyfed_rates.csv"), "w").write(r.text); print("  NY Fed rates", r.text.count("\n"), "rows")
    # repo operations (SRF usage): try the documented endpoints
    for u in ["https://markets.newyorkfed.org/api/rp/repo/all/results/search.csv?startDate=2021-07-01&endDate=2030-01-01",
              "https://markets.newyorkfed.org/api/rp/all/all/results/search.csv?startDate=2021-07-01&endDate=2030-01-01"]:
        r = get(u)
        if r and r.text.strip():
            open(os.path.join(out, "nyfed_repo_ops.csv"), "w").write(r.text); print("  NY Fed repo ops", r.text.count("\n"), "rows"); break


def fomc_statements(out):
    os.makedirs(os.path.join(out, "fomc"), exist_ok=True)
    links = set()
    for u in ["https://www.federalreserve.gov/monetarypolicy/fomccalendars.htm"] + [f"https://www.federalreserve.gov/monetarypolicy/fomchistorical{y}.htm" for y in range(2015, 2021)]:
        r = get(u)
        if r: links |= set(re.findall(r'/newsevents/pressreleases/monetary(\d{8})a\.htm', r.text))
    n = 0
    for d in sorted(links):
        p = os.path.join(out, "fomc", f"{d}.txt")
        if os.path.exists(p):
            n += 1; continue
        r = get(f"https://www.federalreserve.gov/newsevents/pressreleases/monetary{d}a.htm")
        if not r: continue
        m = re.search(r'<div class="col-xs-12 col-sm-8 col-md-8">(.*?)</div>\s*<div class="col-xs-12 col-sm-4', r.text, re.S)
        body = m.group(1) if m else r.text
        text = re.sub(r"<[^>]+>", " ", body); text = re.sub(r"&nbsp;|&#160;", " ", text); text = re.sub(r"\s+", " ", text).strip()
        open(p, "w", encoding="utf-8").write(text); n += 1; time.sleep(0.5)
    print("  FOMC statements", n)


def ecb_statements(out):
    os.makedirs(os.path.join(out, "ecb"), exist_ok=True)
    n = 0
    for y in range(2015, 2027):
        for u in [f"https://www.ecb.europa.eu/press/press_conference/monetary-policy-statement/{y}/html/index_include.en.html",
                  f"https://www.ecb.europa.eu/press/pressconf/{y}/html/index_include.en.html"]:
            r = get(u)
            if not r: continue
            for href in sorted(set(re.findall(r'href="(/press/(?:press_conference/monetary-policy-statement|pressconf)/\d{4}/html/[^"]+\.en\.html)"', r.text))):
                m = re.search(r"(?:is|ecb\.is)(\d{6})", href)
                if not m: continue
                d = "20" + m.group(1)
                p = os.path.join(out, "ecb", f"{d}.txt")
                if os.path.exists(p): n += 1; continue
                rr = get("https://www.ecb.europa.eu" + href)
                if not rr: continue
                text = re.sub(r"<[^>]+>", " ", rr.text); text = re.sub(r"\s+", " ", text).strip()
                open(p, "w", encoding="utf-8").write(text); n += 1; time.sleep(0.5)
            break
    print("  ECB statements", n)


def sffed_usmpd(out):
    r = get("https://www.frbsf.org/wp-content/uploads/USMPD.xlsx")
    if r: open(os.path.join(out, "USMPD.xlsx"), "wb").write(r.content); print("  USMPD.xlsx", len(r.content))


def cftc_tff(out, years):
    os.makedirs(os.path.join(out, "cftc"), exist_ok=True)
    for y in years:
        p = os.path.join(out, "cftc", f"fut_fin_{y}.txt")
        if os.path.exists(p): continue
        r = get(f"https://www.cftc.gov/files/dea/history/fut_fin_txt_{y}.zip")
        if not r: continue
        z = zipfile.ZipFile(io.BytesIO(r.content)); name = z.namelist()[0]
        open(p, "wb").write(z.read(name)); print("  CFTC TFF", y)


def treasury_auctions(out):
    rows = []; page = 1
    while True:
        r = get(f"https://api.fiscaldata.treasury.gov/services/api/fiscal_service/v1/accounting/od/auctions_query?page%5Bsize%5D=1000&page%5Bnumber%5D={page}&sort=-auction_date&filter=auction_date:gte:2015-01-01")
        if not r: break
        d = r.json(); rows += d["data"]
        if page >= d["meta"]["total-pages"]: break
        page += 1
    keep = ["auction_date", "security_type", "security_term", "cusip", "maturity_date", "high_yield", "bid_to_cover_ratio", "offering_amt", "primary_dealer_accepted", "indirect_bidder_accepted", "direct_bidder_accepted", "allocation_pctage", "avg_med_yield", "int_rate"]
    with open(os.path.join(out, "treasury_auctions.csv"), "w") as f:
        f.write(",".join(keep) + "\n")
        for r in rows: f.write(",".join(str(r.get(k, "") or "") for k in keep) + "\n")
    print("  auctions", len(rows))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "data/raw"; os.makedirs(out, exist_ok=True)
    years = range(2015, dt.date.today().year + 1)
    for name, fn in [("Treasury CMT", lambda: treasury_cmt(out, years)), ("NY Fed", lambda: nyfed_rates(out)), ("FOMC", lambda: fomc_statements(out)),
                     ("ECB", lambda: ecb_statements(out)), ("SF Fed", lambda: sffed_usmpd(out)), ("CFTC", lambda: cftc_tff(out, years)), ("auctions", lambda: treasury_auctions(out))]:
        print(name); fn()


if __name__ == "__main__":
    main()
