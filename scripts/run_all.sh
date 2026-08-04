#!/usr/bin/env bash
# Full pipeline: data -> curve (daily builds, LOO, hedge test) -> factors -> closed-form frontier -> text/NLP/event features
# -> real-day replay -> figures -> summary -> report.   Usage: scripts/run_all.sh [--skip-download] [--skip-text]
set -euo pipefail
B=./build/rcmm; [[ -x $B.exe ]] && B=$B.exe
export HF_HOME="${HF_HOME:-$PWD/storage/hf}"
mkdir -p results/figures data/derived
[[ " $* " == *" --skip-download "* ]] || python tools/download.py
$B curve --date 2025-04-14
$B curvehist --from 2015-01-01 --to 2026-12-31 > results/curvehist.log 2>&1
python tools/hedge_test.py
$B factors
$B mm --days 40 --seeds 8 --gammas 5e-7,1e-6,2e-6,5e-6,1e-5,2e-5,5e-5 > results/mm.log 2>&1
[[ " $* " == *" --skip-text "* ]] || python tools/text_features.py
[[ " $* " == *" --skip-text "* ]] || python tools/nlp_features.py train_end=2023-12-31
python tools/event_features.py scorer=fomc_rob train_end=2023-12-31
EV="$B events --features data/derived/events.csv --from 2024-01-01 --to 2026-12-31 --train 2015-01-01:2023-12-31 --seeds 3"
$EV > results/events.log 2>&1
$EV --exclude-2020-22 --out results/events_ex2020.json > results/events_ex2020.log 2>&1
./build/rcmm_tests | tee results/tests.txt
python scripts/plots.py results results/figures
python scripts/summarize.py results
python scripts/report.py results report.pdf
echo done
