# Results summary

## Curve construction (2025-04-14)

Bootstrap (Hagan-West monotone convex, global iteration): max repricing error 2.84e-10 bp in 5 sweeps. Kernel ridge (lambda = 1e-04): max repricing error 0.092 bp. Forward roughness int f'(t)^2 dt: bootstrap 2.13e-04, kernel ridge 1.65e-04.

Leave-one-out pricing error, interior tenors 6m-10y (mean |err|): bootstrap 10.85 bp, kernel ridge 9.08 bp. Endpoints are extrapolation and excluded (20y: 224 / 347 bp; 30y: -361 / -1511 bp).

| tenor | par % | LOO boot (bp) | LOO KR (bp) |
|---|---|---|---|
| 0.0833333y | 4.340 | -0.48 | -0.40 |
| 0.166667y | 4.370 | +0.47 | +0.43 |
| 0.25y | 4.330 | -0.24 | -0.36 |
| 0.5y | 4.210 | +4.15 | +1.05 |
| 1y | 3.990 | -2.36 | -5.87 |
| 2y | 3.840 | -2.14 | +0.08 |
| 3y | 3.870 | -0.35 | +0.46 |
| 5y | 4.020 | -8.81 | -7.28 |
| 7y | 4.200 | +14.18 | +17.61 |
| 10y | 4.380 | -54.56 | -39.92 |
| 20y | 4.840 | +224.49 | +347.12 |
| 30y | 4.800 | -360.87 | -1511.17 |

Ladder key-rate DV01s (per 1 notional per bp; row = bond, columns = bumped input):

| bond | coupon % | DV01 | 0.0833333y | 0.166667y | 0.25y | 0.5y | 1y | 2y | 3y | 5y | 7y | 10y | 20y | 30y |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2y | 3.840 | 1.91e-04 | -2.2e-16 | -2.2e-16 | -2.2e-16 | 2.9e-15 | 0.0e+00 | 1.9e-04 | -8.9e-16 | -2.2e-16 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 |
| 5y | 4.020 | 4.50e-04 | 0.0e+00 | 1.1e-16 | 0.0e+00 | 1.1e-16 | 0.0e+00 | 0.0e+00 | -1.3e-15 | 4.5e-04 | -3.1e-14 | 1.1e-14 | -6.7e-16 | -4.4e-16 |
| 10y | 4.380 | 8.10e-04 | 0.0e+00 | -2.2e-16 | -2.2e-16 | 0.0e+00 | 0.0e+00 | -2.2e-16 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 8.1e-04 | 0.0e+00 | -2.2e-16 |
| 30y | 4.800 | 1.59e-03 | 2.2e-16 | 0.0e+00 | 0.0e+00 | 2.2e-16 | 5.6e-16 | 0.0e+00 | 0.0e+00 | 4.4e-16 | 4.4e-16 | -2.2e-16 | 0.0e+00 | 1.6e-03 |

Meeting-date step front end (bills + fixing (proxy)): levels 4.34%, 4.30%, 4.10%, 3.93%, 3.79%, 3.68%, 3.61%, 3.57%, 3.57% on 0.12y, 0.25y, 0.37y, 0.50y, 0.62y, 0.75y, 0.87y, 1.00y.

## Daily curve builds, 2015-01-02 .. 2026-09-15 (2927 days)

Selection rule fixed in advance: the smoothest kernel ridge (largest lambda) whose max repricing error over all days is within 1 bp; that is lambda = 1e-05. Leave-one-out = each interior tenor (6m-10y) dropped and repriced off the curve built from the rest; hedge test = a 3y (7y) par bond priced off the curve built without that point and hedged with the remaining inputs using that curve's key-rate DV01s, realised next-day hedged P&L in yield bp.

| method | lambda | max reprice bp | LOO mean |err| yield bp | LOO median | 10y DV01 day-to-day std | forward roughness | hedged 3y std bp | hedged 7y std bp |
|---|---|---|---|---|---|---|---|---|
| bootstrap | - | 2.8e-07 | 3.686 | 3.457 | 2.248 | 3.44e-04 | 1.822 | 0.965 |
| kernel_ridge | 1e-02 | 10 | 3.655 | 3.089 | 2.251 | 1.52e-04 | 1.331 | 0.827 |
| kernel_ridge | 1e-03 | 5.3 | 3.503 | 2.996 | 2.252 | 2.07e-04 | 1.424 | 0.866 |
| kernel_ridge | 1e-04 | 2.2 | 3.514 | 3.066 | 2.252 | 2.64e-04 | 1.446 | 0.872 |
| kernel_ridge | 1e-05 | 0.32 | 3.553 | 3.097 | 2.252 | 3.00e-04 | 1.451 | 0.873 |
| kernel_ridge | 1e-06 | 0.034 | 3.560 | 3.105 | 2.252 | 3.07e-04 | 1.451 | 0.873 |

Selected kernel ridge vs bootstrap: LOO error -3.6% (median -10.4%), roughness -12.8%, hedged 3y / 7y std -20.4% / -9.6%, 10y DV01 day-to-day std +0.2%. Unhedged 7y daily std 5.55 bp.

| year | 2015 | 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 | 2025 | 2026 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bootstrap | 3.42 | 3.27 | 2.59 | 2.08 | 2.13 | 2.43 | 3.71 | 6.01 | 5.66 | 4.38 | 3.70 | 5.36 |
| KR 1e-02 | 2.87 | 2.92 | 2.30 | 2.15 | 2.07 | 2.57 | 3.06 | 7.61 | 6.02 | 4.65 | 3.41 | 4.50 |
| KR 1e-03 | 3.00 | 3.00 | 2.27 | 2.12 | 2.07 | 2.57 | 2.99 | 6.82 | 5.57 | 4.27 | 3.10 | 4.57 |
| KR 1e-04 | 3.08 | 3.07 | 2.27 | 2.11 | 2.14 | 2.61 | 3.01 | 6.70 | 5.53 | 4.26 | 3.15 | 4.57 |
| KR 1e-05 | 3.09 | 3.08 | 2.27 | 2.12 | 2.18 | 2.64 | 3.03 | 6.80 | 5.62 | 4.30 | 3.23 | 4.60 |
| KR 1e-06 | 3.09 | 3.09 | 2.27 | 2.12 | 2.18 | 2.64 | 3.04 | 6.82 | 5.64 | 4.31 | 3.24 | 4.61 |

## Learned vs model hedge ratios (walk-forward from 2016-01-04, trailing window 250 days)

| bond | hedge ratios | days | hedged P&L std (yield bp) | reduction vs unhedged |
|---|---|---|---|---|
| 3y | unhedged | 2675 | 5.473 | 0.0% |
| 3y | 50/50 neighbours | 2675 | 1.159 | 78.8% |
| 3y | bootstrap key-rate DV01s | 2675 | 1.868 | 65.9% |
| 3y | kernel-ridge key-rate DV01s | 2675 | 1.456 | 73.4% |
| 3y | PCA level+slope neutral (trailing window) | 2675 | 1.230 | 77.5% |
| 3y | ridge regression on neighbours (trailing window) | 2675 | 1.158 | 78.8% |
| 7y | unhedged | 2675 | 5.555 | 0.0% |
| 7y | 50/50 neighbours | 2675 | 0.745 | 86.6% |
| 7y | bootstrap key-rate DV01s | 2675 | 0.974 | 82.5% |
| 7y | kernel-ridge key-rate DV01s | 2675 | 0.873 | 84.3% |
| 7y | PCA level+slope neutral (trailing window) | 2675 | 0.731 | 86.8% |
| 7y | ridge regression on neighbours (trailing window) | 2675 | 0.728 | 86.9% |
| 20y | unhedged | 2676 | 5.084 | 0.0% |
| 20y | 50/50 neighbours | 2676 | 0.935 | 81.6% |
| 20y | PCA level+slope neutral (trailing window) | 2676 | 0.885 | 82.6% |
| 20y | ridge regression on neighbours (trailing window) | 2676 | 0.868 | 82.9% |

3y: learned / bootstrap-key-rate std ratio 0.620 (bootstrap 95% CI [0.552, 0.692]).
7y: learned / bootstrap-key-rate std ratio 0.747 (bootstrap 95% CI [0.699, 0.793]).

## Factor model (PCA of daily par-yield changes, 12 tenors)

| window | days | explained 1/2/3 | factor sd bp/day (level, slope, curvature) |
|---|---|---|---|
| 2015-01-01:2026-12-31 | 1977 | 0.708 / 0.817 / 0.908 | 15.0, 5.9, 5.4 |
| 2023-01-01:2026-12-31 | 926 | 0.713 / 0.830 / 0.924 | 16.2, 6.6, 5.9 |
| 2022-01-01:2022-12-31 | 249 | 0.723 / 0.823 / 0.902 | 20.8, 7.7, 6.8 |
| 2020-01-01:2022-12-31 | 751 | 0.702 / 0.831 / 0.908 | 15.0, 6.4, 5.0 |

## Closed-form quoters on synthetic days (40 days x 8 seeds, common random numbers)

P&L per day with 95% bootstrap CI over seeds; 5-minute P&L std; decomposition (spread + inventory - hedge cost = total, residual is machine zero); hit ratios by tier (informed / real money / other).

| quoter | gamma | P&L/day $k | CI | 5-min std $k | P&L/var x1e-3 | spread | inventory | hedge | max |resid| | hit ratios |
|---|---|---|---|---|---|---|---|---|---|---|
| Symmetric | 5e-07 | 122.0 | [107.1, 135.8] | 11.9 | 0.86 | 170.4 | -29.3 | -19.1 | 0.0e+00 | 0.34/0.36/0.35 |
| Symmetric | 1e-06 | 122.0 | [107.1, 135.8] | 11.9 | 0.86 | 170.4 | -29.3 | -19.1 | 0.0e+00 | 0.34/0.36/0.35 |
| Symmetric | 2e-06 | 122.0 | [107.1, 135.8] | 11.9 | 0.86 | 170.4 | -29.3 | -19.1 | 0.0e+00 | 0.34/0.36/0.35 |
| Symmetric | 5e-06 | 122.0 | [107.1, 135.8] | 11.9 | 0.86 | 170.4 | -29.3 | -19.1 | 0.0e+00 | 0.34/0.36/0.35 |
| Symmetric | 1e-05 | 122.0 | [107.1, 135.8] | 11.9 | 0.86 | 170.4 | -29.3 | -19.1 | 0.0e+00 | 0.34/0.36/0.35 |
| Symmetric | 2e-05 | 122.0 | [107.1, 135.8] | 11.9 | 0.86 | 170.4 | -29.3 | -19.1 | 0.0e+00 | 0.34/0.36/0.35 |
| Symmetric | 5e-05 | 122.0 | [107.1, 135.8] | 11.9 | 0.86 | 170.4 | -29.3 | -19.1 | 0.0e+00 | 0.34/0.36/0.35 |
| Bergault | 5e-07 | 124.7 | [113.1, 135.2] | 10.0 | 1.26 | 166.6 | -25.9 | -16.0 | 0.0e+00 | 0.37/0.35/0.32 |
| Bergault | 1e-06 | 125.6 | [114.8, 134.9] | 9.4 | 1.43 | 163.7 | -23.2 | -14.9 | 0.0e+00 | 0.36/0.35/0.32 |
| Bergault | 2e-06 | 127.8 | [118.7, 136.6] | 8.9 | 1.61 | 161.5 | -19.8 | -13.9 | 0.0e+00 | 0.35/0.35/0.31 |
| Bergault | 5e-06 | 126.6 | [119.4, 132.7] | 8.2 | 1.87 | 152.8 | -14.5 | -11.7 | 0.0e+00 | 0.34/0.33/0.30 |
| Bergault | 1e-05 | 122.3 | [117.5, 127.5] | 7.5 | 2.17 | 144.8 | -12.7 | -9.8 | 0.0e+00 | 0.33/0.32/0.29 |
| Bergault | 2e-05 | 114.2 | [108.4, 120.8] | 6.8 | 2.45 | 131.8 | -9.6 | -8.0 | 0.0e+00 | 0.30/0.30/0.28 |
| Bergault | 5e-05 | 99.5 | [93.8, 104.7] | 5.9 | 2.82 | 111.4 | -6.4 | -5.5 | 0.0e+00 | 0.26/0.27/0.26 |
| BarzykinCiceri | 5e-07 | 119.9 | [107.1, 130.1] | 9.9 | 1.21 | 165.4 | -27.7 | -17.8 | 0.0e+00 | 0.35/0.36/0.35 |
| BarzykinCiceri | 1e-06 | 122.0 | [112.5, 130.4] | 9.4 | 1.37 | 163.8 | -24.4 | -17.4 | 0.0e+00 | 0.35/0.36/0.35 |
| BarzykinCiceri | 2e-06 | 124.0 | [115.4, 131.1] | 8.9 | 1.55 | 161.5 | -20.5 | -16.9 | 0.0e+00 | 0.35/0.36/0.35 |
| BarzykinCiceri | 5e-06 | 123.4 | [114.3, 131.0] | 8.2 | 1.82 | 155.6 | -16.0 | -16.2 | 0.0e+00 | 0.36/0.36/0.36 |
| BarzykinCiceri | 1e-05 | 118.6 | [110.3, 125.8] | 7.8 | 1.96 | 148.3 | -14.0 | -15.6 | 0.0e+00 | 0.36/0.36/0.36 |
| BarzykinCiceri | 2e-05 | 108.1 | [101.3, 115.3] | 7.3 | 2.03 | 136.8 | -13.7 | -14.9 | 0.0e+00 | 0.36/0.36/0.36 |
| BarzykinCiceri | 5e-05 | 89.2 | [81.0, 97.2] | 6.8 | 1.92 | 118.5 | -14.8 | -14.5 | 0.0e+00 | 0.37/0.37/0.36 |
| BarzykinCiceri_QAHR | 5e-07 | 125.4 | [118.1, 132.7] | 9.8 | 1.30 | 158.6 | -19.2 | -14.0 | 0.0e+00 | 0.26/0.28/0.28 |
| BarzykinCiceri_QAHR | 1e-06 | 128.2 | [120.1, 136.0] | 9.4 | 1.46 | 157.1 | -15.6 | -13.4 | 0.0e+00 | 0.26/0.28/0.28 |
| BarzykinCiceri_QAHR | 2e-06 | 124.0 | [117.5, 130.4] | 8.8 | 1.59 | 155.8 | -18.7 | -13.2 | 0.0e+00 | 0.26/0.28/0.28 |
| BarzykinCiceri_QAHR | 5e-06 | 122.7 | [115.7, 129.5] | 8.1 | 1.86 | 150.5 | -15.3 | -12.6 | 0.0e+00 | 0.26/0.28/0.28 |
| BarzykinCiceri_QAHR | 1e-05 | 121.9 | [114.8, 129.7] | 7.7 | 2.06 | 146.2 | -12.2 | -12.1 | 0.0e+00 | 0.26/0.28/0.28 |
| BarzykinCiceri_QAHR | 2e-05 | 117.3 | [110.2, 124.5] | 7.2 | 2.28 | 138.2 | -9.2 | -11.7 | 0.0e+00 | 0.26/0.28/0.28 |
| BarzykinCiceri_QAHR | 5e-05 | 104.3 | [98.4, 108.6] | 6.7 | 2.34 | 125.4 | -10.0 | -11.1 | 0.0e+00 | 0.27/0.28/0.29 |
| CarteaWang | 5e-07 | 126.3 | [113.9, 138.3] | 9.8 | 1.31 | 164.0 | -21.7 | -16.0 | 0.0e+00 | 0.35/0.35/0.32 |
| CarteaWang | 1e-06 | 127.8 | [116.6, 138.5] | 9.3 | 1.48 | 161.4 | -18.6 | -15.0 | 0.0e+00 | 0.35/0.35/0.32 |
| CarteaWang | 2e-06 | 129.0 | [120.7, 137.8] | 8.9 | 1.64 | 158.8 | -15.9 | -13.8 | 0.0e+00 | 0.34/0.34/0.31 |
| CarteaWang | 5e-06 | 130.0 | [123.7, 137.4] | 8.2 | 1.93 | 150.7 | -9.0 | -11.8 | 0.0e+00 | 0.32/0.33/0.30 |
| CarteaWang | 1e-05 | 123.7 | [117.7, 130.9] | 7.5 | 2.20 | 141.8 | -8.2 | -9.9 | 0.0e+00 | 0.31/0.32/0.29 |
| CarteaWang | 2e-05 | 115.4 | [109.2, 123.3] | 6.9 | 2.43 | 129.3 | -6.0 | -7.9 | 0.0e+00 | 0.29/0.30/0.28 |
| CarteaWang | 5e-05 | 101.1 | [96.8, 104.3] | 6.0 | 2.80 | 110.2 | -3.6 | -5.5 | 0.0e+00 | 0.26/0.27/0.26 |
| Bergault_oracle | 5e-07 | 126.3 | [116.9, 136.3] | 9.6 | 1.38 | 155.9 | -16.4 | -13.3 | 0.0e+00 | 0.22/0.31/0.32 |
| Bergault_oracle | 1e-06 | 129.2 | [119.9, 138.0] | 9.0 | 1.59 | 153.9 | -12.3 | -12.4 | 0.0e+00 | 0.21/0.31/0.32 |
| Bergault_oracle | 2e-06 | 128.7 | [124.0, 133.2] | 8.5 | 1.77 | 150.2 | -10.1 | -11.5 | 0.0e+00 | 0.21/0.30/0.31 |
| Bergault_oracle | 5e-06 | 127.0 | [121.6, 131.8] | 7.8 | 2.09 | 143.7 | -6.9 | -9.8 | 0.0e+00 | 0.21/0.29/0.30 |
| Bergault_oracle | 1e-05 | 116.6 | [112.8, 120.5] | 7.2 | 2.26 | 134.4 | -9.5 | -8.3 | 0.0e+00 | 0.20/0.28/0.29 |
| Bergault_oracle | 2e-05 | 110.5 | [106.1, 116.1] | 6.6 | 2.51 | 123.9 | -6.6 | -6.8 | 0.0e+00 | 0.19/0.26/0.28 |
| Bergault_oracle | 5e-05 | 93.8 | [88.3, 99.2] | 5.8 | 2.78 | 103.9 | -5.5 | -4.6 | 0.0e+00 | 0.18/0.24/0.25 |

Mark-outs at 30 minutes (bp against the dealer), gamma = 5e-06:

| quoter | tier | window | fills | half-spread bp | 5m | 30m | 60m |
|---|---|---|---|---|---|---|---|
| Symmetric | informed | normal | 278 | 0.45 | +0.35 | +0.34 | +0.31 |
| Symmetric | informed | pre-event | 14 | 0.45 | +2.94 | +2.94 | +3.10 |
| Symmetric | real_money | normal | 575 | 0.54 | +0.04 | +0.03 | +0.02 |
| Symmetric | real_money | pre-event | 39 | 0.54 | +0.14 | +0.13 | +0.21 |
| Symmetric | other | normal | 814 | 0.57 | +0.00 | +0.00 | +0.03 |
| Symmetric | other | pre-event | 16 | 0.57 | -1.01 | -1.01 | -0.90 |
| Bergault | informed | normal | 270 | 0.40 | +0.28 | +0.28 | +0.30 |
| Bergault | informed | pre-event | 11 | 0.37 | +1.90 | +1.90 | +2.19 |
| Bergault | real_money | normal | 519 | 0.53 | +0.01 | +0.01 | +0.04 |
| Bergault | real_money | pre-event | 39 | 0.50 | -0.29 | -0.29 | -0.31 |
| Bergault | other | normal | 700 | 0.61 | -0.06 | -0.06 | -0.04 |
| Bergault | other | pre-event | 14 | 0.55 | -0.50 | -0.39 | -0.33 |
| BarzykinCiceri | informed | normal | 282 | 0.39 | +0.35 | +0.34 | +0.28 |
| BarzykinCiceri | informed | pre-event | 15 | 0.38 | +2.15 | +2.15 | +2.32 |
| BarzykinCiceri | real_money | normal | 569 | 0.49 | -0.04 | -0.04 | -0.03 |
| BarzykinCiceri | real_money | pre-event | 49 | 0.40 | -0.81 | -0.81 | -0.80 |
| BarzykinCiceri | other | normal | 830 | 0.48 | -0.04 | -0.04 | -0.03 |
| BarzykinCiceri | other | pre-event | 17 | 0.49 | -0.59 | -0.57 | -0.42 |
| BarzykinCiceri_QAHR | informed | normal | 198 | 0.52 | +0.42 | +0.41 | +0.35 |
| BarzykinCiceri_QAHR | informed | pre-event | 13 | 0.41 | +4.12 | +4.12 | +4.47 |
| BarzykinCiceri_QAHR | real_money | normal | 429 | 0.62 | +0.09 | +0.09 | -0.00 |
| BarzykinCiceri_QAHR | real_money | pre-event | 42 | 0.47 | -1.04 | -1.03 | -1.01 |
| BarzykinCiceri_QAHR | other | normal | 666 | 0.61 | -0.08 | -0.08 | -0.07 |
| BarzykinCiceri_QAHR | other | pre-event | 13 | 0.57 | -1.97 | -1.97 | -1.83 |
| CarteaWang | informed | normal | 277 | 0.39 | +0.33 | +0.33 | +0.34 |
| CarteaWang | informed | pre-event | 4 | 0.80 | +1.30 | +1.30 | +1.47 |
| CarteaWang | real_money | normal | 518 | 0.52 | -0.03 | -0.02 | -0.02 |
| CarteaWang | real_money | pre-event | 31 | 0.62 | -0.63 | -0.62 | -0.58 |
| CarteaWang | other | normal | 715 | 0.59 | -0.09 | -0.09 | -0.06 |
| CarteaWang | other | pre-event | 9 | 0.70 | -0.80 | -0.80 | -0.89 |
| Bergault_oracle | informed | normal | 190 | 0.53 | +0.39 | +0.39 | +0.36 |
| Bergault_oracle | informed | pre-event | 0 | 0.00 | +0.00 | +0.00 | +0.00 |
| Bergault_oracle | real_money | normal | 469 | 0.57 | -0.02 | -0.01 | -0.02 |
| Bergault_oracle | real_money | pre-event | 8 | 1.29 | +0.18 | +0.02 | +0.10 |
| Bergault_oracle | other | normal | 699 | 0.61 | -0.04 | -0.03 | -0.04 |
| Bergault_oracle | other | pre-event | 14 | 0.55 | -0.52 | -0.40 | -0.48 |

## Text signal (fomc_rob), training window to 2023-12-31

Regression of the 10y reaction on the change in statement hawkishness, pooled FOMC + ECB: b1 = 4.4 bp per unit (bootstrap 95% CI [-0.1, 9.0]), n = 147, in-sample corr 0.150, R2 0.023. Out of sample (n = 44): corr 0.2407691934597632, R2 0.06749162310294932, sign hit rate 0.5227272727272727.

By bank (training window): FOMC: n 76, corr(dh) 0.125, corr(level) 0.117; ECB: n 71, corr(dh) 0.174, corr(level) -0.088.

Look-ahead placebo (FOMC-RoBERTa checkpoint 2023-09-26): corr(dh, reaction) inside the model's training period 0.1647988045987716 (n = 66) vs after the checkpoint 0.4481200049629283 (n = 24).

## NLP: statement text -> event-window move (train to 2023-12-31, n = 147; test n = 44)

Ridge on hand-scored features (lexicon and FOMC-RoBERTa, whole document and redline = sentences added minus removed vs the previous statement, novelty, dispersion) or on the 384-d MiniLM embedding of the redline; lambda by leave-one-out on the training set; placebo = permutation test of the out-of-sample correlation (95th percentile of |corr| under permutation, and the two-sided p-value).

| features | target | lambda | LOO corr (train) | OOS corr [CI] | OOS R2 | sign hit [CI] | placebo corr p95 | placebo p |
|---|---|---|---|---|---|---|---|---|
| lexicon | 2y | 100 | -0.33 | 0.13 [-0.08, 0.35] | 0.008 | 0.50 [0.34, 0.64] | 0.29 | 0.392 |
| lexicon | 5y | 10 | -0.28 | 0.02 [-0.29, 0.33] | -0.007 | 0.57 [0.43, 0.70] | 0.30 | 0.871 |
| lexicon | 10y | 10 | -0.12 | -0.02 [-0.34, 0.30] | -0.025 | 0.50 [0.34, 0.64] | 0.30 | 0.898 |
| lexicon | 30y | 1 | -0.00 | -0.05 [-0.31, 0.25] | -0.061 | 0.50 [0.36, 0.64] | 0.30 | 0.750 |
| roberta | 2y | 100 | 0.07 | 0.34 [0.11, 0.53] | 0.075 | 0.48 [0.32, 0.64] | 0.30 | 0.023 |
| roberta | 5y | 0.1 | 0.08 | 0.29 [0.02, 0.50] | 0.076 | 0.57 [0.43, 0.70] | 0.29 | 0.052 |
| roberta | 10y | 0.1 | 0.00 | 0.18 [-0.07, 0.41] | 0.015 | 0.48 [0.34, 0.61] | 0.30 | 0.235 |
| roberta | 30y | 30 | -0.12 | 0.01 [-0.26, 0.26] | -0.036 | 0.45 [0.32, 0.61] | 0.30 | 0.941 |
| all_hand | 2y | 100 | 0.00 | 0.31 [0.08, 0.51] | 0.070 | 0.48 [0.34, 0.64] | 0.30 | 0.037 |
| all_hand | 5y | 30 | 0.02 | 0.25 [-0.02, 0.45] | 0.060 | 0.59 [0.45, 0.73] | 0.29 | 0.100 |
| all_hand | 10y | 10 | -0.02 | 0.15 [-0.12, 0.40] | 0.011 | 0.50 [0.36, 0.64] | 0.30 | 0.325 |
| all_hand | 30y | 30 | -0.12 | 0.01 [-0.24, 0.26] | -0.039 | 0.43 [0.30, 0.57] | 0.30 | 0.935 |
| embedding | 2y | 1 | 0.09 | -0.03 [-0.36, 0.28] | -4.243 | 0.52 [0.39, 0.66] | 0.30 | 0.843 |
| embedding | 5y | 1 | -0.06 | -0.02 [-0.34, 0.27] | -7.494 | 0.55 [0.39, 0.68] | 0.30 | 0.905 |
| embedding | 10y | 300 | -0.12 | 0.05 [-0.26, 0.29] | -0.223 | 0.57 [0.41, 0.70] | 0.29 | 0.737 |
| embedding | 30y | 1000 | -0.14 | -0.12 [-0.47, 0.19] | -0.147 | 0.45 [0.30, 0.59] | 0.29 | 0.440 |
| roberta_redline_only | 2y | 0.1 | 0.12 | 0.23 [0.05, 0.41] | 0.040 | 0.55 [0.39, 0.70] | 0.30 | 0.131 |
| roberta_redline_only | 5y | 0.1 | 0.11 | 0.19 [0.04, 0.36] | 0.031 | 0.55 [0.39, 0.70] | 0.30 | 0.209 |
| roberta_redline_only | 10y | 0.1 | 0.03 | 0.12 [-0.04, 0.27] | 0.011 | 0.50 [0.34, 0.64] | 0.30 | 0.415 |
| roberta_redline_only | 30y | 0.1 | -0.37 | 0.03 [-0.17, 0.18] | -0.004 | 0.57 [0.43, 0.73] | 0.29 | 0.866 |
| roberta_doc_only | 2y | 0.1 | 0.17 | 0.34 [0.09, 0.54] | 0.091 | 0.57 [0.43, 0.73] | 0.29 | 0.022 |
| roberta_doc_only | 5y | 0.1 | 0.17 | 0.27 [0.03, 0.46] | 0.074 | 0.57 [0.43, 0.73] | 0.30 | 0.068 |
| roberta_doc_only | 10y | 0.1 | 0.07 | 0.16 [-0.09, 0.37] | 0.025 | 0.52 [0.39, 0.68] | 0.29 | 0.282 |
| roberta_doc_only | 30y | 0.1 | -0.53 | -0.01 [-0.27, 0.24] | -0.005 | 0.59 [0.45, 0.73] | 0.29 | 0.947 |

Dated-checkpoint arm (FOMC-RoBERTa checkpoint 2023-09-26, training data to 2022-10): corr(redline score, 2y reaction) inside the training period 0.23 (n = 66) vs after the checkpoint 0.37 (n = 24); document-level dh 0.27 vs 0.48.

Size of the move (|2y reaction|), OOS corr: text_only: -0.31 [-0.60, 0.10]; vol_only: 0.21 [-0.06, 0.49]; text_plus_vol: -0.00 [-0.25, 0.27]. Replay alphas come from the roberta_doc_only feature set.

## Real-day replay with event features: 2024-01-01..2026-12-31, 676 days, 186 event days (21 with measured USMPD jumps), gamma 2e-06

Events by type: AUCTION 111, ECB 22, FOMC 22, NFP 31.

| quoter | variant | P&L/day $k | CI | 5-min std $k | P&L/var x1e-3 | spread | inventory | hedge | hit ratios |
|---|---|---|---|---|---|---|---|---|---|
| Symmetric | text+calendar | 152.8 | [147.7, 157.2] | 12.5 | 0.97 | 185.4 | -10.1 | -22.4 | 0.35/0.35/0.35 |
| Bergault | text+calendar | 152.1 | [149.1, 153.8] | 9.7 | 1.62 | 173.5 | -6.2 | -15.1 | 0.35/0.33/0.31 |
| BarzykinCiceri_QAHR | text+calendar | 148.3 | [145.8, 151.9] | 9.7 | 1.59 | 169.4 | -6.0 | -15.1 | 0.26/0.27/0.28 |
| CarteaWang | text+calendar | 153.4 | [149.9, 155.5] | 9.6 | 1.65 | 172.3 | -4.0 | -14.9 | 0.34/0.33/0.31 |
| CarteaWang | placebo_text+calendar | 153.5 | [149.7, 155.6] | 9.7 | 1.64 | 172.2 | -3.8 | -14.9 | 0.34/0.33/0.31 |
| CarteaWang | calendar_only | 153.6 | [150.3, 155.4] | 9.6 | 1.65 | 172.4 | -3.9 | -14.9 | 0.34/0.33/0.31 |
| CarteaWang | none | 152.1 | [149.1, 153.8] | 9.7 | 1.62 | 173.5 | -6.2 | -15.1 | 0.35/0.33/0.31 |
| Bergault_oracle | text+calendar | 146.9 | [144.7, 148.2] | 9.3 | 1.70 | 160.6 | -1.3 | -12.4 | 0.21/0.29/0.31 |

| quoter | tier | window | fills | half-spread bp | 30m mark-out bp |
|---|---|---|---|---|---|
| Symmetric | informed | normal | 4843 | 0.45 | +0.14 |
| Symmetric | informed | pre-event | 282 | 0.45 | +1.99 |
| Symmetric | real_money | normal | 9630 | 0.54 | +0.04 |
| Symmetric | real_money | pre-event | 564 | 0.54 | +0.06 |
| Symmetric | other | normal | 14390 | 0.57 | -0.00 |
| Symmetric | other | pre-event | 294 | 0.57 | +0.21 |
| Bergault | informed | normal | 4956 | 0.41 | +0.13 |
| Bergault | informed | pre-event | 293 | 0.44 | +1.93 |
| Bergault | real_money | normal | 9167 | 0.54 | +0.05 |
| Bergault | real_money | pre-event | 541 | 0.54 | +0.03 |
| Bergault | other | normal | 12723 | 0.62 | -0.01 |
| Bergault | other | pre-event | 252 | 0.60 | -0.02 |
| BarzykinCiceri_QAHR | informed | normal | 3617 | 0.53 | +0.15 |
| BarzykinCiceri_QAHR | informed | pre-event | 209 | 0.54 | +2.15 |
| BarzykinCiceri_QAHR | real_money | normal | 7487 | 0.63 | +0.05 |
| BarzykinCiceri_QAHR | real_money | pre-event | 474 | 0.59 | +0.02 |
| BarzykinCiceri_QAHR | other | normal | 11547 | 0.67 | -0.01 |
| BarzykinCiceri_QAHR | other | pre-event | 243 | 0.65 | -0.07 |
| CarteaWang | informed | normal | 4954 | 0.41 | +0.14 |
| CarteaWang | informed | pre-event | 50 | 0.96 | +1.73 |
| CarteaWang | real_money | normal | 9181 | 0.54 | +0.05 |
| CarteaWang | real_money | pre-event | 415 | 0.66 | +0.07 |
| CarteaWang | other | normal | 12726 | 0.62 | -0.00 |
| CarteaWang | other | pre-event | 205 | 0.73 | -0.10 |
| Bergault_oracle | informed | normal | 3176 | 0.55 | +0.06 |
| Bergault_oracle | informed | pre-event | 0 | 0.00 | +0.00 |
| Bergault_oracle | real_money | normal | 8364 | 0.58 | +0.05 |
| Bergault_oracle | real_money | pre-event | 81 | 1.26 | -0.01 |
| Bergault_oracle | other | normal | 12742 | 0.62 | -0.02 |
| Bergault_oracle | other | pre-event | 253 | 0.61 | -0.02 |

## Real-day replay, factor model fitted without 2020-22: 2024-01-01..2026-12-31, 676 days, 186 event days (21 with measured USMPD jumps), gamma 2e-06

Events by type: AUCTION 111, ECB 22, FOMC 22, NFP 31.

| quoter | variant | P&L/day $k | CI | 5-min std $k | P&L/var x1e-3 | spread | inventory | hedge | hit ratios |
|---|---|---|---|---|---|---|---|---|---|
| Symmetric | text+calendar | 153.7 | [148.2, 157.0] | 20.1 | 0.38 | 185.4 | -9.2 | -22.4 | 0.35/0.35/0.35 |
| Bergault | text+calendar | 147.5 | [144.9, 150.0] | 12.1 | 1.01 | 166.5 | -5.0 | -14.1 | 0.34/0.32/0.30 |
| BarzykinCiceri_QAHR | text+calendar | 144.4 | [142.9, 146.4] | 12.0 | 1.01 | 163.6 | -4.8 | -14.4 | 0.25/0.26/0.27 |
| CarteaWang | text+calendar | 149.2 | [146.1, 152.5] | 12.0 | 1.04 | 165.3 | -2.3 | -13.9 | 0.33/0.32/0.30 |
| CarteaWang | placebo_text+calendar | 149.7 | [145.7, 154.2] | 12.0 | 1.04 | 165.2 | -1.7 | -13.8 | 0.33/0.32/0.30 |
| CarteaWang | calendar_only | 149.3 | [146.2, 153.3] | 12.0 | 1.04 | 165.4 | -2.3 | -13.9 | 0.33/0.32/0.30 |
| CarteaWang | none | 147.5 | [144.9, 150.0] | 12.1 | 1.01 | 166.5 | -5.0 | -14.1 | 0.34/0.32/0.30 |
| Bergault_oracle | text+calendar | 142.3 | [138.8, 144.4] | 11.5 | 1.08 | 154.0 | -0.0 | -11.6 | 0.21/0.28/0.30 |

| quoter | tier | window | fills | half-spread bp | 30m mark-out bp |
|---|---|---|---|---|---|
| Symmetric | informed | normal | 4843 | 0.45 | +0.10 |
| Symmetric | informed | pre-event | 282 | 0.45 | +2.14 |
| Symmetric | real_money | normal | 9630 | 0.54 | +0.01 |
| Symmetric | real_money | pre-event | 564 | 0.54 | -0.01 |
| Symmetric | other | normal | 14390 | 0.57 | -0.02 |
| Symmetric | other | pre-event | 294 | 0.57 | +0.34 |
| Bergault | informed | normal | 4809 | 0.39 | +0.09 |
| Bergault | informed | pre-event | 260 | 0.42 | +2.38 |
| Bergault | real_money | normal | 8930 | 0.52 | +0.05 |
| Bergault | real_money | pre-event | 526 | 0.53 | -0.18 |
| Bergault | other | normal | 12455 | 0.61 | +0.01 |
| Bergault | other | pre-event | 253 | 0.59 | +0.19 |
| BarzykinCiceri_QAHR | informed | normal | 3512 | 0.50 | +0.16 |
| BarzykinCiceri_QAHR | informed | pre-event | 212 | 0.50 | +2.66 |
| BarzykinCiceri_QAHR | real_money | normal | 7232 | 0.62 | +0.04 |
| BarzykinCiceri_QAHR | real_money | pre-event | 453 | 0.58 | +0.02 |
| BarzykinCiceri_QAHR | other | normal | 11073 | 0.67 | -0.02 |
| BarzykinCiceri_QAHR | other | pre-event | 229 | 0.66 | +0.17 |
| CarteaWang | informed | normal | 4806 | 0.39 | +0.08 |
| CarteaWang | informed | pre-event | 51 | 0.88 | +3.26 |
| CarteaWang | real_money | normal | 8957 | 0.52 | +0.08 |
| CarteaWang | real_money | pre-event | 399 | 0.65 | +0.03 |
| CarteaWang | other | normal | 12464 | 0.60 | +0.01 |
| CarteaWang | other | pre-event | 208 | 0.73 | +0.16 |
| Bergault_oracle | informed | normal | 3155 | 0.52 | +0.03 |
| Bergault_oracle | informed | pre-event | 0 | 0.00 | +0.00 |
| Bergault_oracle | real_money | normal | 8142 | 0.56 | +0.08 |
| Bergault_oracle | real_money | pre-event | 15 | 1.62 | +1.49 |
| Bergault_oracle | other | normal | 12447 | 0.61 | -0.00 |
| Bergault_oracle | other | pre-event | 258 | 0.60 | +0.08 |

## Tests

1831 / 1831 checks passed
