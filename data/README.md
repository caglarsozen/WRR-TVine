# Input data

The experiment uses daily adjusted closing prices for six exchange-traded funds and the closing level of the Cboe Volatility Index (VIX). The requested period is 1 January 2010 to 31 December 2025.

The default local input is `data/prices.csv` with one row per common trading date:

| Column | Exposure |
|---|---|
| `date` | Calendar date in `YYYY-MM-DD` format |
| `XLE` | Conventional energy |
| `ICLN` | Clean energy |
| `TAN` | Solar energy |
| `XLF` | Financial sector |
| `GLD` | Gold |
| `UUP` | US dollar |
| `VIX` | Implied volatility index level |

Each price must be numeric. Dates must be valid and unique. The engine sorts dates, retains the common observations with finite positive prices, records excluded rows, and computes log returns. It does not interpolate missing prices.

If the local CSV is absent, the run scripts request the six ETF series and `^VIX` through Yahoo Finance using `quantmod`; the resulting price file is retained for subsequent runs. To use a frozen input instead, place it at the default location or pass its path as the first argument to the run script.

The raw price series is not distributed in this repository. The authors retain the study's frozen aligned input. Its fingerprint is recorded in [environment/input_snapshot.sha256](../environment/input_snapshot.sha256). A fresh provider download may differ because prices or historical adjustments can change. Exact replication requires the matching input and computation environment; consult the authors about access to the study snapshot.
