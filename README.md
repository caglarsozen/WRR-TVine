# WRR-TVine

Research and replication code for **Forecasting Energy–Finance Tail Risk with Weighted Rolling Regime-Adaptive Vine Copulas**, by Mervenur Sözen and Çağlar Sözen.

The weighted rolling regime-adaptive tail-vine copula (WRR-TVine) combines chronological marginal modelling, observable volatility states, historical relevance weights and regular vine dependence. The experiment examines whether relevance weighting adds predictive value beyond matched rolling estimators. It also distinguishes weighted-sample marginal distortion, estimation uncertainty and scenario simulation error.

## Research design

The common risk target is a simultaneous lower-decile standardized shortfall in at least three of six funds. The principal dependence universe contains XLE, ICLN, TAN, XLF, GLD and UUP. The extended universe additionally models the log change in the Cboe Volatility Index (VIX). Both use lagged VIX levels to classify volatility states, and both retain the same six-fund event and equal-weight portfolio.

| Setting | Reference experiment |
|---|---|
| Requested price period | 1 January 2010 – 31 December 2025 |
| Aligned returns | 4,023 common return vectors |
| After initial marginal history | 3,273 observations |
| Training / validation / test | 2,270 / 501 / 502 (through 2021 / 2022–2023 / 2024–2025) |
| Marginal candidates | ARMA(1,1) with standard, GJR or exponential GARCH and skewed Student-t innovations |
| Rolling candidates | Windows of 500, 750 or 1,000 observations; refit intervals of 10 or 20; two family pools |
| Volatility states | 60th and 85th training percentiles of lagged VIX |
| Relevance weights | Recency decay 2, same-state factor 3, VIX proximity decay 1.25, historical tail bonus 2 |
| Effective sample-size floor | 15% of the rolling window |
| Weighted resampling | 750 rows per component; mixtures of 1, 5 and 20 fitted copulas |
| Daily predictive scenarios | 10,000 |
| Local estimation uncertainty | 199 circular-block refits at each of nine origins per universe |
| Paired loss uncertainty | 1,999 circular-block resamples, block lengths 10, 20 and 40 |
| Robustness | Component removal, family pools, alternative state cutoffs, evaluation-event thresholds, empirical margins and training-only clipping |

Model and rolling settings are selected chronologically. Forecasts use observed histories, and the dependence estimator is updated on its selected refit schedule. The historical test assessment is retrospective; the results are not a claim of universal forecasting superiority.

## Repository structure

| Path | Purpose |
|---|---|
| `R/wrr_tvine.R` | Complete estimation, simulation, inference, checkpointing and export engine |
| `00_INSTALL_PACKAGES.R` | Install missing R dependencies |
| `01_RUN_SMOKE.R` | Short execution check with reduced computational settings |
| `02_RUN_FULL.R` | Full empirical experiment |
| `03_PACKAGE_RESULTS.R` | Package an already completed full run |
| `04_VERIFY_RESULTS.R` | Independent base R checks and supplementary summaries |
| `R/verify_results.R` | Read-only verification functions |
| `tests/` | Numerical invariants and synthetic interface/integration tests |
| `analysis/` | Audit, supplementary diagnostics, table/figure generation and reference comparison |
| `reference_results/` | Observed aggregate numerical reference outputs |
| `data/README.md` | Input-data specification and retrieval |
| `environment/` | Recorded computation environment and frozen-input fingerprint |
| `CITATION.cff` | Citation metadata |

## Quick start

Clone or download the repository and open its root directory in a terminal or RStudio. An R installation and the listed dependencies are required for model estimation.

~~~bash
git clone https://github.com/caglarsozen/WRR-TVine.git
cd WRR-TVine
Rscript 00_INSTALL_PACKAGES.R
Rscript 01_RUN_SMOKE.R
~~~

In RStudio, the equivalent commands from the repository root are:

~~~r
source("00_INSTALL_PACKAGES.R")
source("01_RUN_SMOKE.R")
~~~

The smoke run uses a shortened evaluation horizon and reduced simulation and bootstrap settings. It is a software check, not an empirical result for the study.

### Full empirical run

~~~bash
Rscript 02_RUN_FULL.R
~~~

Or in RStudio:

~~~r
source("02_RUN_FULL.R")
~~~

The full experiment is computationally intensive. It runs both dependence universes, the complete comparison set, sensitivity exercises and local uncertainty analysis. The reference calculation required approximately 235 hours of recorded elapsed computation on its original hardware, including recovery and robustness calculations.

If a local `data/prices.csv` is available, it is reused. Otherwise the engine requests the daily prices through Yahoo Finance and retains the downloaded input. A different frozen input may be supplied explicitly:

~~~bash
Rscript 02_RUN_FULL.R "/path/to/prices.csv"
~~~

See [Data](data/README.md) and [Computational environment](environment/README.md) for the exact schema, reference environment and implications of provider revisions.

### Output, progress and restarting

- Smoke results: `results/smoke/`.
- Full results: `results/full/`.
- Full-run console: `FULL_RUN_CONSOLE.log`; latest progress: `results/full/PROGRESS.txt`.
- Successful full completion: a `COMPLETED full` record in `results/full/RUN_STATUS.txt`.
- The completed full run automatically creates `WRR_TVINE_FULL_RESULTS.zip`. To package existing completed results again without estimation, use `Rscript 03_PACKAGE_RESULTS.R`.

Running the same configuration again reuses compatible saved marginal, copula and daily-forecast checkpoints. A stable local hashing directory under the output cache avoids reliance on a temporary folder that may be removed during a long run. An interrupted run can still require recomputation of local uncertainty calculations. Do not delete the cache when resuming. Only if the previous R process has stopped and a stale `RUN_LOCK` remains, remove that lock directory before restarting. Do not run two writers against the same output directory.

### Read-only verification and reporting

The following commands operate on completed outputs; they do not re-estimate the models:

~~~bash
Rscript 04_VERIFY_RESULTS.R results/full analysis_output/verification
~~~

The base R verification also accepts a completed ZIP, detecting its results directory from the completion record. The separate Python audit reconstructs the reported losses, diagnostics, weight identities and inferential summaries from the archived predictions:

~~~bash
python -m pip install -r analysis/requirements.txt
python analysis/audit_results.py results/full analysis_output/derived
python analysis/extended_diagnostics.py results/full analysis_output/derived
python analysis/make_paper_assets.py results/full analysis_output/derived analysis_output/paper
python analysis/compare_reference.py analysis_output/derived
~~~

The independent audit checks the specific empirical sample and design described above. A different price history can change the sample or results. The audit also handles undefined zero-variance marginal hit tests as missing significance values, retaining the multiplicity family size. Raw diagnostics and the documented reporting diagnostics are therefore distinguished.

Additional synthetic tests exercise weighted vine estimation, marginal filtering and the end-to-end pipeline after the model dependencies are installed:

~~~bash
Rscript tests/test_vine_interface.R
Rscript tests/test_garch.R
Rscript tests/test_integration.R
~~~

## Reference results and interpretation

The [reference result tables](reference_results/README.md) contain the aggregate numerical outputs of the completed research experiment, including rolling selections, paired loss inference, marginal calibration, state-conditional differences and uncertainty diagnostics. Larger full-density log scores are preferable within a fixed dimension and marginal specification. Lower losses are preferable for the scenario, probability and portfolio scores.

The direct weighted estimator does not outperform the matched rolling estimator on mean aggregate log score in either universe. Removing the tail bonus improves common variogram loss, and the 20-component resampling mixture improves log scores relative to direct weighting. These differences survive the reported multiplicity adjustment in both universes. The conditional results depend on the marginal model and the risk target. They should not be interpreted as validation of a trading or allocation policy.

## Citation and license

Please acknowledge the authors and cite the associated research when using this software. Machine-readable citation information is supplied in [CITATION.cff](CITATION.cff).

The research software is distributed under the [MIT License](LICENSE). The raw third-party market-price data is not redistributed; see the input-data instructions.
