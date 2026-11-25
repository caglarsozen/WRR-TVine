# Aggregate reference results

These CSV files are derived from the completed empirical experiment. They retain the observed numerical results rather than substituting results from a smoke or synthetic run. The two universes use the same six-ETF risk target.

| File | Content |
|---|---|
| `score_summary.csv` | Model-level density, scenario, event and portfolio losses |
| `validation_selection.csv` | Selected rolling settings and validation scores |
| `paired_block20_summary.csv` | Paired loss differences, pointwise block intervals and adjusted tests |
| `pit_calibration_checked.csv` | Marginal hit frequencies, KS and residual-dependence diagnostics |
| `regime_logscore_contrasts.csv` | Direct-minus-matched mean log score by state |
| `exploratory_regime_inference.csv` | Exploratory state-conditional contrasts and uncertainty |
| `tail_grid_summary.csv` | Alternative finite-threshold joint-event scores |
| `state_event_forecasts.csv` | State-conditioned event and density summaries |
| `forecast_refit_state_alignment.csv` | Daily state versus latest scheduled refit state, aggregated |
| `training_cohit_directed.csv` | Historical conditional co-hit counts and frequencies |
| `joint6_vix_upper_event.csv` | Joint fund distress and upper-decile VIX-change probabilities |
| `uncertainty_summary.csv` | Local bootstrap and numerical simulation summaries |
| `local_logscore_uncertainty.csv` | Local log-density and probability uncertainty |
| `tail_tilt_mass.csv` | Historical weight-mass and effective-sample-size diagnostics at refits |

The fitted models, raw prices, daily PIT series and daily forecasts are not included here. The R engine produces those files locally. The read-only audit and supplementary analysis scripts recreate these reporting summaries from a completed results directory. Use `analysis/compare_reference.py` to compare the regenerated tables with these files within an explicit numerical tolerance.

The reference results are from a retrospective evaluation of the historical period. They do not demonstrate universal superiority of relevance weighting. Local estimation intervals condition on selected marginal models and settings; the nine local assessment origins in each universe do not include a stress date.
