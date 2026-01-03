# Results snapshot (2026-01-03)

This folder is a lightweight snapshot of the main analysis outputs for RQ1–RQ3.

Notes:
- **Held-out trial prediction uses plug-in posterior means** (not full posterior predictive integration). See `R/12_trial_holdout_predict.R` and `R/14_holdout_eval_homo.R`.
- **Stroop homoskedastic holdout** is evaluated via `R/14_holdout_eval_homo.R` and saved as `holdout_homo_stroop_ml3_holdout_homo_summary.csv`.

## RQ1: shift table (homo vs location–scale)

| dataset | beta2_homo | beta2_hetero | delta_beta2 | tau_site2_homo | tau_site2_hetero | delta_tau_site2 |
|---|---|---|---|---|---|---|
| Stroop | -0.085 | -0.080 | 0.005 | 0.010 | 0.010 | 0.000 |
| PSA001_Attractive | -0.602 | -0.594 | 0.008 | 0.072 | 0.077 | 0.005 |
| PSA001_Dominant | 0.063 | 0.075 | 0.012 | 0.164 | 0.196 | 0.032 |

Source: `rq1_shift_table.csv`.

## RQ2: prevalence categories (percent of participants)

| dataset | n | responders | near_zero | opposite | direction_uncertain | responders_pct | near_zero_pct | opposite_pct | direction_uncertain_pct |
|---|---|---|---|---|---|---|---|---|---|
| PSA001_Attractive | 279 | 206 | 0 | 12 | 61 | 73.8 | 0.0 | 4.3 | 21.9 |
| PSA001_Dominant | 306 | 73 | 0 | 69 | 164 | 23.9 | 0.0 | 22.5 | 53.6 |
| Stroop | 3337 | 3286 | 0 | 0 | 51 | 98.5 | 0.0 | 0.0 | 1.5 |

Source: `person_prevalence_summary.csv` (details in `person_prevalence_detail.csv`).

## RQ3: held-out trial prediction (plug-in posterior means)

| dataset | n_holdout | mean_loglik_hetero | mean_loglik_homo | delta_mean_loglik | rmse_hetero | rmse_homo | mae_hetero | mae_homo | corr_hetero | corr_homo |
|---|---|---|---|---|---|---|---|---|---|---|
| PSA001_Attractive | 10044 | -1.605 | -2.092 | 0.488 | 1.328 | 1.771 | 1.009 | 1.458 | 0.673 | 0.167 |
| PSA001_Dominant | 11016 | -1.641 | -2.003 | 0.361 | 1.393 | 1.700 | 1.087 | 1.387 | 0.573 | 0.011 |
| Stroop | 63395 | -0.277 | -0.498 | 0.221 | 0.341 | 0.387 | 0.229 | 0.269 | 0.482 | 0.103 |

Sources:
- Location–scale (hetero): `holdout_psa001_attractive_holdout_summary.csv`, `holdout_psa001_dominant_holdout_summary.csv`, `holdout_stroop_ml3_holdout_summary.csv`
- Homoskedastic baseline (homo): `holdout_homo_psa001_attractive_holdout_homo_summary.csv`, `holdout_homo_psa001_dominant_holdout_homo_summary.csv`, `holdout_homo_stroop_ml3_holdout_homo_summary.csv`

