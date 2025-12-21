# Draft Results and Interpretation (RQ1–RQ4)

## RQ2: What participant-slope distributions imply

| dataset | n | mean_slope | median_slope | p05 | p95 | responders_pct | opposite_pct | near_zero_pct |
|---|---|---|---|---|---|---|---|---|
| PSA001_Attractive | 279.000 | -0.588 | -0.537 | -1.411 | 0.156 | 73.835 | 4.301 | 0.000 |
| PSA001_Dominant | 306.000 | 0.041 | -0.013 | -0.666 | 0.829 | 23.856 | 22.549 | 0.000 |
| Stroop | 3337.000 | -0.079 | -0.080 | -0.107 | -0.051 | 98.472 | 0.000 | 0.000 |

Interpretation:
- **Stroop**: slopes are tightly negative; almost all participants are in the expected direction (robust effect, heterogeneity mostly magnitude).
- **PSA001 Attractive**: average effect negative with meaningful spread; most are in the expected direction, with a small opposite-direction minority.
- **PSA001 Dominant**: mean near zero with sign-mixing (responders and opposites are comparable), consistent with mixed individual trajectories.

## RQ1: Does modeling precision shift effects and site heterogeneity?
| dataset | beta2_homo_mean | beta2_hetero_mean | delta_beta2 | tau_site2_homo_mean | tau_site2_hetero_mean | delta_tau_site2 |
|---|---|---|---|---|---|---|
| Stroop | -0.085 | -0.080 | 0.005 | 0.010 | 0.010 | 0.000 |
| PSA001_Attractive | -0.602 | -0.594 | 0.008 | 0.072 | 0.077 | 0.005 |
| PSA001_Dominant | 0.070 | 0.075 | 0.005 | 0.165 | 0.196 | 0.031 |

Interpretation: effects are similar across homo vs hetero fits; site heterogeneity does **not** shrink under the heteroskedastic model (often slightly larger).

## RQ3: Are site differences explained by mix/precision?
| dataset | baseline_site_sd | residual_sd | variance_explained_pct | n_sites |
|---|---|---|---|---|
| PSA001_Attractive | 0.271 | 0.197 | 46.899 | 34.000 |
| PSA001_Dominant | 0.410 | 0.193 | 77.889 | 37.000 |
| Stroop | 0.008 | 0.005 | 70.102 | 22.000 |

Interpretation: a substantial portion of site-level variance is explained by composition/precision (strong for Stroop and Dominant; partial for Attractive).

## RQ4: Predictive comparison (site K-fold)
| dataset | comparison_method | comparison_unit | k_folds | elpd_diff | elpd_diff_se | delta_looic |
|---|---|---|---|---|---|---|
| Stroop | kfold | site | 5.000 | 3606.230 | 278.009 | 7212.459 |
| PSA001_Attractive | kfold | site | 5.000 | -352.306 | 15.519 | -704.612 |
| PSA001_Dominant | kfold | site | 5.000 | -434.005 | 17.674 | -868.009 |

Interpretation: heteroskedastic model is strongly favored for Stroop but disfavored for PSA001 tasks under site-held-out prediction; benefits are task-dependent.
