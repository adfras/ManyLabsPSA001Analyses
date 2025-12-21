# Critique resolution report

## 1) Direction coding check
- Stroop: sign_match=true
- PSA001_Attractive: sign_match=true
- PSA001_Dominant: sign_match=true

## 2) Uncertainty beyond posterior means
| dataset | p_in_mean | p_in_median | p_in_gt_0_5 | p_in_gt_0_9 | p_op_gt_0_9 | p_near_mean |
|---|---|---|---|---|---|---|
| PSA001_Attractive | 0.868 | 0.998 | 0.882 | 0.738 | 0.043 | 0.018 |
| PSA001_Dominant | 0.500 | 0.477 | 0.484 | 0.239 | 0.225 | 0.038 |
| Stroop | 0.984 | 0.992 | 1.000 | 0.985 | 0.000 | 0.051 |

## 3) Near-zero sensitivity (ROPE)
| dataset | rope | p_near_mean | near_zero_pct |
|---|---|---|---|
| PSA001_Attractive | 0.050 | 0.045 | 0.000 |
| PSA001_Dominant | 0.050 | 0.094 | 0.003 |
| Stroop | 0.050 | 0.214 | 0.048 |

## 4) Outliers file
- reports/participant_slope_outliers.csv

## 5) Site-size sensitivity (variance decomposition)
| dataset | min_n | n_sites | baseline_sd | residual_sd | variance_explained |
|---|---|---|---|---|---|
| PSA001_Attractive | 5.000 | 16.000 | 0.159 | 0.131 | 0.319 |
| PSA001_Attractive | 10.000 | 7.000 | 0.160 | 0.059 | 0.863 |
| PSA001_Dominant | 5.000 | 13.000 | 0.167 | 0.125 | 0.442 |
| PSA001_Dominant | 10.000 | 8.000 | 0.072 | 0.139 | -2.729 |
| Stroop | 5.000 | 22.000 | 0.008 | 0.005 | 0.701 |
| Stroop | 10.000 | 22.000 | 0.008 | 0.005 | 0.701 |

## 6) Stimulus coverage
| dataset | n_stimuli | median_stim_per_person | min_stim_per_person | max_stim_per_person |
|---|---|---|---|---|
| PSA001_Attractive | 120 | 120 | 120 | 120 |
| PSA001_Dominant | 120 | 120 | 120 | 120 |

## 7) Diagnostics (core parameters)
| summary_file | max_rhat | n_gt_1_01 | n_gt_1_05 | n_params |
|---|---|---|---|---|
| location_scale_homo_psa001_attractive_gender_site_ad99_summary.csv | 1.031 | 5.000 | 0.000 | 12.000 |
| location_scale_homo_psa001_dominant_gender_site_ad99_summary.csv | 1.015 | 2.000 | 0.000 | 12.000 |
| location_scale_homo_stroop_ml3_site_summary.csv | 1.010 | 0.000 | 0.000 | 12.000 |
| location_scale_psa001_attractive_gender_site_ad99_summary.csv | 1.031 | 8.000 | 0.000 | 292.000 |
| location_scale_psa001_attractive_gender_summary.csv | 1.053 | 8.000 | 2.000 | 292.000 |
| location_scale_psa001_dominant_gender_site_ad99_summary.csv | 1.041 | 9.000 | 0.000 | 319.000 |
| location_scale_stroop_ml3_site_sub10_summary.csv | 1.028 | 6.000 | 0.000 | 233.000 |
| location_scale_stroop_ml3_site_sub30_mm_summary.csv | 1.026 | 8.000 | 0.000 | 673.000 |
| location_scale_stroop_ml3_site_sub_summary.csv | 1.021 | 10.000 | 0.000 | 453.000 |
| location_scale_stroop_ml3_site_summary.csv | 1.026 | 4.000 | 0.000 | 3350.000 |
