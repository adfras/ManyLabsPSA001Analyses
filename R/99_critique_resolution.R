#!/usr/bin/env Rscript
# Resolve and document core critiques for the participant-slope results.
# Produces:
# - results_checks.md
# - results_draft.md
# - reports/participant_slope_outliers.csv
# - reports/figs/hist_*_p_in_direction.png

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(file.path("R", "lib", "cli_utils.R"))

md_table <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return("(no data)")
  fmt <- function(x) {
    if (is.numeric(x)) return(sprintf(paste0("%.", digits, "f"), x))
    as.character(x)
  }
  df2 <- as.data.frame(lapply(df, fmt), stringsAsFactors = FALSE)
  header <- paste0("| ", paste(names(df2), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df2)), collapse = "|"), "|")
  rows <- apply(df2, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

# Inputs
prevalence_detail <- parse_flag("detail", "reports/person_prevalence_detail.csv")
prevalence_summary <- parse_flag("summary", "reports/person_prevalence_summary.csv")
rq1_shift <- parse_flag("rq1", "reports/rq1_shift_table.csv")
rq3_path <- parse_flag("rq3", "reports/site_variance_decomposition.csv")
site_kfold_path <- parse_flag("site_kfold", "reports/site_model_comparisons.csv")
site_mix <- parse_flag("site_mix", "reports/site_level_mix_precision.csv")
rope <- parse_num(parse_flag("rope", "0.05"), 0.05)
p_dir_thresh <- parse_num(parse_flag("p_dir_thresh", "0.90"), 0.90)
p_rope_thresh <- parse_num(parse_flag("p_rope_thresh", "0.90"), 0.90)
refit_if_bad <- parse_bool(parse_flag("refit_if_bad", "false"), FALSE)
rhat_threshold <- parse_num(parse_flag("rhat_threshold", "1.05"), 1.05)

refit_data <- parse_flag("refit_data", "data/processed/psa001_dominant_gender.csv")
refit_tag <- parse_flag("refit_tag", "psa001_dominant_gender_site_ad99")
refit_chains <- parse_int(parse_flag("refit_chains", "4"), 4)
refit_iter <- parse_int(parse_flag("refit_iter", "2000"), 2000)
refit_adapt <- parse_flag("refit_adapt_delta", "0.995")
refit_treedepth <- parse_flag("refit_max_treedepth", "20")
refit_init <- parse_flag("refit_init", "0.1")

out_md <- parse_flag("out_md", "results_checks.md")
out_draft <- parse_flag("out_draft", "results_draft.md")
out_outliers <- parse_flag("out_outliers", "reports/participant_slope_outliers.csv")
fig_dir <- parse_flag("fig_dir", "reports/figs")

# Helper: check R-hat on core params
rhat_ok <- function(path, threshold = 1.05) {
  if (!file.exists(path)) return(FALSE)
  df <- read_csv(path, show_col_types = FALSE)
  core <- df %>% filter(grepl("^beta\\[|^tau_site|^tau_person|^sigma|^tau_sigma|^lp__", variable))
  if (!nrow(core)) return(FALSE)
  rh <- core$rhat
  rh <- rh[is.finite(rh)]
  if (!length(rh)) return(FALSE)
  max_rhat <- max(rh, na.rm = TRUE)
  is.finite(max_rhat) && max_rhat <= threshold
}

# Optional refit if diagnostics are bad
homo_summary <- file.path("reports", paste0("location_scale_homo_", refit_tag, "_summary.csv"))
if (refit_if_bad && !rhat_ok(homo_summary, rhat_threshold)) {
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  system2(rscript, args = c(
    "R/04_fit_stroop_location_scale.R",
    refit_data,
    "--tag", refit_tag,
    "--homo_only", "true", "--compare_homo", "true", "--loo", "false",
    "--chains", as.character(refit_chains), "--parallel_chains", as.character(refit_chains),
    "--iter", as.character(refit_iter), "--init", refit_init,
    "--adapt_delta", refit_adapt, "--max_treedepth", refit_treedepth,
    "--refresh", "1"
  ))
}

# Load data
if (!file.exists(prevalence_detail) || !file.exists(prevalence_summary) || !file.exists(rq1_shift)) {
  stop("Missing required inputs. Expected: ", prevalence_detail, ", ", prevalence_summary, ", ", rq1_shift)
}

detail <- read_csv(prevalence_detail, show_col_types = FALSE)
summary <- read_csv(prevalence_summary, show_col_types = FALSE)
rq1 <- read_csv(rq1_shift, show_col_types = FALSE)
rq3 <- if (file.exists(rq3_path)) read_csv(rq3_path, show_col_types = FALSE) else NULL
site_kfold <- if (file.exists(site_kfold_path)) read_csv(site_kfold_path, show_col_types = FALSE) else NULL

# 1) Direction coding check
sign_check <- rq1 %>%
  select(dataset, beta2_hetero_mean) %>%
  left_join(summary %>% select(dataset, effect_direction), by = "dataset") %>%
  mutate(sign_match = sign(beta2_hetero_mean) == sign(effect_direction))

# 2) Uncertainty beyond posterior means
pdir <- detail %>%
  group_by(dataset) %>%
  summarise(
    p_in_mean = mean(p_in_direction),
    p_in_median = median(p_in_direction),
    p_in_gt_0_5 = mean(p_in_direction > 0.5),
    p_in_gt_0_9 = mean(p_in_direction >= p_dir_thresh),
    p_op_gt_0_9 = mean(p_opposite_direction >= p_dir_thresh),
    p_near_mean = mean(p_near),
    .groups = "drop"
  )

# 3) Near-zero sensitivity (ROPE)
pnear <- detail %>%
  mutate(sd_eff = pmax(beta_sd, 1e-6)) %>%
  mutate(p_near_rope = pnorm(rope, beta_mean, sd_eff) - pnorm(-rope, beta_mean, sd_eff)) %>%
  group_by(dataset) %>%
  summarise(
    rope = rope,
    p_near_mean = mean(p_near_rope),
    near_zero_pct = mean(p_near_rope >= p_rope_thresh),
    .groups = "drop"
  )

# 4) Outliers (top 1% absolute slopes)
outliers <- detail %>%
  group_by(dataset) %>%
  mutate(abs_beta = abs(beta_mean), thr = quantile(abs_beta, 0.99)) %>%
  filter(abs_beta >= thr) %>%
  select(dataset, person, site, beta_mean, beta_sd, p_in_direction)
write_csv(outliers, out_outliers)

# 5) Site-size sensitivity (variance decomposition)
site_sens <- NULL
if (file.exists(site_mix)) {
  smp <- read_csv(site_mix, show_col_types = FALSE)
  site_sens <- smp %>%
    group_by(dataset) %>%
    group_modify(~{
      bind_rows(lapply(c(5, 10), function(min_n) {
        s <- .x %>% filter(n >= min_n)
        baseline_sd <- sd(s$mean_beta)
        residual_sd <- sd(s$residual)
        tibble(
          min_n = min_n,
          n_sites = nrow(s),
          baseline_sd = baseline_sd,
          residual_sd = residual_sd,
          variance_explained = 1 - (residual_sd^2 / baseline_sd^2)
        )
      }))
    })
}

# 6) Stimulus coverage (faces)
stim_cov <- NULL
if (file.exists("data/processed/psa001_attractive_gender.csv")) {
  d <- read_csv("data/processed/psa001_attractive_gender.csv", show_col_types = FALSE)
  stim_counts <- d %>% distinct(person, stim_id) %>% group_by(person) %>% summarise(n = n(), .groups = "drop")
  stim_cov <- bind_rows(stim_cov, tibble(
    dataset = "PSA001_Attractive",
    n_stimuli = n_distinct(d$stim_id),
    median_stim_per_person = median(stim_counts$n),
    min_stim_per_person = min(stim_counts$n),
    max_stim_per_person = max(stim_counts$n)
  ))
}
if (file.exists("data/processed/psa001_dominant_gender.csv")) {
  d <- read_csv("data/processed/psa001_dominant_gender.csv", show_col_types = FALSE)
  stim_counts <- d %>% distinct(person, stim_id) %>% group_by(person) %>% summarise(n = n(), .groups = "drop")
  stim_cov <- bind_rows(stim_cov, tibble(
    dataset = "PSA001_Dominant",
    n_stimuli = n_distinct(d$stim_id),
    median_stim_per_person = median(stim_counts$n),
    min_stim_per_person = min(stim_counts$n),
    max_stim_per_person = max(stim_counts$n)
  ))
}

# 7) Diagnostics snapshot (core params)
summary_paths <- list.files("reports", pattern = "^location_scale_.*_summary\\.csv$", full.names = TRUE)
core_diag <- bind_rows(lapply(summary_paths, function(p) {
  df <- read_csv(p, show_col_types = FALSE)
  if (!"rhat" %in% names(df)) return(NULL)
  core <- df %>% filter(grepl("^beta\\[|^tau_site|^tau_person|^sigma|^tau_sigma|^lp__", variable))
  if (!nrow(core)) return(NULL)
  rh <- core$rhat
  rh <- rh[is.finite(rh)]
  if (!length(rh)) return(NULL)
  tibble(
    summary_file = basename(p),
    max_rhat = max(rh, na.rm = TRUE),
    n_gt_1_01 = sum(rh > 1.01, na.rm = TRUE),
    n_gt_1_05 = sum(rh > 1.05, na.rm = TRUE),
    n_params = length(rh)
  )
}))

# 8) p(direction) histograms
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
detail %>%
  group_by(dataset) %>%
  group_walk(~{
    slug <- tolower(gsub("[^a-z0-9]+", "_", .y$dataset))
    p <- ggplot(.x, aes(p_in_direction)) +
      geom_histogram(bins = 40, fill = "#2b6cb0", color = "white", alpha = 0.85) +
      labs(title = paste0(.y$dataset, ": p(in direction)"), x = "p(in direction)", y = "Count") +
      theme_minimal(base_size = 12)
    ggsave(file.path(fig_dir, paste0("hist_", slug, "_p_in_direction.png")), p, width = 6, height = 4, dpi = 300)
  })

# Additional RQ2 summary table (for results draft)
rq2_stats <- detail %>%
  group_by(dataset) %>%
  summarise(
    n = n(),
    mean_slope = mean(beta_mean),
    median_slope = median(beta_mean),
    p05 = quantile(beta_mean, 0.05),
    p95 = quantile(beta_mean, 0.95),
    .groups = "drop"
  ) %>%
  left_join(
    summary %>% select(dataset, responders_pct, opposite_pct, near_zero_pct),
    by = "dataset"
  )

# Write critique-resolution markdown report
lines <- c(
  "# Critique resolution report",
  "",
  "## 1) Direction coding check",
  paste0("- ", sign_check$dataset, ": sign_match=", tolower(sign_check$sign_match)),
  "",
  "## 2) Uncertainty beyond posterior means",
  md_table(pdir, digits = 3),
  "",
  "## 3) Near-zero sensitivity (ROPE)",
  md_table(pnear, digits = 3),
  "",
  "## 4) Outliers file",
  paste0("- ", out_outliers),
  "",
  "## 5) Site-size sensitivity (variance decomposition)",
  if (is.null(site_sens)) "(site mix file not found)" else md_table(site_sens, digits = 3),
  "",
  "## 6) Stimulus coverage",
  if (is.null(stim_cov)) "(PSA001 processed files not found)" else md_table(stim_cov, digits = 0),
  "",
  "## 7) Diagnostics (core parameters)",
  md_table(core_diag, digits = 3)
)

writeLines(lines, out_md)
cat("Wrote:", out_md, "\n")
cat("Wrote:", out_outliers, "\n")

# Write results draft summary (RQ1–RQ3 + supplemental holdouts)
draft <- c(
  "# Draft Results and Interpretation (RQ1–RQ3)",
  "",
  "## RQ1/RQ2: Task comparison and direction prevalence (participant slopes)",
  "",
  md_table(
    rq2_stats %>% mutate(
      responders_pct = responders_pct * 100,
      opposite_pct = opposite_pct * 100,
      near_zero_pct = near_zero_pct * 100
    ) %>%
      select(dataset, n, mean_slope, median_slope, p05, p95, responders_pct, opposite_pct, near_zero_pct),
    digits = 3
  ),
  "",
  "Interpretation:",
  {
    ds_present <- unique(rq2_stats$dataset)
    bullets <- character(0)
    if ("Stroop" %in% ds_present) {
      bullets <- c(bullets, "- **Stroop**: slopes are tightly negative; almost all participants are in the expected direction (robust effect, heterogeneity mostly magnitude).")
    }
    if ("PSA001_Attractive" %in% ds_present) {
      bullets <- c(bullets, "- **PSA001 Attractive**: average effect negative with meaningful spread; most are in the expected direction, with a small opposite-direction minority.")
    }
    if ("PSA001_Dominant" %in% ds_present) {
      bullets <- c(bullets, "- **PSA001 Dominant**: mean near zero with sign-mixing (responders and opposites are comparable), consistent with mixed individual trajectories.")
    }
    if (!length(bullets)) bullets <- "- (No datasets found in prevalence detail.)"
    bullets
  },
  "",
  "## Supplemental: Homoskedastic vs heteroskedastic shift",
  if (is.null(rq1)) "(rq1_shift_table.csv not found)" else md_table(
    rq1 %>% select(dataset, beta2_homo_mean, beta2_hetero_mean, delta_beta2, tau_site2_homo_mean, tau_site2_hetero_mean, delta_tau_site2),
    digits = 3
  ),
  "",
  "Interpretation: effects are similar across homo vs hetero fits; site heterogeneity does **not** shrink under the heteroskedastic model (often slightly larger).",
  "",
  "## Supplemental: Site variance decomposition (mix/precision)",
  if (is.null(rq3)) "(site_variance_decomposition.csv not found)" else md_table(
    rq3 %>% select(dataset, baseline_site_sd, residual_sd, variance_explained_pct, n_sites),
    digits = 3
  ),
  "",
  {
    if (is.null(rq3) || !nrow(rq3)) {
      "Interpretation: (no site variance decomposition results found)."
    } else {
      ds <- rq3$dataset
      msg <- "Interpretation: a substantial portion of site-level variance is explained by composition/precision"
      if ("PSA001_Attractive" %in% ds) {
        paste0(msg, " (partial for Attractive; stronger for Stroop/Dominant).")
      } else {
        paste0(msg, ".")
      }
    }
  },
  "",
  "## Supplemental: Site K-fold predictive comparison",
  if (is.null(site_kfold)) "(site_model_comparisons.csv not found)" else {
    cols <- c("dataset", "comparison_method", "comparison_unit", "k_folds", "elpd_diff", "elpd_diff_se", "delta_looic")
    md_table(site_kfold %>% select(any_of(cols)), digits = 3)
  },
  "",
  "Interpretation: heteroskedastic model is strongly favored for Stroop but disfavored for PSA001 tasks under site-held-out prediction; benefits are task-dependent."
)

writeLines(draft, out_draft)
cat("Wrote:", out_draft, "\n")
