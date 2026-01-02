#!/usr/bin/env Rscript
# Decompose site-level variation using participant mix/precision summaries.
# Prefers the person_prevalence_detail.csv produced by summarise_participants.R;
# will fall back to raw participant tables if that file is absent.

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path("R", "lib", "cli_utils.R"))
source(file.path("R", "lib", "file_utils.R"))
source(file.path("R", "lib", "prevalence_utils.R"))

calc_probs <- function(df, dataset, mean_col, sd_col, sigma_col = NULL, eps = 0.02) {
  if (!mean_col %in% names(df) || !sd_col %in% names(df)) {
    stop("Columns ", mean_col, " and ", sd_col, " must exist for ", dataset)
  }
  sd_eff <- pmax(df[[sd_col]], 1e-6)
  sigma_use <- if (!is.null(sigma_col) && sigma_col %in% names(df)) dplyr::coalesce(df[[sigma_col]], sd_eff) else sd_eff
  direction <- detect_direction(df[[mean_col]], dataset)
  tibble(
    dataset = dataset,
    person = df[["person"]],
    site = if ("site" %in% names(df)) df[["site"]] else NA_character_,
    beta_mean = df[[mean_col]],
    beta_sd = sd_eff,
    sigma = sigma_use
  ) %>%
    mutate(
      p_gt0 = 1 - pnorm(0, beta_mean, beta_sd),
      p_in_direction = if (direction > 0) p_gt0 else 1 - p_gt0,
      effect_direction = direction,
      responder = p_in_direction > 0.9,
      high_precision = sigma < median(sigma, na.rm = TRUE)
    )
}

summarise_site <- function(df) {
  df %>%
    group_by(dataset, site) %>%
    summarise(
      n = n(),
      mean_beta = mean(beta_mean, na.rm = TRUE),
      sd_beta = sd(beta_mean, na.rm = TRUE),
      mean_sigma = mean(sigma, na.rm = TRUE),
      frac_high_precision = mean(high_precision, na.rm = TRUE),
      frac_responder = mean(responder, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(mean_beta))
}

fit_decomposition <- function(site_df) {
  baseline_sd <- sd(site_df$mean_beta, na.rm = TRUE)
  if (!is.finite(baseline_sd) || baseline_sd == 0 || nrow(site_df) < 3) {
    site_df$predicted <- NA_real_
    site_df$residual <- NA_real_
    return(list(
      summary = tibble(baseline_site_sd = baseline_sd, residual_sd = NA_real_, variance_explained_pct = NA_real_, r2 = NA_real_, n_sites = nrow(site_df)),
      site = site_df
    ))
  }
  mod <- lm(mean_beta ~ frac_high_precision + frac_responder + sd_beta + mean_sigma, data = site_df)
  # Use newdata to get predictions aligned to the full input (lm() drops rows with NA).
  pred <- as.numeric(predict(mod, newdata = site_df))
  site_df <- site_df %>% mutate(predicted = pred, residual = mean_beta - predicted)
  resid_sd <- sd(site_df$residual, na.rm = TRUE)
  var_expl <- 1 - (resid_sd^2 / baseline_sd^2)
  r2 <- summary(mod)$r.squared
  list(
    summary = tibble(baseline_site_sd = baseline_sd, residual_sd = resid_sd, variance_explained_pct = 100 * var_expl, r2 = r2, n_sites = nrow(site_df)),
    site = site_df
  )
}

defaults <- list(
  detail_file = "reports/person_prevalence_detail.csv",
  stroop_file = "reports/location_scale_stroop_ml3_site_participants.csv",
  stroop_alt = "reports/location_scale_stroop_ml3_site_sub_participants.csv",
  stroop_trials = "data/processed/trials_stroop_ml3_with_site.csv",
  psa001_attractive_file = "reports/location_scale_psa001_attractive_gender_site_ad99_participants.csv",
  psa001_attractive_alt = "reports/location_scale_psa001_attractive_gender_participants.csv",
  psa001_attractive_trials = "data/processed/psa001_attractive_gender.csv",
  psa001_dominant_file = "reports/location_scale_psa001_dominant_gender_site_ad99_participants.csv",
  psa001_dominant_alt = "reports/location_scale_psa001_dominant_gender_participants.csv",
  psa001_dominant_trials = "data/processed/psa001_dominant_gender.csv",
  include_psa001_attractive = "false",
  include_psa001_dominant = "true",
  eps = "0.02",
  out_summary = "reports/site_variance_decomposition.csv",
  out_site = "reports/site_level_mix_precision.csv",
  covars_file = ""
)

opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

eps <- as.numeric(opts$eps)
include_psa001_attractive <- parse_bool(opts$include_psa001_attractive, default = FALSE)
include_psa001_dominant <- parse_bool(opts$include_psa001_dominant, default = TRUE)

stroop_site_map <- load_person_site_map(opts$stroop_trials)
psa001_attractive_site_map <- load_person_site_map(opts$psa001_attractive_trials)
psa001_dominant_site_map <- load_person_site_map(opts$psa001_dominant_trials)

stroop_file <- resolve_report_path(opts$stroop_file)
stroop_alt <- resolve_report_path(opts$stroop_alt)
psa001_attractive_file <- resolve_report_path(opts$psa001_attractive_file)
psa001_attractive_alt <- resolve_report_path(opts$psa001_attractive_alt)
psa001_dominant_file <- resolve_report_path(opts$psa001_dominant_file)
psa001_dominant_alt <- resolve_report_path(opts$psa001_dominant_alt)

if (file.exists(opts$detail_file)) {
  detail <- readr::read_csv(opts$detail_file, show_col_types = FALSE)
  if (!"site" %in% names(detail)) detail$site <- NA_character_
  keep_ds <- c("Stroop")
  if (include_psa001_attractive) keep_ds <- c(keep_ds, "PSA001_Attractive")
  if (include_psa001_dominant) keep_ds <- c(keep_ds, "PSA001_Dominant")
  detail <- detail %>% filter(.data$dataset %in% keep_ds)
  if (any(detail$dataset == "Stroop")) {
    detail <- bind_rows(
      detail %>% filter(dataset != "Stroop"),
      detail %>% filter(dataset == "Stroop") %>% attach_site(stroop_site_map)
    )
  }
  if (include_psa001_attractive && any(detail$dataset == "PSA001_Attractive")) {
    detail <- bind_rows(
      detail %>% filter(dataset != "PSA001_Attractive"),
      detail %>% filter(dataset == "PSA001_Attractive") %>% attach_site(psa001_attractive_site_map)
    )
  }
  if (include_psa001_dominant && any(detail$dataset == "PSA001_Dominant")) {
    detail <- bind_rows(
      detail %>% filter(dataset != "PSA001_Dominant"),
      detail %>% filter(dataset == "PSA001_Dominant") %>% attach_site(psa001_dominant_site_map)
    )
  }
} else {
  rows <- list()
  if (file.exists(stroop_file)) {
    df_raw <- readr::read_csv(stroop_file, show_col_types = FALSE)
    df_raw <- attach_site(df_raw, stroop_site_map)
    rows[[length(rows) + 1]] <- calc_probs(df_raw, "Stroop", "beta_mean_X_congruent", "beta_sd_X_congruent", "sigma_mean", eps)
  } else if (file.exists(stroop_alt)) {
    df_raw <- readr::read_csv(stroop_alt, show_col_types = FALSE)
    df_raw <- attach_site(df_raw, stroop_site_map)
    rows[[length(rows) + 1]] <- calc_probs(df_raw, "Stroop", "beta_mean_X_congruent", "beta_sd_X_congruent", "sigma_mean", eps)
  }
  if (include_psa001_attractive) {
    psa001_attractive_part_path <- dplyr::coalesce(
      if (file.exists(psa001_attractive_file)) psa001_attractive_file else NA_character_,
      if (file.exists(psa001_attractive_alt)) psa001_attractive_alt else NA_character_
    )
    if (!is.na(psa001_attractive_part_path)) {
      df_raw <- readr::read_csv(psa001_attractive_part_path, show_col_types = FALSE)
      df_raw <- attach_site(df_raw, psa001_attractive_site_map)
      rows[[length(rows) + 1]] <- calc_probs(df_raw, "PSA001_Attractive", "beta_mean_X_male", "beta_sd_X_male", "sigma_mean", eps)
    }
  }
  if (include_psa001_dominant) {
    psa001_dominant_part_path <- dplyr::coalesce(
      if (file.exists(psa001_dominant_file)) psa001_dominant_file else NA_character_,
      if (file.exists(psa001_dominant_alt)) psa001_dominant_alt else NA_character_
    )
    if (!is.na(psa001_dominant_part_path)) {
      df_raw <- readr::read_csv(psa001_dominant_part_path, show_col_types = FALSE)
      df_raw <- attach_site(df_raw, psa001_dominant_site_map)
      rows[[length(rows) + 1]] <- calc_probs(df_raw, "PSA001_Dominant", "beta_mean_X_male", "beta_sd_X_male", "sigma_mean", eps)
    }
  }
  if (!length(rows)) stop("No participant data found for site decomposition.")
  detail <- bind_rows(rows)
}

# Optional covariates file: must contain a 'site' column; otherwise ignored.
if (nzchar(opts$covars_file) && file.exists(opts$covars_file)) {
  covars <- readr::read_csv(opts$covars_file, show_col_types = FALSE)
  if (!"site" %in% names(covars)) {
    warning("Covariates file lacks a 'site' column; skipping: ", opts$covars_file)
  } else {
    cov_keep <- covars %>% select(where(~ !is.list(.)))
    detail <- detail %>% left_join(cov_keep, by = "site")
  }
}

site_tbl <- summarise_site(detail)

site_results <- site_tbl %>% group_split(dataset) %>% lapply(function(df) {
  res <- fit_decomposition(df)
  list(
    summary = res$summary %>% mutate(dataset = df$dataset[1]),
    site = res$site %>% mutate(dataset = df$dataset[1])
  )
})

summary_out <- bind_rows(lapply(site_results, function(x) x$summary)) %>%
  select(dataset, everything())
site_out <- bind_rows(lapply(site_results, function(x) x$site)) %>%
  select(dataset, site, n, mean_beta, sd_beta, mean_sigma, frac_high_precision, frac_responder, predicted, residual)

dir.create(dirname(opts$out_summary), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(summary_out, opts$out_summary)
readr::write_csv(site_out, opts$out_site)

message("Wrote ", opts$out_summary, " and ", opts$out_site)
