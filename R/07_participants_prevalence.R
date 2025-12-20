#!/usr/bin/env Rscript
# Compute participant-level probabilities and prevalence summaries for the main datasets.
# Uses saved participant summary tables; no model refits required.

suppressPackageStartupMessages({
  library(tidyverse)
})

parse_args <- function(args, defaults) {
  out <- defaults
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    key <- sub("^--", "", kv[1])
    out[[key]] <- kv[2]
  }
  out
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "yes", "y", "t")
}

load_person_site_map <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    warning("Site map input file not found: ", path)
    return(tibble(person = character(), site = character()))
  }
  # Read only the mapping columns (works for any processed CSV with person + site columns).
  m <- readr::read_csv(
    path,
    show_col_types = FALSE,
    col_types = readr::cols_only(
      person = readr::col_character(),
      site = readr::col_character()
    )
  ) %>%
    transmute(person = as.character(person), site = as.character(site)) %>%
    distinct()
  # Defensive: if a person appears in multiple sites (shouldn't), keep the first.
  dup <- m %>% count(person, name = "n_sites") %>% filter(n_sites > 1)
  if (nrow(dup)) {
    warning("Multiple sites per person in ", path, "; taking the first site per person.")
    m <- m %>% group_by(person) %>% slice_head(n = 1) %>% ungroup()
  }
  m
}

attach_site <- function(df, site_map) {
  if ("site" %in% names(df) && any(!is.na(df$site))) return(df)
  if (!nrow(site_map)) {
    if (!"site" %in% names(df)) df$site <- NA_character_
    return(df)
  }
  df %>%
    mutate(person = as.character(.data$person)) %>%
    select(-any_of("site")) %>%
    left_join(site_map, by = "person")
}

detect_direction <- function(beta_mean, dataset, min_abs_mean = 1e-4) {
  m <- mean(beta_mean, na.rm = TRUE)
  if (!is.finite(m)) {
    warning("Non-finite mean beta for ", dataset, "; using direction=+1.")
    return(1L)
  }
  if (abs(m) < min_abs_mean) {
    warning("Mean beta near 0 for ", dataset, " (mean=", signif(m, 3), "); using direction=+1.")
    return(1L)
  }
  if (m >= 0) 1L else -1L
}

calc_probs <- function(df, dataset, mean_col, sd_col, sigma_col = NULL, eps = 0.02) {
  if (!nrow(df)) stop("Empty participant table for ", dataset)
  if (!mean_col %in% names(df) || !sd_col %in% names(df)) {
    stop("Columns ", mean_col, " and ", sd_col, " must exist for ", dataset)
  }
  sd_eff <- pmax(df[[sd_col]], 1e-6)
  sigma_use <- if (!is.null(sigma_col) && sigma_col %in% names(df)) dplyr::coalesce(df[[sigma_col]], sd_eff) else sd_eff
  direction <- detect_direction(df[[mean_col]], dataset)
  tibble(
    dataset = dataset,
    person = as.character(df[["person"]]),
    site = if ("site" %in% names(df)) as.character(df[["site"]]) else NA_character_,
    beta_mean = df[[mean_col]],
    beta_sd = sd_eff,
    sigma = sigma_use,
    p_outlier = if ("p_outlier_valence" %in% names(df)) df[["p_outlier_valence"]] else NA_real_
  ) %>%
    mutate(
      p_gt0 = 1 - pnorm(0, beta_mean, beta_sd),
      p_lt0 = pnorm(0, beta_mean, beta_sd),
      p_near = pnorm(eps, beta_mean, beta_sd) - pnorm(-eps, beta_mean, beta_sd),
      p_in_direction = if (direction > 0) p_gt0 else p_lt0,
      p_opposite_direction = if (direction > 0) p_lt0 else p_gt0,
      effect_direction = direction,
      responder = p_in_direction > 0.9,
      near_zero = p_near > 0.5,
      opposite = p_opposite_direction > 0.9,
      high_precision = sigma < median(sigma, na.rm = TRUE)
    )
}

summarise_prevalence <- function(df) {
  df %>%
    summarise(
      n = n(),
      responders = sum(responder, na.rm = TRUE),
      near_zero = sum(near_zero, na.rm = TRUE),
      opposite = sum(opposite, na.rm = TRUE),
      high_precision = sum(high_precision, na.rm = TRUE),
      responders_pct = responders / n,
      near_zero_pct = near_zero / n,
      opposite_pct = opposite / n,
      high_precision_pct = high_precision / n,
      effect_direction = dplyr::coalesce(first(effect_direction), NA_integer_),
      mean_sigma = mean(sigma, na.rm = TRUE),
      mean_beta = mean(beta_mean, na.rm = TRUE)
    )
}

make_hist <- function(df, dataset, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  p <- ggplot(df, aes(beta_mean)) +
    geom_histogram(bins = 40, fill = "#2b6cb0", color = "white", alpha = 0.85) +
    labs(title = paste0(dataset, " participant slopes"), x = "Posterior mean slope", y = "Count") +
    theme_minimal(base_size = 12)
  out_path <- file.path(out_dir, paste0("hist_", tolower(dataset), "_beta_mean.png"))
  ggsave(out_path, p, width = 6, height = 4, dpi = 300)
  out_path
}

defaults <- list(
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
  out_detail = "reports/person_prevalence_detail.csv",
  out_summary = "reports/person_prevalence_summary.csv",
  fig_dir = "reports/figs"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

eps <- as.numeric(opts$eps)
include_psa001_attractive <- parse_bool(opts$include_psa001_attractive, default = FALSE)
include_psa001_dominant <- parse_bool(opts$include_psa001_dominant, default = TRUE)
rows <- list()

stroop_site_map <- load_person_site_map(opts$stroop_trials)
psa001_attractive_site_map <- load_person_site_map(opts$psa001_attractive_trials)
psa001_dominant_site_map <- load_person_site_map(opts$psa001_dominant_trials)

if (file.exists(opts$stroop_file)) {
  df_raw <- readr::read_csv(opts$stroop_file, show_col_types = FALSE)
  df_raw <- attach_site(df_raw, stroop_site_map)
  rows[[length(rows) + 1]] <- calc_probs(df_raw, "Stroop", "beta_mean_X_congruent", "beta_sd_X_congruent", "sigma_mean", eps)
} else if (file.exists(opts$stroop_alt)) {
  df_raw <- readr::read_csv(opts$stroop_alt, show_col_types = FALSE)
  df_raw <- attach_site(df_raw, stroop_site_map)
  rows[[length(rows) + 1]] <- calc_probs(df_raw, "Stroop", "beta_mean_X_congruent", "beta_sd_X_congruent", "sigma_mean", eps)
} else {
  warning("Missing Stroop participant file: ", opts$stroop_file, " or alt: ", opts$stroop_alt)
}

if (include_psa001_attractive) {
  psa001_attractive_part_path <- dplyr::coalesce(
    if (file.exists(opts$psa001_attractive_file)) opts$psa001_attractive_file else NA_character_,
    if (file.exists(opts$psa001_attractive_alt)) opts$psa001_attractive_alt else NA_character_
  )

  if (!is.na(psa001_attractive_part_path)) {
    df_raw <- readr::read_csv(psa001_attractive_part_path, show_col_types = FALSE)
    df_raw <- attach_site(df_raw, psa001_attractive_site_map)
    rows[[length(rows) + 1]] <- calc_probs(df_raw, "PSA001_Attractive", "beta_mean_X_male", "beta_sd_X_male", "sigma_mean", eps)
  } else {
    warning("Missing PSA001 attractive participant file. Expected one of: ", opts$psa001_attractive_file, ", ", opts$psa001_attractive_alt)
  }
}

if (include_psa001_dominant) {
  psa001_dominant_part_path <- dplyr::coalesce(
    if (file.exists(opts$psa001_dominant_file)) opts$psa001_dominant_file else NA_character_,
    if (file.exists(opts$psa001_dominant_alt)) opts$psa001_dominant_alt else NA_character_
  )

  if (!is.na(psa001_dominant_part_path)) {
    df_raw <- readr::read_csv(psa001_dominant_part_path, show_col_types = FALSE)
    df_raw <- attach_site(df_raw, psa001_dominant_site_map)
    rows[[length(rows) + 1]] <- calc_probs(df_raw, "PSA001_Dominant", "beta_mean_X_male", "beta_sd_X_male", "sigma_mean", eps)
  } else {
    warning("Missing PSA001 dominant participant file. Expected one of: ", opts$psa001_dominant_file, ", ", opts$psa001_dominant_alt)
  }
}

if (!length(rows)) stop("No participant tables found.")

dir.create(dirname(opts$out_detail), recursive = TRUE, showWarnings = FALSE)

detail <- bind_rows(rows)
readr::write_csv(detail, opts$out_detail)

summary_tbl <- detail %>% group_by(dataset) %>% summarise_prevalence()
readr::write_csv(summary_tbl, opts$out_summary)

hist_paths <- purrr::imap(detail %>% split(.$dataset), ~ make_hist(.x, .y, out_dir = opts$fig_dir))

message("Wrote ", opts$out_detail, "; ", opts$out_summary)
message("Figures: ", paste(unlist(hist_paths), collapse = ", "))
