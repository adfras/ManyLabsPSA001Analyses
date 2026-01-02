#!/usr/bin/env Rscript
# Quick per-site counts table to sanity-check what the data support.
#
# For each site:
# - participants
# - trials total
# - median/mean/min/max trials per participant
# Optionally also by a condition column (e.g., CONGRUENT or X_smile).
#
# Usage:
#   Rscript R/utils/00_counts_table.R --in data/processed/trials_stroop_ml3_with_site.csv --out reports/counts_stroop_by_site.csv
#   Rscript R/utils/00_counts_table.R --in data/processed/psa001_dominant_gender.csv --out reports/counts_psa001_dominant_by_site.csv --condition_col X_male

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path("R", "lib", "cli_utils.R"))

in_path <- parse_flag("in", NA)
out_path <- parse_flag("out", "reports/counts_by_site.csv")
out_condition_path <- parse_flag("out_condition", NA)
person_col <- parse_flag("person_col", "person")
site_col <- parse_flag("site_col", "site")
condition_col <- parse_flag("condition_col", NA)

if (is.na(in_path) || !nzchar(in_path)) stop("--in is required")
if (!file.exists(in_path)) stop("Input file not found: ", in_path)

df <- readr::read_csv(in_path, show_col_types = FALSE)
if (!all(c(person_col, site_col) %in% names(df))) {
  stop("Input must have columns ", person_col, " and ", site_col)
}

detect_condition_col <- function(df) {
  # Prefer within-person binary predictors when present.
  if ("X_congruent" %in% names(df)) return(if ("CONGRUENT" %in% names(df)) "CONGRUENT" else "X_congruent")
  if ("X_smile" %in% names(df)) return("X_smile")
  x_cols <- grep("^X_", names(df), value = TRUE)
  for (c in x_cols) {
    if (length(unique(df[[c]])) == 2) return(c)
  }
  if ("condition" %in% names(df)) return("condition")
  NULL
}

cond_col <- if (!is.na(condition_col) && nzchar(condition_col)) condition_col else detect_condition_col(df)
if (!is.null(cond_col) && !(cond_col %in% names(df))) stop("condition_col not found: ", cond_col)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

per_person <- df %>%
  count(.data[[site_col]], .data[[person_col]], name = "trials_per_person")

site_totals <- df %>%
  count(.data[[site_col]], name = "trials_total")

site_summary <- per_person %>%
  group_by(.data[[site_col]]) %>%
  summarise(
    participants = n(),
    median_trials_per_participant = median(trials_per_person),
    mean_trials_per_participant = mean(trials_per_person),
    min_trials_per_participant = min(trials_per_person),
    max_trials_per_participant = max(trials_per_person),
    .groups = "drop"
  ) %>%
  left_join(site_totals, by = site_col) %>%
  arrange(desc(participants), .data[[site_col]])

readr::write_csv(site_summary, out_path)
message("Wrote site counts to ", out_path)

if (!is.null(cond_col)) {
  cond_label <- function(x) {
    if (cond_col == "X_smile") return(ifelse(x > 0, "smile", "neutral"))
    if (cond_col == "X_congruent") return(ifelse(x > 0, "congruent", "incongruent"))
    x
  }
  df2 <- df %>%
    mutate(.cond = cond_label(.data[[cond_col]]))

  per_person_cond <- df2 %>%
    count(.data[[site_col]], .cond, .data[[person_col]], name = "trials_per_person")
  totals_cond <- df2 %>%
    count(.data[[site_col]], .cond, name = "trials_total")

  site_cond_summary <- per_person_cond %>%
    group_by(.data[[site_col]], .cond) %>%
    summarise(
      participants = n(),
      median_trials_per_participant = median(trials_per_person),
      mean_trials_per_participant = mean(trials_per_person),
      min_trials_per_participant = min(trials_per_person),
      max_trials_per_participant = max(trials_per_person),
      .groups = "drop"
    ) %>%
    left_join(totals_cond, by = c(site_col, ".cond")) %>%
    rename(condition = .cond) %>%
    arrange(.data[[site_col]], condition)

  if (is.na(out_condition_path) || !nzchar(out_condition_path)) {
    out_condition_path <- sub("\\.csv$", "_by_site_condition.csv", out_path)
  }
  dir.create(dirname(out_condition_path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(site_cond_summary, out_condition_path)
  message("Wrote site-by-condition counts to ", out_condition_path, " (condition_col=", cond_col, ")")
} else {
  message("No condition column detected; wrote site-only table.")
}
