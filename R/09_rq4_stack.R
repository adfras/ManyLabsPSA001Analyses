#!/usr/bin/env Rscript
# Stack dataset-level comparison results (RQ4).

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

defaults <- list(
  site_model_file = "reports/site_model_comparisons.csv",
  site_var_file = "reports/site_variance_decomposition.csv",
  out = "reports/rq4_comparison_stack.csv"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop("None of these paths exist: ", paste(paths, collapse = ", "))
  hit[[1]]
}

site_model_path <- first_existing(c(opts$site_model_file, "stroop_ml3_bundle/reports/site_model_comparisons.csv"))
site_var_path   <- first_existing(c(opts$site_var_file, "stroop_ml3_bundle/reports/site_variance_decomposition.csv"))

site_mod <- readr::read_csv(site_model_path, show_col_types = FALSE)
site_var <- readr::read_csv(site_var_path, show_col_types = FALSE)

out_tbl <- site_mod %>%
  left_join(site_var %>% select(dataset, residual_sd, baseline_site_sd, variance_explained_pct, r2), by = "dataset") %>%
  mutate(
    evidence = case_when(
      bf_hetero_vs_homo == Inf ~ ">100",
      bf_hetero_vs_homo >= 30 ~ "30-100",
      bf_hetero_vs_homo >= 10 ~ "10-30",
      bf_hetero_vs_homo >= 3 ~ "3-10",
      TRUE ~ "<3"
    )
  )

readr::write_csv(out_tbl, opts$out)
message("Wrote ", opts$out)
