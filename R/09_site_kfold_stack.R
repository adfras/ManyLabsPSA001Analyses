#!/usr/bin/env Rscript
# Stack dataset-level comparison results (site K-fold holdout).

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path("R", "lib", "cli_utils.R"))
source(file.path("R", "lib", "file_utils.R"))

defaults <- list(
  site_model_file = "reports/site_model_comparisons.csv",
  site_var_file = "reports/site_variance_decomposition.csv",
  out = "reports/site_kfold_comparison_stack.csv"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

site_model_path <- first_existing(c(opts$site_model_file, "stroop_ml3_bundle/reports/site_model_comparisons.csv"))
site_var_path   <- first_existing(c(opts$site_var_file, "stroop_ml3_bundle/reports/site_variance_decomposition.csv"))
if (is.na(site_model_path) || !nzchar(site_model_path)) stop("Missing site_model_file: ", opts$site_model_file)
if (is.na(site_var_path) || !nzchar(site_var_path)) stop("Missing site_var_file: ", opts$site_var_file)

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
