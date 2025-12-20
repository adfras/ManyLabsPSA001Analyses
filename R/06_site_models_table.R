#!/usr/bin/env Rscript
# Build compact hetero vs homo comparison tables for the main datasets.
# Default paths assume working directory is the project root (e.g., D:/ManyLabsAnalyses).
#
# This script intentionally separates:
# - parameter summaries (from full-data fits), and
# - predictive comparison (from PSIS-LOO on a subsample or site K-fold CV).

suppressPackageStartupMessages({
  library(tidyverse)
})

# Simple key=value arg parser (e.g., --out=reports/site_model_comparisons.csv)
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

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop("None of these paths exist: ", paste(paths, collapse = ", "))
  hit[[1]]
}

pick_path <- function(paths, label = "", hint = NULL) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) return(hit[[1]])
  msg <- paste0("Missing ", label, " file. Tried: ", paste(paths, collapse = ", "))
  if (!is.null(hint)) msg <- paste0(msg, "\n\n", hint)
  stop(msg)
}

pull_stat <- function(df, vars, col = "mean", path = "") {
  # vars can be a single name or a vector of fallbacks.
  if (length(vars) == 1) vars <- c(vars)
  hit <- vars[vars %in% df$variable]
  if (!length(hit)) stop("None of vars {", paste(vars, collapse = ", "), "} found in ", path)
  row <- df[df$variable == hit[[1]], col, drop = TRUE]
  row[[1]]
}

extract_effects <- function(summary_path, effect_vars, tau_vars) {
  df <- readr::read_csv(summary_path, show_col_types = FALSE)
  tibble(
    effect_mean = pull_stat(df, effect_vars, "mean", summary_path),
    effect_sd = pull_stat(df, effect_vars, "sd", summary_path),
    tau_site_mean = pull_stat(df, tau_vars, "mean", summary_path),
    tau_site_sd = pull_stat(df, tau_vars, "sd", summary_path)
  )
}

extract_fit <- function(bf_path) {
  bf <- readr::read_csv(bf_path, show_col_types = FALSE)
  hetero <- bf %>% filter(model == "heterogeneous") %>% slice_head(n = 1)
  homo   <- bf %>% filter(model == "homogeneous") %>% slice_head(n = 1)
  if (!nrow(hetero) || !nrow(homo)) stop("Expected hetero and homo rows in ", bf_path)
  log_bf <- hetero$elpd - homo$elpd
  bf_num <- exp(pmin(log_bf, 700)) # cap to avoid overflow while keeping Inf when truly huge

  method <- if ("method" %in% names(bf)) hetero$method[[1]] else "psis_loo"
  unit <- if ("unit" %in% names(bf)) hetero$unit[[1]] else "obs"
  k_folds <- if ("k_folds" %in% names(bf)) hetero$k_folds[[1]] else NA_integer_
  elpd_diff_se <- if ("delta_elpd_se" %in% names(bf)) hetero$delta_elpd_se[[1]] else NA_real_
  tibble(
    comparison_method = method,
    comparison_unit = unit,
    k_folds = k_folds,
    elpd_hetero = hetero$elpd,
    elpd_homo = homo$elpd,
    elpd_se_hetero = hetero$elpd_se,
    elpd_se_homo = homo$elpd_se,
    looic_hetero = hetero$looic,
    looic_homo = homo$looic,
    elpd_diff = hetero$elpd - homo$elpd,
    elpd_diff_se = elpd_diff_se,
    bf_hetero_vs_homo = bf_num,
    log_bf = log_bf
  )
}

build_row <- function(dataset, summary_path, bf_path, effect_vars, tau_vars) {
  effects <- extract_effects(summary_path, effect_vars, tau_vars)
  fit <- extract_fit(bf_path)
  bind_cols(tibble(dataset = dataset, summary_file = summary_path, bf_file = bf_path), effects, fit)
}

defaults <- list(
  out = "reports/site_model_comparisons.csv",
  stroop_tag = "stroop_ml3_site",
  stroop_cv_tag = "stroop_ml3_site_sub30_mm",
  psa001_attractive_tag = "none",
  psa001_attractive_cv_tag = "",
  psa001_dominant_tag = "psa001_dominant_gender_site_ad99",
  psa001_dominant_cv_tag = ""
)

opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
stroop_tag <- opts$stroop_tag
stroop_cv_tag <- opts$stroop_cv_tag
psa001_attractive_tag <- opts$psa001_attractive_tag
psa001_attractive_cv_tag <- opts$psa001_attractive_cv_tag
psa001_dominant_tag <- opts$psa001_dominant_tag
psa001_dominant_cv_tag <- opts$psa001_dominant_cv_tag

stroop_summary <- pick_path(
  c(
    file.path("reports", paste0("location_scale_", stroop_tag, "_summary.csv")),
    file.path("stroop_ml3_bundle", "reports", paste0("location_scale_", stroop_tag, "_summary.csv"))
  ),
  label = "Stroop summary",
  hint = "Run: Rscript R/04_fit_stroop_location_scale.R data/processed/trials_stroop_ml3_with_site.csv --tag stroop_ml3_site --loo false --compare_homo false"
)
stroop_bf <- pick_path(
  c(
    file.path("reports", "kfold", paste0("location_scale_", stroop_cv_tag, "_kfold_site_bf.csv")),
    file.path("reports", paste0("location_scale_", stroop_cv_tag, "_bf.csv")),
    file.path("stroop_ml3_bundle", "reports", paste0("location_scale_", stroop_cv_tag, "_bf.csv"))
  ),
  label = "Stroop predictive comparison (hetero vs homo)",
  hint = paste0(
    "Preferred (site K-fold, on Stroop subsample):\n",
    "  Rscript R/11_kfold_location_scale.R data/processed/trials_stroop_ml3_with_site_sub30.csv \\\n",
    "    --tag ", stroop_cv_tag, " --unit site --k_site 10 --chains 2 --iter 600 --init 0 --adapt_delta 0.99 --max_treedepth 15\n\n",
    "Fallback (PSIS-LOO on Stroop subsample):\n",
    "  Rscript R/04_fit_stroop_location_scale.R data/processed/trials_stroop_ml3_with_site_sub30.csv \\\n",
    "    --tag ", stroop_cv_tag, " --loo true --compare_homo true --moment_match true --init 0 --adapt_delta 0.95"
  )
)

psa001_attractive_summary <- if (nzchar(psa001_attractive_tag) && psa001_attractive_tag != "none") {
  pick_path(
    c(
      file.path("reports", paste0("location_scale_", psa001_attractive_tag, "_summary.csv")),
      file.path("reports", "location_scale_psa001_attractive_gender_summary.csv")
    ),
    label = "PSA001 attractive summary",
    hint = paste0(
      "Run: Rscript R/04_fit_stroop_location_scale.R data/processed/psa001_attractive_gender.csv ",
      "--tag ", psa001_attractive_tag, " --loo false --compare_homo false --init 0 --adapt_delta 0.99 --max_treedepth 15"
    )
  )
} else NA_character_
psa001_attractive_bf <- if (!is.na(psa001_attractive_summary)) {
  if (!nzchar(psa001_attractive_cv_tag)) psa001_attractive_cv_tag <- psa001_attractive_tag
  pick_path(
    c(
      file.path("reports", "kfold", paste0("location_scale_", psa001_attractive_cv_tag, "_kfold_site_bf.csv")),
      file.path("reports", paste0("location_scale_", psa001_attractive_cv_tag, "_bf.csv"))
    ),
    label = "PSA001 attractive predictive comparison (hetero vs homo)",
    hint = paste0(
      "Run site K-fold CV:\n",
      "  Rscript R/11_kfold_location_scale.R data/processed/psa001_attractive_gender.csv \\\n",
      "    --tag ", psa001_attractive_tag, " --unit site --k_site 10 --chains 2 --iter 600 --init 0 --adapt_delta 0.99 --max_treedepth 15"
    )
  )
} else NA_character_

psa001_dominant_summary <- if (nzchar(psa001_dominant_tag) && psa001_dominant_tag != "none") {
  pick_path(
    c(
      file.path("reports", paste0("location_scale_", psa001_dominant_tag, "_summary.csv")),
      file.path("reports", "location_scale_psa001_dominant_gender_summary.csv")
    ),
    label = "PSA001 dominant summary",
    hint = paste0(
      "Run: Rscript R/04_fit_stroop_location_scale.R data/processed/psa001_dominant_gender.csv ",
      "--tag ", psa001_dominant_tag, " --loo false --compare_homo false --init 0 --adapt_delta 0.99 --max_treedepth 15"
    )
  )
} else NA_character_
psa001_dominant_bf <- if (!is.na(psa001_dominant_summary)) {
  if (!nzchar(psa001_dominant_cv_tag)) psa001_dominant_cv_tag <- psa001_dominant_tag
  pick_path(
    c(
      file.path("reports", "kfold", paste0("location_scale_", psa001_dominant_cv_tag, "_kfold_site_bf.csv")),
      file.path("reports", paste0("location_scale_", psa001_dominant_cv_tag, "_bf.csv"))
    ),
    label = "PSA001 dominant predictive comparison (hetero vs homo)",
    hint = paste0(
      "Run site K-fold CV:\n",
      "  Rscript R/11_kfold_location_scale.R data/processed/psa001_dominant_gender.csv \\\n",
      "    --tag ", psa001_dominant_tag, " --unit site --k_site 10 --chains 2 --iter 600 --init 0 --adapt_delta 0.99 --max_treedepth 15"
    )
  )
} else NA_character_

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)

rows <- list()
rows[[length(rows) + 1]] <- build_row("Stroop", stroop_summary, stroop_bf, "beta[2]", "tau_site[2]")
if (!is.na(psa001_attractive_summary)) {
  rows[[length(rows) + 1]] <- build_row("PSA001_Attractive", psa001_attractive_summary, psa001_attractive_bf, "beta[2]", "tau_site[2]")
}
if (!is.na(psa001_dominant_summary)) {
  rows[[length(rows) + 1]] <- build_row("PSA001_Dominant", psa001_dominant_summary, psa001_dominant_bf, "beta[2]", "tau_site[2]")
}

out_tbl <- bind_rows(rows) %>%
  mutate(
    delta_elpd = elpd_hetero - elpd_homo,
    delta_looic = looic_homo - looic_hetero
  )

readr::write_csv(out_tbl, opts$out)
message("Wrote ", opts$out)
