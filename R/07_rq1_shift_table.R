#!/usr/bin/env Rscript
# RQ1 helper: quantify how the location–scale model changes the story
# (effect size + between-site heterogeneity) relative to the homoskedastic baseline.
#
# Output: reports/rq1_shift_table.csv
#
# Notes:
# - Heteroskedastic (location–scale) summaries are read from reports/location_scale_<tag>_summary.csv
# - Homoskedastic summaries are computed from CmdStan CSV outputs under models/cmdstan_homo_<tag>/

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
})

source(file.path("R", "lib", "cli_utils.R"))
source(file.path("R", "lib", "file_utils.R"))

extract_from_summary_csv <- function(path, variables) {
  path <- resolve_report_path(path, reports_dir = dirname(path))
  if (!file.exists(path)) stop("Missing summary file: ", path)
  df <- readr::read_csv(path, show_col_types = FALSE)
  df <- df %>% mutate(variable = trimws(as.character(.data$variable)))
  variables <- trimws(variables)
  need <- c("variable", "mean", "sd")
  if (!all(need %in% names(df))) stop("Summary file missing columns ", paste(need, collapse = ", "), ": ", path)
  out <- df %>%
    filter(.data$variable %in% variables) %>%
    select(variable, mean, sd)
  miss <- setdiff(variables, out$variable)
  if (length(miss)) stop("Summary file missing variables: ", paste(miss, collapse = ", "), " in ", path)
  out
}

to_dot_vars <- function(vars) {
  gsub("\\]", "", gsub("\\[", ".", vars))
}

normalize_cmdstan_vars <- function(vars) {
  # Convert simple one-index dot notation to bracket form (e.g., beta.2 -> beta[2]).
  sub("^(.*)\\.([0-9]+)$", "\\1[\\2]", vars)
}

extract_from_cmdstan_csv <- function(tag, variables) {
  dir <- file.path("models", paste0("cmdstan_homo_", tag))
  if (!dir.exists(dir)) stop("Missing CmdStan homo output dir: ", dir)
  pat <- paste0("^location_scale_homo_", tag, "-[0-9]+\\.csv$")
  files <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!length(files)) stop("No CmdStan CSV files found in ", dir, " (pattern=", pat, ")")
  vars_dot <- unique(to_dot_vars(variables))
  x <- tryCatch(
    cmdstanr::read_cmdstan_csv(files, variables = variables),
    error = function(e) {
      cmdstanr::read_cmdstan_csv(files, variables = vars_dot)
    }
  )
  post <- x$post_warmup_draws
  if (is.null(post)) stop("No post-warmup draws loaded for ", tag, " from ", dir)
  posterior::summarise_draws(post, mean, sd) %>%
    select(variable, mean, sd) %>%
    mutate(variable = trimws(as.character(.data$variable))) %>%
    mutate(variable = normalize_cmdstan_vars(.data$variable)) %>%
    filter(.data$variable %in% variables)
}

stroop_tag <- trimws(parse_flag("stroop_tag", "stroop_ml3_site"))
psa001_attractive_tag <- trimws(parse_flag("psa001_attractive_tag", "psa001_attractive_gender_site_ad99"))
psa001_dominant_tag <- trimws(parse_flag("psa001_dominant_tag", "psa001_dominant_gender_site_ad99"))
out_path <- parse_flag("out", "reports/rq1_shift_table.csv")

vars <- c("beta[2]", "tau_site[2]")

datasets <- tibble(
  dataset = c("Stroop", "PSA001_Attractive", "PSA001_Dominant"),
  tag = c(stroop_tag, psa001_attractive_tag, psa001_dominant_tag)
)

results <- datasets %>%
  filter(!is.na(.data$tag), nzchar(.data$tag), .data$tag != "none") %>%
  mutate(
    hetero_summary_file = file.path("reports", paste0("location_scale_", tag, "_summary.csv"))
  ) %>%
  pmap_dfr(function(dataset, tag, hetero_summary_file) {
    debug <- tolower(Sys.getenv("RQ1_DEBUG", "false")) %in% c("1", "true", "yes", "y", "t")
    hetero <- extract_from_summary_csv(hetero_summary_file, vars) %>%
      mutate(model = "heterogeneous")
    homo_summary_file <- file.path("reports", paste0("location_scale_homo_", tag, "_summary.csv"))
    homo_summary_file <- resolve_report_path(homo_summary_file, reports_dir = dirname(homo_summary_file))
    homo <- if (file.exists(homo_summary_file)) {
      extract_from_summary_csv(homo_summary_file, vars) %>% mutate(model = "homogeneous")
    } else {
      extract_from_cmdstan_csv(tag, vars) %>% mutate(model = "homogeneous")
    }

    both <- bind_rows(hetero, homo)
    if (isTRUE(debug)) {
      message("RQ1_DEBUG tag=", tag, " homo_summary_file=", homo_summary_file, " exists=", file.exists(homo_summary_file))
      message("RQ1_DEBUG homo vars: ", paste(homo$variable, collapse = ", "))
      message("RQ1_DEBUG hetero vars: ", paste(hetero$variable, collapse = ", "))
    }

    pull_one <- function(model_name, variable_name, col) {
      if (isTRUE(debug)) {
        message("RQ1_DEBUG models: ", paste(unique(both$model), collapse = ", "))
        message("RQ1_DEBUG variables: ", paste(unique(both$variable), collapse = ", "))
      }
      rows <- both %>% filter(.data$model == model_name, .data$variable == variable_name)
      if (!nrow(rows)) stop("Missing ", variable_name, " for model=", model_name, " (tag=", tag, ")")
      rows[[col]][[1]]
    }

    beta_homo <- pull_one("homogeneous", "beta[2]", "mean")
    beta_hetero <- pull_one("heterogeneous", "beta[2]", "mean")
    tau_homo <- pull_one("homogeneous", "tau_site[2]", "mean")
    tau_hetero <- pull_one("heterogeneous", "tau_site[2]", "mean")

    tibble(
      dataset = dataset,
      tag = tag,
      beta2_homo_mean = beta_homo,
      beta2_homo_sd = pull_one("homogeneous", "beta[2]", "sd"),
      beta2_hetero_mean = beta_hetero,
      beta2_hetero_sd = pull_one("heterogeneous", "beta[2]", "sd"),
      tau_site2_homo_mean = tau_homo,
      tau_site2_homo_sd = pull_one("homogeneous", "tau_site[2]", "sd"),
      tau_site2_hetero_mean = tau_hetero,
      tau_site2_hetero_sd = pull_one("heterogeneous", "tau_site[2]", "sd"),
      delta_beta2 = beta_hetero - beta_homo,
      delta_tau_site2 = tau_hetero - tau_homo,
      tau_site2_reduction_pct = ifelse(is.finite(tau_homo) && tau_homo > 0, 100 * (tau_homo - tau_hetero) / tau_homo, NA_real_)
    )
  })

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(results, out_path)
message("Wrote ", out_path)
