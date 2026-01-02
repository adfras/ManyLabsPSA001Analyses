#!/usr/bin/env Rscript
# One-shot holdout prediction pipeline (make holdout -> fit -> evaluate).

suppressPackageStartupMessages({
  library(dplyr)
})

source(file.path("R", "lib", "cli_utils.R"))

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript.exe")
}

holdout_frac <- parse_num(parse_flag("holdout_frac", "0.30"), 0.30)
chains <- as.integer(parse_flag("chains", 2))
iter <- as.integer(parse_flag("iter", 800))
parallel_chains <- as.integer(parse_flag("parallel_chains", chains))
adapt_delta <- parse_num(parse_flag("adapt_delta", "0.99"), 0.99)
max_treedepth <- as.integer(parse_flag("max_treedepth", 15))
refresh <- as.integer(parse_flag("refresh", 1))

include_stroop <- parse_bool(parse_flag("include_stroop", "false"), FALSE)
include_attractive <- parse_bool(parse_flag("include_attractive", "true"), TRUE)
include_dominant <- parse_bool(parse_flag("include_dominant", "true"), TRUE)

run_one <- function(tag, in_path, effect_direction) {
  holdout_path <- file.path("data", "processed", paste0(tag, "_holdout.csv"))

  # 1) Make holdout file
  system2(rscript, args = c(
    "R/12_trial_holdout_predict.R",
    "--in", in_path,
    "--tag", tag,
    "--holdout_frac", holdout_frac,
    "--holdout_out", holdout_path,
    "--eval", "false",
    "--make_holdout", "true"
  ))

  # 2) Fit model on w==1
  system2(rscript, args = c(
    "R/05_fit_location_scale.R",
    holdout_path,
    "--tag", tag,
    "--loo", "false",
    "--compare_homo", "false",
    "--chains", chains,
    "--parallel_chains", parallel_chains,
    "--iter", iter,
    "--init", "0",
    "--adapt_delta", adapt_delta,
    "--max_treedepth", max_treedepth,
    "--refresh", refresh
  ))

  # 3) Evaluate holdout predictions
  system2(rscript, args = c(
    "R/12_trial_holdout_predict.R",
    "--in", holdout_path,
    "--tag", tag,
    "--make_holdout", "false",
    "--eval", "true",
    "--effect_direction", effect_direction
  ))
}

if (include_attractive) {
  run_one(
    tag = "psa001_attractive_holdout",
    in_path = "data/processed/psa001_attractive_gender.csv",
    effect_direction = -1
  )
}

if (include_dominant) {
  run_one(
    tag = "psa001_dominant_holdout",
    in_path = "data/processed/psa001_dominant_gender.csv",
    effect_direction = 1
  )
}

if (include_stroop) {
  run_one(
    tag = "stroop_ml3_holdout",
    in_path = "data/processed/trials_stroop_ml3_with_site.csv",
    effect_direction = -1
  )
}

message("Holdout pipeline complete.")
