#!/usr/bin/env Rscript
# Update only missing outputs for revised RQs and holdout comparisons.
# This script tries to rebuild lightweight outputs from existing fits/CmdStan CSVs.
# It will only refit models if explicitly allowed via --allow_refit true and/or
# --allow_kfold true.

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path("R", "lib", "cli_utils.R"))

run_rscript <- function(args, desc) {
  cat("\n==>", desc, "\n")
  status <- system2(rscript, args = args)
  if (status != 0) stop("Step failed: ", desc, " (exit status ", status, ")")
}

has_cmdstan_csv <- function(dir, pattern) {
  dir.exists(dir) && length(list.files(dir, pattern = pattern, full.names = TRUE)) > 0
}

rscript <- file.path(R.home("bin"), "Rscript")

# Tags + data paths
stroop_tag <- parse_flag("stroop_tag", "stroop_ml3_site")
stroop_data <- parse_flag("stroop_data", "data/processed/trials_stroop_ml3_with_site.csv")
stroop_cv_tag <- parse_flag("stroop_cv_tag", "stroop_ml3_site_sub30")
stroop_sub_data <- parse_flag("stroop_sub_data", "data/processed/trials_stroop_ml3_with_site_sub30.csv")

psa_attr_tag <- parse_flag("psa001_attractive_tag", "psa001_attractive_gender_site_ad99")
psa_attr_data <- parse_flag("psa001_attractive_data", "data/processed/psa001_attractive_gender.csv")
psa_attr_cv_tag <- parse_flag("psa001_attractive_cv_tag", "psa001_attractive_gender_sub30")
psa_attr_sub_data <- parse_flag("psa001_attractive_sub_data", "data/processed/psa001_attractive_gender_sub30.csv")

psa_dom_tag <- parse_flag("psa001_dominant_tag", "psa001_dominant_gender_site_ad99")
psa_dom_data <- parse_flag("psa001_dominant_data", "data/processed/psa001_dominant_gender.csv")
psa_dom_cv_tag <- parse_flag("psa001_dominant_cv_tag", "psa001_dominant_gender_sub30")
psa_dom_sub_data <- parse_flag("psa001_dominant_sub_data", "data/processed/psa001_dominant_gender_sub30.csv")

# Controls
allow_refit <- parse_bool(parse_flag("allow_refit", "false"), default = FALSE)
allow_kfold <- parse_bool(parse_flag("allow_kfold", "false"), default = FALSE)
per_site <- parse_int(parse_flag("per_site", "30"), default = 30)
seed <- parse_int(parse_flag("seed", "2027"), default = 2027)
k_site <- parse_int(parse_flag("k_site", "5"), default = 5)

# Quick refit settings (used only if allow_refit/allow_kfold true)
refit_chains <- parse_int(parse_flag("refit_chains", "2"), default = 2)
refit_iter <- parse_int(parse_flag("refit_iter", "400"), default = 400)
refit_adapt <- parse_flag("refit_adapt_delta", "0.99")
refit_treedepth <- parse_flag("refit_max_treedepth", "15")

kfold_chains <- parse_int(parse_flag("kfold_chains", "1"), default = 1)
kfold_iter <- parse_int(parse_flag("kfold_iter", "400"), default = 400)
kfold_adapt <- parse_flag("kfold_adapt_delta", "0.99")
kfold_treedepth <- parse_flag("kfold_max_treedepth", "15")

made_inputs <- FALSE

ensure_subsample <- function(in_path, out_path, label) {
  if (file.exists(out_path)) return()
  if (!file.exists(in_path)) {
    warning("Missing input for subsample: ", in_path, " (", label, ")")
    return()
  }
  run_rscript(
    args = c("R/05_make_stroop_subsample.R",
             "--in", in_path,
             "--out", out_path,
             "--per_site", as.character(per_site),
             "--seed", as.character(seed)),
    desc = paste0("Make subsample (", label, ")")
  )
  made_inputs <<- TRUE
}

ensure_hetero_reports <- function(tag, data_path, label) {
  summ <- file.path("reports", paste0("location_scale_", tag, "_summary.csv"))
  part <- file.path("reports", paste0("location_scale_", tag, "_participants.csv"))
  if (file.exists(summ) && file.exists(part)) return(TRUE)

  rds <- file.path("models", paste0("location_scale_", tag, ".rds"))
  csv_dir <- file.path("models", paste0("cmdstan_", tag))
  csv_pat <- paste0("^location_scale_", tag, "-[0-9]+\\.csv$")
  can_reuse <- file.exists(rds) || has_cmdstan_csv(csv_dir, csv_pat)

  if (can_reuse) {
    run_rscript(
      args = c("R/04_fit_stroop_location_scale.R", data_path,
               "--tag", tag,
               "--reuse_hetero", "true",
               "--loo", "false",
               "--compare_homo", "false"),
      desc = paste0("Rebuild hetero reports (", label, ")")
    )
    made_inputs <<- TRUE
    return(TRUE)
  }

  if (!allow_refit) {
    warning("Missing hetero outputs for ", label, " (tag=", tag, "). Set --allow_refit true to fit.")
    return(FALSE)
  }

  run_rscript(
    args = c("R/04_fit_stroop_location_scale.R", data_path,
             "--tag", tag,
             "--loo", "false",
             "--compare_homo", "false",
             "--chains", as.character(refit_chains),
             "--parallel_chains", as.character(refit_chains),
             "--iter", as.character(refit_iter),
             "--init", "0",
             "--adapt_delta", refit_adapt,
             "--max_treedepth", refit_treedepth),
    desc = paste0("Fit hetero (", label, ")")
  )
  made_inputs <<- TRUE
  TRUE
}

ensure_homo_summary <- function(tag, data_path, label) {
  summ <- file.path("reports", paste0("location_scale_homo_", tag, "_summary.csv"))
  if (file.exists(summ)) return(TRUE)

  rds <- file.path("models", paste0("location_scale_homo_", tag, ".rds"))
  csv_dir <- file.path("models", paste0("cmdstan_homo_", tag))
  csv_pat <- paste0("^location_scale_homo_", tag, "-[0-9]+\\.csv$")
  can_reuse <- file.exists(rds) || has_cmdstan_csv(csv_dir, csv_pat)

  if (can_reuse) {
    run_rscript(
      args = c("R/04_fit_stroop_location_scale.R", data_path,
               "--tag", tag,
               "--homo_only", "true",
               "--reuse_homo", "true",
               "--loo", "false"),
      desc = paste0("Rebuild homo summary (", label, ")")
    )
    made_inputs <<- TRUE
    return(TRUE)
  }

  if (!allow_refit) {
    warning("Missing homo summary for ", label, " (tag=", tag, "). Set --allow_refit true to fit.")
    return(FALSE)
  }

  run_rscript(
    args = c("R/04_fit_stroop_location_scale.R", data_path,
             "--tag", tag,
             "--homo_only", "true",
             "--compare_homo", "true",
             "--loo", "false",
             "--chains", as.character(refit_chains),
             "--parallel_chains", as.character(refit_chains),
             "--iter", as.character(refit_iter),
             "--init", "0",
             "--adapt_delta", refit_adapt,
             "--max_treedepth", refit_treedepth),
    desc = paste0("Fit homo baseline (", label, ")")
  )
  made_inputs <<- TRUE
  TRUE
}

ensure_kfold <- function(data_path, tag, label) {
  bf <- file.path("reports", "kfold", paste0("location_scale_", tag, "_kfold_site_bf.csv"))
  if (file.exists(bf)) return(TRUE)
  if (!allow_kfold) {
    warning("Missing k-fold outputs for ", label, " (tag=", tag, "). Set --allow_kfold true to run.")
    return(FALSE)
  }
  run_rscript(
    args = c("R/11_kfold_location_scale.R", data_path,
             "--tag", tag,
             "--unit", "site",
             "--k_site", as.character(k_site),
             "--chains", as.character(kfold_chains),
             "--parallel_chains", as.character(kfold_chains),
             "--iter", as.character(kfold_iter),
             "--init", "0",
             "--adapt_delta", kfold_adapt,
             "--max_treedepth", kfold_treedepth,
             "--resume", "true"),
    desc = paste0("Site K-fold (", label, ")")
  )
  made_inputs <<- TRUE
  TRUE
}

# 1) Ensure subsamples for K-fold (if needed)
ensure_subsample(stroop_data, stroop_sub_data, "Stroop")
ensure_subsample(psa_attr_data, psa_attr_sub_data, "PSA001_Attractive")
ensure_subsample(psa_dom_data, psa_dom_sub_data, "PSA001_Dominant")

# 2) Ensure hetero summaries/participants
stroop_ok <- ensure_hetero_reports(stroop_tag, stroop_data, "Stroop")
psa_attr_ok <- ensure_hetero_reports(psa_attr_tag, psa_attr_data, "PSA001_Attractive")
psa_dom_ok <- ensure_hetero_reports(psa_dom_tag, psa_dom_data, "PSA001_Dominant")

# 3) Ensure homo summaries (needed for shift table)
stroop_homo_ok <- ensure_homo_summary(stroop_tag, stroop_data, "Stroop")
psa_attr_homo_ok <- ensure_homo_summary(psa_attr_tag, psa_attr_data, "PSA001_Attractive")
psa_dom_homo_ok <- ensure_homo_summary(psa_dom_tag, psa_dom_data, "PSA001_Dominant")

# 4) Ensure K-fold outputs
stroop_kfold_ok <- ensure_kfold(stroop_sub_data, stroop_cv_tag, "Stroop")
psa_attr_kfold_ok <- ensure_kfold(psa_attr_sub_data, psa_attr_cv_tag, "PSA001_Attractive")
psa_dom_kfold_ok <- ensure_kfold(psa_dom_sub_data, psa_dom_cv_tag, "PSA001_Dominant")

# 5) Build RQ tables (only if inputs exist)

# Shift table (homo vs hetero)
rq1_stroop_tag <- if (stroop_ok && stroop_homo_ok) stroop_tag else "none"
rq1_psa_attr_tag <- if (psa_attr_ok && psa_attr_homo_ok) psa_attr_tag else "none"
rq1_psa_dom_tag <- if (psa_dom_ok && psa_dom_homo_ok) psa_dom_tag else "none"

if (rq1_stroop_tag != "none" || rq1_psa_attr_tag != "none" || rq1_psa_dom_tag != "none") {
  run_rscript(
    args = c("R/12_rq1_shift_table.R",
             paste0("--stroop_tag=", rq1_stroop_tag),
             paste0("--psa001_attractive_tag=", rq1_psa_attr_tag),
             paste0("--psa001_dominant_tag=", rq1_psa_dom_tag),
             "--out=reports/rq1_shift_table.csv"),
    desc = "Shift table (homo vs hetero)"
  )
}

# Participant prevalence
if (stroop_ok || psa_attr_ok || psa_dom_ok) {
  run_rscript(
    args = c("R/07_participants_prevalence.R",
             paste0("--include_psa001_attractive=", ifelse(psa_attr_ok, "true", "false")),
             paste0("--include_psa001_dominant=", ifelse(psa_dom_ok, "true", "false")),
             "--out_detail=reports/person_prevalence_detail.csv",
             "--out_summary=reports/person_prevalence_summary.csv"),
    desc = "Participant prevalence summary"
  )
}

# Site variance decomposition (supplemental)
if (stroop_ok || psa_attr_ok || psa_dom_ok) {
  run_rscript(
    args = c("R/08_site_variance_decomp.R",
             paste0("--include_psa001_attractive=", ifelse(psa_attr_ok, "true", "false")),
             paste0("--include_psa001_dominant=", ifelse(psa_dom_ok, "true", "false")),
             "--out_summary=reports/site_variance_decomposition.csv",
             "--out_site=reports/site_level_mix_precision.csv"),
    desc = "Site variance decomposition (supplemental)"
  )
}

# Holdout comparison tables (requires summaries + k-fold)
kfold_stroop_tag <- if (stroop_ok && stroop_kfold_ok) stroop_tag else "none"
kfold_psa_attr_tag <- if (psa_attr_ok && psa_attr_kfold_ok) psa_attr_tag else "none"
kfold_psa_dom_tag <- if (psa_dom_ok && psa_dom_kfold_ok) psa_dom_tag else "none"

if (kfold_stroop_tag != "none" || kfold_psa_attr_tag != "none" || kfold_psa_dom_tag != "none") {
  run_rscript(
    args = c("R/06_site_models_table.R",
             "--out=reports/site_model_comparisons.csv",
             paste0("--stroop_tag=", kfold_stroop_tag),
             paste0("--stroop_cv_tag=", stroop_cv_tag),
             paste0("--psa001_attractive_tag=", kfold_psa_attr_tag),
             paste0("--psa001_attractive_cv_tag=", psa_attr_cv_tag),
             paste0("--psa001_dominant_tag=", kfold_psa_dom_tag),
             paste0("--psa001_dominant_cv_tag=", psa_dom_cv_tag)),
    desc = "Site model comparisons"
  )

  run_rscript(
    args = c("R/09_site_kfold_stack.R"),
    desc = "Site K-fold comparison stack"
  )
}

cat("\nDone. If anything was skipped, rerun with --allow_refit true and/or --allow_kfold true.\n")
