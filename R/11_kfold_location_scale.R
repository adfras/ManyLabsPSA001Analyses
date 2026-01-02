#!/usr/bin/env Rscript
# K-fold cross-validation for the location–scale model (heteroskedastic vs homoskedastic).
#
# Why this exists:
# - Grouped PSIS-LOO can be unreliable for hierarchical random-effects models (high Pareto-k).
# - K-fold CV avoids the PSIS importance-weight failure mode, at the cost of refits.
#
# Key idea:
# - We add observation weights `w` to the Stan programs.
# - For fold f, set w[n]=0 for held-out observations and w[n]=1 for training observations.
# - This lets us evaluate held-out log predictive density using generated-quantities `log_lik`
#   while fitting on training data only, without changing parameter dimensions.

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
})
if (file.exists("R/00_run_manifest.R")) source("R/00_run_manifest.R")

override <- Sys.getenv("CMDSTAN_OVERRIDE", unset = "")
if (nzchar(override)) Sys.setenv(CMDSTAN = override)
cp <- Sys.getenv("CMDSTAN", unset = "")
if (nzchar(cp)) cmdstanr::set_cmdstan_path(cp)

source(file.path("R", "lib", "cli_utils.R"))

log_mean_exp_cols <- function(mat) {
  # mat: draws x n_points
  m <- apply(mat, 2, max)
  m + log(colMeans(exp(sweep(mat, 2, m, "-"))))
}

se_sum <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sqrt(length(x) * stats::var(x))
}

make_balanced_folds <- function(n_groups, k, seed) {
  k <- as.integer(k)
  if (!is.finite(k) || k < 2) stop("k must be >= 2")
  if (k > n_groups) {
    warning("Requested k=", k, " folds but only n_groups=", n_groups, "; using k=n_groups.")
    k <- n_groups
  }
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n_groups))
}

safe_unlink_tmp <- function(path, tmp_root) {
  # Avoid accidental deletion outside the expected tmp root.
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(tmp_root, winslash = "/", mustWork = TRUE)
  if (!startsWith(path_norm, root_norm)) {
    stop("Refusing to delete non-tmp path: ", path_norm, " (tmp root: ", root_norm, ")")
  }
  unlink(path, recursive = TRUE, force = TRUE)
}

pos_args <- commandArgs(trailingOnly = TRUE)
in_path <- if (length(pos_args) >= 1 && !startsWith(pos_args[[1]], "--")) pos_args[[1]] else parse_flag("data", "data/processed/trials_stroop_ml3_with_site_sub30.csv")
tag <- parse_flag("tag", "stroop_ml3_site_sub30_mm")
unit_raw <- tolower(parse_flag("unit", "both")) # person|site|both or comma-separated list

chains <- as.integer(parse_flag("chains", 2))
iter <- as.integer(parse_flag("iter", 600))
thin <- as.integer(parse_flag("thin", 1))
parallel_chains <- as.integer(parse_flag("parallel_chains", chains))

k_global <- parse_num(parse_flag("k", NA), default = NA)
k_person <- as.integer(parse_flag("k_person", 10))
k_site <- as.integer(parse_flag("k_site", 10))
if (!is.na(k_global)) {
  k_person <- as.integer(k_global)
  k_site <- as.integer(k_global)
}

seed <- as.integer(parse_flag("seed", 2027))
adapt_delta <- parse_num(parse_flag("adapt_delta", 0.99), default = 0.99)
max_treedepth <- as.integer(parse_flag("max_treedepth", 15))
stepsize <- parse_num(parse_flag("stepsize", NA), default = NA)
init_val <- parse_flag("init", NA)
refresh <- as.integer(parse_flag("refresh", 50))
save_warmup <- parse_bool(parse_flag("save_warmup", "false"), default = FALSE)

homo_person_slope <- parse_bool(parse_flag("homo_person_slope", "false"), default = FALSE)
keep_draws <- parse_bool(parse_flag("keep_draws", "false"), default = FALSE)
resume <- parse_bool(parse_flag("resume", "true"), default = TRUE)
write_pointwise <- parse_bool(parse_flag("write_pointwise", "true"), default = TRUE)

units <- strsplit(unit_raw, ",", fixed = TRUE)[[1]] |> trimws()
if (identical(units, "both") || all(units == "")) units <- c("person", "site")
units <- unique(units)
if (!all(units %in% c("person", "site"))) {
  stop("--unit must be one of: person, site, both (or comma-separated person,site)")
}

if (!file.exists(in_path)) stop("Input file not found: ", in_path)
df <- readr::read_csv(in_path, show_col_types = FALSE)
if (!all(c("y", "person", "site") %in% names(df))) {
  stop("Input must have columns y, person, and site (got: ", paste(names(df), collapse = ", "), ")")
}

preds <- grep("^X_", names(df), value = TRUE)
X <- as.matrix(cbind(Intercept = 1, if (length(preds)) df[, preds, drop = FALSE] else NULL))
P <- ncol(X)

person_factor <- factor(df$person)
site_factor <- factor(df$site)
J <- nlevels(person_factor)
K <- nlevels(site_factor)
person_int <- as.integer(person_factor)
site_int <- as.integer(site_factor)
person_site <- as.integer(tapply(site_int, person_factor, function(x) x[[1]]))

stan_data_base <- list(
  N = nrow(df),
  J = J,
  K = K,
  P = P,
  X = X,
  person = person_int,
  site = site_int,
  person_site = person_site,
  y = df$y,
  w = rep(1L, nrow(df))
)

# Record the K-fold run configuration.
if (exists("append_run_manifest")) {
  append_run_manifest(list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"),
    script = "R/11_kfold_location_scale.R",
    action = "kfold_run",
    model = "kfold",
    tag = tag,
    data = in_path,
    N = stan_data_base$N,
    J = stan_data_base$J,
    K = stan_data_base$K,
    P = stan_data_base$P,
    unit_raw = unit_raw,
    units = paste(units, collapse = ";"),
    k_person = k_person,
    k_site = k_site,
    k_global = k_global,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter,
    iter_sampling = iter,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    stepsize = ifelse(is.na(stepsize), NA, stepsize),
    init = init_val,
    refresh = refresh,
    save_warmup = save_warmup,
    homo_person_slope = homo_person_slope,
    keep_draws = keep_draws,
    resume = resume,
    write_pointwise = write_pointwise,
    cmd = paste(commandArgs(), collapse = " ")
  ))
}

dir.create("models", showWarnings = FALSE, recursive = TRUE)
dir.create("reports", showWarnings = FALSE, recursive = TRUE)
out_dir <- file.path("reports", "kfold")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

compile_dir <- file.path("models", "cmdstan_compiled")
dir.create(compile_dir, showWarnings = FALSE, recursive = TRUE)

hetero_mod <- cmdstan_model(
  # K-fold uses holdout-only log-lik to keep CmdStan CSV outputs small.
  file.path("stan", "location_scale_kfold_holdout_lik.stan"),
  dir = compile_dir,
  compile = TRUE,
  cpp_options = list(STAN_NO_PCH = TRUE)
)

homo_stan_file <- if (homo_person_slope) {
  "location_scale_homoskedastic_kfold_holdout_lik.stan"
} else {
  "location_scale_homoskedastic_no_person_slope_kfold_holdout_lik.stan"
}
homo_mod <- cmdstan_model(
  file.path("stan", homo_stan_file),
  dir = compile_dir,
  compile = TRUE,
  cpp_options = list(STAN_NO_PCH = TRUE)
)

run_one_unit <- function(unit, k_folds) {
  group_levels <- if (unit == "person") levels(person_factor) else levels(site_factor)
  n_groups <- length(group_levels)
  if (n_groups < 2) stop("Need at least 2 groups for unit=", unit)
  k_folds <- as.integer(min(k_folds, n_groups))

  assign_path <- file.path(out_dir, paste0("location_scale_", tag, "_kfold_", unit, "_assignments.csv"))
  progress_path <- file.path(out_dir, paste0("location_scale_", tag, "_kfold_", unit, "_progress.rds"))
  folds_path <- file.path(out_dir, paste0("location_scale_", tag, "_kfold_", unit, "_folds.csv"))
  bf_path <- file.path(out_dir, paste0("location_scale_", tag, "_kfold_", unit, "_bf.csv"))
  pointwise_path <- file.path(out_dir, paste0("location_scale_", tag, "_kfold_", unit, "_pointwise.csv"))

  folds_by_group <- NULL
  if (resume && file.exists(assign_path)) {
    a <- readr::read_csv(assign_path, show_col_types = FALSE)
    if (!all(c("group", "fold") %in% names(a))) stop("Invalid assignments file: ", assign_path)
    a <- a %>% mutate(group = as.character(group), fold = as.integer(fold))
    a <- a %>% filter(group %in% group_levels)
    folds_by_group <- a$fold[match(group_levels, a$group)]
    if (anyNA(folds_by_group)) stop("Assignments file does not cover all groups for unit=", unit)
    k_folds <- max(folds_by_group, na.rm = TRUE)
  } else {
    folds_by_group <- make_balanced_folds(n_groups, k_folds, seed = seed)
    readr::write_csv(tibble(group = group_levels, fold = folds_by_group), assign_path)
  }

  obs_fold <- if (unit == "person") folds_by_group[person_int] else folds_by_group[site_int]

  progress <- NULL
  if (resume && file.exists(progress_path)) {
    progress <- readRDS(progress_path)
    ok <- is.list(progress) &&
      length(progress$elpd_hetero) == stan_data_base$N &&
      length(progress$elpd_homo) == stan_data_base$N &&
      length(progress$fold_complete) == k_folds
    if (!ok) {
      warning("Ignoring incompatible progress file: ", progress_path)
      progress <- NULL
    }
  }
  if (is.null(progress)) {
    progress <- list(
      elpd_hetero = rep(NA_real_, stan_data_base$N),
      elpd_homo = rep(NA_real_, stan_data_base$N),
      fold_complete = rep(FALSE, k_folds),
      fold_summary = vector("list", k_folds)
    )
  }

  tmp_root <- file.path("models", "kfold_tmp")
  dir.create(tmp_root, showWarnings = FALSE, recursive = TRUE)

  for (f in seq_len(k_folds)) {
    if (isTRUE(progress$fold_complete[[f]])) {
      message("K-fold (", unit, "): fold ", f, "/", k_folds, " already complete; skipping.")
      next
    }

    hold_idx <- which(obs_fold == f)
    if (!length(hold_idx)) {
      warning("Fold ", f, " has 0 held-out observations for unit=", unit, "; skipping.")
      progress$fold_complete[[f]] <- TRUE
      saveRDS(progress, progress_path)
      next
    }

    w <- as.integer(obs_fold != f)
    stan_data <- stan_data_base
    stan_data$w <- w
    stan_data$N_holdout <- length(hold_idx)
    stan_data$holdout_idx <- as.integer(hold_idx)

    message(
      "K-fold (", unit, "): starting fold ", f, "/", k_folds,
      " (holdout obs=", length(hold_idx),
      ", train obs=", sum(w),
      ")"
    )

    base_dir <- if (keep_draws) file.path("models", "kfold", tag, unit, paste0("fold_", f)) else file.path(tmp_root, tag, unit, paste0("fold_", f))
    hetero_dir <- file.path(base_dir, "hetero")
    homo_dir <- file.path(base_dir, "homo")
    dir.create(hetero_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(homo_dir, showWarnings = FALSE, recursive = TRUE)

    t_start <- Sys.time()

    hetero_args <- list(
      data = stan_data,
      seed = seed + f,
      chains = chains,
      parallel_chains = min(parallel_chains, chains),
      iter_warmup = iter,
      iter_sampling = iter,
      thin = thin,
      output_dir = hetero_dir,
      output_basename = paste0("location_scale_", tag, "_kfold_", unit, "_fold", f),
      refresh = refresh
    )
    if (!is.na(adapt_delta)) hetero_args$adapt_delta <- adapt_delta
    if (!is.na(max_treedepth)) hetero_args$max_treedepth <- max_treedepth
    if (!is.na(stepsize)) hetero_args$step_size <- stepsize
    if (!is.na(init_val)) hetero_args$init <- eval(parse(text = init_val))
    if (isTRUE(save_warmup)) hetero_args$save_warmup <- TRUE

    hetero_fit <- do.call(hetero_mod$sample, hetero_args)
    t_after_hetero <- Sys.time()
    message(
      "K-fold (", unit, "): fold ", f, "/", k_folds, " hetero done (",
      round(as.numeric(difftime(t_after_hetero, t_start, units = "mins")), 2),
      " min). Running homo..."
    )
    ll_hetero <- hetero_fit$draws("log_lik_holdout", format = "matrix")
    elpd_hetero_i <- log_mean_exp_cols(ll_hetero)
    progress$elpd_hetero[hold_idx] <- elpd_hetero_i

    homo_args <- list(
      data = stan_data,
      seed = seed + 10000 + f,
      chains = chains,
      parallel_chains = min(parallel_chains, chains),
      iter_warmup = iter,
      iter_sampling = iter,
      thin = thin,
      output_dir = homo_dir,
      output_basename = paste0("location_scale_homo_", tag, "_kfold_", unit, "_fold", f),
      refresh = refresh
    )
    if (!is.na(adapt_delta)) homo_args$adapt_delta <- adapt_delta
    if (!is.na(max_treedepth)) homo_args$max_treedepth <- max_treedepth
    if (!is.na(stepsize)) homo_args$step_size <- stepsize
    if (!is.na(init_val)) homo_args$init <- eval(parse(text = init_val))
    if (isTRUE(save_warmup)) homo_args$save_warmup <- TRUE

    homo_fit <- do.call(homo_mod$sample, homo_args)
    t_after_homo <- Sys.time()
    message(
      "K-fold (", unit, "): fold ", f, "/", k_folds, " homo done (",
      round(as.numeric(difftime(t_after_homo, t_after_hetero, units = "mins")), 2),
      " min). Computing ELPD..."
    )
    ll_homo <- homo_fit$draws("log_lik_holdout", format = "matrix")
    elpd_homo_i <- log_mean_exp_cols(ll_homo)
    progress$elpd_homo[hold_idx] <- elpd_homo_i

    t_end <- Sys.time()
    elapsed_s <- as.numeric(difftime(t_end, t_start, units = "secs"))
    message(
      "K-fold (", unit, "): fold ", f, "/", k_folds, " complete (",
      round(elapsed_s / 60, 2),
      " min)."
    )

    progress$fold_summary[[f]] <- tibble(
      unit = unit,
      fold = f,
      n_holdout = length(hold_idx),
      elpd_hetero = sum(elpd_hetero_i),
      elpd_homo = sum(elpd_homo_i),
      delta_elpd = sum(elpd_hetero_i - elpd_homo_i),
      elapsed_seconds = elapsed_s,
      chains = chains,
      iter = iter
    )

    progress$fold_complete[[f]] <- TRUE
    saveRDS(progress, progress_path)

    if (!keep_draws) {
      safe_unlink_tmp(base_dir, tmp_root = tmp_root)
    }
  }

  if (anyNA(progress$elpd_hetero) || anyNA(progress$elpd_homo)) {
    stop("K-fold incomplete for unit=", unit, ". Missing ELPD values remain; see: ", progress_path)
  }

  diff_i <- progress$elpd_hetero - progress$elpd_homo
  elpd_hetero <- sum(progress$elpd_hetero)
  elpd_homo <- sum(progress$elpd_homo)
  elpd_diff <- sum(diff_i)

  se_hetero <- se_sum(progress$elpd_hetero)
  se_homo <- se_sum(progress$elpd_homo)
  se_diff <- se_sum(diff_i)

  bf_hetero_vs_homo <- if (elpd_diff > 700) Inf else exp(elpd_diff)
  bf_homo_vs_hetero <- if (elpd_diff < -700) Inf else exp(-elpd_diff)

  bf_tbl <- tibble(
    model = c("heterogeneous", "homogeneous"),
    method = "kfold",
    unit = unit,
    k_folds = k_folds,
    chains = chains,
    iter_warmup = iter,
    iter_sampling = iter,
    elpd = c(elpd_hetero, elpd_homo),
    elpd_se = c(se_hetero, se_homo),
    looic = c(-2 * elpd_hetero, -2 * elpd_homo),
    looic_se = c(2 * se_hetero, 2 * se_homo),
    log_bayes_factor_hetero_vs_homo = c(elpd_diff, -elpd_diff),
    bayes_factor_hetero_vs_homo = c(bf_hetero_vs_homo, bf_homo_vs_hetero),
    delta_elpd = c(elpd_diff, -elpd_diff),
    delta_elpd_se = c(se_diff, se_diff)
  )

  fold_tbl <- bind_rows(progress$fold_summary) %>%
    arrange(fold)

  readr::write_csv(fold_tbl, folds_path)
  readr::write_csv(bf_tbl, bf_path)

  if (write_pointwise) {
    point_tbl <- tibble(
      row = seq_len(stan_data_base$N),
      person = as.character(df$person),
      site = as.character(df$site),
      fold = obs_fold,
      elpd_hetero = progress$elpd_hetero,
      elpd_homo = progress$elpd_homo,
      delta_elpd = diff_i
    )
    readr::write_csv(point_tbl, pointwise_path)
  }

  message("K-fold (", unit, ") complete. Wrote: ", bf_path)
  bf_tbl
}

bf_results <- list()
if ("person" %in% units) bf_results[["person"]] <- run_one_unit("person", k_person)
if ("site" %in% units) bf_results[["site"]] <- run_one_unit("site", k_site)

summary_path <- file.path(out_dir, paste0("location_scale_", tag, "_kfold_summary.csv"))
summary_tbl <- bind_rows(bf_results) %>%
  group_by(unit) %>%
  summarise(
    k_folds = first(k_folds),
    elpd_hetero = elpd[model == "heterogeneous"][1],
    elpd_homo = elpd[model == "homogeneous"][1],
    delta_elpd = elpd_hetero - elpd_homo,
    delta_elpd_se = delta_elpd_se[model == "heterogeneous"][1],
    looic_hetero = looic[model == "heterogeneous"][1],
    looic_homo = looic[model == "homogeneous"][1],
    .groups = "drop"
  )
readr::write_csv(summary_tbl, summary_path)
message("Wrote combined summary: ", summary_path)
