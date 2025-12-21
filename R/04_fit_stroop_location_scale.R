#!/usr/bin/env Rscript
# Hierarchical location-scale model with participant-level variance
# Usage: Rscript R/model_location_scale.R data.csv [--tag myrun] [--chains 4] [--iter 1500]
# Expects columns: y (numeric), person (id), optional predictors prefixed X_ (will auto add intercept).

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(loo)
})
if (file.exists("R/00_run_manifest.R")) source("R/00_run_manifest.R")

# Honor CMDSTAN env var if provided (useful when the default install is read-only)
override <- Sys.getenv("CMDSTAN_OVERRIDE", unset = "")
if (nzchar(override)) Sys.setenv(CMDSTAN = override)
cp <- Sys.getenv("CMDSTAN", unset = "")
if (nzchar(cp)) cmdstanr::set_cmdstan_path(cp)

parse_flag <- function(name, default = NULL) {
  args <- paste(commandArgs(), collapse = " ")
  m <- regexpr(paste0("--", name, "(=| )([^ ]+)"), args)
  if (m[1] == -1) return(default)
  sub(".*--", "", regmatches(args, m)) |>
    sub(paste0(name, "(=| )"), "", x = _) |>
    trimws()
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "yes", "y", "t")
}

parse_num <- function(x, default = NULL) {
  if (is.null(x)) return(default)
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out)) return(default)
  out
}

pos_args <- commandArgs(trailingOnly = TRUE)
in_path <- if (length(pos_args) >= 1 && !startsWith(pos_args[[1]], "--")) pos_args[[1]] else parse_flag("data", "data/processed/trials.csv")
tag <- parse_flag("tag", tools::file_path_sans_ext(basename(in_path)))
chains <- as.integer(parse_flag("chains", 4))
iter <- as.integer(parse_flag("iter", 1500))
compare_homo <- parse_bool(parse_flag("compare_homo", "true"), default = TRUE)
do_loo <- parse_bool(parse_flag("loo", "true"), default = TRUE)
loo_unit <- tolower(parse_flag("loo_unit", "obs"))
moment_match <- parse_bool(parse_flag("moment_match", "false"), default = FALSE)
force_loo <- parse_bool(parse_flag("force_loo", "false"), default = FALSE)
loo_max_n <- as.integer(parse_flag("loo_max_n", "50000"))
sigma_outlier_prob <- parse_bool(parse_flag("sigma_outlier_prob", "false"), default = FALSE)
reuse_hetero <- parse_bool(parse_flag("reuse_hetero", "false"), default = FALSE)
reuse_homo <- parse_bool(parse_flag("reuse_homo", "false"), default = FALSE)
homo_person_slope <- parse_bool(parse_flag("homo_person_slope", "true"), default = TRUE)
homo_only <- parse_bool(parse_flag("homo_only", "false"), default = FALSE)
adapt_delta <- parse_num(parse_flag("adapt_delta", NA), default = NA)
max_treedepth <- parse_num(parse_flag("max_treedepth", NA), default = NA)
stepsize <- parse_num(parse_flag("stepsize", NA), default = NA)
thin <- as.integer(parse_flag("thin", 1))
parallel_chains <- as.integer(parse_flag("parallel_chains", chains))
likelihood <- parse_flag("likelihood", "normal")
init_val <- parse_flag("init", NA) # allow --init 0 for safer starts
use_site_effects <- parse_bool(parse_flag("site_effects", "true"), default = TRUE)
refresh <- as.integer(parse_flag("refresh", 10))
save_warmup <- parse_bool(parse_flag("save_warmup", "false"), default = FALSE)

# Capture requested flags before any automatic adjustments.
req_compare_homo <- compare_homo
req_do_loo <- do_loo
req_use_site_effects <- use_site_effects
req_loo_unit <- loo_unit
req_likelihood <- likelihood
cmd_args <- paste(commandArgs(), collapse = " ")
loo_disabled_reason <- NA_character_
# The Student-t variant was experimentally explored for robustness but showed
# poor geometry and unstable diagnostics for the ManyLabs data. To avoid
# accidentally using it in production runs, we currently disable it here.
if (likelihood == "student_t") {
  stop(
    "The student_t likelihood variant is currently disabled for this script ",
    "because it produced unstable fits for these data. ",
    "Use --likelihood normal (the default) or edit R/model_location_scale.R ",
    "and stan/location_scale_student_t.stan if you intentionally want to ",
    "experiment with the Student-t model."
  )
}
if (!identical(likelihood, "normal")) {
  stop("Unknown likelihood: ", likelihood, " (expected 'normal')")
}

if (!file.exists(in_path)) stop("Input file not found: ", in_path)
df <- readr::read_csv(in_path, show_col_types = FALSE)
if (!all(c("y", "person") %in% names(df))) stop("Input must have columns y and person")

if (homo_only) {
  compare_homo <- TRUE
  do_loo <- FALSE
  loo_disabled_reason <- "homo_only"
}

if (!loo_unit %in% c("obs", "person", "site")) {
  stop("Unknown --loo_unit: ", loo_unit, " (expected one of: obs, person, site)")
}
loo_var <- switch(loo_unit, obs = "log_lik", person = "log_lik_person", site = "log_lik_site")

if (do_loo && loo_unit != "obs") {
  message(
    "Note: --loo_unit=", loo_unit, " requests PSIS-LOO over entire groups. ",
    "For hierarchical random-effects models, PSIS can be unreliable when ",
    "leaving out whole persons/sites (high Pareto-k). ",
    "For defensible group-level generalization, prefer K-fold CV via ",
    "R/11_kfold_location_scale.R."
  )
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
# map each person to their first observed site (aligned with factor levels)
person_site <- as.integer(tapply(site_int, person_factor, function(x) x[[1]]))

if (!use_site_effects) {
  K <- 1L
  site_int <- rep(1L, nrow(df))
  person_site <- rep(1L, J)
  compare_homo <- FALSE  # site comparison not meaningful without site effects
  if (is.na(loo_disabled_reason)) loo_disabled_reason <- "site_effects_disabled"
}

stan_data <- list(
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

# Output paths (used by both sampling + reuse)
dir.create("models", showWarnings = FALSE, recursive = TRUE)
dir.create("reports", showWarnings = FALSE, recursive = TRUE)
model_path <- file.path("models", paste0("location_scale_", tag, ".rds"))
homo_model_path <- file.path("models", paste0("location_scale_homo_", tag, ".rds"))

# Keep compiled CmdStan binaries out of `stan/` to reduce clutter (and to avoid
# generating platform-specific executables next to the Stan source files).
compile_dir <- file.path("models", "cmdstan_compiled")
dir.create(compile_dir, showWarnings = FALSE, recursive = TRUE)

# Full-data PSIS-LOO is often infeasible at the trial level for very large N
# because fit$loo() must first read all draws into memory (including log_lik).
# Default behavior: automatically disable --loo/--compare_homo for large N unless
# the user explicitly opts in via --force_loo true.
if (do_loo && loo_unit == "obs" && !force_loo && !is.na(loo_max_n) && stan_data$N > loo_max_n) {
  message(
    "N=", stan_data$N, " is larger than loo_max_n=", loo_max_n, ". ",
    "Disabling --loo/--compare_homo to avoid huge log_lik memory use. ",
    "Use subsample fits for LOO/BF, or re-run with --force_loo true if you ",
    "really want full-data LOO and have sufficient RAM."
  )
  do_loo <- FALSE
  compare_homo <- FALSE
  loo_disabled_reason <- "loo_max_n"
}

# If --loo=false, we can still fit the homoskedastic baseline for RQ1-style
# comparisons, but we will not compute LOO/BF outputs.
if (!do_loo && compare_homo) {
  message("Fitting homoskedastic baseline without LOO/BF (--loo=false).")
}
mod <- NULL
if (!homo_only) {
  stan_file <- if (!do_loo) {
    file.path("stan", "location_scale_light.stan")
  } else if (loo_unit == "obs") {
    file.path("stan", "location_scale.stan")
  } else {
    # Grouped PSIS-LOO over persons/sites only needs a length-J/length-K log-lik,
    # not a length-N log_lik vector. This keeps CmdStan CSVs small and makes full
    # multi-site datasets tractable on normal machines.
    file.path("stan", "location_scale_group_lik.stan")
  }
  mod <- cmdstan_model(
    stan_file,
    dir = compile_dir,
    compile = TRUE,
    cpp_options = list(STAN_NO_PCH = TRUE),
    force_recompile = moment_match,
    compile_model_methods = moment_match
  )
}

output_dir <- file.path("models", paste0("cmdstan_", tag))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
sample_args <- list(
  data = stan_data,
  seed = 2027,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter,
  iter_sampling = iter,
  thin = thin,
  output_dir = output_dir,
  output_basename = paste0("location_scale_", tag),
  refresh = refresh
)
if (!is.na(adapt_delta)) sample_args$adapt_delta <- adapt_delta
if (!is.na(max_treedepth)) sample_args$max_treedepth <- max_treedepth
if (!is.na(stepsize)) sample_args$step_size <- stepsize
if (!is.na(init_val)) sample_args$init <- eval(parse(text = init_val))
if (isTRUE(save_warmup)) sample_args$save_warmup <- TRUE

append_manifest_row <- function(model_type, action, seed_val) {
  if (!exists("append_run_manifest")) return(invisible(FALSE))
  append_run_manifest(list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"),
    script = "R/04_fit_stroop_location_scale.R",
    action = action,
    model = model_type,
    tag = tag,
    data = in_path,
    N = stan_data$N,
    J = stan_data$J,
    K = stan_data$K,
    P = stan_data$P,
    likelihood = req_likelihood,
    use_site_effects = use_site_effects,
    homo_only = homo_only,
    compare_homo_requested = req_compare_homo,
    compare_homo_effective = compare_homo,
    loo_requested = req_do_loo,
    loo_effective = do_loo,
    loo_unit = req_loo_unit,
    loo_disabled_reason = loo_disabled_reason,
    moment_match = moment_match,
    force_loo = force_loo,
    loo_max_n = loo_max_n,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter,
    iter_sampling = iter,
    thin = thin,
    adapt_delta = ifelse(is.na(adapt_delta), NA, adapt_delta),
    max_treedepth = ifelse(is.na(max_treedepth), NA, max_treedepth),
    stepsize = ifelse(is.na(stepsize), NA, stepsize),
    init = init_val,
    seed = seed_val,
    refresh = refresh,
    save_warmup = save_warmup,
    reuse_hetero = reuse_hetero,
    reuse_homo = reuse_homo,
    cmd = cmd_args
  ))
}

fit <- NULL
if (!homo_only) {
  if (reuse_hetero && file.exists(model_path)) {
    message("Loading existing heteroskedastic fit from ", model_path)
    append_manifest_row("hetero", "reuse_rds", 2027)
    fit <- tryCatch(readRDS(model_path), error = function(e) {
      message(
        "Failed to read ", model_path, " (", conditionMessage(e), "). ",
        "Falling back to CmdStan CSV outputs under ", output_dir, "."
      )
      NULL
    })
  }
  if (is.null(fit) && reuse_hetero) {
    pat <- paste0("^location_scale_", tag, "-[0-9]+\\.csv$")
    csv_files <- list.files(output_dir, pattern = pat, full.names = TRUE)
    if (length(csv_files)) {
      fit <- cmdstanr::as_cmdstan_fit(csv_files, check_diagnostics = TRUE)
    }
  }
  if (is.null(fit)) {
    append_manifest_row("hetero", "fit_cmdstan", 2027)
    fit <- do.call(mod$sample, sample_args)
  }
}

summ <- NULL
person_df <- NULL
if (!homo_only) {
  # Summaries
  summ <- as.data.frame(fit$summary())

  # Participant-level summaries (memory-light via summaries + limited draws)
  person_ids <- levels(person_factor)
  Pnames <- colnames(X)
  vars <- fit$metadata()$model_params
  has_beta_person <- any(grepl("^beta_person", vars))
  has_sigma_person <- any(grepl("^sigma_person\\[", vars))
  has_z_sigma <- any(grepl("^z_sigma\\[", vars))

  # Beta means/SDs
  if (has_beta_person) {
    beta_summ <- fit$summary(variables = "beta_person", c("mean", "sd"))
    beta_idx <- stringr::str_match(beta_summ$variable, "beta_person\\[(\\d+),(\\d+)\\]")
    beta_tbl <- tibble(
      person_idx = as.integer(beta_idx[, 2]),
      p_idx = as.integer(beta_idx[, 3]),
      mean = beta_summ$mean,
      sd = beta_summ$sd
    ) %>%
      mutate(
        person = person_ids[person_idx],
        p_name = Pnames[p_idx]
      )
  } else {
    beta_summ <- fit$summary(variables = "beta", c("mean", "sd"))
    beta_tbl <- tidyr::crossing(
      person_idx = seq_along(person_ids),
      p_idx = seq_along(Pnames)
    ) %>%
      mutate(
        person = person_ids[person_idx],
        p_name = Pnames[p_idx],
        mean = beta_summ$mean[p_idx],
        sd = beta_summ$sd[p_idx]
      )
  }

  beta_means_w <- beta_tbl %>%
    select(person, p_name, mean) %>%
    pivot_wider(names_from = p_name, values_from = mean, names_prefix = "beta_mean_")
  beta_sds_w <- beta_tbl %>%
    select(person, p_name, sd) %>%
    pivot_wider(names_from = p_name, values_from = sd, names_prefix = "beta_sd_")

  # sigma summaries
  if (has_sigma_person) {
    sigma_summ <- fit$summary(variables = "sigma_person", c("mean", "sd"))
    sigma_idx <- as.integer(stringr::str_match(sigma_summ$variable, "sigma_person\\[(\\d+)\\]")[, 2])
    sigma_df <- tibble(
      person = person_ids[sigma_idx],
      sigma_mean = sigma_summ$mean,
      sigma_sd = sigma_summ$sd
    )
  } else if ("sigma" %in% vars) {
    sigma_summ <- fit$summary(variables = "sigma", c("mean", "sd"))
    sigma_df <- tibble(
      person = person_ids,
      sigma_mean = sigma_summ$mean[1],
      sigma_sd = sigma_summ$sd[1]
    )
  } else {
    sigma_df <- tibble(person = person_ids, sigma_mean = NA_real_, sigma_sd = NA_real_)
  }

  # z_sigma outlier probability (optional).
  # NOTE: Any call to fit$draws() loads *all* draws into memory. For large models
  # this can be extremely slow and can balloon the size of saved .rds files.
  if (sigma_outlier_prob && has_z_sigma) {
    if (stan_data$N > 50000) {
      message(
        "Computing p_outlier_sigma requires fit$draws() and may be very slow ",
        "and memory-intensive for large N. Consider omitting --sigma_outlier_prob."
      )
    }
    z_mat <- fit$draws("z_sigma", format = "matrix")
    colnames(z_mat) <- person_ids[as.integer(stringr::str_match(colnames(z_mat), "z_sigma\\[(\\d+)\\]")[, 2])]
    po <- colMeans(abs(z_mat) > 1.96)
    po_df <- tibble(person = names(po), p_outlier_sigma = as.numeric(po))
  } else {
    po_df <- tibble(person = person_ids, p_outlier_sigma = NA_real_)
  }

  person_df <- beta_means_w %>%
    left_join(beta_sds_w, by = "person") %>%
    left_join(sigma_df, by = "person") %>%
    left_join(po_df, by = "person") %>%
    mutate(variance_flag = p_outlier_sigma > 0.5)

  # EARLY SAVE: persist main fit outputs before LOO/compare_homo to avoid
  # losing long runs if diagnostics fail or are interrupted.
  tryCatch(saveRDS(fit, model_path), error = function(e) {
    warning("Failed to save heteroskedastic fit to ", model_path, ": ", conditionMessage(e))
  })
  readr::write_csv(summ, file.path("reports", paste0("location_scale_", tag, "_summary.csv")))
  readr::write_csv(person_df, file.path("reports", paste0("location_scale_", tag, "_participants.csv")))
  message("Saved main model outputs to ", model_path)
}

# LOO (optional to avoid heavy log_lik memory usage)
loo_res <- NULL
if (!homo_only && do_loo) {
  loo_err <- NULL
  loo_res <- tryCatch({
    fit$loo(
      variables = loo_var,
      r_eff = TRUE,
      moment_match = moment_match,
      k_threshold = 0.7,
      cores = parallel_chains
    )
  }, error = function(e) {
    loo_err <<- e
    NULL
  })
  if (is.null(loo_res) && isTRUE(moment_match)) {
    message(
      "fit$loo failed with moment_match=TRUE (", conditionMessage(loo_err), "). ",
      "Retrying with moment_match=FALSE."
    )
    loo_err <- NULL
    loo_res <- tryCatch({
      fit$loo(
        variables = loo_var,
        r_eff = TRUE,
        moment_match = FALSE,
        k_threshold = 0.7,
        cores = parallel_chains
      )
    }, error = function(e) {
      loo_err <<- e
      NULL
    })
  }

  # Fallback: if we have grouped `log_lik_*`, `loo::loo()` is cheap. If we are
  # on an older fit without group log-lik, fall back to grouping the obs-level
  # log_lik vector in R (memory-heavy, but keeps older runs usable).
  if (is.null(loo_res)) {
    if (!is.null(loo_err)) message("fit$loo failed: ", conditionMessage(loo_err))

    vars <- fit$metadata()$model_params
    has_group_ll <- any(grepl(paste0("^", loo_var, "\\\\["), vars))

    if (has_group_ll) {
      loo_res <- tryCatch({
        ll_array <- fit$draws(loo_var, format = "matrix")
        loo::loo(ll_array, k_threshold = 0.7)
      }, error = function(e) NULL)
    } else if (loo_unit != "obs") {
      loo_res <- tryCatch({
        ll_array <- fit$draws("log_lik", format = "matrix")
        iters_per_chain <- nrow(ll_array) / chains
        if (!isTRUE(all.equal(iters_per_chain, round(iters_per_chain)))) {
          stop("Unexpected draws shape; cannot infer per-chain iterations for r_eff.")
        }
        chain_id <- rep(seq_len(chains), each = iters_per_chain)
        group <- if (loo_unit == "person") factor(stan_data$person, levels = seq_len(stan_data$J)) else factor(stan_data$site, levels = seq_len(stan_data$K))
        ll_group <- t(rowsum(t(ll_array), group = group))
        r_eff <- loo::relative_eff(exp(ll_group), chain_id = chain_id)
        loo::loo(ll_group, r_eff = r_eff, cores = 1, k_threshold = 0.7)
      }, error = function(e) NULL)
    } else if (stan_data$N <= 50000) {
      loo_res <- tryCatch({
        ll_array <- fit$draws("log_lik", format = "matrix")
        loo::loo(ll_array, k_threshold = 0.7)
      }, error = function(e) NULL)
    } else {
      message("Skipping loo::loo fallback because N=", stan_data$N, " is large.")
    }
  }
}

# Optional homogeneous-variance comparison for Bayes factor / membership test
homo_fit <- NULL
homo_loo <- NULL
bf_row <- NULL
if (compare_homo) {
  # Choose a homoskedastic Stan file that matches the chosen LOO mode.
  # - For --loo=true, use log_lik-enabled programs (obs-level or grouped).
  # - For --loo=false, use light programs with no log_lik to keep CSV outputs small.
  homo_stan_file <- if (homo_person_slope) {
    if (do_loo) {
      if (loo_unit == "obs") "location_scale_homoskedastic.stan" else "location_scale_homoskedastic_group_lik.stan"
    } else {
      "location_scale_homoskedastic_light.stan"
    }
  } else {
    if (do_loo) {
      if (loo_unit == "obs") "location_scale_homoskedastic_no_person_slope.stan" else "location_scale_homoskedastic_no_person_slope_group_lik.stan"
    } else {
      "location_scale_homoskedastic_no_person_slope_light.stan"
    }
  }
  homo_mod <- cmdstan_model(
    file.path("stan", homo_stan_file),
    dir = compile_dir,
    compile = TRUE,
    cpp_options = list(STAN_NO_PCH = TRUE),
    force_recompile = moment_match,
    compile_model_methods = moment_match
  )
  homo_output_dir <- file.path("models", paste0("cmdstan_homo_", tag))
  dir.create(homo_output_dir, showWarnings = FALSE, recursive = TRUE)
  if (reuse_homo && file.exists(homo_model_path)) {
    message("Loading existing homoskedastic fit from ", homo_model_path)
    append_manifest_row("homo", "reuse_rds", 2029)
    homo_fit <- tryCatch(readRDS(homo_model_path), error = function(e) {
      message(
        "Failed to read ", homo_model_path, " (", conditionMessage(e), "). ",
        "Falling back to CmdStan CSV outputs under ", homo_output_dir, "."
      )
      NULL
    })
    if (is.null(homo_fit)) {
      pat <- paste0("^location_scale_homo_", tag, "-[0-9]+\\.csv$")
      csv_files <- list.files(homo_output_dir, pattern = pat, full.names = TRUE)
      if (length(csv_files)) {
        homo_fit <- cmdstanr::as_cmdstan_fit(csv_files, check_diagnostics = TRUE)
      }
    }
  }
  if (is.null(homo_fit)) {
    homo_sample_args <- list(
      data = stan_data,
      seed = 2029,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter,
      iter_sampling = iter,
      thin = thin,
      output_dir = homo_output_dir,
      output_basename = paste0("location_scale_homo_", tag),
      refresh = refresh
    )
    if (!is.na(adapt_delta)) homo_sample_args$adapt_delta <- adapt_delta
    if (!is.na(max_treedepth)) homo_sample_args$max_treedepth <- max_treedepth
    if (!is.na(stepsize)) homo_sample_args$step_size <- stepsize
    if (!is.na(init_val)) homo_sample_args$init <- eval(parse(text = init_val))
    if (isTRUE(save_warmup)) homo_sample_args$save_warmup <- TRUE
    append_manifest_row("homo", "fit_cmdstan", 2029)
    homo_fit <- do.call(homo_mod$sample, homo_sample_args)
    # Save homoskedastic fit immediately in case LOO/BF fails later
    tryCatch(saveRDS(homo_fit, homo_model_path), error = function(e) {
      warning("Failed to save homoskedastic fit to ", homo_model_path, ": ", conditionMessage(e))
    })
  }

  # Always write a homoskedastic summary so RQ1 shift tables can be built without
  # requiring PSIS-LOO runs.
  readr::write_csv(
    as.data.frame(homo_fit$summary()),
    file.path("reports", paste0("location_scale_homo_", tag, "_summary.csv"))
  )

  if (do_loo && !is.null(loo_res)) {
    homo_loo <- NULL
    homo_loo_err <- NULL
    homo_loo <- tryCatch({
      homo_fit$loo(
        variables = loo_var,
        r_eff = TRUE,
        moment_match = moment_match,
        k_threshold = 0.7,
        cores = parallel_chains
      )
    }, error = function(e) {
      homo_loo_err <<- e
      NULL
    })
    if (is.null(homo_loo) && isTRUE(moment_match)) {
      message(
        "homo_fit$loo failed with moment_match=TRUE (", conditionMessage(homo_loo_err), "). ",
        "Retrying with moment_match=FALSE."
      )
      homo_loo_err <- NULL
      homo_loo <- tryCatch({
        homo_fit$loo(
          variables = loo_var,
          r_eff = TRUE,
          moment_match = FALSE,
          k_threshold = 0.7,
          cores = parallel_chains
        )
      }, error = function(e) {
        homo_loo_err <<- e
        NULL
      })
    }
    if (is.null(homo_loo)) {
      if (!is.null(homo_loo_err)) message("homo_fit$loo failed: ", conditionMessage(homo_loo_err))
      vars <- homo_fit$metadata()$model_params
      has_group_ll <- any(grepl(paste0("^", loo_var, "\\\\["), vars))
      if (has_group_ll) {
        homo_loo <- tryCatch({
          ll_array <- homo_fit$draws(loo_var, format = "matrix")
          loo::loo(ll_array, k_threshold = 0.7)
        }, error = function(e) NULL)
      } else if (loo_unit != "obs") {
        homo_loo <- tryCatch({
          ll_array <- homo_fit$draws("log_lik", format = "matrix")
          iters_per_chain <- nrow(ll_array) / chains
          if (!isTRUE(all.equal(iters_per_chain, round(iters_per_chain)))) {
            stop("Unexpected draws shape; cannot infer per-chain iterations for r_eff.")
          }
          chain_id <- rep(seq_len(chains), each = iters_per_chain)
          group <- if (loo_unit == "person") factor(stan_data$person, levels = seq_len(stan_data$J)) else factor(stan_data$site, levels = seq_len(stan_data$K))
          ll_group <- t(rowsum(t(ll_array), group = group))
          r_eff <- loo::relative_eff(exp(ll_group), chain_id = chain_id)
          loo::loo(ll_group, r_eff = r_eff, cores = 1, k_threshold = 0.7)
        }, error = function(e) NULL)
      } else if (stan_data$N <= 50000) {
        homo_loo <- tryCatch({
          ll_array <- homo_fit$draws("log_lik", format = "matrix")
          loo::loo(ll_array, k_threshold = 0.7)
        }, error = function(e) NULL)
      } else {
        message("Skipping loo::loo fallback for homo model because N=", stan_data$N, " is large.")
      }
    }
    if (!is.null(loo_res) && !is.null(homo_loo)) {
      elpd_hetero <- loo_res$estimates["elpd_loo", "Estimate"]
      elpd_homo <- homo_loo$estimates["elpd_loo", "Estimate"]
      elpd_diff <- elpd_hetero - elpd_homo
      # This is not a marginal-likelihood Bayes factor; it's the script's
      # long-standing exp(ΔELPD) heuristic for model comparison.
      bf_hetero_vs_homo <- ifelse(elpd_diff > 700, Inf, exp(elpd_diff))
      bf_homo_vs_hetero <- ifelse(elpd_diff < -700, Inf, exp(-elpd_diff))
      bf_row <- tibble(
        model = c("heterogeneous", "homogeneous"),
        elpd = c(elpd_hetero, elpd_homo),
        elpd_se = c(loo_res$estimates["elpd_loo", "SE"], homo_loo$estimates["elpd_loo", "SE"]),
        looic = c(loo_res$estimates["looic", "Estimate"], homo_loo$estimates["looic", "Estimate"]),
        log_bayes_factor_hetero_vs_homo = c(elpd_diff, -elpd_diff),
        bayes_factor_hetero_vs_homo = c(bf_hetero_vs_homo, bf_homo_vs_hetero)
      )
    }
  } else if (do_loo && is.null(loo_res)) {
    message("Skipping homoskedastic PSIS-LOO/BF because heteroskedastic LOO did not compute successfully.")
  }
}

if (!is.null(bf_row)) {
  if (!is.null(person_df)) {
    bf_val <- bf_row$bayes_factor_hetero_vs_homo[bf_row$model == "heterogeneous"][1]
    person_df <- person_df %>%
      mutate(
        variance_bayes_factor = bf_val,
        variance_membership = ifelse(bf_val > 3 | dplyr::coalesce(variance_flag, FALSE), "heterogeneous", "homogeneous")
      )
  }
}

# Save outputs (fits were saved before LOO/compare_homo to avoid bloating .rds files)
dir.create("reports", showWarnings = FALSE, recursive = TRUE)
if (!is.null(summ)) readr::write_csv(summ, file.path("reports", paste0("location_scale_", tag, "_summary.csv")))
if (!is.null(person_df)) readr::write_csv(person_df, file.path("reports", paste0("location_scale_", tag, "_participants.csv")))
if (!is.null(loo_res)) capture.output(print(loo_res), file = file.path("reports", paste0("location_scale_", tag, "_loo.txt")))
if (!is.null(homo_loo)) capture.output(print(homo_loo), file = file.path("reports", paste0("location_scale_homo_", tag, "_loo.txt")))
if (!is.null(bf_row)) readr::write_csv(bf_row, file.path("reports", paste0("location_scale_", tag, "_bf.csv")))

if (homo_only) {
  message("Homoskedastic model complete. Saved to ", homo_model_path)
} else {
  message("Location-scale model complete. Saved to ", model_path)
}
