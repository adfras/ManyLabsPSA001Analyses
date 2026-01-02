#!/usr/bin/env Rscript
# Comprehensive covariate sweep for PSA001 divergence.
# - Builds a rich participant covariate table (demographics + rating style + RT + stimulus exposure).
# - Tests all covariates against slope magnitude and cluster membership (Responder/Opposite/Partial).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(forcats)
  library(broom)
})

source(file.path("R", "lib", "cli_utils.R"))

build_covariates <- function(trait_name) {
  ind <- read_csv("data/psa001_ind.csv", show_col_types = FALSE) %>%
    mutate(user_id = as.character(user_id)) %>%
    filter(tolower(.data$trait) == trait_name)

  demo <- ind %>%
    group_by(user_id) %>%
    summarise(
      age = first(na.omit(age)),
      sex = first(na.omit(sex)),
      ethnicity = first(na.omit(ethnicity)),
      language = first(na.omit(language)),
      country = first(na.omit(country)),
      region = first(na.omit(region)),
      lab = first(na.omit(lab)),
      n_trials = n(),
      mean_rating = mean(rating, na.rm = TRUE),
      sd_rating = sd(rating, na.rm = TRUE),
      median_rating = median(rating, na.rm = TRUE),
      iqr_rating = IQR(rating, na.rm = TRUE),
      min_rating = suppressWarnings(min(rating, na.rm = TRUE)),
      max_rating = suppressWarnings(max(rating, na.rm = TRUE)),
      range_rating = max_rating - min_rating,
      prop_low = mean(rating <= 2, na.rm = TRUE),
      prop_high = mean(rating >= 8, na.rm = TRUE),
      prop_mid = mean(rating >= 4 & rating <= 6, na.rm = TRUE),
      rating_skew = {
        m <- mean(rating, na.rm = TRUE)
        s <- sd(rating, na.rm = TRUE)
        if (is.finite(s) && s > 0) mean((rating - m)^3, na.rm = TRUE) / (s^3) else NA_real_
      },
      unique_ratings = n_distinct(rating),
      prop_integer = mean(abs(rating - round(rating)) < 1e-8, na.rm = TRUE),
      n_faces = n_distinct(stim_id),
      trials_per_face = n_trials / n_faces,
      .groups = "drop"
    )

  faces <- read_csv("data/psa001_cfd_faces.csv", show_col_types = FALSE)
  stim <- ind %>%
    left_join(faces, by = c("stim_id" = "Target")) %>%
    group_by(user_id) %>%
    summarise(
      prop_male_face = mean(Gender == "M", na.rm = TRUE),
      mean_face_age = mean(Age, na.rm = TRUE),
      sd_face_age = sd(Age, na.rm = TRUE),
      prop_race_A = mean(Race == "A", na.rm = TRUE),
      prop_race_B = mean(Race == "B", na.rm = TRUE),
      prop_race_L = mean(Race == "L", na.rm = TRUE),
      prop_race_W = mean(Race == "W", na.rm = TRUE),
      .groups = "drop"
    )

  exp <- read_csv("data/psa001_exp_data.csv", show_col_types = FALSE) %>%
    mutate(user_id = as.character(user_id))
  exp_pat <- if (trait_name == "dominant") "_Dom_" else if (trait_name == "attractive") "_Att_" else NULL
  exp_t <- if (!is.null(exp_pat)) exp %>% filter(grepl(exp_pat, exp_name, fixed = TRUE)) else exp[0, ]

  rt_feat <- exp_t %>%
    mutate(rt = as.numeric(rt)) %>%
    group_by(user_id) %>%
    summarise(
      rt_n = sum(!is.na(rt)),
      rt_missing_prop = mean(is.na(rt)),
      rt_neg_prop = mean(rt <= 0, na.rm = TRUE),
      rt_fast_prop = mean(rt > 0 & rt < 200, na.rm = TRUE),
      rt_slow_prop = mean(rt > 10000, na.rm = TRUE),
      rt_valid_n = sum(rt >= 200 & rt <= 10000, na.rm = TRUE),
      rt_mean = mean(rt[rt >= 200 & rt <= 10000], na.rm = TRUE),
      rt_median = median(rt[rt >= 200 & rt <= 10000], na.rm = TRUE),
      rt_sd = sd(rt[rt >= 200 & rt <= 10000], na.rm = TRUE),
      rt_iqr = IQR(rt[rt >= 200 & rt <= 10000], na.rm = TRUE),
      log_rt_mean = mean(log(rt[rt >= 200 & rt <= 10000]), na.rm = TRUE),
      .groups = "drop"
    )

  sess <- NULL
  if (file.exists("data/psa001_session.csv")) {
    sess <- read_csv("data/psa001_session.csv", show_col_types = FALSE) %>%
      mutate(user_id = as.character(user_id)) %>%
      group_by(user_id) %>%
      summarise(user_status = first(na.omit(user_status)), .groups = "drop")
  }

  demo %>%
    left_join(stim, by = "user_id") %>%
    left_join(rt_feat, by = "user_id") %>%
    left_join(sess, by = "user_id")
}

prep_predictors <- function(df) {
  fac_vars <- c("sex", "ethnicity", "language", "country", "region", "lab", "user_status")
  drop_vars <- c(
    "user_id", "dataset", "person", "cluster", "cluster_label", "site",
    "beta_mean", "beta_sd", "sigma",
    "p_in_direction", "p_opposite_direction", "p_near",
    "p_gt0", "p_lt0", "p_outlier",
    "effect_direction", "responder", "near_zero", "opposite",
    "high_precision"
  )
  num_vars <- setdiff(names(df), drop_vars)

  df <- df %>%
    mutate(across(all_of(fac_vars), ~ fct_na_value_to_level(as.factor(.x), "Unknown"))) %>%
    mutate(
      sex = fct_lump_min(sex, min = 10),
      ethnicity = fct_lump_min(ethnicity, min = 10),
      language = fct_lump_min(language, min = 10),
      country = fct_lump_min(country, min = 10),
      region = fct_lump_min(region, min = 10),
      lab = fct_lump_min(lab, min = 10),
      user_status = fct_lump_min(user_status, min = 10)
    )

  miss_indicators <- character(0)
  for (v in num_vars) {
    if (!is.numeric(df[[v]])) next
    miss_col <- paste0(v, "_missing")
    df[[miss_col]] <- ifelse(is.na(df[[v]]), 1L, 0L)
    miss_indicators <- c(miss_indicators, miss_col)
    med <- suppressWarnings(median(df[[v]], na.rm = TRUE))
    if (is.finite(med)) df[[v]][is.na(df[[v]])] <- med
  }

  # scale numeric predictors
  for (v in num_vars) {
    if (!is.numeric(df[[v]])) next
    s <- suppressWarnings(sd(df[[v]], na.rm = TRUE))
    if (is.finite(s) && s > 0) df[[v]] <- as.numeric(scale(df[[v]]))
  }

  list(df = df, fac_vars = fac_vars, num_vars = num_vars, miss_vars = miss_indicators)
}

univariate_glm <- function(df, outcome, predictors) {
  out_rows <- list()
  for (v in predictors) {
    if (!v %in% names(df)) next
    d <- df %>% select(outcome = all_of(outcome), x = all_of(v)) %>%
      filter(!is.na(outcome), !is.na(x))
    if (nrow(d) < 20) next
    if (is.factor(d$x) && length(unique(d$x)) < 2) next
    if (is.numeric(d$x) && sd(d$x, na.rm = TRUE) == 0) next
    fit <- glm(outcome ~ x, data = d, family = binomial())
    tt <- tidy(fit) %>%
      mutate(predictor = v, n = nrow(d), levels = if (is.factor(d$x)) length(unique(d$x)) else NA_integer_)
    out_rows[[length(out_rows) + 1]] <- tt
  }
  if (!length(out_rows)) return(NULL)
  bind_rows(out_rows)
}

univariate_lm <- function(df, outcome, predictors, weights = NULL) {
  out_rows <- list()
  for (v in predictors) {
    if (!v %in% names(df)) next
    d <- df %>% select(outcome = all_of(outcome), x = all_of(v)) %>%
      filter(!is.na(outcome), !is.na(x))
    if (nrow(d) < 20) next
    if (is.factor(d$x) && length(unique(d$x)) < 2) next
    if (is.numeric(d$x) && sd(d$x, na.rm = TRUE) == 0) next
    if (is.null(weights)) {
      fit <- lm(outcome ~ x, data = d)
    } else {
      fit <- lm(outcome ~ x, data = d, weights = weights[dplyr::row_number()])
    }
    tt <- tidy(fit) %>%
      mutate(predictor = v, n = nrow(d), levels = if (is.factor(d$x)) length(unique(d$x)) else NA_integer_)
    out_rows[[length(out_rows) + 1]] <- tt
  }
  if (!length(out_rows)) return(NULL)
  bind_rows(out_rows)
}

run_trait <- function(dataset_label, trait_name) {
  cluster_path <- paste0("reports/psa001_", trait_name, "_clusters.csv")
  if (!file.exists(cluster_path)) {
    stop("Missing cluster file: ", cluster_path, ". Run R/15_psa001_cluster_analysis.R first.")
  }

  clusters <- read_csv(cluster_path, show_col_types = FALSE)
  covs <- build_covariates(trait_name)
  df <- clusters %>%
    mutate(user_id = as.character(user_id)) %>%
    left_join(covs, by = "user_id")

  prep <- prep_predictors(df)
  df <- prep$df
  predictors <- unique(c(prep$fac_vars, prep$num_vars, prep$miss_vars))

  # Save covariate table
  write_csv(df, paste0("reports/psa001_", trait_name, "_covariates.csv"))

  # Magnitude sweep (beta_mean)
  mag <- univariate_lm(df, "beta_mean", predictors)
  if (!is.null(mag)) {
    mag <- mag %>% mutate(outcome = "beta_mean")
    write_csv(mag, paste0("reports/psa001_", trait_name, "_covariate_sweep_magnitude.csv"))
  }

  # Cluster label sweeps
  if ("cluster_label" %in% names(df)) {
    for (lab in c("Responder", "Opposite", "Partial")) {
      df$target <- df$cluster_label == lab
      if (length(unique(df$target[!is.na(df$target)])) < 2) next
      out <- univariate_glm(df, "target", predictors)
      if (!is.null(out)) {
        out <- out %>% mutate(outcome = lab)
        write_csv(out, paste0("reports/psa001_", trait_name, "_covariate_sweep_", tolower(lab), ".csv"))
      }
    }
  }
}

run_all <- function() {
  run_trait("PSA001_Dominant", "dominant")
  run_trait("PSA001_Attractive", "attractive")
}

run_all()
