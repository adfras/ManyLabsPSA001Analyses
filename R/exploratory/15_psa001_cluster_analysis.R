#!/usr/bin/env Rscript
# Cluster PSA001 participants by posterior slope profiles and test predictors.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(forcats)
  library(broom)
})

source(file.path("R", "lib", "cli_utils.R"))

as_int <- function(x, default) {
  out <- suppressWarnings(as.integer(x))
  if (is.na(out)) default else out
}

build_covariates <- function(trait_name) {
  ind <- read_csv("data/psa001_ind.csv", show_col_types = FALSE) %>%
    mutate(user_id = as.character(user_id)) %>%
    filter(tolower(.data$trait) == trait_name)

  demo <- ind %>%
    group_by(user_id) %>%
    summarise(
      age = first(na.omit(age)),
      sex = first(na.omit(sex)),
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
      prop_white_face = mean(Race == "W", na.rm = TRUE),
      .groups = "drop"
    )

  # RT features from exp_data (per trait).
  exp <- read_csv("data/psa001_exp_data.csv", show_col_types = FALSE) %>%
    mutate(user_id = as.character(user_id))
  exp_pat <- if (trait_name == "dominant") "_Dom_" else if (trait_name == "attractive") "_Att_" else NULL
  exp_t <- if (!is.null(exp_pat)) exp %>% filter(grepl(exp_pat, exp_name, fixed = TRUE)) else exp[0, ]
  rt_feat <- exp_t %>%
    mutate(rt = as.numeric(rt)) %>%
    group_by(user_id) %>%
    summarise(
      rt_n = sum(!is.na(rt)),
      rt_neg_prop = mean(rt <= 0, na.rm = TRUE),
      rt_fast_prop = mean(rt > 0 & rt < 200, na.rm = TRUE),
      rt_slow_prop = mean(rt > 10000, na.rm = TRUE),
      rt_valid_n = sum(rt >= 200 & rt <= 10000, na.rm = TRUE),
      rt_mean = mean(rt[rt >= 200 & rt <= 10000], na.rm = TRUE),
      rt_median = median(rt[rt >= 200 & rt <= 10000], na.rm = TRUE),
      rt_sd = sd(rt[rt >= 200 & rt <= 10000], na.rm = TRUE),
      log_rt_mean = mean(log(rt[rt >= 200 & rt <= 10000]), na.rm = TRUE),
      .groups = "drop"
    )

  demo %>%
    left_join(stim, by = "user_id") %>%
    left_join(rt_feat, by = "user_id")
}

prep_predictors <- function(df) {
  fac_vars <- c("sex", "region", "country", "language")
  num_vars <- c(
    "age", "n_trials", "mean_rating", "sd_rating", "median_rating", "iqr_rating",
    "range_rating", "prop_low", "prop_high", "prop_mid", "rating_skew",
    "n_faces", "trials_per_face",
    "prop_male_face", "mean_face_age", "prop_white_face",
    "rt_n", "rt_neg_prop", "rt_fast_prop", "rt_slow_prop", "rt_valid_n",
    "rt_mean", "rt_median", "rt_sd", "log_rt_mean"
  )

  df <- df %>%
    mutate(across(all_of(fac_vars), ~ fct_na_value_to_level(as.factor(.x), "Unknown"))) %>%
    mutate(
      sex = fct_lump_min(sex, min = 10),
      region = fct_lump_min(region, min = 10),
      country = fct_lump_min(country, min = 10),
      language = fct_lump_min(language, min = 10)
    )

  miss_indicators <- character(0)
  for (v in num_vars) {
    miss_col <- paste0(v, "_missing")
    df[[miss_col]] <- ifelse(is.na(df[[v]]), 1L, 0L)
    miss_indicators <- c(miss_indicators, miss_col)
    med <- suppressWarnings(median(df[[v]], na.rm = TRUE))
    if (is.finite(med)) df[[v]][is.na(df[[v]])] <- med
  }

  num_keep <- num_vars[sapply(num_vars, function(v) {
    s <- suppressWarnings(sd(df[[v]], na.rm = TRUE))
    !is.na(s) && s > 0
  })]
  df <- df %>%
    mutate(across(all_of(num_keep), ~ as.numeric(scale(.x))))

  miss_keep <- miss_indicators[sapply(miss_indicators, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) > 1
  })]

  list(df = df, predictors = c(fac_vars, num_keep, miss_keep))
}

label_clusters <- function(stats) {
  stats %>%
    mutate(
      cluster_label = case_when(
        mean_p_in >= 0.6 ~ "Responder",
        mean_p_opp >= 0.6 ~ "Opposite",
        TRUE ~ "Partial"
      )
    )
}

run_cluster <- function(dataset_label, trait_name, k = 3) {
  prev <- read_csv("reports/person_prevalence_detail.csv", show_col_types = FALSE) %>%
    filter(dataset == dataset_label) %>%
    mutate(user_id = as.character(person))

  feats <- c("beta_mean", "beta_sd", "sigma", "p_in_direction", "p_opposite_direction", "p_near")
  X <- prev[, feats]
  for (v in feats) {
    med <- suppressWarnings(median(X[[v]], na.rm = TRUE))
    X[[v]][is.na(X[[v]])] <- med
  }
  Xs <- scale(as.matrix(X))
  set.seed(2027)
  km <- kmeans(Xs, centers = k, nstart = 50)
  prev$cluster <- as.integer(km$cluster)

  cluster_stats <- prev %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      mean_beta = mean(beta_mean, na.rm = TRUE),
      mean_p_in = mean(p_in_direction, na.rm = TRUE),
      mean_p_opp = mean(p_opposite_direction, na.rm = TRUE),
      mean_p_near = mean(p_near, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    label_clusters()

  prev <- prev %>%
    left_join(cluster_stats %>% select(cluster, cluster_label), by = "cluster")

  write_csv(prev, paste0("reports/psa001_", trait_name, "_clusters.csv"))
  write_csv(cluster_stats, paste0("reports/psa001_", trait_name, "_cluster_summary.csv"))

  covs <- build_covariates(trait_name)
  df <- prev %>% left_join(covs, by = "user_id")

  prep <- prep_predictors(df)
  df <- prep$df
  predictors <- prep$predictors

  run_logit <- function(target_label, suffix) {
    if (!any(df$cluster_label == target_label, na.rm = TRUE)) return(invisible(NULL))
    df$target <- df$cluster_label == target_label
    if (length(unique(df$target[!is.na(df$target)])) < 2) return(invisible(NULL))
    form <- as.formula(paste("target ~", paste(predictors, collapse = " + ")))
    fit <- glm(form, data = df, family = binomial())
    out <- tidy(fit)
    write_csv(out, paste0("reports/psa001_", trait_name, "_cluster_predict_", suffix, ".csv"))
  }

  run_logit("Responder", "responder")
  run_logit("Opposite", "opposite")
  run_logit("Partial", "partial")

  message("Wrote PSA001 cluster outputs for ", trait_name)
}

run_all <- function() {
  k <- as_int(parse_flag("k", "3"), 3)
  run_cluster("PSA001_Dominant", "dominant", k = k)
  run_cluster("PSA001_Attractive", "attractive", k = k)
}

run_all()
