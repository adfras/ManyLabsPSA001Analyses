#!/usr/bin/env Rscript
# Identify participant-level predictors of divergence in PSA001 (Dominant/Attractive).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(forcats)
  library(broom)
})

make_psa_analysis <- function(trait, part_path, proc_path) {
  if (!file.exists(part_path)) stop("Missing participant file: ", part_path)
  if (!file.exists(proc_path)) stop("Missing processed file: ", proc_path)
  if (!file.exists("data/psa001_ind.csv")) stop("Missing data/psa001_ind.csv")
  if (!file.exists("data/psa001_cfd_faces.csv")) stop("Missing data/psa001_cfd_faces.csv")

  trait_name <- tolower(trait)

  # The participant file already stores person as the original user_id.
  part <- read_csv(part_path, show_col_types = FALSE) %>%
    mutate(user_id = as.character(person))

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

  part %>%
    left_join(demo, by = "user_id") %>%
    left_join(stim, by = "user_id") %>%
    mutate(across(where(is.character), as.factor))
}

prep_predictors <- function(df) {
  fac_vars <- c("sex", "ethnicity", "language", "region")
  num_vars <- c(
    "age", "n_trials", "mean_rating", "sd_rating",
    "prop_male_face", "mean_face_age", "prop_white_face"
  )

  # Treat missing categories explicitly so we don't drop nearly all rows.
  df <- df %>%
    mutate(across(all_of(fac_vars), ~ fct_explicit_na(as.factor(.x), na_level = "Unknown"))) %>%
    mutate(
      sex = fct_lump_min(sex, min = 10),
      ethnicity = fct_lump_min(ethnicity, min = 10),
      language = fct_lump_min(language, min = 10),
      region = fct_lump_min(region, min = 10)
    ) %>%
    mutate(across(all_of(fac_vars), droplevels))

  # Impute numeric missings with median + missingness indicator.
  miss_indicators <- character(0)
  for (v in num_vars) {
    miss_col <- paste0(v, "_missing")
    df[[miss_col]] <- ifelse(is.na(df[[v]]), 1L, 0L)
    miss_indicators <- c(miss_indicators, miss_col)
    med <- suppressWarnings(median(df[[v]], na.rm = TRUE))
    if (is.finite(med)) df[[v]][is.na(df[[v]])] <- med
  }

  # Standardize numeric predictors that have variance after imputation.
  num_keep <- num_vars[sapply(num_vars, function(v) {
    s <- suppressWarnings(sd(df[[v]], na.rm = TRUE))
    !is.na(s) && s > 0
  })]
  df <- df %>%
    mutate(across(all_of(num_keep), ~ as.numeric(scale(.x))))

  fac_keep <- fac_vars[sapply(fac_vars, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) > 1
  })]

  # Keep missingness indicators that vary.
  miss_keep <- miss_indicators[sapply(miss_indicators, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) > 1
  })]

  list(df = df, predictors = c(fac_keep, num_keep, miss_keep))
}

fit_divergence_models <- function(df, dataset_label) {
  prep <- prep_predictors(df)
  df <- prep$df
  predictors <- prep$predictors

  if (!all(c("beta_mean_X_male", "beta_sd_X_male") %in% names(df))) {
    stop("Missing beta_mean_X_male / beta_sd_X_male columns.")
  }

  # Magnitude model (weighted by slope precision).
  df_mag <- df %>%
    filter(is.finite(beta_mean_X_male), is.finite(beta_sd_X_male)) %>%
    filter(if_all(all_of(predictors), ~ is.finite(.x) | is.na(.x)))

  if (length(predictors) == 0) {
    stop("No usable predictors after filtering (all single-level or zero-variance).")
  }

  form_mag <- as.formula(paste("beta_mean_X_male ~", paste(predictors, collapse = " + ")))
  fit_mag <- lm(
    form_mag,
    data = df_mag,
    weights = 1 / (beta_sd_X_male^2)
  )
  write_csv(tidy(fit_mag), paste0("reports/psa001_", dataset_label, "_divergence_magnitude.csv"))
  message("Wrote reports/psa001_", dataset_label, "_divergence_magnitude.csv")

  # Direction model (responder vs not).
  prev <- read_csv("reports/person_prevalence_detail.csv", show_col_types = FALSE) %>%
    filter(dataset == paste0("PSA001_", dataset_label)) %>%
    mutate(person = as.character(person))

  df_dir <- df %>%
    mutate(person = as.character(person)) %>%
    left_join(prev %>% select(person, responder), by = "person") %>%
    mutate(responder = as.integer(responder))

  if (length(unique(df_dir$responder[!is.na(df_dir$responder)])) < 2) {
    message("Responder has <2 levels for PSA001_", dataset_label, "; skipping direction model.")
    return(invisible(NULL))
  }

  form_dir <- as.formula(paste("responder ~", paste(predictors, collapse = " + ")))
  fit_dir <- glm(form_dir, data = df_dir, family = binomial())
  write_csv(tidy(fit_dir), paste0("reports/psa001_", dataset_label, "_divergence_direction.csv"))
  message("Wrote reports/psa001_", dataset_label, "_divergence_direction.csv")
}

run_univariate <- function(df, dataset_label) {
  fac_vars <- c("sex", "ethnicity", "language", "region")
  num_vars <- c(
    "age", "n_trials", "mean_rating", "sd_rating",
    "prop_male_face", "mean_face_age", "prop_white_face"
  )

  df <- df %>%
    mutate(across(all_of(fac_vars), ~ fct_explicit_na(as.factor(.x), na_level = "Unknown")))

  # Magnitude: weighted LM per predictor
  mag_rows <- list()
  for (v in c(fac_vars, num_vars)) {
    d <- df %>% select(beta_mean_X_male, beta_sd_X_male, all_of(v)) %>% filter(!is.na(beta_sd_X_male))
    if (is.numeric(d[[v]])) {
      d <- d %>% mutate(x = as.numeric(scale(.data[[v]])))
    } else {
      d <- d %>% mutate(x = .data[[v]])
    }
    d <- d %>% filter(!is.na(x))
    if (nrow(d) < 20) next
    if (is.factor(d$x) && length(unique(d$x)) < 2) next
    fit <- lm(beta_mean_X_male ~ x, data = d, weights = 1 / (beta_sd_X_male^2))
    tt <- tidy(fit) %>% mutate(predictor = v, n = nrow(d), outcome = "magnitude")
    mag_rows[[length(mag_rows) + 1]] <- tt
  }
  if (length(mag_rows)) {
    write_csv(bind_rows(mag_rows),
              paste0("reports/psa001_", dataset_label, "_divergence_magnitude_univariate.csv"))
    message("Wrote reports/psa001_", dataset_label, "_divergence_magnitude_univariate.csv")
  }

  # Direction: logistic per predictor
  prev <- read_csv("reports/person_prevalence_detail.csv", show_col_types = FALSE) %>%
    filter(dataset == paste0("PSA001_", dataset_label)) %>%
    mutate(person = as.character(person))

  df_dir <- df %>%
    mutate(person = as.character(person)) %>%
    left_join(prev %>% select(person, responder), by = "person") %>%
    mutate(responder = as.integer(responder))

  if (length(unique(df_dir$responder[!is.na(df_dir$responder)])) >= 2) {
    dir_rows <- list()
    for (v in c(fac_vars, num_vars)) {
      d <- df_dir %>% select(responder, all_of(v)) %>% filter(!is.na(responder))
      if (is.numeric(d[[v]])) {
        d <- d %>% mutate(x = as.numeric(scale(.data[[v]])))
      } else {
        d <- d %>% mutate(x = .data[[v]])
      }
      d <- d %>% filter(!is.na(x))
      if (nrow(d) < 20) next
      if (is.factor(d$x) && length(unique(d$x)) < 2) next
      fit <- glm(responder ~ x, data = d, family = binomial())
      tt <- tidy(fit) %>% mutate(predictor = v, n = nrow(d), outcome = "direction")
      dir_rows[[length(dir_rows) + 1]] <- tt
    }
    if (length(dir_rows)) {
      write_csv(bind_rows(dir_rows),
                paste0("reports/psa001_", dataset_label, "_divergence_direction_univariate.csv"))
      message("Wrote reports/psa001_", dataset_label, "_divergence_direction_univariate.csv")
    }
  } else {
    message("Responder has <2 levels for PSA001_", dataset_label, "; skipping univariate direction models.")
  }
}

run_psa_divergence <- function() {
  df_dom <- make_psa_analysis(
    "dominant",
    "reports/location_scale_psa001_dominant_gender_site_ad99_participants.csv",
    "data/processed/psa001_dominant_gender.csv"
  )
  df_att <- make_psa_analysis(
    "attractive",
    "reports/location_scale_psa001_attractive_gender_site_ad99_participants.csv",
    "data/processed/psa001_attractive_gender.csv"
  )

  # Multivariate (may be singular with sparse covariates) + robust univariate runs.
  try(fit_divergence_models(df_dom, "Dominant"), silent = TRUE)
  try(fit_divergence_models(df_att, "Attractive"), silent = TRUE)
  run_univariate(df_dom, "Dominant")
  run_univariate(df_att, "Attractive")
}

run_psa_divergence()
