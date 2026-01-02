#!/usr/bin/env Rscript
# Hold out trials per participant, fit on w==1, and evaluate prediction on w==0.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source(file.path("R", "lib", "cli_utils.R"))

in_path <- parse_flag("in", "data/processed/psa001_attractive_gender.csv")
tag <- parse_flag("tag", tools::file_path_sans_ext(basename(in_path)))
holdout_frac <- parse_num(parse_flag("holdout_frac", "0.3"), 0.3)
seed <- as.integer(parse_flag("seed", "2027"))
make_holdout <- as_bool(parse_flag("make_holdout", "true"), TRUE)
effect_direction <- parse_num(parse_flag("effect_direction", NA), NA)
do_eval <- as_bool(parse_flag("eval", "true"), TRUE)

if (!file.exists(in_path)) stop("Missing input: ", in_path)

df <- readr::read_csv(in_path, show_col_types = FALSE)
if (!all(c("y", "person") %in% names(df))) stop("Input must have columns y and person")

holdout_path <- parse_flag("holdout_out", file.path("data", "processed", paste0(tag, "_holdout.csv")))

if (make_holdout || !"w" %in% names(df)) {
  set.seed(seed)
  df <- df %>%
    group_by(person) %>%
    mutate(
      w = if_else(row_number() %in% sample.int(n(), size = ceiling(n() * holdout_frac)), 0L, 1L)
    ) %>%
    ungroup()
  dir.create(dirname(holdout_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, holdout_path)
  message("Wrote holdout dataset: ", holdout_path)
} else {
  message("Using existing w column in input.")
}

# Allow creating the holdout file without evaluation.
if (!do_eval) {
  message("Skipping evaluation (--eval false).")
  quit(status = 0)
}

# Evaluate predictions using participant posterior means.
participants_path <- file.path("reports", paste0("location_scale_", tag, "_participants.csv"))
if (!file.exists(participants_path)) {
  stop("Missing participants file: ", participants_path, "\nRun the fit on the holdout data first.")
}

part <- readr::read_csv(participants_path, show_col_types = FALSE)
if (!"person" %in% names(part)) stop("Participants file missing 'person' column.")

preds <- grep("^X_", names(df), value = TRUE)
beta_cols <- c("beta_mean_Intercept", paste0("beta_mean_", preds))
missing_cols <- setdiff(beta_cols, names(part))
if (length(missing_cols)) {
  stop("Participants file missing columns: ", paste(missing_cols, collapse = ", "))
}

holdout <- df %>% filter(w == 0)
if (!nrow(holdout)) stop("No holdout rows (w==0). Increase holdout_frac or check w.")

extra_cols <- intersect(c("beta_mean_X_male", "beta_sd_X_male"), names(part))
tmp <- holdout %>%
  left_join(part %>% select(person, all_of(beta_cols), sigma_mean, all_of(extra_cols)), by = "person")

mu <- tmp$beta_mean_Intercept
for (p in preds) {
  mu <- mu + tmp[[paste0("beta_mean_", p)]] * tmp[[p]]
}
tmp$mu <- mu
tmp$abs_err <- abs(tmp$y - tmp$mu)
tmp$sq_err <- (tmp$y - tmp$mu)^2

rmse <- sqrt(mean(tmp$sq_err, na.rm = TRUE))
mae <- mean(tmp$abs_err, na.rm = TRUE)
cor_y <- suppressWarnings(cor(tmp$y, tmp$mu, use = "complete.obs"))

loglik <- if ("sigma_mean" %in% names(tmp)) {
  dnorm(tmp$y, tmp$mu, tmp$sigma_mean, log = TRUE)
} else {
  rep(NA_real_, nrow(tmp))
}

summary <- tibble(
  tag = tag,
  n_holdout = nrow(tmp),
  rmse = rmse,
  mae = mae,
  corr = cor_y,
  mean_loglik = mean(loglik, na.rm = TRUE)
)

summary_path <- file.path("reports", paste0("holdout_", tag, "_summary.csv"))
detail_path <- file.path("reports", paste0("holdout_", tag, "_detail.csv"))
readr::write_csv(summary, summary_path)
readr::write_csv(tmp, detail_path)
message("Wrote: ", summary_path)
message("Wrote: ", detail_path)

# Optional: responder prediction accuracy vs full-data labels.
if ("beta_mean_X_male" %in% names(part) && "beta_sd_X_male" %in% names(part)) {
  if (is.na(effect_direction)) {
    if (grepl("dominant", tag, ignore.case = TRUE)) effect_direction <- 1
    if (grepl("attractive|stroop", tag, ignore.case = TRUE)) effect_direction <- -1
  }
  if (!is.na(effect_direction) && file.exists("reports/person_prevalence_detail.csv")) {
    full <- readr::read_csv("reports/person_prevalence_detail.csv", show_col_types = FALSE)
    # match dataset by max overlap
    overlap <- full %>%
      mutate(person = as.character(person)) %>%
      group_by(dataset) %>%
      summarise(overlap = sum(person %in% as.character(part$person)), .groups = "drop") %>%
      arrange(desc(overlap))
    if (nrow(overlap)) {
      ds <- overlap$dataset[1]
      truth <- full %>% filter(dataset == ds) %>% mutate(person = as.character(person))
      pred <- part %>%
        mutate(
          person = as.character(person),
          p_in = pnorm(effect_direction * beta_mean_X_male / beta_sd_X_male),
          responder_pred = p_in >= 0.9
        ) %>%
        select(person, p_in, responder_pred)
      acc <- pred %>%
        left_join(truth %>% select(person, responder), by = "person") %>%
        filter(!is.na(responder)) %>%
        summarise(
          dataset = ds,
          n = n(),
          accuracy = mean(responder_pred == responder),
          responder_rate = mean(responder_pred),
          true_rate = mean(responder),
          .groups = "drop"
        )
      acc_path <- file.path("reports", paste0("holdout_", tag, "_responder_accuracy.csv"))
      readr::write_csv(acc, acc_path)
      message("Wrote: ", acc_path)
    }
  }
}
