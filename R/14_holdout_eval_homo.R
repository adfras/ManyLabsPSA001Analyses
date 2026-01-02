#!/usr/bin/env Rscript
# Evaluate holdout predictions for a homoskedastic fit using summary means.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

source(file.path("R", "lib", "cli_utils.R"))

in_path <- parse_flag("in", "data/processed/psa001_attractive_holdout.csv")
tag <- parse_flag("tag", tools::file_path_sans_ext(basename(in_path)))
summary_path <- parse_flag("summary", file.path("reports", paste0("location_scale_homo_", tag, "_summary.csv")))

if (!file.exists(in_path)) stop("Missing input: ", in_path)

df <- read_csv(in_path, show_col_types = FALSE)
if (!"w" %in% names(df)) stop("Input missing w column (holdout indicator).")

preds <- grep("^X_", names(df), value = TRUE)
X <- as.matrix(cbind(Intercept = 1, if (length(preds)) df[, preds, drop = FALSE] else NULL))
P <- ncol(X)

person_factor <- factor(df$person)
site_factor <- factor(df$site)
J <- nlevels(person_factor)
K <- nlevels(site_factor)
person_idx <- as.integer(person_factor)
site_idx <- as.integer(site_factor)

summ <- NULL
if (file.exists(summary_path)) {
  summ <- read_csv(summary_path, show_col_types = FALSE)
}
if (is.null(summ) || !all(c("variable", "mean") %in% names(summ))) {
  # Build summary from CmdStan CSV if needed.
  suppressPackageStartupMessages(library(cmdstanr))
  csv_dir <- file.path("models", paste0("cmdstan_homo_", tag))
  pat <- paste0("^location_scale_homo_", tag, "-[0-9]+\\.csv$")
  csv_files <- list.files(csv_dir, pattern = pat, full.names = TRUE)
  if (!length(csv_files)) stop("No CmdStan CSVs found at ", csv_dir)
  fit <- cmdstanr::as_cmdstan_fit(csv_files, check_diagnostics = FALSE)
  summ <- as.data.frame(fit$summary())
  readr::write_csv(summ, summary_path)
}

get_vec <- function(name, len) {
  out <- rep(0, len)
  for (i in seq_len(len)) {
    v <- summ$mean[summ$variable == sprintf("%s[%d]", name, i)]
    if (length(v)) out[i] <- v[1]
  }
  out
}

get_mat <- function(name, nrow, ncol) {
  out <- matrix(0, nrow = nrow, ncol = ncol)
  pat <- paste0("^", name, "\\\\[(\\d+),(\\d+)\\\\]$")
  hits <- summ %>% filter(str_detect(variable, pat))
  if (nrow(hits)) {
    idx <- str_match(hits$variable, pat)
    i <- as.integer(idx[, 2])
    j <- as.integer(idx[, 3])
    out[cbind(i, j)] <- hits$mean
  }
  out
}

beta <- get_vec("beta", P)
tau_site <- get_vec("tau_site", P)
tau_person <- get_vec("tau_person", P)
site_raw <- get_mat("site_raw", K, P)
person_raw <- get_mat("person_raw", J, P)
sigma <- summ$mean[summ$variable == "sigma"][1]
if (is.na(sigma)) sigma <- mean(summ$mean[summ$variable == "sigma[1]"], na.rm = TRUE)

site_eff <- sweep(site_raw, 2, tau_site, "*")
person_eff <- sweep(person_raw, 2, tau_person, "*")

mu <- as.numeric(X %*% beta)
mu <- mu + rowSums(X * site_eff[site_idx, , drop = FALSE])
mu <- mu + rowSums(X * person_eff[person_idx, , drop = FALSE])

holdout <- df$w == 0
if (!any(holdout)) stop("No holdout rows (w==0).")

mu_h <- mu[holdout]
y_h <- df$y[holdout]

rmse <- sqrt(mean((y_h - mu_h)^2, na.rm = TRUE))
mae <- mean(abs(y_h - mu_h), na.rm = TRUE)
cor_y <- suppressWarnings(cor(y_h, mu_h, use = "complete.obs"))
loglik <- if (is.finite(sigma)) dnorm(y_h, mu_h, sigma, log = TRUE) else rep(NA_real_, length(y_h))

summary <- tibble(
  tag = tag,
  n_holdout = sum(holdout),
  rmse = rmse,
  mae = mae,
  corr = cor_y,
  mean_loglik = mean(loglik, na.rm = TRUE)
)

summary_path_out <- file.path("reports", paste0("holdout_homo_", tag, "_summary.csv"))
detail_path_out <- file.path("reports", paste0("holdout_homo_", tag, "_detail.csv"))

write_csv(summary, summary_path_out)
write_csv(tibble(mu = mu_h, y = y_h), detail_path_out)

message("Wrote: ", summary_path_out)
message("Wrote: ", detail_path_out)
