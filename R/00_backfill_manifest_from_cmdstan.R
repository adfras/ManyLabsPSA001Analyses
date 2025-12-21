#!/usr/bin/env Rscript
# Backfill run_manifest.csv from existing CmdStan CSV headers.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

if (file.exists("R/00_run_manifest.R")) source("R/00_run_manifest.R")

parse_header <- function(path, max_lines = 200) {
  lines <- readLines(path, n = max_lines, warn = FALSE)
  header <- lines[grepl("^#", lines)]
  header <- sub("^#\\s*", "", header)
  kv <- header[grepl("=", header)]
  out <- list()
  for (line in kv) {
    parts <- strsplit(line, "=", fixed = TRUE)[[1]]
    key <- trimws(parts[1])
    val <- trimws(paste(parts[-1], collapse = "="))
    out[[key]] <- val
  }
  out
}

infer_tag <- function(fname) {
  tag <- fname
  tag <- sub("^-*", "", tag)
  tag <- sub("\\.csv$", "", tag)
  tag <- sub("-[0-9]+$", "", tag)
  if (startsWith(tag, "location_scale_homo_")) {
    tag <- sub("^location_scale_homo_", "", tag)
  } else if (startsWith(tag, "location_scale_")) {
    tag <- sub("^location_scale_", "", tag)
  }
  tag
}

infer_model <- function(fname) {
  if (startsWith(fname, "location_scale_homo_")) return("homo")
  if (startsWith(fname, "location_scale_")) return("hetero")
  return("unknown")
}

find_chain1 <- function(dir_path) {
  list.files(dir_path, pattern = "-1\\.csv$", full.names = TRUE)
}

csv_dirs <- list.dirs("models", recursive = TRUE, full.names = TRUE)
csv_dirs <- csv_dirs[grepl("cmdstan", csv_dirs)]

chain1_files <- unlist(lapply(csv_dirs, find_chain1))
chain1_files <- chain1_files[file.exists(chain1_files)]
if (!length(chain1_files)) {
  message("No CmdStan chain-1 CSVs found. Nothing to backfill.")
  quit(status = 0)
}

# De-duplicate by output_basename (in case multiple chain-1 files exist).
basenames <- basename(chain1_files)
output_basename <- sub("-1\\.csv$", "", basenames)
chain1_files <- chain1_files[!duplicated(output_basename)]

manifest_path <- "reports/run_manifest.csv"
existing <- NULL
if (file.exists(manifest_path)) {
  existing <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
}

for (path in chain1_files) {
  fname <- basename(path)
  out_base <- sub("-1\\.csv$", "", fname)
  out_dir <- dirname(path)

  # Skip if already recorded (by output_basename + output_dir if available).
  if (!is.null(existing) && nrow(existing) > 0) {
    if ("output_basename" %in% names(existing) && "output_dir" %in% names(existing)) {
      if (any(existing$output_basename == out_base & existing$output_dir == out_dir, na.rm = TRUE)) next
    }
  }

  h <- parse_header(path)
  model <- infer_model(fname)
  tag <- infer_tag(fname)

  # Count chains for this output basename.
  chain_pattern <- paste0("^", gsub("\\.", "\\\\.", out_base), "-[0-9]+\\.csv$")
  chains <- length(list.files(out_dir, pattern = chain_pattern, full.names = FALSE))

  append_run_manifest(list(
    timestamp = h[["start_datetime"]],
    script = "cmdstan_header",
    action = "backfill_cmdstan",
    model = model,
    tag = tag,
    data = h[["file"]],
    output_dir = out_dir,
    output_basename = out_base,
    cmdstan_model = h[["model"]],
    N = NA,
    J = NA,
    K = NA,
    P = NA,
    likelihood = NA,
    use_site_effects = NA,
    homo_only = NA,
    compare_homo_requested = NA,
    compare_homo_effective = NA,
    loo_requested = NA,
    loo_effective = NA,
    loo_unit = NA,
    loo_disabled_reason = NA,
    moment_match = NA,
    force_loo = NA,
    loo_max_n = NA,
    chains = chains,
    parallel_chains = NA,
    iter_warmup = as.integer(h[["num_warmup"]]),
    iter_sampling = as.integer(h[["num_samples"]]),
    thin = as.integer(h[["thin"]]),
    adapt_delta = as.numeric(h[["delta"]]),
    max_treedepth = as.integer(h[["max_depth"]]),
    stepsize = as.numeric(h[["stepsize"]]),
    init = h[["init"]],
    seed = as.integer(h[["seed"]]),
    refresh = as.integer(h[["refresh"]]),
    save_warmup = h[["save_warmup"]],
    reuse_hetero = NA,
    reuse_homo = NA,
    cmd = NA
  ), path = manifest_path)
}

message("Backfill complete: ", manifest_path)
