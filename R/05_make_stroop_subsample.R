#!/usr/bin/env Rscript
# Create a smaller, site-balanced Stroop dataset for fast LOO / hetero-vs-homo comparisons.
#
# The full Stroop trials file (~210k rows) makes PSIS-LOO and homoskedastic comparisons
# expensive. This script subsamples participants within each site (stratified by site),
# then keeps all trials for the sampled participants.
#
# Usage:
#   Rscript R/05_make_stroop_subsample.R \
#     --in data/processed/trials_stroop_ml3_with_site.csv \
#     --out data/processed/trials_stroop_ml3_with_site_sub.csv \
#     --per_site 20 \
#     --seed 2027

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path("R", "lib", "cli_utils.R"))

in_path <- parse_flag("in", "data/processed/trials_stroop_ml3_with_site.csv")
out_path <- parse_flag("out", "data/processed/trials_stroop_ml3_with_site_sub.csv")
per_site <- as.integer(parse_flag("per_site", 20))
seed <- as.integer(parse_flag("seed", 2027))

if (!file.exists(in_path)) stop("Input file not found: ", in_path)
if (is.na(per_site) || per_site < 1) stop("--per_site must be >= 1")
if (is.na(seed)) stop("--seed must be an integer")

df <- readr::read_csv(in_path, show_col_types = FALSE)
if (!all(c("person", "site") %in% names(df))) stop("Input must have columns person and site")

person_site <- df %>% distinct(person, site)
sites <- sort(unique(person_site$site))

set.seed(seed)
sampled_people <- split(person_site, person_site$site) |>
  lapply(function(d) {
    people <- d$person
    sample(people, size = min(per_site, length(people)), replace = FALSE)
  }) |>
  unlist(use.names = FALSE) |>
  unique()

out_df <- df %>% filter(person %in% sampled_people)

site_counts <- out_df %>% distinct(person, site) %>% count(site, sort = TRUE)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(out_df, out_path)

message("Wrote ", out_path)
message("Sites in input: ", length(sites), "; sites in output: ", nrow(site_counts))
message("Sampled participants: ", length(unique(out_df$person)), " (target per site: ", per_site, ")")
message("Rows: ", nrow(out_df))
message("Participants per site:")
print(site_counts, n = Inf)
