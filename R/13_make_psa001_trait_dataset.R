#!/usr/bin/env Rscript
# Prepare a PSA001 (Social Faces) trait dataset in the "ManyLabsAnalyses" format.
#
# Input (placed by user in ./data):
# - data/psa001_ind.csv         (trial-level ratings: user_id x stim_id)
# - data/psa001_cfd_faces.csv   (stimulus metadata: Target, Race, Gender, Age)
#
# Output (default):
# - data/processed/psa001_<trait>_gender.csv
#
# The output format is compatible with:
# - R/00_counts_table.R
# - R/04_fit_stroop_location_scale.R
#
# Usage:
#   Rscript R/13_make_psa001_trait_dataset.R --trait dominant
#   Rscript R/13_make_psa001_trait_dataset.R --trait trustworthy --out data/processed/psa001_trustworthy_gender.csv

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path("R", "lib", "cli_utils.R"))

ind_path <- parse_flag("ind", "data/psa001_ind.csv")
faces_path <- parse_flag("faces", "data/psa001_cfd_faces.csv")
trait_name <- tolower(parse_flag("trait", "dominant"))
site_var <- tolower(parse_flag("site_var", "lab")) # lab|country|region
out_path <- parse_flag("out", file.path("data", "processed", paste0("psa001_", trait_name, "_gender.csv")))

if (!file.exists(ind_path)) stop("Missing ind file: ", ind_path)
if (!file.exists(faces_path)) stop("Missing faces file: ", faces_path)

ind <- readr::read_csv(ind_path, show_col_types = FALSE)
need <- c("user_id", "lab", "country", "region", "trait", "stim_id", "rating")
miss <- setdiff(need, names(ind))
if (length(miss)) stop("psa001_ind is missing columns: ", paste(miss, collapse = ", "))

faces <- readr::read_csv(faces_path, show_col_types = FALSE)
need_faces <- c("Target", "Race", "Gender", "Age")
miss_faces <- setdiff(need_faces, names(faces))
if (length(miss_faces)) stop("psa001_cfd_faces is missing columns: ", paste(miss_faces, collapse = ", "))

traits_avail <- sort(unique(tolower(ind$trait)))
if (!trait_name %in% traits_avail) {
  stop("Unknown --trait '", trait_name, "'. Available: ", paste(traits_avail, collapse = ", "))
}

site_col <- dplyr::case_when(
  site_var == "lab" ~ "lab",
  site_var == "country" ~ "country",
  site_var == "region" ~ "region",
  TRUE ~ NA_character_
)
if (is.na(site_col)) stop("--site_var must be one of: lab, country, region")

df <- ind %>%
  mutate(trait = tolower(.data$trait)) %>%
  filter(.data$trait == .env$trait_name) %>%
  left_join(faces, by = c("stim_id" = "Target"))

if (anyNA(df$Gender)) stop("Join failed: missing Gender after joining faces metadata.")

out <- df %>%
  transmute(
    y = as.numeric(.data$rating),
    person = as.character(.data$user_id),
    site = as.character(.data[[site_col]]),
    X_male = if_else(.data$Gender == "M", 0.5, -0.5),
    stim_id = as.character(.data$stim_id),
    face_race = as.character(.data$Race),
    face_age = as.numeric(.data$Age)
  ) %>%
  filter(is.finite(.data$y), !is.na(.data$person), !is.na(.data$site))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(out, out_path)

message(
  "Wrote ", out_path, "\n",
  "  trait=", trait_name, "\n",
  "  N=", nrow(out), " rows\n",
  "  J=", n_distinct(out$person), " participants\n",
  "  K=", n_distinct(out$site), " sites (", site_var, ")\n",
  "  trials/person (median)=", median(table(out$person)), "\n"
)
