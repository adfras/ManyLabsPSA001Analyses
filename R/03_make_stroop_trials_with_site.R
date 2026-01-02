#!/usr/bin/env Rscript
# Build Stroop (ML3) trial-level dataset with site labels.
# Inputs:
#   data/raw/ml3/StroopCleanSet.csv            (trial-level Stroop data)
#   data/raw/ml3/ML3AllSitesandmTurk.csv       (session_id + Site)
# Output:
#   data/processed/trials_stroop_ml3_with_site.csv

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path("R", "lib", "cli_utils.R"))
source(file.path("R", "lib", "file_utils.R"))

out_path <- parse_flag("out", "data/processed/trials_stroop_ml3_with_site.csv")
rt_path  <- parse_flag("stroop_rt", NA)
site_path<- parse_flag("site_map", NA)

rt_path <- rt_path %||% first_existing(c(
  "data/raw/ml3/StroopCleanSet.csv",
  "data/raw/ml3_full/osfstorage/StroopCleanSet.csv"
))
site_path <- site_path %||% first_existing(c(
  "data/raw/ml3/ML3AllSitesandmTurk.csv",
  "data/raw/ml3/ML3 PPool and MTurk Data.csv",
  "data/raw/ml3/ML3_PPool_and_MTurk_Data.csv"
))

if (!file.exists(rt_path)) {
  stop("Missing StroopCleanSet.csv. Expected at data/raw/ml3/StroopCleanSet.csv (or provide --stroop_rt).")
}
if (!file.exists(site_path)) {
  stop("Missing ML3 site map (ML3AllSitesandmTurk.csv). Expected at data/raw/ml3/ML3AllSitesandmTurk.csv (or provide --site_map).")
}

rt <- readr::read_csv(rt_path, show_col_types = FALSE)

required <- c("session_id", "trial_latency", "trial_error", "congruent")
if (!all(required %in% names(rt))) {
  stop("StroopCleanSet.csv must contain columns: ", paste(required, collapse = ", "))
}

# Clean participants following ML3 conventions (mirrors R/24_make_trial_dataset_stroop_ml3.R)
Latency <- rt %>%
  rename(
    SESSION_ID   = session_id,
    TRIAL_LATENCY = trial_latency,
    TRIAL_ERROR   = trial_error,
    CONGRUENT     = congruent
  )

myTbl <- dplyr::tbl_df(Latency) %>% dplyr::group_by(SESSION_ID)
myTbl$SUBEXCL <- 0L
myTblNoLong <- dplyr::filter(myTbl, TRIAL_LATENCY < 10000, TRIAL_LATENCY >= 0)

myFastTbl <- myTbl %>%
  summarise(FASTM = sum(TRIAL_LATENCY < 300) / dplyr::n())
isTooFast <- dplyr::filter(myFastTbl, FASTM > 0.10) %>%
  dplyr::select(SESSION_ID)
if (nrow(isTooFast) > 0) {
  myTbl[myTbl$SESSION_ID %in% isTooFast$SESSION_ID, ]$SUBEXCL <- 1L
}

Latency_filt <- myTblNoLong %>%
  dplyr::filter(SUBEXCL == 0L)

myTblNotFast <- dplyr::group_by(Latency_filt, SESSION_ID, CONGRUENT)

meanReplace <- myTblNotFast %>%
  dplyr::filter(TRIAL_ERROR == 1) %>%
  summarise(blockMean = mean(TRIAL_LATENCY) + 600)

mergeTbl <- dplyr::left_join(
  myTblNotFast,
  meanReplace,
  by = c("SESSION_ID", "CONGRUENT")
)

Correct <- dplyr::filter(mergeTbl, TRIAL_ERROR == 1)
Incorrect <- dplyr::filter(mergeTbl, TRIAL_ERROR == 0)
Incorrect$TRIAL_LATENCY <- Incorrect$blockMean

Corrected <- dplyr::bind_rows(Correct, Incorrect)
Corrected <- dplyr::filter(Corrected, TRIAL_LATENCY > 0)

trial_df <- Corrected %>%
  dplyr::transmute(
    SESSION_ID,
    CONGRUENT,
    y = log(TRIAL_LATENCY),
    person = as.character(SESSION_ID),
    X_congruent = dplyr::case_when(
      CONGRUENT == "Congruent"   ~  0.5,
      CONGRUENT == "Incongruent" ~ -0.5,
      TRUE                       ~ NA_real_
    )
  ) %>%
  tidyr::drop_na(y, person, X_congruent)

# Site mapping
site_map <- readr::read_csv(site_path, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  transmute(SESSION_ID = as.integer(session_id), site = Site)

trial_df <- trial_df %>%
  mutate(SESSION_ID = as.integer(SESSION_ID)) %>%
  left_join(site_map, by = "SESSION_ID") %>%
  mutate(site = coalesce(site, "UNKNOWN")) %>%
  select(SESSION_ID, CONGRUENT, y, person, X_congruent, site)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(trial_df, out_path)
message("Wrote ", out_path, " (", nrow(trial_df), " rows; ",
        length(unique(trial_df$person)), " participants)")
