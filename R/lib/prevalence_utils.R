# Shared helper functions for participant/site mapping and effect direction.
# Requires tidyverse/readr in the calling script.

load_person_site_map <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    warning("Site map input file not found: ", path)
    return(tibble(person = character(), site = character()))
  }
  # Read only the mapping columns (works for any processed CSV with person + site columns).
  m <- readr::read_csv(
    path,
    show_col_types = FALSE,
    col_types = readr::cols_only(
      person = readr::col_character(),
      site = readr::col_character()
    )
  ) %>%
    transmute(person = as.character(person), site = as.character(site)) %>%
    distinct()
  # Defensive: if a person appears in multiple sites (shouldn't), keep the first.
  dup <- m %>% count(person, name = "n_sites") %>% filter(n_sites > 1)
  if (nrow(dup)) {
    warning("Multiple sites per person in ", path, "; taking the first site per person.")
    m <- m %>% group_by(person) %>% slice_head(n = 1) %>% ungroup()
  }
  m
}

attach_site <- function(df, site_map) {
  if ("site" %in% names(df) && any(!is.na(df$site))) return(df)
  if (!nrow(site_map)) {
    if (!"site" %in% names(df)) df$site <- NA_character_
    return(df)
  }
  df %>%
    mutate(person = as.character(.data$person)) %>%
    select(-any_of("site")) %>%
    left_join(site_map, by = "person")
}

detect_direction <- function(beta_mean, dataset, min_abs_mean = 1e-4) {
  m <- mean(beta_mean, na.rm = TRUE)
  if (!is.finite(m)) {
    warning("Non-finite mean beta for ", dataset, "; using direction=+1.")
    return(1L)
  }
  if (abs(m) < min_abs_mean) {
    warning("Mean beta near 0 for ", dataset, " (mean=", signif(m, 3), "); using direction=+1.")
    return(1L)
  }
  if (m >= 0) 1L else -1L
}

