#!/usr/bin/env Rscript
# Helper to append a single-row run manifest (CSV).

append_run_manifest <- function(row, path = "reports/run_manifest.csv") {
  if (is.null(row) || !length(row)) return(invisible(FALSE))
  # Coerce to 1-row data.frame, flattening any vector values.
  flat <- lapply(row, function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) > 1) return(paste(x, collapse = ";"))
    x
  })
  df <- as.data.frame(flat, stringsAsFactors = FALSE, check.names = FALSE)

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    write.table(df, path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
    return(invisible(TRUE))
  }

  # Align columns with existing manifest; if new columns appear, rewrite with union header.
  existing_header <- names(read.csv(path, nrows = 0, check.names = FALSE))
  all_cols <- union(existing_header, names(df))

  if (!setequal(existing_header, all_cols)) {
    existing <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    for (col in setdiff(all_cols, names(existing))) existing[[col]] <- NA
    for (col in setdiff(all_cols, names(df))) df[[col]] <- NA
    existing <- existing[, all_cols, drop = FALSE]
    df <- df[, all_cols, drop = FALSE]
    out <- rbind(existing, df)
    write.table(out, path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
  } else {
    for (col in setdiff(existing_header, names(df))) df[[col]] <- NA
    df <- df[, existing_header, drop = FALSE]
    write.table(df, path, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE, append = TRUE)
  }
  invisible(TRUE)
}
