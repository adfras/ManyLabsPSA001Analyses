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
  df <- as.data.frame(flat, stringsAsFactors = FALSE)

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    write.table(df, path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
  } else {
    write.table(df, path, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE, append = TRUE)
  }
  invisible(TRUE)
}
