# Shared filesystem/path helpers (base-R only).

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) && !is.na(a) && nzchar(as.character(a))) return(a)
  b
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

resolve_report_path <- function(path, reports_dir = "reports") {
  if (file.exists(path)) return(path)
  search_dir <- reports_dir
  if (!dir.exists(search_dir)) {
    # Walk up to the nearest existing parent directory (e.g., reports/kfold -> reports).
    cur <- search_dir
    while (!dir.exists(cur) && cur != "." && cur != dirname(cur)) {
      cur <- dirname(cur)
    }
    if (dir.exists(cur)) search_dir <- cur else return(path)
  }
  base <- basename(path)
  hits <- list.files(search_dir, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == base]
  if (length(hits)) {
    message("Using archived report file: ", hits[[1]])
    return(hits[[1]])
  }
  path
}

pick_path <- function(paths, label = "", hint = NULL) {
  for (p in paths) {
    if (!nzchar(p)) next
    if (file.exists(p)) return(p)
    resolved <- resolve_report_path(p, reports_dir = dirname(p))
    if (file.exists(resolved)) return(resolved)
  }
  msg <- paste0("Missing ", label, " file. Tried: ", paste(paths, collapse = ", "))
  if (!is.null(hint)) msg <- paste0(msg, "\n\n", hint)
  stop(msg)
}
