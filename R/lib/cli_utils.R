# Shared CLI helpers for this repo's Rscript entrypoints.
# Base-R only: keep these usable before loading packages.

parse_flag <- function(name, default = NULL, args = commandArgs(trailingOnly = FALSE)) {
  if (is.null(name) || !nzchar(name)) stop("parse_flag: name must be non-empty")
  if (!length(args)) return(default)

  key <- paste0("--", name)

  # --key value
  idx <- which(args == key)
  if (length(idx)) {
    i <- idx[[1]]
    if (i < length(args)) return(args[[i + 1]])
    return(default)
  }

  # --key=value
  prefix <- paste0(key, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit)) return(sub(prefix, "", hit[[1]], fixed = TRUE))

  default
}

# Parse trailing args like: --a=1 --b=foo
parse_args <- function(args, defaults = list()) {
  out <- defaults
  if (!length(args)) return(out)
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    key <- sub("^--", "", kv[1])
    out[[key]] <- kv[2]
  }
  out
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

as_bool <- parse_bool

parse_num <- function(x, default = NULL) {
  if (is.null(x) || !nzchar(as.character(x))) return(default)
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out) || !is.finite(out)) return(default)
  out
}

parse_int <- function(x, default = NULL) {
  if (is.null(x) || !nzchar(as.character(x))) return(default)
  out <- suppressWarnings(as.integer(x))
  if (is.na(out) || !is.finite(out)) return(default)
  out
}

