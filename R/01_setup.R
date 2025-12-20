#!/usr/bin/env Rscript
options(warn = 1)

## Use a user-writable library to avoid admin permissions on Windows.
user_lib <- Sys.getenv(
  "R_LIBS_USER",
  unset = file.path(Sys.getenv("HOME"), "R", "win-library",
                    paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
)
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

need <- c(
  "tidyverse", "cmdstanr", "posterior", "jsonlite", "glue", "loo", "janitor",
  "logspline", "metafor"
)

ensure <- function(pkgs) {
  inst <- rownames(installed.packages())
  to_install <- pkgs[!pkgs %in% inst]
  if (length(to_install)) install.packages(to_install, lib = user_lib, repos = "https://cloud.r-project.org")
}

ensure(need)

library(cmdstanr)
ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
if (is.null(ver)) {
  message("CmdStan not found; attempting installation (may take several minutes)...")
  try(cmdstanr::install_cmdstan(), silent = TRUE)
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver)) stop("CmdStan installation failed. Install toolchain and re-run.")
}
message("CmdStan version: ", ver)

dir.create("models", showWarnings = FALSE, recursive = TRUE)
dir.create("reports", showWarnings = FALSE, recursive = TRUE)
message("Setup complete.")
