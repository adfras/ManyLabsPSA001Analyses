#!/usr/bin/env Rscript
# One-shot runner to rebuild summary outputs from existing fits/CV results.
# Run from repo root: Rscript run_all.R
#
# This does NOT download data and does NOT run long model fits by default.
# See README.md for end-to-end commands (including K-fold CV).

steps <- list(
  list(cmd = "Rscript R/07_rq1_shift_table.R", desc = "Shift table (effect + tau_site)"),
  list(cmd = "Rscript R/08_participants_prevalence.R --include_psa001_attractive=true --include_psa001_dominant=true", desc = "Participant prevalence summary"),
  list(cmd = "Rscript R/09_site_variance_decomp.R --include_psa001_attractive=true --include_psa001_dominant=true", desc = "Site variance decomposition (supplemental)"),
  list(cmd = "Rscript R/10_site_models_table.R", desc = "Site model comparison table (uses K-fold outputs when present)"),
  list(cmd = "Rscript R/11_site_kfold_stack.R", desc = "Site K-fold comparison stack"),
  list(cmd = "Rscript R/99_critique_resolution.R", desc = "Critique-resolution report and checks")
)

for (s in steps) {
  cat("\n==>", s$desc, "\n")
  status <- system(s$cmd)
  if (status != 0) {
    stop("Step failed: ", s$desc, " (command: ", s$cmd, ")")
  }
}

cat("\nAll steps completed successfully.\n")
