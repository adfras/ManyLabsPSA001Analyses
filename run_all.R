#!/usr/bin/env Rscript
# One-shot runner to rebuild the RQ outputs from existing fits/CV results.
# Run from repo root: Rscript run_all.R
#
# This does NOT download data and does NOT run long model fits by default.
# See `run_steps.md` for the full end-to-end commands (including K-fold CV).

steps <- list(
  list(cmd = "Rscript R/12_rq1_shift_table.R", desc = "RQ1 shift table (effect + tau_site)"),
  list(cmd = "Rscript R/07_participants_prevalence.R", desc = "RQ2 participant prevalence"),
  list(cmd = "Rscript R/08_site_variance_decomp.R", desc = "RQ3 site variance decomposition"),
  list(cmd = "Rscript R/06_site_models_table.R", desc = "Model comparison table (uses K-fold outputs when present)"),
  list(cmd = "Rscript R/09_rq4_stack.R", desc = "RQ4 stacked table")
)

for (s in steps) {
  cat("\n==>", s$desc, "\n")
  status <- system(s$cmd)
  if (status != 0) {
    stop("Step failed: ", s$desc, " (command: ", s$cmd, ")")
  }
}

cat("\nAll steps completed successfully.\n")
