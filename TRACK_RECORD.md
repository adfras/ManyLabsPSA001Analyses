# Track record: setbacks, decisions, and fixes

This file captures the practical issues we hit while running the ManyLabs (Stroop ML3 + PSA001) location–scale pipeline,
plus the concrete decisions and fixes we used to keep the workflow reproducible and fast.

## Key decisions
- Many Smiles is excluded from the main RQs because it has too few repeated measures per person to support precision/heterogeneity claims.
- RQ4 uses site K-fold CV rather than PSIS-LOO (Pareto-k issues for hierarchical random-effects models).
- Stroop K-fold runs on a site-balanced subsample (default 30 per site) to keep runtime manageable.

## Recurring setbacks and solutions
- Runs appear "stuck" at iteration 1 or 0%.
  - Cause: old CmdStan/Rscript processes still running or multiple overlapping jobs.
  - Fix: check active processes, kill stale ones, then rerun.

- K-fold CV takes forever.
  - Cause: each fold refits two models (heteroskedastic + homoskedastic).
  - Fix: reduce folds (e.g., 5), use fewer chains/iters for CV, and run on the Stroop subsample.

- Post-fit steps are unexpectedly slow.
  - Cause: summarizing all draws/parameters (large log_lik and random-effect blocks).
  - Fix: summarize only core parameters (beta, tau_site, tau_person, sigma, etc.).

- Too many reruns / "start over" behavior.
  - Fix: added `run_update_missing.R` so only missing outputs are built; optional refits via `--allow_refit` and `--allow_kfold`.

- Divergences in homoskedastic fits (esp. PSA001).
  - Fix: increase `adapt_delta` and `max_treedepth` (e.g., 0.995 / 20) for the baseline model.

## Practical runtime guidance
- Prefer running long fits from PowerShell/Terminal instead of the RStudio console.
- Match `--parallel_chains` to `--chains` to use available CPU cores.
- Use `--resume true` for K-fold jobs to avoid losing progress.

## Outputs used in writeups
- RQ1: `reports/rq1_shift_table.csv`
- RQ2: `reports/person_prevalence_summary.csv`
- RQ3: `reports/site_variance_decomposition.csv` and `reports/site_level_mix_precision.csv`
- RQ4: `reports/site_model_comparisons.csv` and `reports/rq4_comparison_stack.csv`
