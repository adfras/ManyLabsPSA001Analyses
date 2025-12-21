# Track record: setbacks, decisions, and fixes

This log records the methodological and operational issues encountered while running the ManyLabs (Stroop ML3 + PSA001)
location–scale pipeline, along with the rationale for the design choices and the practical fixes adopted.

## Decisions and rationale
- **Exclude Many Smiles from main RQs.** The dataset has too few repeated measures per participant to support the
  person-level precision/heterogeneity claims required by the location–scale model.
- **Use site K-fold CV for RQ4 (generalization to unseen sites).** PSIS-LOO produced unreliable Pareto-k diagnostics in
  hierarchical random-effects models and does not match the *site-level* generalization target. Site K-fold directly
  tests prediction on held-out sites and avoids LOO instability.
- **Run Stroop K-fold on a site-balanced subsample (default 30/site).** Full trial-level CV is computationally
  prohibitive; a balanced subsample preserves site structure while keeping runtime tractable.

## Model-fitting and pipeline setbacks (with resolutions)
- **MCMC runs appear stalled at iteration 1/0%.**
  - *Cause:* orphaned CmdStan/Rscript processes or overlapping jobs holding resources.
  - *Resolution:* enumerate active Rscript/CmdStan processes, terminate stale ones, then rerun.
- **K-fold CV runtime is excessive.**
  - *Cause:* each fold fits two models (heteroskedastic + homoskedastic), multiplying total fits.
  - *Resolution:* reduce folds when needed (e.g., 5), lower chains/iters for CV-only runs, and use the Stroop subsample.
- **Post-fit summarization is unexpectedly slow.**
  - *Cause:* summarizing all draws/parameters (large `log_lik` vectors and random-effect blocks).
  - *Resolution:* summarize only core parameters needed for RQs (e.g., `beta`, `tau_site`, `tau_person`, `sigma`).
- **Frequent reruns due to missing outputs.**
  - *Resolution:* added `run_update_missing.R` to build only missing artifacts; optional refits via
    `--allow_refit` and `--allow_kfold`.
- **Divergences in homoskedastic baseline fits (notably PSA001).**
  - *Resolution:* increase `adapt_delta` and `max_treedepth` (e.g., 0.995 / 20) for the baseline model.

## Operational guidance
- Run long fits from PowerShell/Terminal rather than the RStudio console to avoid UI interruptions.
- Match `--parallel_chains` to `--chains` to utilize available cores.
- Use `--resume true` for K-fold jobs to preserve progress after interruptions.

## Outputs used in writeups
- RQ1: `reports/rq1_shift_table.csv`
- RQ2: `reports/person_prevalence_summary.csv`
- RQ3: `reports/site_variance_decomposition.csv` and `reports/site_level_mix_precision.csv`
- RQ4: `reports/site_model_comparisons.csv` and `reports/rq4_comparison_stack.csv`
