# ManyLabsPSA001Analyses

Code-only repo for the Stroop ML3 + PSA001 (Social Faces) location–scale analyses.
Raw data, processed data, model draws, and reports are intentionally *not* tracked in git.

## Purpose
This project examines how trial-level variability and sampling across random facets (participants, stimuli, and sites) shape the interpretation of psychological effects. Using hierarchical location-scale models fit at the trial level, we estimate mean effects and within-person volatility jointly while allowing effects to vary across people and, where relevant, across labs and stimuli. The objective is diagnostic rather than corrective: to determine whether apparent instability reflects mean shifts, variance shifts, or differences in precision across facets, and to evaluate models using held-out trial prediction against a homoskedastic baseline.

## Scope and decisions
- Many Smiles is excluded from the main RQs (too few repeated measures per person).
- Holdout evaluation for Stroop uses a **site-balanced subsample** for speed.
- Site K-fold comparisons are **supplemental** to RQ3 (generalization to unseen sites), not a core RQ.

## Research questions
- **RQ1 (Task comparison)**: How do within-person volatility and person-level effect heterogeneity differ between Stroop and PSA001 under the same hierarchical location-scale framework?
- **RQ2 (Direction / prevalence)**: For each task, what is the prevalence of positive, near-zero, and negative person-level effects once trial noise and person-specific variance are modeled?
- **RQ3 (Held-out prediction)**: Does modeling person-specific residual variance improve out-of-sample predictive performance (and calibration) relative to a homoskedastic mixed-effects baseline, particularly in the more variable task?

## Quick start (reproducible)
1) Install packages + CmdStan:
```bash
Rscript R/01_setup.R
```

2) Download raw data (see **Data download** below).

3) Build processed datasets:
```bash
Rscript R/03_make_stroop_trials_with_site.R
Rscript R/13_make_psa001_trait_dataset.R --trait attractive
Rscript R/13_make_psa001_trait_dataset.R --trait dominant
```

4) (Optional) Build the Stroop subsample for site K-fold holdout (supplemental):
```bash
Rscript R/05_make_stroop_subsample.R \
  --in data/processed/trials_stroop_ml3_with_site.csv \
  --out data/processed/trials_stroop_ml3_with_site_sub30.csv \
  --per_site 30 --seed 2027
```

5) Update RQ outputs (only what’s missing):
```bash
Rscript run_update_missing.R
```
To let the script refit missing models or run missing K-fold CV:
```bash
Rscript run_update_missing.R --allow_refit true --allow_kfold true
```

6) Full rebuild (assumes fits/CV already exist):
```bash
Rscript run_all.R
```

## Data download (not included in repo)
### Stroop (Many Labs 3) — OSF
Download and place under `data/raw/ml3/`:
- `https://osf.io/n8xa7/download` → `data/raw/ml3/StroopCleanSet.csv`
- `https://osf.io/bxw8j/download` → `data/raw/ml3/ML3AllSitesandmTurk.csv`

### PSA001 (PSA Social Faces) — OSF
Use the PSA001 project page and download the **Full data** files:
- OSF project: `https://osf.io/f7v3n/`
- Files needed: `psa001_ind.csv` and `psa001_cfd_faces.csv`

Place them at:
- `data/psa001_ind.csv`
- `data/psa001_cfd_faces.csv`

If you only have the smaller exploratory subset, run:
```bash
Rscript R/13_make_psa001_trait_dataset.R --trait dominant --ind path/to/subset_ind.csv --faces path/to/subset_faces.csv
```

## Key outputs (generated locally)
- RQ1: `reports/rq1_shift_table.csv`
- RQ2: `reports/person_prevalence_summary.csv` (figs under `reports/figs/`)
- RQ3: `reports/holdout_*_holdout_summary.csv` and `reports/holdout_homo_*_holdout_homo_summary.csv`
- Supplemental (site K-fold holdout): `reports/site_model_comparisons.csv` and `reports/site_kfold_comparison_stack.csv`

## Results snapshot (tracked)
A lightweight snapshot of the core RQ outputs is in:
- `results/rq_results_2026-01-02/` (latest)
- Older snapshots are kept under `results/` for reference.
