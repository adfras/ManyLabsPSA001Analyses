# ManyLabsPSA001Analyses

Code-only repo for the Stroop ML3 + PSA001 (Social Faces) location–scale analyses.
Raw data, processed data, model draws, and reports are intentionally *not* tracked in git.

## Scope and decisions
- Many Smiles is excluded from the main RQs (too few repeated measures per person).
- RQ4 uses **site K-fold CV** (PSIS-LOO is unreliable for these hierarchical models).
- Stroop K-fold runs on a **site-balanced subsample** for speed.

## Research questions
- **RQ1**: Does modeling trial-level precision change effect size and between-site heterogeneity?
- **RQ2**: How much person-level heterogeneity remains after modeling precision?
- **RQ3**: Are site differences explained by the mix of individual trajectories and their precision?
- **RQ4**: Does the location–scale model predict better than a homoskedastic mixed model for *unseen sites*?

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

4) (Recommended) Build the Stroop subsample for K-fold CV:
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
- RQ3: `reports/site_variance_decomposition.csv` and `reports/site_level_mix_precision.csv`
- RQ4: `reports/site_model_comparisons.csv` and `reports/rq4_comparison_stack.csv`

## Results snapshot (tracked)
A lightweight snapshot of the core RQ outputs is in:
- `results/rq_results_2025-12-20/`

## Track record
See `TRACK_RECORD.md` for setbacks, decisions, and fixes.
