# ManyLabsPSA001Analyses

Code-only repo for the Stroop ML3 + PSA001 location–scale analyses.
Raw data, processed data, model draws, and reports are intentionally *not* tracked in git.

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

4) Update RQ outputs (only what’s missing):
```bash
Rscript run_update_missing.R
```
If you want the script to refit missing models or run missing k-fold CV:
```bash
Rscript run_update_missing.R --allow_refit true --allow_kfold true
```

5) Full rebuild (assumes fits/CV already exist):
```bash
Rscript run_all.R
```

Refer to `analysis_plan.md` for RQ definitions and recommended settings.

## Data download (not included in repo)
See `data/README.md` for the exact file names/paths.

### Stroop (Many Labs 3) — OSF
Download and place under `data/raw/ml3/`:
- `https://osf.io/n8xa7/download` → `data/raw/ml3/StroopCleanSet.csv`
- `https://osf.io/bxw8j/download` → `data/raw/ml3/ML3AllSitesandmTurk.csv`

### PSA001 (PSA Social Faces) — OSF
Use the PSA001 project page and download the **Full data** files:
- OSF project: `https://osf.io/f7v3n/`
- Files needed: `psa001_ind.csv` and `psa001_cfd_faces.csv`

If you only have the smaller exploratory subset, see `data/README.md` for how to point
the trait builder at the subset files.

## Output files (generated)
- `reports/` contains RQ tables and figures (generated).
- `models/` contains CmdStan outputs and `.rds` fits (generated).

Note: Many Smiles artifacts are archived under `archive/datasets/many_smiles/` and are not part of the main RQ pipeline.
