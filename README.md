# Tirocinio — Shannon Entropy & Stemness in Meningioma Spatial Transcriptomics

## Overview
This project investigates whether Shannon entropy, computed over spatial
transcriptomics data, can identify high-stemness regions in meningioma
tumors. ~20 spatial samples are processed and analyzed using Seurat and
related R tooling.

For the reasoning behind specific methodological choices (normalization,
entropy formulation, stemness reference, statistical tests), see
`docs/DECISIONS.md` — that file is the source of truth for *why* the
pipeline works the way it does, not just what it does.

## Repository structure

```
Tirocinio/
├── data/
│   └── sampleN/                        # one folder per sample
│       └── raw_data/
│           ├── spatial/                # .tiff, .csv, .png, .json from sample
│           └── GSM..._raw_feature_bc_matrix.h5
├── results/
│   ├── seurat_objects/                 # normalized counts + entropy values (Seurat objects)
│   └── <other-folders-from-each-script-output>
├── docs/
│   └── DECISIONS.md                    # methodology decisions + rationale, updated as the project evolves
├── AGENTS.md                           # rules/context for AI coding agents (Antigravity, Claude Code)
├── environment.yml                     # conda environment
├── shannon_entropy.R                   # function
├── entropy_correlation.R               # function
├── analyze_entropy.R
└── analyze_stemness.R

```

## Pipeline order

1. `analyze_entropy.R` — uses `shannon_entropy.R` and `entropy_correlation.R` to calculate entropy on SpaNorm normalized data and tests the results
2. `analyze_stemness.R` — uses `AddModuleScore()` to map where entropy markers are and if they match high entropy regions

## Environment setup

```bash
conda env create -f environment.yml
conda activate <env-name>
```

## Status

Thesis project, in progress. See `docs/DECISIONS.md` for the current state
of each pipeline stage (draft / implemented / under review / finalized).

```
