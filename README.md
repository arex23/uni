# Shannon Entropy & Stemness in Meningioma Spatial Transcriptomics

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
│   ├── seurat_objects/                 # normalized counts + entropy values + stemness scores (Seurat objects)
│   ├── analyze_entropy/                # spatial entropy plots and QC metrics tables
│   ├── statistical_tests/              # sequencing depth vs entropy correlations (raw & normalized)
│   ├── stemness_analysis/              # stemness score vs entropy correlations, spatial plots, and marker QC
│   └── <other-folders-from-each-script-output>
├── docs/
    ├── PROPOSED_FIXES.md
│   └── DECISIONS.md                    # methodology decisions + rationale, updated as the project evolves
├── R/
│   ├── shannon_entropy.R               # Shannon entropy calculation function
│   └── entropy_correlation.R           # Generalized correlation & publication scatterplot function
├── AGENTS.md                           # rules/context for AI coding agents (Antigravity, Claude Code)
├── environment.yml                     # conda environment
├── analyze_entropy.R                   # Stage 1: preprocessing, SpaNorm normalization, entropy calculation
├── analyze_stemness.R                  # Stage 2: stemness module scoring, correlation, and comparison
├── check_cohort_chemistry.R
└── diagnose_entropy_scaling.R
```

## Pipeline order

1. `analyze_entropy.R` — creates a Seurat object with the sample data; normalizes it with SpaNorm, calculates per-spot entropy using `R/shannon_entropy.R`; uses `R/entropy_correlation.R` to test whether there's still correlation between entropy and sequencing depth.
2. `analyze_stemness.R` — computes stemness module score using `AddModuleScore()`, evaluates entropy-stemness correlation using `R/entropy_correlation.R`, generates spatial comparison plots, and updates the Seurat object.
3. `find_entropy_markers.R` — downstream differential expression and marker discovery.

## Environment setup

```bash
conda env create -f environment.yml
conda activate stemness
```

## Status

Thesis project, in progress. See `docs/DECISIONS.md` for the current state
of each pipeline stage (draft / implemented / under review / finalized).
