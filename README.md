# Shannon Entropy & Stemness in Meningioma Spatial Transcriptomics

## Overview
This project investigates whether Shannon entropy, computed over spatial
transcriptomics data, can identify high-stemness regions in meningioma
tumors. ~20 spatial samples are processed and analyzed using Seurat and
related R tooling.

For the reasoning behind specific methodological choices (normalization,
entropy formulation, gene reference sets, statistical tests), see
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
│   ├── seurat_objects/                 # normalized counts + rarefied entropy + stemness scores (Seurat objects)
│   ├── analyze_entropy/                # spatial rarefied entropy plots and QC metrics tables
│   ├── statistical_tests/              # sequencing depth vs entropy correlations, depth sweeps, diagnostics
│   ├── stemness_analysis/              # stemness score vs entropy correlations, spatial plots, marker QC
│   ├── cohort_qc/                      # Stage 0 chemistry check and cohort gene universe
│   └── <other-folders-from-each-script-output>
├── docs/
│   ├── PROPOSED_FIXES.md
│   ├── DECISIONS.md                    # methodology decisions + rationale, updated as the project evolves
│   └── reference_unoptimized/          # Pre-optimization reference copies of pipeline scripts
├── R/
│   ├── shannon_entropy.R               # Rarefied & plug-in Shannon entropy calculation (column-blocked sparse)
│   ├── gene_universe.R                 # Cohort gene universe loading and filtering helpers
│   ├── spanorm_lowmem.R                # Gene-blocked low-memory SpaNorm logpac adjustment kernel
│   └── entropy_correlation.R           # Generalized correlation & publication scatterplot function
├── AGENTS.md                           # rules/context for AI coding agents (Antigravity, Claude Code)
├── environment.yml                     # conda environment
├── check_cohort_chemistry.R            # Stage 0: Chemistry and feature grid verification across 22 samples
├── build_gene_universe.R               # Stage 0/3: Cohort gene universe generation across all 22 samples
├── sweep_rarefaction_depth.R           # Stage 3: Downsampling depth sweep diagnostic (D in {1000..8000})
├── analyze_entropy.R                   # Stage 3: QC (nCount >= 3000), rarefied entropy, SpaNorm normalization
├── diagnose_entropy_scaling.R          # Stage 2/3: Diagnostic comparing raw, logpac, linear, and rarefied entropy
└── analyze_stemness.R                  # Stage 3/5: Stemness module scoring, correlation, and spatial comparison
```

## Pipeline order

1. `check_cohort_chemistry.R` — verifies cohort chemistry, feature list identity, and on-tissue metrics across all 22 samples.
2. `build_gene_universe.R` — generates the deterministic cohort gene universe (`results/cohort_qc/gene_universe.csv`, 18,535 expressed non-deprecated features across 22 samples).
3. `sweep_rarefaction_depth.R` — evaluates candidate downsampling depths $D \in \{1000, 2000, 3000, 5000, 8000\}$ for spot retention and depth decoupling. Reads raw counts from the sample `.h5`, not from the Seurat object, so it runs before `analyze_entropy.R` and its retention column is not circular.
4. `analyze_entropy.R` — filters on-tissue spots ($nCount \ge 3000, nFeature \ge 250$) and gene universe, then re-gates spots on the post-exclusion entropy-matrix depth ($\ge D$); computes primary rarefied Shannon entropy ($D = 3000, n_{\text{draws}} = 5$) and baseline raw plug-in entropy; normalizes via SpaNorm; evaluates entropy vs depth correlations.
5. `diagnose_entropy_scaling.R` — diagnostic script comparing raw counts, log-scale SpaNorm, linear back-transformed ($2^x - 1$) SpaNorm, and rarefied entropy scaling against sequencing depth and feature support.
6. `analyze_stemness.R` — computes stemness module score using `AddModuleScore()`, evaluates per-gene marker detection and entropy-stemness correlation using `R/entropy_correlation.R`, generates spatial comparison plots, and updates the Seurat object.
7. `find_entropy_markers.R` — downstream differential expression and marker discovery.

## Environment setup

```bash
conda env create -f environment.yml
conda activate stemness
```

## Status

Thesis project, in progress. See `docs/DECISIONS.md` for the current state
of each pipeline stage (draft / implemented / under review / finalized).
