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
│   ├── seurat_objects/                 # normalized counts + entropy columns + stemness scores (Seurat objects)
│   ├── analyze_entropy/                # spatial entropy plots and QC metrics tables
│   ├── statistical_tests/              # sequencing depth vs entropy correlations, scaling diagnostics
│   ├── rarefaction/                    # outputs of the off-pipeline rarefaction_entropy.R sweep
│   ├── stemness_analysis/              # stemness score vs entropy correlations, spatial plots, marker QC
│   ├── cohort_qc/                      # Stage 0 chemistry check and cohort gene universe
│   └── <other-folders-from-each-script-output>
├── docs/
│   ├── PROPOSED_FIXES.md
│   ├── STAGE_3_ANALYSIS.md             # superseded: review of the abandoned rarefaction run
│   ├── DECISIONS.md                    # methodology decisions + rationale, updated as the project evolves
│   └── reference_unoptimized/          # Pre-optimization reference copies of pipeline scripts
├── R/
│   ├── shannon_entropy.R               # Plug-in Shannon entropy calculation (column-blocked sparse)
│   ├── gene_universe.R                 # Cohort gene universe loading and filtering helpers
│   ├── spanorm_lowmem.R                # Gene-blocked low-memory SpaNorm logpac adjustment kernel
│   ├── entropy_correlation.R           # Generalized correlation & publication scatterplot function
│   ├── quality_confound.R              # Spot-quality confound control: partial & subset correlations
│   └── cohort.R                        # Frozen 16-sample analysis cohort (D6) and retention threshold
├── AGENTS.md                           # rules/context for AI coding agents (Antigravity, Claude Code)
├── environment.yml                     # conda environment
├── check_cohort_chemistry.R            # Stage 0: Chemistry and feature grid verification across 22 samples
├── build_gene_universe.R               # Stage 0/3: Cohort gene universe generation across all 22 samples
├── analyze_entropy.R                   # Stage 3: QC (nCount >= 500), plug-in entropy, SpaNorm normalization
├── diagnose_entropy_scaling.R          # Stage 2: Diagnostic comparing raw, logpac and linear back-transformed entropy
├── rarefaction_entropy.R               # Off-pipeline: the abandoned rarefied estimator + its depth sweep
├── analyze_stemness.R                  # Stage 3/5: Stemness module scoring, correlation, and spatial comparison
├── run_cohort.sh                       # Cohort driver: one Rscript process per sample per stage
└── check_cohort_retention.R            # D6 retention gate, run after the cohort completes
```

## Pipeline order

1. `check_cohort_chemistry.R` — verifies cohort chemistry, feature list identity, and on-tissue metrics across all 22 samples.
2. `build_gene_universe.R` — generates the deterministic cohort gene universe (`results/cohort_qc/gene_universe.csv`, 18,535 expressed non-deprecated features across 22 samples).
3. `analyze_entropy.R` — runs on the 16-sample cohort only (D6). Filters on-tissue spots ($nCount \ge 500, nFeature \ge 250$) and the gene universe as-is; computes plug-in Shannon entropy on raw counts (`entropy_raw_plugin`); normalizes via SpaNorm; computes plug-in entropy on the normalized `data` layer (`entropy_spanorm_plugin`); evaluates entropy vs depth correlations. Both entropy columns are documented baselines with known biases (D2), not committed metrics.
4. `diagnose_entropy_scaling.R` — diagnostic script comparing raw counts, log-scale SpaNorm and linear back-transformed ($2^x - 1$) SpaNorm entropy scaling against sequencing depth and feature support.
5. `analyze_stemness.R` — runs on the 16-sample cohort only (D6). Computes stemness module score using `AddModuleScore()`, evaluates per-gene marker detection and entropy-stemness correlation using `R/entropy_correlation.R`, runs the spot-quality confound control (`R/quality_confound.R`: score-vs-quality correlations, entropy-stemness partial correlations given `percent.mt` and $\log(nCount)$, and `percent.mt` subset/quartile sensitivity), generates spatial comparison plots, and updates the Seurat object.
6. `check_cohort_retention.R` — applies the frozen D6 retention floor (≥ 60% of on-tissue spots kept through the QC gates) to the per-sample QC tables and writes `results/cohort_qc/retention_gate.csv`. Run after the cohort completes, before any cross-sample inference.
7. `find_entropy_markers.R` — downstream differential expression and marker discovery. *(Not yet written.)*

Off the pipeline: `rarefaction_entropy.R` holds the rarefied entropy estimator and its downsampling-depth sweep. Rarefaction was implemented, evaluated and set aside (D2); the script is standalone, nothing in the pipeline sources it, and it writes to `results/rarefaction/`. Run it as `Rscript rarefaction_entropy.R <sample> [D ...]` if the sweep needs reproducing.

### Running the cohort

`analyze_entropy.R` and `analyze_stemness.R` each process one sample per invocation. Use the driver rather than an in-process loop — SpaNorm peaks near 12.5 GB and only process exit reliably returns it:

```bash
./run_cohort.sh                 # entropy then stemness, over all 16 cohort samples
./run_cohort.sh entropy         # one stage only
./run_cohort.sh entropy sample4 sample21
Rscript check_cohort_retention.R
```

Per-sample logs land in `results/cohort_qc/logs/`. A failing sample is skipped rather than aborting the run, and the exit status is non-zero if anything failed.

## Environment setup

```bash
conda env create -f environment.yml
conda activate stemness
```

## Status

Thesis project, in progress. See `docs/DECISIONS.md` for the current state
of each pipeline stage (draft / implemented / under review / finalized).
