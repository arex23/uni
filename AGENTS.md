# AGENTS.md — meningioma spatial entropy project

## Context
Spatial transcriptomics project (~20 meningioma samples). Goal: test
whether Shannon entropy identifies high-stemness regions. See
`docs/DECISIONS.md` for the full methodology and rationale — this file
only holds operational rules, not the science.

## Before implementing anything
Read `docs/DECISIONS.md` for the decision matching the pipeline stage
you're touching. Do not invent methodological choices (normalization,
entropy formulation, gene reference sets, statistical tests) — if a
decision isn't documented yet, ask rather than assume, and propose an
entry for `docs/DECISIONS.md` rather than silently picking an approach.

## Before reviewing or modifying existing code
Check the implementation against the matching entry in `docs/DECISIONS.md`.
Flag mismatches between the stated method and the actual code instead of
silently rewriting toward whatever looks more "standard." Note findings
under that entry's "Review notes," don't just fix and move on silently.

## Environment
Use the conda environment defined in `environment.yml`. Don't add,
remove, or upgrade packages without logging it in `docs/DECISIONS.md`.

## Pipeline order
`check_cohort_chemistry.R` (Stage 0) → `build_gene_universe.R` (calls `R/gene_universe.R`) → `sweep_rarefaction_depth.R` (Stage 3 depth sweep) → `analyze_entropy.R` (calls `R/spanorm_lowmem.R`, `R/gene_universe.R`, `R/shannon_entropy.R`, `R/entropy_correlation.R`, `R/cohort.R`) → `diagnose_entropy_scaling.R` (Stage 2/3 scaling diagnostic) → `analyze_stemness.R` (calls `R/entropy_correlation.R`, `R/quality_confound.R`, `R/cohort.R`) → `check_cohort_retention.R` (D6 retention gate) → `find_entropy_markers.R` (not yet written)

Cohort runs go through `./run_cohort.sh`, one `Rscript` process per sample per stage — SpaNorm peaks near 12.5 GB on a 15 GB machine and an in-process loop gets OOM-killed partway through.

## Cohort scope
`R/cohort.R` holds the frozen 16-sample analysis cohort (D6). Any script producing a *biological* result takes its sample set from `cohort_samples()` and resolves its CLI argument through `resolve_target_sample()`; do not call `list.dirs("data")` directly in such a script. Only `check_cohort_chemistry.R` (which produces the exclusion criterion) and `build_gene_universe.R` (reference data, D5) run over all 22 via `all_samples()`. The exclusion list and the retention floor are frozen upfront on technical criteria and must never be revised after seeing entropy or stemness results.

## Data handling
- `data/sampleN/raw_data/` — raw per-sample Visium output. Read-only,
  never modify in place.
- All generated output goes to `results/`, mirroring the script that
  produced it:
  - `results/cohort_qc/` ← `check_cohort_chemistry.R`, `build_gene_universe.R`
  - `results/seurat_objects/` ← normalized counts + rarefied entropy values + stemness scores (Seurat objects from `analyze_entropy.R` / `analyze_stemness.R`)
  - `results/analyze_entropy/` ← `analyze_entropy.R` (QC metrics, spatial rarefied entropy plots)
  - `results/statistical_tests/` ← `entropy_correlation.R`, `sweep_rarefaction_depth.R`, `diagnose_entropy_scaling.R`
  - `results/stemness_analysis/` ← `analyze_stemness.R` (marker QC, stemness correlations, quality confound tables, spatial comparisons)
  - `results/entropy_deg_plots/` and `results/entropy_deg/` ← `find_entropy_markers.R`
