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
`check_cohort_chemistry.R` (Stage 0) → `build_gene_universe.R` (calls `R/gene_universe.R`) → `analyze_entropy.R` (calls `R/load_sample.R`, `R/spanorm_lowmem.R`, `R/gene_universe.R`, `R/shannon_entropy.R`, `R/entropy_correlation.R`, `R/cohort.R`) → `analyze_stemness.R` (calls `R/entropy_correlation.R`, `R/quality_confound.R`, `R/cohort.R`) → `check_cohort_retention.R` (D6 retention gate) → `find_entropy_markers.R` (not yet written)

`R/load_sample.R` holds D1's load-and-QC (`load_qc_sample()`): in-tissue subset, array lattice coordinates, `percent.mt`/`percent.ribo`, the `nCount ≥ 500` / `nFeature ≥ 250` floors, coordinate validation, gene-universe filter. It is the single implementation of those steps — `analyze_entropy.R` calls it rather than carrying its own copy, so a diagnostic can reach the identical spot set without also fitting SpaNorm. Anything needing a QC'd object should call it instead of re-deriving the QC.

**Off-pipeline diagnostics** (they produce no Seurat object and nothing downstream reads their output):

- `diagnose_entropy_scaling.R` — Stage A3 estimator bake-off. Loads through `load_qc_sample()`, computes the four D2 estimators, and evaluates the A0 criteria checkable without downstream input.
- `diagnose_subsampling_stability.R` — Stage A4 depth-invariance ladder. Common-proportion binomial thinning, all four estimators at each level. This is a *measuring instrument*, not an estimator: it perturbs depth to test invariance and is not rarefaction (D2). Also emits `<sample>_subsampling_reliability.csv`, the depth-free split-half reliability that replaced the retired rank half of A0 criterion 3 (D2).
- `diagnose_correction_routes.R` — Stage A5a. Compares the two routes to a depth-corrected entropy: the corrected estimator itself against the plug-in residualised on `log(nCount)`. Also measures the residual depth structure a *linear* adjustment leaves behind, via a spline of `log(nCount)`. Its decision rule was frozen in D2 before the script was first run.

All three take `--all` to sweep the frozen 16-sample cohort (D6) instead of one sample, writing `cohort_*.csv` alongside the per-sample files. The driver is `R/diagnostic_sweep.R`: it resolves the sample set through `cohort_samples()`, reseeds the RNG per sample so results do not depend on loop position, and records a failing sample and continues rather than losing a 45-minute unattended run. None of them runs SpaNorm — that is what makes a cohort sweep affordable.

`R/spatial_neighbors.R` holds the Visium hex-lattice adjacency and Moran's I with permutation inference (D7), built for Stage B2 and reused by the spatial inference. No `spdep`/`sf` dependency — see D7.

`rarefaction_entropy.R` is deliberately outside this chain. Rarefaction was tried and set aside (D2); the script keeps the estimator and its depth sweep runnable on their own. The dependency runs one way only — it may `source()` files in `R/`, but nothing in `R/` and neither `analyze_entropy.R` nor `analyze_stemness.R` may reference it. Keep it that way: if a new entropy approach needs something from it, copy the piece rather than wiring the pipeline back to it.

## Tests
`tests/` holds standalone check scripts, run from the repository root with
`Rscript tests/<file>.R`. They use base R plus the packages the pipeline
already requires — no test framework, and no test-only dependency in
`environment.yml`. A script exits non-zero on failure.

`tests/test_entropy_estimators.R` validates the D2 diversity estimators in
`R/shannon_entropy.R`. Its expected values were verified once against the CRAN
`entropy` package and are frozen as literals so the dependency could be
dropped; if you change the estimator arithmetic, re-verify against a throwaway
install of `entropy` and regenerate the fixtures rather than editing the
numbers to match the new output.

Cohort runs go through `./run_cohort.sh`, one `Rscript` process per sample per stage — SpaNorm peaks near 12.5 GB on a 15 GB machine and an in-process loop gets OOM-killed partway through.

## Cohort scope
`R/cohort.R` holds the frozen 16-sample analysis cohort (D6). Any script producing a *biological* result takes its sample set from `cohort_samples()` and resolves its CLI argument through `resolve_target_sample()`; do not call `list.dirs("data")` directly in such a script. Only `check_cohort_chemistry.R` (which produces the exclusion criterion) and `build_gene_universe.R` (reference data, D5) run over all 22 via `all_samples()`. The exclusion list and the retention floor are frozen upfront on technical criteria and must never be revised after seeing entropy or stemness results.

## Data handling
- `data/sampleN/raw_data/` — raw per-sample Visium output. Read-only,
  never modify in place.
- All generated output goes to `results/`, mirroring the script that
  produced it:
  - `results/cohort_qc/` ← `check_cohort_chemistry.R`, `build_gene_universe.R`
  - `results/seurat_objects/` ← normalized counts + entropy columns + stemness scores (Seurat objects from `analyze_entropy.R` / `analyze_stemness.R`)
  - `results/analyze_entropy/` ← `analyze_entropy.R` (QC metrics, spatial entropy plots)
  - `results/statistical_tests/` ← `entropy_correlation.R`, `diagnose_entropy_scaling.R`
  - `results/rarefaction/` ← `rarefaction_entropy.R` (off-pipeline, D2)
  - `results/stemness_analysis/` ← `analyze_stemness.R` (marker QC, stemness correlations, quality confound tables, spatial comparisons)
  - `results/entropy_deg_plots/` and `results/entropy_deg/` ← `find_entropy_markers.R`
