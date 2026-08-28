# Decisions log

This file is the methodology record for the project: what was decided,
why, what alternatives were considered, and what happened when the
implementation was reviewed. Treat it as a running draft of the thesis
Methods section — write entries in full sentences, not just keywords.

Update `Status` as things progress: `draft` → `implemented` →
`under review` → `finalized`.

---

## D1 — Normalization & preprocessing
**Script:** `analyze_entropy.R`
**Decision:** Raw Visium spatial transcriptomics data are loaded strictly from `data/<sample_name>/raw_data/` requiring exactly one `.h5` file. Mitochondrial (`percent.mt` from `^MT-`) and ribosomal (`percent.ribo` from `^RP[SL]`) transcript percentages are calculated for each spot. Spots are filtered using minimum thresholds of $nCount \ge 500$ and $nFeature \ge 250$, and aligned to imaged tissue coordinates via `align_spots_to_coords()`. Gene expression filtering is computed *after* spot subsetting, retaining genes detected in at least $\max(20, \lceil 0.02 \times N_{\text{valid\_spots}} \rceil)$ spots. Normalization is performed directly via `SpaNorm::SpaNorm()` using negative binomial modeling (`gene.model = "nb"`), `adj.method = "logpac"`, `df.tps = 6`, `lambda.a = 1e-04`, a fixed random seed (`set.seed(23)`), and default sampling fraction (`sample.p = 0.25`). An S4 method for `SpaNorm::SpaNorm` is registered on `Seurat` objects using public Bioconductor APIs (`SpatialExperiment` and `SummarizedExperiment`), completely eliminating private `:::` function calls. Spot and gene retention counts, along with average mitochondrial and ribosomal percentages, are logged per sample to `results/analyze_entropy/<sample_name>_qc_metrics.csv`.
**Rationale:** Computing gene filtering thresholds after spot subsetting prevents dropped low-quality or off-tissue spots from falsely inflating gene retention. SpaNorm removes technical spatial noise and sequencing depth confounders while preserving true spatial tissue architecture. A sampling fraction of `sample.p = 0.25` matches the package default, reducing variance in the fitted spatial spline surface. Interfacing through the standard `SpaNorm::SpaNorm()` generic maintains full compliance with package interfaces and avoids depending on unexported internal routines.
**Alternatives considered / rejected:** QC thresholds were raised from 100/50 to 500/250 for computational and memory reasons; however, this aggressive cutoff carries the biological consequence of preferentially removing low-cellularity spots (e.g., necrotic centers, sample edges, and sparse stroma). SCTransform and standard log-normalization were rejected because they do not explicitly model spatial coordinates or adjust for spatial technical autocorrelation.
**Status:** implemented
**Review notes:** Replaced internal `:::` fitting architecture with direct `SpaNorm::SpaNorm()` method dispatch on `Seurat` objects. Pinned `bioconductor-spanorm=1.4.0=r45hdfd78af_0` and `bioconductor-spatialexperiment=1.20.0=r45hdfd78af_0` in `environment.yml`. Log spot/gene pre- and post-filtering counts and QC metrics per sample.

---

## D2 — Shannon entropy formulation
**Script:** `R/shannon_entropy.R` (called from `analyze_entropy.R`)
**Decision:** Shannon entropy is computed per spot as $H = -\sum_{i} p_i \log_2(p_i)$, where $p_i = y_i / \sum_j y_j$ represents the relative expression proportion of gene $i$ within that spot. Mitochondrial (`^MT-`) and ribosomal (`^RP[SL]`) gene families are strictly excluded from the expression matrix prior to entropy calculation. In `analyze_entropy.R`, entropy is computed twice: once on raw counts (`layer = "counts"`, stored as `shannon_entropy_raw`) before normalization, and once on normalized expression (`layer = "data"`, stored as `shannon_entropy`) after SpaNorm. Calculations are executed via sparse matrix operations (`dgCMatrix`) to optimize memory efficiency.
**Rationale:** In Visium spatial transcriptomics, mitochondrial and ribosomal transcripts can constitute 10–30% of total spot counts and predominantly reflect technical tissue quality, hypoxia, necrosis, or translational machinery rather than true phenotypic plasticity. Excluding MT and RP genes ensures that Shannon entropy reflects genuine biological transcriptomic diversity. Computing entropy on both raw and normalized counts allows an integrated within-pipeline QC assessment of depth decoupling.
**Alternatives considered / rejected:** Including MT and RP transcripts in entropy calculation was rejected because spot-level entropy would be heavily dominated by non-specific metabolic and translational variations. Top highly variable genes (HVGs) only was rejected because transcriptional plasticity spans broader gene sets.
**Status:** implemented
**Review notes:** Added `exclude_pattern = "^(MT-|RP[SL])"` and `col.name` parameter to `calculate_shannon_entropy()` to support dual calculation of raw and normalized entropy layers without MT/RP contamination.

---

## D3 — Entropy–depth & covariate QC and stemness correlation tests
**Script:** `R/entropy_correlation.R` (called from `analyze_entropy.R` and `analyze_stemness.R`)
**Decision:** Both parametric Pearson correlation ($r$) and non-parametric Spearman rank correlation ($\rho$) are computed alongside two-sided $p$-values:
1. In `analyze_entropy.R`: evaluates correlations for both `shannon_entropy_raw` and `shannon_entropy` against sequencing depth (`nCount`, `nFeature`) and quality/metabolic covariates (`percent.mt`, `percent.ribo`), outputting results to `results/statistical_tests/<sample_name>_entropy_correlations.csv` and an 8-panel scatterplot figure.
2. In `analyze_stemness.R`: evaluates biological correlation between normalized `shannon_entropy` and `Stemness_Score1`, outputting results to `results/stemness_analysis/<sample_name>_stemness_correlations.csv` and `_plot.png`.
All outputs contain a leading `Sample` column for automated cohort-level multi-sample concatenation (`rbind`).
**Rationale:** Pearson correlation assesses linear relationships, whereas Spearman correlation captures monotonic non-linear dependencies without normality assumptions. Testing entropy against `percent.mt` and `percent.ribo` confirms whether spatial entropy patterns are confounded by tissue degradation or necrosis. Comparing raw vs normalized entropy against sequencing depth quantitatively proves that SpaNorm successfully decoupled entropy from library size.
**Alternatives considered / rejected:** Unstructured console text logging via `sink()` was rejected because it is not machine-readable for cross-sample cohort comparison.
**Status:** implemented
**Review notes:** Refactored `calculate_entropy_correlations()` in `R/entropy_correlation.R` to accept arbitrary named target columns (`targets`), multiple entropy columns (`entropy_cols`), and auto-detect `percent.mt` and `percent.ribo`.

---

## D4 — Stemness gene reference / scoring
**Script:** `analyze_stemness.R`
**Decision:** Spot-level stemness is quantified using Seurat's `AddModuleScore()` with a curated panel of meningioma and neural stemness markers (`PROM1`, `SOX2`, `POU5F1`, `NANOG`, `NES`, `CD44`, `MYC`). Prior to scoring, per-gene spot detection rates and mean expression counts are calculated and saved to `results/stemness_analysis/<sample_name>_stemness_marker_detection.csv`. Variable marker panels are accepted across samples due to Visium dropouts, and the exact subset of detected vs missing markers is logged in `results/stemness_analysis/<sample_name>_stemness_qc.csv`. Random seed 23 (`set.seed(23)`) is initialized before scoring to ensure deterministic control feature binning and sampling. Updated Seurat objects containing `Stemness_Score1` are saved back to `results/seurat_objects/<sample_name>_spatial_obj.rds` for downstream marker discovery (`find_entropy_markers.R`).
**Rationale:** Evaluating expression of established stem cell transcription factors and surface markers provides a validated baseline reference for stemness in meningioma. Pre-evaluating marker detection ensures that researchers can verify whether module scores reflect multi-gene programs or are predominantly driven by single robust markers. Persisting the resulting scores in the Seurat object prevents redundant re-computation in downstream pipeline steps.
**Alternatives considered / rejected:** Pluripotency transcription factors including Oct-4 (`POU5F1`) and `SOX2` are known to be present in non-stem tumor cells and normal meningeal cells, which limits their specificity for cleanly distinguishing stem from non-stem populations. Furthermore, `POU5F1`, `NANOG`, and `PROM1` are low-abundance transcripts frequently detected in <2% of Visium spots; when dropped or near zero, `AddModuleScore` is effectively driven by robustly detected markers like `CD44`, `NES`, and `MYC`. Candidate alternative or supplementary markers for future panels include Frizzled 9 (`FZD9`), `GFAP`, Vimentin (`VIM`), SSEA-4, `CD44`, `CD73` (`NT5E`), and `CD105` (`ENG`), of which `VIM` and `CD44` exhibit especially high detection sensitivity in Visium.
**Status:** implemented
**Review notes:** Added per-gene detection evaluation (`_stemness_marker_detection.csv`), per-sample marker presence tracking (`_stemness_qc.csv`), deterministic control gene selection (`set.seed(23)`), and unified correlation reporting via `R/entropy_correlation.R`.

---

## Template for new entries

```
## D<N> — <short title>
**Script:** <file>
**Decision:**
**Rationale:**
**Alternatives considered / rejected:**
**Status:** draft
**Review notes:**
```
