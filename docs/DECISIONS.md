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
**Decision:** Raw Visium spatial transcriptomics data are loaded from `data/<sample_name>/`. Spots are explicitly subsetted to in-tissue spots (`in_tissue == 1`) using `spatial/tissue_positions.csv` *before* applying depth quality control. Spots are then filtered using minimum thresholds of $nCount \ge 500$ and $nFeature \ge 250$, followed by coordinate alignment verification via `validate_imaged_coordinates()`. Gene expression filtering is computed *after* spot subsetting, retaining genes detected in at least $\max(20, \lceil 0.02 \times N_{\text{valid\_spots}} \rceil)$ spots. Normalization is performed directly via `SpaNorm::SpaNorm()` using negative binomial modeling (`gene.model = "nb"`), `adj.method = "logpac"`, `df.tps = 6`, `lambda.a = 1e-04`, a fixed random seed (`set.seed(23)`), and default sampling fraction (`sample.p = 0.25`). An S4 method for `SpaNorm::SpaNorm` is registered on `Seurat` objects using public Bioconductor APIs (`SpatialExperiment` and `SummarizedExperiment`), eliminating private `:::` function calls. Honest QC metrics with explicit denominators (`Raw_Spots_Grid`, `Spots_On_Tissue`, `Spots_Post_Depth_QC`, `Pct_OnTissue_Retained`, `Raw_Genes`, `Genes_Post_Filter`, `Pct_Genes_Retained`, `Mean_Percent_MT`, `Mean_Percent_Ribo`) are logged per sample to `results/analyze_entropy/<sample_name>_qc_metrics.csv`.
**Rationale:** `Load10X_Spatial(filter.matrix = TRUE)` loads all 14,336 grid barcodes when reading raw feature matrices. Explicitly subsetting on `in_tissue == 1` before depth QC ensures that QC denominators reflect true tissue area rather than the unutilized slide background (which contains ~2,600 empty spots with ambient RNA). Computing gene filtering thresholds after spot subsetting prevents dropped low-quality spots from falsely inflating gene retention. SpaNorm removes technical spatial noise and sequencing depth confounders while preserving spatial tissue architecture.
**Alternatives considered / rejected:** SCTransform and standard log-normalization were rejected because they do not explicitly model spatial coordinates or adjust for spatial technical autocorrelation. QC thresholds of $nCount \ge 500$ & $nFeature \ge 250$ retain 99.97% of on-tissue spots in sample1 (11,709 of 11,713 spots, dropping only 4 marginal spots) and 91.2%–99.9% across the 22-sample cohort, proving that this threshold removes low-quality artifacts without discarding anatomical zones or necrotic centers.
**Status:** implemented
**Review notes:** Replaced internal `:::` fitting architecture with direct `SpaNorm::SpaNorm()` method dispatch on `Seurat` objects. Pinned `bioconductor-spanorm=1.4.0=r45hdfd78af_0` and `bioconductor-spatialexperiment=1.20.0=r45hdfd78af_0` in `environment.yml`. Logged explicit on-tissue QC metrics (`Raw_Spots_Grid`, `Spots_On_Tissue`, `Spots_Post_Depth_QC`, `Pct_OnTissue_Retained`) per sample.

---

## D2 — Shannon entropy formulation
**Script:** `R/shannon_entropy.R` (called from `analyze_entropy.R`)
**Decision:** Shannon entropy is computed per spot as $H = -\sum_{i} p_i \log_2(p_i)$, where $p_i = y_i / \sum_j y_j$ represents the relative expression proportion of gene $i$ within that spot. Mitochondrial (`^MT-`) and ribosomal (`^RP[SL]`) gene families are filtered from the expression matrix prior to entropy calculation via `exclude_pattern = "^(MT-|RP[SL])"`. In `analyze_entropy.R`, entropy is computed twice: once on raw counts (`layer = "counts"`, stored as `shannon_entropy_raw`) before normalization, and once on normalized expression (`layer = "data"`, stored as `shannon_entropy`) after SpaNorm. Calculations are executed via sparse matrix operations (`dgCMatrix`) to optimize memory efficiency.
**Rationale:** In the Visium CytAssist Human Transcriptome probe set (Visium V5 Slide chemistry), cytoplasmic ribosomal protein probes are deliberately excluded from probe panel design to prevent translational machinery from overwhelming sequencing bandwidth. Stage 0 cohort analysis confirmed that all 103 `^RP[SL]` rows have exactly 0 counts across all 22 samples (`n_RPSL_nonzero = 0`). Non-zero mitochondrial transcript content was measured at 2.46%–10.31% of total counts (mean spot MT ~3.4%–16.8% across samples). Excluding MT transcripts removes metabolic and hypoxia-driven artifacts, while `exclude_pattern = "^(MT-|RP[SL])"` acts in practice as an MT-only filter given zero ribosomal probe representation.
**Alternatives considered / rejected:** Including MT transcripts in entropy calculation was rejected because spot-level entropy would be distorted by variable mitochondrial transcription. Top highly variable genes (HVGs) only was rejected because transcriptional plasticity spans broader gene sets.
**Status:** implemented
**Review notes:** Confirmed ribosomal probe absence across cohort chemistry (`n_RPSL_nonzero = 0`). Documented empirical mitochondrial fractions (2.5%–10.5%). Shannon entropy filtering excludes MT features to capture phenotypic biological diversity.

---

## D3 — Entropy–depth & covariate QC and stemness correlation tests
**Script:** `R/entropy_correlation.R` (called from `analyze_entropy.R` and `analyze_stemness.R`)
**Decision:** Both parametric Pearson correlation ($r$) and non-parametric Spearman rank correlation ($\rho$) are computed alongside two-sided $p$-values, log10 $p$-values, 95% confidence intervals, and spot sample size ($N$):
1. In `analyze_entropy.R`: evaluates correlations for both `shannon_entropy_raw` and `shannon_entropy` against sequencing depth (`nCount`, `nFeature`) and quality/metabolic covariates (`percent.mt`, `percent.ribo`), outputting results to `results/statistical_tests/<sample_name>_entropy_correlations.csv` and an 8-panel scatterplot figure.
2. In `analyze_stemness.R`: evaluates biological correlation between normalized `shannon_entropy` and `Stemness_Score1`, outputting results to `results/stemness_analysis/<sample_name>_stemness_correlations.csv` and `_plot.png`.
All summary CSV outputs store exact unfloored double-precision floats for `Pearson_p` and `Spearman_p`, include `Pearson_log10_p` computed via Student's $t$ cumulative distribution in log space ($\log_{10}(p) = \log_{10}(2) + \frac{\text{pt}(-|t|, \text{df}, \text{log.p = TRUE})}{\ln(10)}$) to eliminate saturation at `0.00e+00`, report Pearson 95% CIs and Fisher $z$-transformed Spearman 95% CIs ($z = \text{atanh}(\rho), \text{SE} = 1/\sqrt{N-3}$), and include spot count $N$ alongside the leading `Sample` column.
**Rationale:** Pearson correlation assesses linear relationships, whereas Spearman correlation captures monotonic non-linear dependencies without normality assumptions. Testing entropy against `percent.mt` and `percent.ribo` confirms whether spatial entropy patterns are confounded by tissue degradation. Unfloored floating-point $p$-values and log10 $p$-values enable downstream meta-analysis and FDR corrections across the cohort without truncation artifacts.
**Alternatives considered / rejected:** Storing formatted character strings in CSVs was rejected because it saturated small $p$-values at `"0.00e+00"`, preventing quantitative ranking and meta-analytic pooling.
**Status:** implemented
**Review notes:** Refactored `calculate_entropy_correlations()` in `R/entropy_correlation.R` to output numeric double columns, log10 $p$-values, Pearson and Spearman Fisher-z 95% CIs, and sample size $N$. Formatted strings are restricted to plot display annotations.

---

## D4 — Stemness gene reference / scoring
**Script:** `analyze_stemness.R`
**Decision:** Spot-level stemness is quantified using Seurat's `AddModuleScore()` with a curated panel of meningioma and neural stemness markers (`PROM1`, `SOX2`, `POU5F1`, `NANOG`, `NES`, `CD44`, `MYC`). Prior to scoring, per-gene spot detection rates, mean expression counts, and detection levels are calculated and saved to `results/stemness_analysis/<sample_name>_stemness_marker_detection.csv`. Variable marker panels are accepted across samples due to Visium dropouts, and the exact subset of detected vs missing markers is logged in `results/stemness_analysis/<sample_name>_stemness_qc.csv`. Random seed 23 (`set.seed(23)`) is initialized before scoring to ensure deterministic control feature binning and sampling. Updated Seurat objects containing `Stemness_Score1` are saved back to `results/seurat_objects/<sample_name>_spatial_obj.rds` for downstream marker discovery (`find_entropy_markers.R`).
**Rationale:** Evaluating expression of established stem cell transcription factors and surface markers provides a validated baseline reference for stemness in meningioma. Pre-evaluating marker detection ensures that researchers can verify whether module scores reflect multi-gene programs or are predominantly driven by single robust markers. Persisting the resulting scores in the Seurat object prevents redundant re-computation in downstream pipeline steps.
**Alternatives considered / rejected:** Pluripotency transcription factors including Oct-4 (`POU5F1`) and `SOX2` are known to be present in non-stem tumor cells and normal meningeal cells, which limits their specificity for cleanly distinguishing stem from non-stem populations. Measured detection rates in meningioma Visium data reveal heterogeneous expression: robust markers include `CD44` (~79%), `MYC` (~74%), and `NES` (~60%); moderate markers include `SOX2` (~18.6%), `POU5F1` (~14.2%), and `PROM1` (~11.7%); and low-abundance markers include `NANOG` (~2.1%). Candidate alternative or supplementary markers for future panels include Frizzled 9 (`FZD9`), `GFAP`, Vimentin (`VIM`), SSEA-4, `CD44`, `CD73` (`NT5E`), and `CD105` (`ENG`).
**Status:** implemented
**Review notes:** Added per-gene detection evaluation (`_stemness_marker_detection.csv`), per-sample marker presence tracking (`_stemness_qc.csv`), deterministic control gene selection (`set.seed(23)`), and unified correlation reporting via `R/entropy_correlation.R`.
