set.seed(23)

library(Seurat)
library(SpatialExperiment)
library(SpaNorm)
library(ggplot2)
library(patchwork)

source("R/shannon_entropy.R")
source("R/entropy_correlation.R")
source("R/spanorm_lowmem.R")
source("R/gene_universe.R")
source("R/cohort.R")
source("R/load_sample.R")

# Replace SpaNorm's logpac adjustment with the gene-blocked kernel. The public
# SpaNorm::SpaNorm() call below is unchanged and the output is bit-identical;
# only the peak memory of the adjustment step changes. Unpatched, that step
# holds ~10 dense 15,556 x 11,709 doubles at once (~13.6 GB) and is what pushes
# a 15.1 GB machine into the OOM killer. See R/spanorm_lowmem.R.
enable_spanorm_lowmem()

# Define S4 method for SpaNorm on Seurat objects using public APIs
setMethod("SpaNorm", signature(spe = "Seurat"), function(spe,
                                                         sample.p = 0.25,
                                                         gene.model = "nb",
                                                         adj.method = "logpac",
                                                         scale.factor = 1,
                                                         df.tps = 6,
                                                         lambda.a = 1e-04,
                                                         batch = NULL,
                                                         tol = 1e-04,
                                                         step.factor = 0.5,
                                                         maxit.nb = 50,
                                                         maxit.psi = 25,
                                                         maxn.psi = 500,
                                                         overwrite = FALSE,
                                                         verbose = TRUE,
                                                         ...) {
  counts_mat <- Seurat::GetAssayData(spe, assay = "Spatial", layer = "counts")
  coords_df <- Seurat::GetTissueCoordinates(spe)
  coords_mat <- as.matrix(coords_df[colnames(counts_mat), 1:2, drop = FALSE])

  spe_exp <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts_mat),
    spatialCoords = coords_mat
  )

  spe_norm <- SpaNorm::SpaNorm(
    spe = spe_exp,
    sample.p = sample.p,
    gene.model = gene.model,
    adj.method = adj.method,
    scale.factor = scale.factor,
    df.tps = df.tps,
    lambda.a = lambda.a,
    batch = batch,
    tol = tol,
    step.factor = step.factor,
    maxit.nb = maxit.nb,
    maxit.psi = maxit.psi,
    maxn.psi = maxn.psi,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )

  norm_mat <- SummarizedExperiment::assay(spe_norm, "logcounts")

  # Release every other reference to the dense normalised matrix before building
  # the sparse copy, so the two representations are never alive simultaneously
  # alongside the SpatialExperiment that also holds the dense one.
  rm(spe_norm, spe_exp, counts_mat, coords_mat, coords_df)
  gc(verbose = FALSE)

  norm_mat <- as(norm_mat, "CsparseMatrix")
  gc(verbose = FALSE)

  spe[["Spatial"]]$data <- norm_mat
  rm(norm_mat)
  gc(verbose = FALSE)
  return(spe)
})

samples <- cohort_samples()

analyze_sample <- function(sample_name) {
  cat("==========================================\n")
  cat("Processing sample:", sample_name, "\n")

  # D1 load + spot/feature QC (steps 1-3), in R/load_sample.R so that the
  # Stage A diagnostics can reach the identical spot set without running SpaNorm.
  entropy_exclude_pattern <- "^(MT-|RP[SL])"
  loaded <- load_qc_sample(sample_name, min_counts = 500, min_features = 250)
  spatial_obj <- loaded$obj

  raw_spots_grid <- loaded$qc$raw_spots_grid
  n_raw_genes <- loaded$qc$n_raw_genes
  n_spots_ontissue <- loaded$qc$n_spots_ontissue
  n_spots_post_depth_qc <- loaded$qc$n_spots_post_depth_qc
  n_spots_post_coord <- loaded$qc$n_spots_post_coord
  n_spots_final <- loaded$qc$n_spots_final
  n_genes_in_universe <- loaded$qc$n_genes_in_universe
  n_genes_detected <- loaded$qc$n_genes_detected
  rm(loaded)
  gc()

  # 4. Entropy on raw counts: the Chao-Shen primary metric and the plug-in
  # baseline it replaces, in ONE blocked pass (D2, Stage A5).
  #
  # `plugin` is numerically identical to the calculate_shannon_entropy() call
  # this replaced -- A2 verified the kernel against CRAN `entropy` to 1.26e-12
  # across 22 fixtures -- so `entropy_raw_plugin` keeps its name, its values and
  # every downstream consumer.
  #
  # Both are carried deliberately. Chao-Shen removes most of the depth artefact
  # and amplifies a percent.mt one (D2); keeping the baseline alongside makes
  # every downstream table a comparison of the two correction routes rather than
  # a single unchecked number, and is what the pre-registered criterion 5(b)
  # gate at Stage C needs to adjudicate them.
  spatial_obj <- calculate_entropy(
    spatial_obj,
    estimator = c("plugin", "chao_shen"),
    col.name = c("entropy_raw_plugin", "entropy_chao_shen"),
    assay = "Spatial",
    layer = "counts",
    exclude_pattern = entropy_exclude_pattern
  )

  # 5. SpaNorm normalization using direct public API (used for downstream DE, scoring, viz)
  spatial_obj <- SpaNorm::SpaNorm(
    spatial_obj,
    sample.p = 0.25,
    gene.model = "nb",
    adj.method = "logpac",
    df.tps = 6,
    lambda.a = 1e-04,
    verbose = TRUE
  )

  # 6. Plug-in Shannon entropy on the SpaNorm logpac `data` layer, as the second
  # baseline. Recorded with its known defect attached: on this layer the metric is
  # near-degenerate against the detected-gene count K (sample1: r = 0.9986 against
  # log2(K), diagnose_entropy_scaling.R), because the log2(qnbinom(.) + 1)
  # compression flattens the relative proportions towards 1/K. It is kept as a
  # documented reference point, not as a candidate metric.
  spatial_obj <- calculate_shannon_entropy(
    spatial_obj,
    assay = "Spatial",
    layer = "data",
    col.name = "entropy_spanorm_plugin",
    exclude_pattern = entropy_exclude_pattern
  )
  gc()

  entropy_dir <- file.path("results", "analyze_entropy")
  stat_dir <- file.path("results", "statistical_tests")
  save_dir <- file.path("results", "seurat_objects")

  if (!dir.exists(entropy_dir)) dir.create(entropy_dir, recursive = TRUE)
  if (!dir.exists(stat_dir)) dir.create(stat_dir, recursive = TRUE)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # Save per-sample QC filtering log table with honest on-tissue accounting
  qc_df <- data.frame(
    Sample = sample_name,
    Raw_Spots_Grid = raw_spots_grid,
    Spots_On_Tissue = n_spots_ontissue,
    Spots_Post_Depth_QC = n_spots_post_depth_qc,
    Spots_Post_Coord_Validation = n_spots_post_coord,
    Spots_Final = n_spots_final,
    Pct_OnTissue_Retained = round((n_spots_final / n_spots_ontissue) * 100, 2),
    Raw_Genes = n_raw_genes,
    Genes_In_Universe = n_genes_in_universe,
    Genes_Detected = n_genes_detected,
    # With the per-sample detection filter dropped, the feature set is the frozen
    # universe for every sample, so a "genes retained" percentage is always 100
    # and its detected-denominator variant exceeds 100 whenever a universe gene
    # goes unobserved here. The informative direction is the other one: how much
    # of the frozen feature space this sample actually observes.
    Pct_Universe_Detected = round((n_genes_detected / n_genes_in_universe) * 100, 2),
    Mean_Percent_MT = round(mean(spatial_obj$percent.mt, na.rm = TRUE), 2),
    Mean_Percent_Ribo = round(mean(spatial_obj$percent.ribo, na.rm = TRUE), 2),
    Mean_Raw_Plugin_Entropy = round(mean(spatial_obj$entropy_raw_plugin, na.rm = TRUE), 4),
    SD_Raw_Plugin_Entropy = round(sd(spatial_obj$entropy_raw_plugin, na.rm = TRUE), 4),
    Median_Raw_Plugin_Entropy = round(median(spatial_obj$entropy_raw_plugin, na.rm = TRUE), 4),
    IQR_Raw_Plugin_Entropy = round(IQR(spatial_obj$entropy_raw_plugin, na.rm = TRUE), 4),
    Min_Raw_Plugin_Entropy = round(min(spatial_obj$entropy_raw_plugin, na.rm = TRUE), 4),
    Max_Raw_Plugin_Entropy = round(max(spatial_obj$entropy_raw_plugin, na.rm = TRUE), 4),
    Range_Raw_Plugin_Entropy = round(diff(range(spatial_obj$entropy_raw_plugin, na.rm = TRUE)), 4),
    Mean_Chao_Shen_Entropy = round(mean(spatial_obj$entropy_chao_shen, na.rm = TRUE), 4),
    SD_Chao_Shen_Entropy = round(sd(spatial_obj$entropy_chao_shen, na.rm = TRUE), 4),
    Median_Chao_Shen_Entropy = round(median(spatial_obj$entropy_chao_shen, na.rm = TRUE), 4),
    IQR_Chao_Shen_Entropy = round(IQR(spatial_obj$entropy_chao_shen, na.rm = TRUE), 4),
    Min_Chao_Shen_Entropy = round(min(spatial_obj$entropy_chao_shen, na.rm = TRUE), 4),
    Max_Chao_Shen_Entropy = round(max(spatial_obj$entropy_chao_shen, na.rm = TRUE), 4),
    Range_Chao_Shen_Entropy = round(diff(range(spatial_obj$entropy_chao_shen, na.rm = TRUE)), 4),
    stringsAsFactors = FALSE
  )
  write.csv(qc_df, file.path(entropy_dir, paste0(sample_name, "_qc_metrics.csv")), row.names = FALSE)
  cat(sprintf("Logged QC metrics to %s\n", file.path(entropy_dir, paste0(sample_name, "_qc_metrics.csv"))))

  # Spatial entropy visualization plot
  entropy_panel <- function(col, label) {
    suppressMessages(
      SpatialFeaturePlot(spatial_obj, features = col) +
        scale_fill_viridis_c(option = "magma", name = paste0("Entropy\n(", label, ")"))
    ) +
      ggtitle(label) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
  }

  p_raw <- (entropy_panel("entropy_chao_shen", "Chao-Shen") |
              entropy_panel("entropy_raw_plugin", "Plug-in")) +
    plot_annotation(
      title = paste("Spatial Shannon Entropy -", sample_name),
      subtitle = "Chao-Shen is the primary metric (D2); the plug-in baseline is carried for comparison",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 10))
    )

  ggsave(file.path(entropy_dir, paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 14, height = 7, dpi = 300)

  # Normalization QC & Covariate checks: evaluate both plug-in baselines
  corr_res <- calculate_entropy_correlations(
    seurat_obj = spatial_obj,
    entropy_cols = c("entropy_chao_shen", "entropy_raw_plugin", "entropy_spanorm_plugin"),
    sample_name = sample_name,
    output_dir = stat_dir
  )

  save_path <- file.path(save_dir, paste0(sample_name, "_spatial_obj.rds"))
  saveRDS(spatial_obj, file = save_path)
  cat(sprintf("Saved spatial object to %s\n", save_path))

  rm(spatial_obj, p_raw)
  gc()

  cat("Done:", sample_name, "\n")
  cat("==========================================\n")
  invisible(corr_res)
}

args <- commandArgs(trailingOnly = TRUE)
target_sample <- resolve_target_sample(args)

analyze_sample(target_sample)
