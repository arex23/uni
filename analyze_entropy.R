set.seed(23)

library(Seurat)
library(SpatialExperiment)
library(SpaNorm)
library(ggplot2)
library(patchwork)

source("R/shannon_entropy.R")
source("R/entropy_correlation.R")

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
  spe[["Spatial"]]$data <- as(norm_mat, "CsparseMatrix")
  return(spe)
})

samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
samples <- samples[samples != ""]

align_spots_to_coords <- function(seurat_obj) {
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  coords_df <- Seurat::GetTissueCoordinates(seurat_obj)
  coords_mat <- as.matrix(coords_df[, 1:2])

  common_spots <- intersect(colnames(counts_mat), rownames(coords_mat))
  if (length(common_spots) > 0 && length(common_spots) < ncol(seurat_obj)) {
    cat(sprintf("Aligning spots: subsetting from %d to %d spots with valid coordinates.\n",
                ncol(seurat_obj), length(common_spots)))
    seurat_obj <- subset(seurat_obj, cells = common_spots)
  }
  return(seurat_obj)
}

analyze_sample <- function(sample_name) {
  cat("==========================================\n")
  cat("Processing sample:", sample_name, "\n")

  data_dir <- file.path("data", sample_name, "raw_data")
  if (!dir.exists(data_dir)) {
    stop(sprintf("Directory '%s' does not exist for sample '%s'", data_dir, sample_name))
  }

  h5_files <- list.files(path = data_dir, pattern = "\\.h5$", full.names = FALSE)
  if (length(h5_files) != 1) {
    stop(sprintf("Expected exactly 1 .h5 file in %s for sample '%s', found %d: %s",
                 data_dir, sample_name, length(h5_files), paste(h5_files, collapse = ", ")))
  }
  h5_file <- h5_files[1]

  spatial_obj <- Load10X_Spatial(
    data.dir = data_dir,
    filename = h5_file,
    assay = "Spatial",
    filter.matrix = TRUE
  )

  # Calculate mitochondrial and ribosomal percentages
  spatial_obj[["percent.mt"]] <- PercentageFeatureSet(spatial_obj, pattern = "^MT-")
  spatial_obj[["percent.ribo"]] <- PercentageFeatureSet(spatial_obj, pattern = "^RP[SL]")

  n_raw_spots <- ncol(spatial_obj)
  n_raw_genes <- nrow(spatial_obj)

  # 1. Spot filtering by sequencing depth and feature count
  valid_spots <- colnames(spatial_obj)[spatial_obj$nCount_Spatial >= 500 & spatial_obj$nFeature_Spatial >= 250]
  spatial_obj <- subset(spatial_obj, cells = valid_spots)
  n_spots_post_qc <- ncol(spatial_obj)

  # 2. Align spots with available tissue coordinates
  spatial_obj <- align_spots_to_coords(spatial_obj)
  n_spots_post_align <- ncol(spatial_obj)

  # 3. Gene filtering computed AFTER subsetting to valid spots
  raw_counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  min_spots_per_gene <- max(20, ceiling(0.02 * ncol(spatial_obj)))
  expressed_genes <- rownames(raw_counts)[Matrix::rowSums(raw_counts > 0) >= min_spots_per_gene]

  spatial_obj <- subset(spatial_obj, features = expressed_genes)
  n_genes_post_filter <- nrow(spatial_obj)
  rm(raw_counts, valid_spots, expressed_genes)
  gc()

  # 4. Calculate raw Shannon entropy (excluding MT and ribosomal genes)
  spatial_obj <- calculate_shannon_entropy(
    spatial_obj,
    assay = "Spatial",
    layer = "counts",
    col.name = "shannon_entropy_raw",
    exclude_pattern = "^(MT-|RP[SL])"
  )

  # 5. SpaNorm normalization using direct public API
  spatial_obj <- SpaNorm::SpaNorm(
    spatial_obj,
    sample.p = 0.25,
    gene.model = "nb",
    adj.method = "logpac",
    df.tps = 6,
    lambda.a = 1e-04,
    verbose = TRUE
  )

  # 6. Calculate normalized Shannon entropy (excluding MT and ribosomal genes)
  spatial_obj <- calculate_shannon_entropy(
    spatial_obj,
    assay = "Spatial",
    layer = "data",
    col.name = "shannon_entropy",
    exclude_pattern = "^(MT-|RP[SL])"
  )

  entropy_dir <- file.path("results", "analyze_entropy")
  stat_dir <- file.path("results", "statistical_tests")
  save_dir <- file.path("results", "seurat_objects")

  if (!dir.exists(entropy_dir)) dir.create(entropy_dir, recursive = TRUE)
  if (!dir.exists(stat_dir)) dir.create(stat_dir, recursive = TRUE)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # Save per-sample QC filtering log table
  qc_df <- data.frame(
    Sample = sample_name,
    Raw_Spots = n_raw_spots,
    Raw_Genes = n_raw_genes,
    Spots_Post_QC = n_spots_post_qc,
    Spots_Post_Coord_Align = n_spots_post_align,
    Genes_Post_Filter = n_genes_post_filter,
    Pct_Spots_Retained = round((n_spots_post_align / n_raw_spots) * 100, 2),
    Pct_Genes_Retained = round((n_genes_post_filter / n_raw_genes) * 100, 2),
    Mean_Percent_MT = round(mean(spatial_obj$percent.mt, na.rm = TRUE), 2),
    Mean_Percent_Ribo = round(mean(spatial_obj$percent.ribo, na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
  write.csv(qc_df, file.path(entropy_dir, paste0(sample_name, "_qc_metrics.csv")), row.names = FALSE)
  cat(sprintf("Logged QC metrics to %s\n", file.path(entropy_dir, paste0(sample_name, "_qc_metrics.csv"))))

  # Spatial entropy visualization plot
  p_raw <- suppressMessages(
    SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
      scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy")
  ) +
    ggtitle(paste("Spatial Distribution of Shannon Entropy -", sample_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  ggsave(file.path(entropy_dir, paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 8, height = 7, dpi = 300)

  # Normalization QC & Covariate checks: correlations against nCount, nFeature, percent.mt, percent.ribo
  corr_res <- calculate_entropy_correlations(
    seurat_obj = spatial_obj,
    entropy_cols = c("shannon_entropy_raw", "shannon_entropy"),
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
if (length(args) > 0) {
  target_sample <- args[1]
  if (!target_sample %in% samples) {
    stop(sprintf("Sample '%s' not found in data/. Available samples: %s", target_sample, paste(samples, collapse = ", ")))
  }
} else {
  target_sample <- samples[1]
}

analyze_sample(target_sample)
