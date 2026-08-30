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
samples <- samples[grepl("sample[0-9]+$", samples)]
sample_nums <- as.numeric(gsub(".*sample", "", samples))
samples <- samples[order(sample_nums)]

validate_imaged_coordinates <- function(seurat_obj) {
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  coords_df <- Seurat::GetTissueCoordinates(seurat_obj)
  coords_mat <- as.matrix(coords_df[, 1:2])

  common_spots <- intersect(colnames(counts_mat), rownames(coords_mat))
  if (length(common_spots) > 0 && length(common_spots) < ncol(seurat_obj)) {
    cat(sprintf("Subsetting from %d to %d spots with valid imaged coordinates.\n",
                ncol(seurat_obj), length(common_spots)))
    seurat_obj <- subset(seurat_obj, cells = common_spots)
  }
  return(seurat_obj)
}

analyze_sample <- function(sample_name) {
  cat("==========================================\n")
  cat("Processing sample:", sample_name, "\n")

  sample_dir <- file.path("data", sample_name)
  if (!dir.exists(sample_dir)) {
    stop(sprintf("Directory '%s' does not exist for sample '%s'", sample_dir, sample_name))
  }

  h5_files <- list.files(path = sample_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  if (length(h5_files) == 0) {
    stop(sprintf("No .h5 file found in %s for sample '%s'", sample_dir, sample_name))
  }
  h5_path <- h5_files[1]

  tp_files <- list.files(path = sample_dir, pattern = "tissue_positions.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(tp_files) == 0) {
    stop(sprintf("No tissue_positions.csv found in %s for sample '%s'", sample_dir, sample_name))
  }
  tp_file <- tp_files[1]

  spatial_dir <- dirname(dirname(tp_file))
  h5_rel <- if (dirname(h5_path) == spatial_dir) basename(h5_path) else file.path("..", basename(h5_path))

  spatial_obj <- Load10X_Spatial(
    data.dir = spatial_dir,
    filename = h5_rel,
    assay = "Spatial",
    filter.matrix = TRUE
  )

  raw_spots_grid <- ncol(spatial_obj)
  n_raw_genes <- nrow(spatial_obj)

  # Load tissue positions and explicitly subset to in_tissue == 1 spots
  tp_df <- read.csv(tp_file)
  ontissue_barcodes <- if ("in_tissue" %in% colnames(tp_df)) {
    tp_df$barcode[tp_df$in_tissue == 1]
  } else {
    tp_df[[1]][tp_df[[2]] == 1]
  }

  spatial_obj <- subset(spatial_obj, cells = intersect(colnames(spatial_obj), ontissue_barcodes))
  n_spots_ontissue <- ncol(spatial_obj)

  # Calculate mitochondrial and ribosomal percentages on-tissue
  spatial_obj[["percent.mt"]] <- PercentageFeatureSet(spatial_obj, pattern = "^MT-")
  spatial_obj[["percent.ribo"]] <- PercentageFeatureSet(spatial_obj, pattern = "^RP[SL]")

  # 1. Spot filtering by sequencing depth and feature count on on-tissue spots
  valid_spots <- colnames(spatial_obj)[spatial_obj$nCount_Spatial >= 500 & spatial_obj$nFeature_Spatial >= 250]
  spatial_obj <- subset(spatial_obj, cells = valid_spots)
  n_spots_post_depth_qc <- ncol(spatial_obj)

  # 2. Validate/align spots with imaged tissue coordinates.
  # After the explicit in_tissue subset above this should be a no-op; it is kept
  # as a safety net and its result is recorded so the QC table cannot drift.
  spatial_obj <- validate_imaged_coordinates(spatial_obj)
  n_spots_final <- ncol(spatial_obj)

  # 3. Gene filtering computed AFTER subsetting to valid on-tissue spots
  raw_counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  min_spots_per_gene <- max(20, ceiling(0.02 * ncol(spatial_obj)))
  gene_totals <- Matrix::rowSums(raw_counts > 0)
  # Genes with at least one count in at least one retained spot. The probe panel
  # emits the full GRCh38 feature list, so ~18,100 of the 37,082 rows are
  # structurally zero (no probe targets them) and are not real gene loss.
  n_genes_detected <- sum(gene_totals > 0)
  expressed_genes <- rownames(raw_counts)[gene_totals >= min_spots_per_gene]

  spatial_obj <- subset(spatial_obj, features = expressed_genes)
  n_genes_post_filter <- nrow(spatial_obj)
  rm(raw_counts, valid_spots, expressed_genes, gene_totals)
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

  # Save per-sample QC filtering log table with honest on-tissue accounting
  qc_df <- data.frame(
    Sample = sample_name,
    Raw_Spots_Grid = raw_spots_grid,
    Spots_On_Tissue = n_spots_ontissue,
    Spots_Post_Depth_QC = n_spots_post_depth_qc,
    Spots_Final = n_spots_final,
    Pct_OnTissue_Retained = round((n_spots_final / n_spots_ontissue) * 100, 2),
    Raw_Genes = n_raw_genes,
    Genes_Detected = n_genes_detected,
    Genes_Post_Filter = n_genes_post_filter,
    Pct_Genes_Retained_OfAll = round((n_genes_post_filter / n_raw_genes) * 100, 2),
    Pct_Genes_Retained_OfDetected = round((n_genes_post_filter / n_genes_detected) * 100, 2),
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

  # Normalization QC & Covariate checks: correlations against nCount, nFeature, percent.mt.
  # percent.ribo is not a target -- it is identically zero cohort-wide (D2); the
  # Mean_Percent_Ribo column above is retained as the record of that zero.
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
