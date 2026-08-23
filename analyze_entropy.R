library(Seurat)
library(ggplot2)
library(patchwork)
library(SpaNorm)

source("shannon_entropy.R")
source("entropy_correlation.R")

samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
samples <- samples[samples != ""]

apply_spanorm_normalization <- function(seurat_obj, sample_p = 0.1, chunk_size = 1000) {
  cat("Applying memory-efficient chunked SpaNorm normalization...\n")

  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  coords_df <- Seurat::GetTissueCoordinates(seurat_obj)
  coords_mat <- as.matrix(coords_df[, 1:2])

  common_spots <- intersect(colnames(counts_mat), rownames(coords_mat))
  if (length(common_spots) > 0) {
    counts_mat <- counts_mat[, common_spots, drop = FALSE]
    coords_mat <- coords_mat[common_spots, , drop = FALSE]
    
    # CRITICAL FIX: Subset the actual Seurat object here so metadata and all downstream 
    # steps (entropy, visualization, DE analysis) operate on this exact consistent spot set.
    seurat_obj <- subset(seurat_obj, cells = common_spots)
  }

  cat("Fitting SpaNorm model...\n")
  fit_spanorm <- SpaNorm:::fitSpaNorm(
    Y = counts_mat,
    coords = coords_mat,
    sample.p = sample_p,
    gene.model = "nb",
    df.tps = 6,
    lambda.a = 1e-04,
    batch = NULL,
    LS = NULL,
    msgfun = message
  )

  cat("Normalizing data in chunks...\n")
  ngenes <- nrow(counts_mat)
  sparse_chunks <- list()

  isbio <- fit_spanorm$wtype %in% "biology"
  W_bio <- fit_spanorm$W[, isbio, drop = FALSE]

  for (start_idx in seq(1, ngenes, by = chunk_size)) {
    end_idx <- min(start_idx + chunk_size - 1, ngenes)
    chunk_idx <- start_idx:end_idx

    Y_chunk <- counts_mat[chunk_idx, , drop = FALSE]
    gmean_chunk <- fit_spanorm$gmean[chunk_idx]
    alpha_chunk <- fit_spanorm$alpha[chunk_idx, , drop = FALSE]
    alpha_bio_chunk <- alpha_chunk[, isbio, drop = FALSE]
    psi_chunk <- fit_spanorm$psi[chunk_idx]

    mu_chunk <- gmean_chunk + Matrix::tcrossprod(alpha_chunk, fit_spanorm$W)
    lmu_max <- matrixStats::rowMedians(as.matrix(mu_chunk)) + 4 * matrixStats::rowMads(as.matrix(mu_chunk))
    mu_chunk <- exp(pmin(as.matrix(mu_chunk), lmu_max))

    mu2_chunk <- gmean_chunk + Matrix::tcrossprod(alpha_bio_chunk, W_bio)
    lmu_max2 <- matrixStats::rowMedians(as.matrix(mu2_chunk)) + 4 * matrixStats::rowMads(as.matrix(mu2_chunk))
    mu2_chunk <- exp(pmin(as.matrix(mu2_chunk), lmu_max2))

    psi_max <- exp(median(log(fit_spanorm$psi)) + 3 * mad(log(fit_spanorm$psi)))
    psi_chunk <- pmin(psi_chunk, psi_max)

    Y_dense <- as.matrix(Y_chunk)
    lb <- pnbinom(Y_dense - 1, mu = mu_chunk, size = 1 / psi_chunk)
    ub <- dnbinom(Y_dense, mu = mu_chunk, size = 1 / psi_chunk) + lb
    p <- (lb + ub) / 2
    p <- pmax(pmin(p, 0.999), 0.001)

    normmat_chunk <- log2(qnbinom(p, mu = mu2_chunk, size = 1 / psi_chunk) + 1)
    rownames(normmat_chunk) <- rownames(Y_chunk)
    colnames(normmat_chunk) <- colnames(Y_chunk)

    sparse_chunks[[length(sparse_chunks) + 1]] <- as(normmat_chunk, "CsparseMatrix")

    cat(sprintf("\rNormalized %d / %d genes", end_idx, ngenes))
    gc()
  }
  cat("\n")

  norm_sparse <- do.call(rbind, sparse_chunks)
  rm(sparse_chunks, counts_mat, coords_mat, fit_spanorm)
  gc()

  seurat_obj[["Spatial"]]$data <- norm_sparse

  rm(norm_sparse)
  gc()

  return(seurat_obj)
}

analyze_sample <- function(sample_name) {
  cat("Processing sample:", sample_name, "\n")

  data_dir <- file.path("data", sample_name, "raw_data")
  if (!dir.exists(data_dir)) data_dir <- file.path("data", sample_name)

  h5_file <- list.files(path = data_dir, pattern = "\\.h5$", full.names = FALSE)[1]
  spatial_obj <- Load10X_Spatial(
    data.dir = data_dir,
    filename = h5_file,
    assay = "Spatial",
    filter.matrix = TRUE
  )

  valid_spots <- colnames(spatial_obj)[spatial_obj$nCount_Spatial >= 500 & spatial_obj$nFeature_Spatial >= 250]
  raw_counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  min_spots_per_gene <- max(20, ceiling(0.02 * ncol(spatial_obj)))
  expressed_genes <- rownames(raw_counts)[Matrix::rowSums(raw_counts > 0) >= min_spots_per_gene]

  spatial_obj <- subset(spatial_obj, cells = valid_spots, features = expressed_genes)
  rm(raw_counts, valid_spots, expressed_genes)
  gc()

  spatial_obj <- apply_spanorm_normalization(spatial_obj)
  spatial_obj <- calculate_shannon_entropy(spatial_obj, assay = "Spatial", layer = "data")

  p_raw <- SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy") +
    ggtitle(paste("Spatial Distribution of Shannon Entropy -", sample_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  he_dir <- file.path("results", "he_graphs")
  stat_dir <- file.path("results", "statistical_tests")
  save_dir <- file.path("results", "seurat_objects")

  if (!dir.exists(he_dir)) dir.create(he_dir, recursive = TRUE)
  if (!dir.exists(stat_dir)) dir.create(stat_dir, recursive = TRUE)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  ggsave(file.path(he_dir, paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 8, height = 7, dpi = 300)

  corr_res <- calculate_entropy_correlations(spatial_obj, sample_name = sample_name, output_dir = stat_dir)

  save_path <- file.path(save_dir, paste0(sample_name, "_spatial_obj.rds"))
  saveRDS(spatial_obj, file = save_path)
  cat(sprintf("Saved spatial object to %s\n", save_path))

  rm(spatial_obj, p_raw)
  gc()

  cat("Done:", sample_name, "\n")
  invisible(corr_res)
}

args <- commandArgs(trailingOnly = TRUE)
target_sample <- if (length(args) > 0 && args[1] %in% samples) args[1] else samples[1]

analyze_sample(target_sample)
