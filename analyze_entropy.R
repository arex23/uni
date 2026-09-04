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
  if ("in_tissue" %in% colnames(tp_df)) {
    ontissue_barcodes <- tp_df$barcode[tp_df$in_tissue == 1]
    bc_col <- "barcode"
    row_col <- "array_row"
    col_col <- "array_col"
  } else {
    ontissue_barcodes <- tp_df[[1]][tp_df[[2]] == 1]
    bc_col <- colnames(tp_df)[1]
    row_col <- colnames(tp_df)[3]
    col_col <- colnames(tp_df)[4]
  }

  spatial_obj <- subset(spatial_obj, cells = intersect(colnames(spatial_obj), ontissue_barcodes))
  n_spots_ontissue <- ncol(spatial_obj)

  # Add Visium hexagonal array lattice coordinates to metadata for spatial modeling (Stage 4)
  matched_idx <- match(colnames(spatial_obj), tp_df[[bc_col]])
  spatial_obj$array_row <- tp_df[[row_col]][matched_idx]
  spatial_obj$array_col <- tp_df[[col_col]][matched_idx]

  # Calculate mitochondrial and ribosomal percentages on-tissue
  spatial_obj[["percent.mt"]] <- PercentageFeatureSet(spatial_obj, pattern = "^MT-")
  spatial_obj[["percent.ribo"]] <- PercentageFeatureSet(spatial_obj, pattern = "^RP[SL]")

  # 1. Spot filtering by sequencing depth (nCount >= 2000, matching downsampling depth D = 2000)
  # and feature count (nFeature >= 250) on on-tissue spots
  target_depth_D <- 2000
  entropy_exclude_pattern <- "^(MT-|RP[SL])"
  valid_spots <- colnames(spatial_obj)[spatial_obj$nCount_Spatial >= target_depth_D & spatial_obj$nFeature_Spatial >= 250]
  spatial_obj <- subset(spatial_obj, cells = valid_spots)
  n_spots_post_depth_qc <- ncol(spatial_obj)

  # 2. Validate/align spots with imaged tissue coordinates.
  spatial_obj <- validate_imaged_coordinates(spatial_obj)
  n_spots_post_coord <- ncol(spatial_obj)

  # 3. Gene universe filtering (frozen cohort gene universe, D1/D5).
  # The per-sample spot detection threshold is dropped; rarefaction standardises
  # depth across spots and rare genes contribute minimally to entropy.
  spatial_obj <- filter_by_gene_universe(spatial_obj)
  n_genes_in_universe <- nrow(spatial_obj)

  raw_counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  gene_totals <- Matrix::rowSums(raw_counts > 0)
  n_genes_detected <- sum(gene_totals > 0)
  rm(raw_counts, valid_spots, gene_totals)
  gc()

  # 3b. Second depth gate, on the matrix the entropy is actually computed from.
  # The nCount >= D gate in step 1 is taken over all 37,082 features, but entropy
  # is computed after the universe filter and the MT/ribo exclusion, so a spot can
  # pass step 1 and still hold fewer than D counts here. Rarefaction cannot standardise
  # such a spot -- downsampleMatrix would cap its sampling proportion at 1.0 and
  # leave it at full depth -- so it is dropped rather than silently left unrarefied.
  entropy_counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  entropy_counts <- entropy_counts[!grepl(entropy_exclude_pattern, rownames(entropy_counts), ignore.case = TRUE), , drop = FALSE]
  entropy_depth <- Matrix::colSums(entropy_counts)
  deep_spots <- colnames(spatial_obj)[entropy_depth >= target_depth_D]
  n_spots_dropped_entropy_depth <- ncol(spatial_obj) - length(deep_spots)
  cat(sprintf("Entropy-matrix depth QC: dropping %d spots below D = %d after gene filtering (min depth %.0f).\n",
              n_spots_dropped_entropy_depth, target_depth_D, min(entropy_depth)))
  spatial_obj <- subset(spatial_obj, cells = deep_spots)
  n_spots_final <- ncol(spatial_obj)
  rm(entropy_counts, entropy_depth, deep_spots)
  gc()

  # 4. Calculate primary Rarefied Shannon entropy (depth D = 2000, n_draws = 5).
  # allow_shallow stays FALSE: after 3b every spot has >= D counts in this matrix,
  # and the function errors rather than silently capping if that ever stops holding.
  spatial_obj <- calculate_rarefied_entropy(
    spatial_obj,
    depth = target_depth_D,
    n_draws = 5,
    seed = 23,
    col.name = "entropy_rarefied",
    exclude_pattern = entropy_exclude_pattern,
    allow_shallow = FALSE
  )

  # 5. Calculate baseline full-depth raw plug-in Shannon entropy (demonstrates depth bias)
  spatial_obj <- calculate_shannon_entropy(
    spatial_obj,
    assay = "Spatial",
    layer = "counts",
    col.name = "entropy_raw_plugin",
    exclude_pattern = entropy_exclude_pattern
  )

  # 6. SpaNorm normalization using direct public API (used for downstream DE, scoring, viz)
  spatial_obj <- SpaNorm::SpaNorm(
    spatial_obj,
    sample.p = 0.25,
    gene.model = "nb",
    adj.method = "logpac",
    df.tps = 6,
    lambda.a = 1e-04,
    verbose = TRUE
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
    Spots_Post_Coord_Validation = n_spots_post_coord,
    Spots_Dropped_Entropy_Depth = n_spots_dropped_entropy_depth,
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
    Mean_Rarefied_Entropy = round(mean(spatial_obj$entropy_rarefied, na.rm = TRUE), 4),
    SD_Rarefied_Entropy = round(sd(spatial_obj$entropy_rarefied, na.rm = TRUE), 4),
    Median_Rarefied_Entropy = round(median(spatial_obj$entropy_rarefied, na.rm = TRUE), 4),
    IQR_Rarefied_Entropy = round(IQR(spatial_obj$entropy_rarefied, na.rm = TRUE), 4),
    Min_Rarefied_Entropy = round(min(spatial_obj$entropy_rarefied, na.rm = TRUE), 4),
    Max_Rarefied_Entropy = round(max(spatial_obj$entropy_rarefied, na.rm = TRUE), 4),
    Range_Rarefied_Entropy = round(diff(range(spatial_obj$entropy_rarefied, na.rm = TRUE)), 4),
    stringsAsFactors = FALSE
  )
  write.csv(qc_df, file.path(entropy_dir, paste0(sample_name, "_qc_metrics.csv")), row.names = FALSE)
  cat(sprintf("Logged QC metrics to %s\n", file.path(entropy_dir, paste0(sample_name, "_qc_metrics.csv"))))

  # Spatial rarefied entropy visualization plot
  p_raw <- suppressMessages(
    SpatialFeaturePlot(spatial_obj, features = "entropy_rarefied") +
      scale_fill_viridis_c(option = "magma", name = "Rarefied\nEntropy\n(D=2000)")
  ) +
    ggtitle(paste("Spatial Distribution of Rarefied Shannon Entropy (D=2000) -", sample_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )

  ggsave(file.path(entropy_dir, paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 8, height = 7, dpi = 300)

  # Normalization QC & Covariate checks: evaluate both entropy_rarefied and entropy_raw_plugin
  corr_res <- calculate_entropy_correlations(
    seurat_obj = spatial_obj,
    entropy_cols = c("entropy_rarefied", "entropy_raw_plugin"),
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
