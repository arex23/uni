# R/load_sample.R
# D1 sample loading and spot/feature quality control, factored out of
# `analyze_entropy.R` so that anything needing a QC'd Seurat object can get one
# without also running SpaNorm.
#
# Why this exists: the Stage A estimator bake-off reads raw counts only, but the
# diagnostic that runs it needs exactly the spot set the pipeline analyses -- the
# same in-tissue subset, the same nCount/nFeature floors, the same gene universe.
# Before this extraction the only way to get that was to run `analyze_entropy.R`,
# which fits SpaNorm at ~7.6 GB peak and takes hours, to produce a normalized
# layer the diagnostic then ignores. The QC steps here are moved verbatim from
# `analyze_entropy.R`, not rewritten; `analyze_entropy.R` calls this module, so
# there is one implementation of D1's QC rather than two that can drift.
#
# Provides:
# 1. validate_imaged_coordinates() -- align spots to the imaged coordinates.
# 2. load_qc_sample()              -- D1 steps 1-3, returning the object plus
#                                     the QC counters the per-sample table logs.

suppressPackageStartupMessages({
  library(Seurat)
})

#' Subset a Seurat object to spots that carry valid imaged coordinates
#'
#' @param seurat_obj Seurat object
#' @param verbose Logical, report the subsetting
#' @return The object restricted to spots present in the tissue coordinates
validate_imaged_coordinates <- function(seurat_obj, verbose = TRUE) {
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  coords_df <- Seurat::GetTissueCoordinates(seurat_obj)
  coords_mat <- as.matrix(coords_df[, 1:2])

  common_spots <- intersect(colnames(counts_mat), rownames(coords_mat))
  if (length(common_spots) > 0 && length(common_spots) < ncol(seurat_obj)) {
    if (verbose) {
      cat(sprintf("Subsetting from %d to %d spots with valid imaged coordinates.\n",
                  ncol(seurat_obj), length(common_spots)))
    }
    seurat_obj <- subset(seurat_obj, cells = common_spots)
  }
  return(seurat_obj)
}

#' Load one sample and apply the D1 spot and feature quality control
#'
#' Steps, in the order D1 fixes them:
#'   0. `Load10X_Spatial(filter.matrix = TRUE)` -- which filters the *image*, not
#'      the counts matrix, so the object still holds every grid barcode.
#'   1. Explicit subset to `in_tissue == 1`, before any depth QC, so the
#'      retention denominators reflect tissue area rather than slide background.
#'      Visium lattice coordinates (`array_row`, `array_col`) are carried into
#'      metadata here for the spatial neighbourhood construction (D1, B1/D7).
#'   2. `percent.mt` / `percent.ribo`, then the `nCount >= 500` /
#'      `nFeature >= 250` floors. These are plain quality floors and are
#'      deliberately not tied to any estimator's target depth, so changing the
#'      entropy formulation cannot silently change which tissue is analysed.
#'   3. Coordinate validation, then the frozen cohort gene universe (D5). No
#'      per-sample detection threshold, so every sample shares one feature space.
#'
#' @param sample_name Sample directory name under `data/`
#' @param min_counts Minimum nCount per spot (D1 default 500)
#' @param min_features Minimum nFeature per spot (D1 default 250)
#' @param data_dir Directory holding one subdirectory per sample
#' @param verbose Logical, print progress
#' @return List with `obj` (the QC'd Seurat object) and `qc` (named list of the
#'   counters `analyze_entropy.R` writes into `<sample>_qc_metrics.csv`)
load_qc_sample <- function(sample_name,
                           min_counts = 500,
                           min_features = 250,
                           data_dir = "data",
                           verbose = TRUE) {
  sample_dir <- file.path(data_dir, sample_name)
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

  # Add Visium hexagonal array lattice coordinates to metadata for spatial modeling
  matched_idx <- match(colnames(spatial_obj), tp_df[[bc_col]])
  spatial_obj$array_row <- tp_df[[row_col]][matched_idx]
  spatial_obj$array_col <- tp_df[[col_col]][matched_idx]

  # Calculate mitochondrial and ribosomal percentages on-tissue
  spatial_obj[["percent.mt"]] <- PercentageFeatureSet(spatial_obj, pattern = "^MT-")
  spatial_obj[["percent.ribo"]] <- PercentageFeatureSet(spatial_obj, pattern = "^RP[SL]")

  # 1. Spot filtering by sequencing depth and feature count on on-tissue spots.
  valid_spots <- colnames(spatial_obj)[spatial_obj$nCount_Spatial >= min_counts &
                                         spatial_obj$nFeature_Spatial >= min_features]
  spatial_obj <- subset(spatial_obj, cells = valid_spots)
  n_spots_post_depth_qc <- ncol(spatial_obj)

  # 2. Validate/align spots with imaged tissue coordinates.
  spatial_obj <- validate_imaged_coordinates(spatial_obj, verbose = verbose)
  n_spots_post_coord <- ncol(spatial_obj)

  # 3. Gene universe filtering (frozen cohort gene universe, D1/D5).
  # The per-sample spot detection threshold is dropped so that every sample is
  # processed on the exact same feature space, which is what makes the cohort
  # comparable; rare genes contribute minimally to entropy either way.
  spatial_obj <- filter_by_gene_universe(spatial_obj, verbose = verbose)
  n_genes_in_universe <- nrow(spatial_obj)
  n_spots_final <- ncol(spatial_obj)

  raw_counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  gene_totals <- Matrix::rowSums(raw_counts > 0)
  n_genes_detected <- sum(gene_totals > 0)
  rm(raw_counts, valid_spots, gene_totals)
  gc(verbose = FALSE)

  list(
    obj = spatial_obj,
    qc = list(
      raw_spots_grid = raw_spots_grid,
      n_raw_genes = n_raw_genes,
      n_spots_ontissue = n_spots_ontissue,
      n_spots_post_depth_qc = n_spots_post_depth_qc,
      n_spots_post_coord = n_spots_post_coord,
      n_spots_final = n_spots_final,
      n_genes_in_universe = n_genes_in_universe,
      n_genes_detected = n_genes_detected
    )
  )
}
