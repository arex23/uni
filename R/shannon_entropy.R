# R/shannon_entropy.R
# Per-spot Shannon entropy computation for spatial transcriptomics.
#
# Provides:
# 1. rarefied_entropy_matrix() / calculate_rarefied_entropy() — Primary estimator:
#    downsamples the raw counts matrix to a uniform target depth D across n_draws
#    iterations, averaging entropy across draws. Solves the depth-driven support
#    ceiling.
# 2. plugin_entropy_matrix() / calculate_shannon_entropy() — Full-depth plug-in
#    estimator on a specified layer (counts or data), used as the baseline that
#    demonstrates the depth bias.
#
# The matrix-level functions are the implementation; the calculate_*() functions
# are thin Seurat wrappers around them. `sweep_rarefaction_depth.R` uses the
# matrix-level entry points so the depth sweep can run on raw counts, before any
# D-dependent spot filtering has been applied.

suppressPackageStartupMessages({
  library(Matrix)
  library(scuttle)
})

#' Drop excluded gene families from an expression matrix
#'
#' @param expr_mat Gene x spot matrix
#' @param exclude_pattern Regex of gene families to drop (NULL/"" disables)
#' @param verbose Logical, report how many rows were dropped
#' @param context Short string naming the caller, used in the message only
#' @return The matrix with matching rows removed
exclude_gene_families <- function(expr_mat,
                                  exclude_pattern = "^(MT-|RP[SL])",
                                  verbose = TRUE,
                                  context = "Shannon entropy calculation") {
  if (is.null(exclude_pattern) || nchar(exclude_pattern) == 0) {
    return(expr_mat)
  }
  keep_genes <- !grepl(exclude_pattern, rownames(expr_mat), ignore.case = TRUE)
  n_excluded <- sum(!keep_genes)
  if (n_excluded > 0 && verbose) {
    cat(sprintf("Excluding %d MT/ribosomal genes matching '%s' from %s.\n",
                n_excluded, exclude_pattern, context))
  }
  expr_mat[keep_genes, , drop = FALSE]
}

#' Per-spot plug-in Shannon entropy of a matrix
#'
#' Column-blocked so that the `denom`, `p` and `-p log2 p` temporaries are bounded
#' by the nnz budget rather than by nnz(expr_mat); the returned vector is identical
#' to the unblocked computation because column sums are taken over the full matrix
#' first and entropy is a per-column quantity.
#'
#' @param expr_mat Gene x spot matrix (dgCMatrix or dense)
#' @return Named numeric vector of per-spot entropy in bits
plugin_entropy_matrix <- function(expr_mat) {
  col_sums <- Matrix::colSums(expr_mat)

  if (!inherits(expr_mat, "dgCMatrix")) {
    return(apply(expr_mat, 2, function(x) {
      x <- x[x > 0]
      if (length(x) == 0) return(0)
      p <- x / sum(x)
      -sum(p * log2(p))
    }))
  }

  ncells <- ncol(expr_mat)
  full_entropy <- numeric(ncells)
  nnz_budget <- getOption("entropy.block.nnz", 1e6)
  nnz_per_cell <- max(1, length(expr_mat@x) / max(1, ncells))
  block_cells <- max(1L, min(ncells, as.integer(floor(nnz_budget / nnz_per_cell))))

  for (start in seq.int(1L, ncells, by = block_cells)) {
    cols <- seq.int(start, min(start + block_cells - 1L, ncells))
    sub <- expr_mat[, cols, drop = FALSE]

    counts_per_cell <- diff(sub@p)
    denom <- col_sums[cols][rep.int(seq_along(cols), counts_per_cell)]

    p <- sub@x / denom
    valid_p <- p > 0 & !is.na(p)

    mat_entropy <- sub
    mat_entropy@x[!valid_p] <- 0
    mat_entropy@x[valid_p] <- -p[valid_p] * log2(p[valid_p])

    full_entropy[cols] <- Matrix::colSums(mat_entropy)
  }

  names(full_entropy) <- colnames(expr_mat)
  full_entropy
}

#' Per-spot rarefied Shannon entropy of a counts matrix
#'
#' Every spot is downsampled to exactly `depth` counts, `n_draws` times, and the
#' entropies are averaged. The matrix passed in must already have had any excluded
#' gene families removed: the rarefaction depth is defined on the same matrix the
#' entropy is computed from, so that every spot really is standardised to D.
#'
#' @param counts_mat Gene x spot raw counts matrix (dgCMatrix)
#' @param depth Target downsampling depth D (default: 2000)
#' @param n_draws Number of downsampling draws to average over (default: 5)
#' @param seed Random seed, set immediately before the draw loop (default: 23)
#' @param allow_shallow Logical. If FALSE (default), spots with fewer than `depth`
#'   counts are an error: `scuttle::downsampleMatrix` would cap their sampling
#'   proportion at 1.0 and leave them at full depth, silently reintroducing exactly
#'   the depth bias rarefaction exists to remove. Set TRUE only for diagnostics that
#'   deliberately want the capped behaviour.
#' @param verbose Logical, print progress information (default: TRUE)
#' @return Named numeric vector of per-spot rarefied entropy in bits
rarefied_entropy_matrix <- function(counts_mat,
                                    depth = 2000,
                                    n_draws = 5,
                                    seed = 23,
                                    allow_shallow = FALSE,
                                    verbose = TRUE) {
  col_totals <- Matrix::colSums(counts_mat)
  n_shallow <- sum(col_totals < depth)
  if (n_shallow > 0) {
    msg <- sprintf(
      paste0("%d of %d spots have post-exclusion depth < %d (min: %.0f). Rarefaction ",
             "cannot standardise them; filter spots on the post-exclusion column sums, ",
             "not on nCount, before calling this function."),
      n_shallow, length(col_totals), depth, min(col_totals))
    if (!allow_shallow) stop(msg)
    warning(paste(msg, "Sampling proportion capped at 1.0 for these spots."))
  }

  props <- pmin(1.0, depth / col_totals)
  ncells <- ncol(counts_mat)

  if (verbose) {
    cat(sprintf("Calculating rarefied Shannon entropy (depth = %d, n_draws = %d, seed = %d) across %d spots...\n",
                depth, n_draws, seed, ncells))
  }

  # Initialize seed immediately prior to downsampling draws
  set.seed(seed)

  draw_entropies <- matrix(0, nrow = ncells, ncol = n_draws)
  for (draw in seq_len(n_draws)) {
    ds_mat <- scuttle::downsampleMatrix(counts_mat, prop = props, bycol = TRUE)
    draw_entropies[, draw] <- plugin_entropy_matrix(ds_mat)
  }

  mean_entropy <- rowMeans(draw_entropies)
  names(mean_entropy) <- colnames(counts_mat)
  mean_entropy
}

#' Calculate Rarefied Shannon Entropy across Spots
#'
#' Seurat wrapper around rarefied_entropy_matrix(): pulls the counts layer, drops
#' the excluded gene families, rarefies to `depth` and stores the mean entropy.
#'
#' @param seurat_obj Seurat object containing spatial data
#' @param depth Target downsampling depth D (default: 2000)
#' @param n_draws Number of downsampling draws to average across (default: 5)
#' @param seed Random seed initialized immediately prior to downsampling draws (default: 23)
#' @param assay Assay name to extract counts from (default: NULL -> DefaultAssay)
#' @param col.name Metadata column name to store result (default: "entropy_rarefied")
#' @param exclude_pattern Regex pattern of gene families to exclude (default: "^(MT-|RP[SL])")
#' @param allow_shallow Logical, see rarefied_entropy_matrix() (default: FALSE)
#' @param verbose Logical, print progress information (default: TRUE)
#' @return Seurat object with rarefied entropy added to metadata
calculate_rarefied_entropy <- function(seurat_obj,
                                       depth = 2000,
                                       n_draws = 5,
                                       seed = 23,
                                       assay = NULL,
                                       col.name = "entropy_rarefied",
                                       exclude_pattern = "^(MT-|RP[SL])",
                                       allow_shallow = FALSE,
                                       verbose = TRUE) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)

  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
  if (is.null(counts_mat) || nrow(counts_mat) == 0) {
    stop(sprintf("No counts found in assay '%s', layer 'counts'.", assay))
  }

  counts_mat <- exclude_gene_families(counts_mat, exclude_pattern, verbose,
                                      context = "rarefied entropy calculation")

  mean_entropy <- rarefied_entropy_matrix(
    counts_mat,
    depth = depth,
    n_draws = n_draws,
    seed = seed,
    allow_shallow = allow_shallow,
    verbose = verbose
  )

  Seurat::AddMetaData(seurat_obj, metadata = mean_entropy, col.name = col.name)
}

#' Calculate Plug-in Shannon Entropy across Spots
#'
#' Computes spot-level Shannon entropy directly on the expression matrix
#' (raw counts or normalized data) without downsampling.
#'
#' @param seurat_obj Seurat object containing spatial data
#' @param assay Assay name (default: NULL -> DefaultAssay)
#' @param layer Layer name to extract expression from (default: "counts")
#' @param col.name Metadata column name to store result (default: "entropy_raw_plugin")
#' @param exclude_pattern Regex pattern of gene families to exclude (default: "^(MT-|RP[SL])")
#' @param verbose Logical, print progress information (default: TRUE)
#' @return Seurat object with plug-in entropy added to metadata
calculate_shannon_entropy <- function(seurat_obj,
                                      assay = NULL,
                                      layer = "counts",
                                      col.name = "entropy_raw_plugin",
                                      exclude_pattern = "^(MT-|RP[SL])",
                                      verbose = TRUE) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)

  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  if (is.null(expr_mat) || nrow(expr_mat) == 0) {
    stop(sprintf("No data found in assay '%s', layer '%s'.", assay, layer))
  }

  expr_mat <- exclude_gene_families(expr_mat, exclude_pattern, verbose,
                                    context = sprintf("Shannon entropy calculation (%s layer)", layer))

  Seurat::AddMetaData(seurat_obj, metadata = plugin_entropy_matrix(expr_mat), col.name = col.name)
}
