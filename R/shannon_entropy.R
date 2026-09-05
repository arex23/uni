# R/shannon_entropy.R
# Per-spot Shannon entropy computation for spatial transcriptomics.
#
# Provides:
# 1. exclude_gene_families() — drops the MT/ribosomal rows before estimation.
# 2. plugin_entropy_matrix() / calculate_shannon_entropy() — plug-in estimator
#    H_j = -sum_i p_ij log2(p_ij) with p_ij = y_ij / sum_i y_ij, evaluated on a
#    specified layer (raw `counts` or the SpaNorm-normalized `data`), in
#    column-blocked sparse form to bound peak memory.
#
# plugin_entropy_matrix() is the implementation; calculate_shannon_entropy() is a
# thin Seurat wrapper around it that pulls a layer and writes a metadata column.

suppressPackageStartupMessages({
  library(Matrix)
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
