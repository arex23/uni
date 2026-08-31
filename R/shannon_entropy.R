calculate_shannon_entropy <- function(seurat_obj,
                                      assay = NULL,
                                      layer = "data",
                                      col.name = "shannon_entropy",
                                      exclude_pattern = "^(MT-|RP[SL])") {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)

  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  if (is.null(expr_mat) || nrow(expr_mat) == 0) {
    stop(sprintf("No data found in assay '%s', layer '%s'.", assay, layer))
  }

  # Exclude mitochondrial and ribosomal gene families from entropy calculation
  if (!is.null(exclude_pattern) && nchar(exclude_pattern) > 0) {
    keep_genes <- !grepl(exclude_pattern, rownames(expr_mat), ignore.case = TRUE)
    n_excluded <- sum(!keep_genes)
    if (n_excluded > 0) {
      cat(sprintf("Excluding %d MT/ribosomal genes matching '%s' from Shannon entropy calculation (%s layer).\n",
                  n_excluded, exclude_pattern, layer))
      expr_mat <- expr_mat[keep_genes, , drop = FALSE]
    }
  }

  col_sums <- Matrix::colSums(expr_mat)

  if (inherits(expr_mat, "dgCMatrix")) {
    # Entropy is a per-spot quantity, so the sparse branch runs over blocks of
    # columns rather than the whole matrix at once. Done unblocked, each of
    # `denom`, `p`, `valid_p`, `!valid_p` and the `-p * log2(p)` temporaries is
    # as long as nnz(expr_mat) (~76M for sample1, ~0.6 GB per double vector) and
    # `mat_entropy` is a second full copy of the matrix, which put the peak for
    # this function near 5 GB above baseline. The arithmetic per column is
    # unchanged -- `col_sums` is still computed once over the full matrix -- so
    # the returned vector is identical, block size only bounds the temporaries.
    ncells <- ncol(expr_mat)
    full_entropy <- numeric(ncells)
    nnz_budget <- getOption("entropy.block.nnz", 1e6)
    nnz_per_cell <- max(1, length(expr_mat@x) / max(1, ncells))
    block.cells <- max(1L, min(ncells, as.integer(floor(nnz_budget / nnz_per_cell))))

    for (start in seq.int(1L, ncells, by = block.cells)) {
      cols <- seq.int(start, min(start + block.cells - 1L, ncells))
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
  } else {
    full_entropy <- apply(expr_mat, 2, function(x) {
      x <- x[x > 0]
      if (length(x) == 0) return(0)
      p <- x / sum(x)
      -sum(p * log2(p))
    })
  }

  Seurat::AddMetaData(seurat_obj, metadata = full_entropy, col.name = col.name)
}
