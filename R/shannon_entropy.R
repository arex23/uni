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
    counts_per_cell <- diff(expr_mat@p)
    denom <- col_sums[rep.int(seq_len(ncol(expr_mat)), counts_per_cell)]

    p <- expr_mat@x / denom
    valid_p <- p > 0 & !is.na(p)

    mat_entropy <- expr_mat
    mat_entropy@x[!valid_p] <- 0
    mat_entropy@x[valid_p] <- -p[valid_p] * log2(p[valid_p])

    full_entropy <- Matrix::colSums(mat_entropy)
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
