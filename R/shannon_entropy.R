# R/shannon_entropy.R
# Per-spot transcriptional-diversity estimation for spatial transcriptomics.
#
# Provides:
# 1. exclude_gene_families() — drops the MT/ribosomal rows before estimation.
# 2. entropy_estimator_matrix() — the single column-blocked sparse kernel. One
#    pass over the matrix yields the per-spot count statistics (N, K, f1, f2)
#    and every requested estimator, because all four estimators are functions of
#    the same four statistics plus the observed proportions.
# 3. plugin_entropy_matrix() / calculate_shannon_entropy() — the plug-in
#    estimator H_j = -sum_i p_ij log2(p_ij), kept as-is so that the existing
#    pipeline columns (`entropy_raw_plugin`, `entropy_spanorm_plugin`, D2) and
#    the off-pipeline scripts that call them do not change.
# 4. calculate_entropy() — Seurat wrapper over the kernel, taking one or more of
#    "plugin", "chao_shen", "chao_shen_cwj", "miller_madow" (D2, Stage A1).
# 5. check_entropy_ordering() — the ordering invariant used as a bug check.
#
# All estimators are in bits (log base 2) and are evaluated on integer counts.
# Pointing them at a normalized layer is a category error, not a variant: f1 and
# f2 are counts of singleton and doubleton *reads*, and are meaningless once the
# matrix has been rescaled. The kernel asserts integrality rather than assuming
# it, because that exact mistake is the origin of the `entropy_spanorm_plugin`
# defect recorded in D2.

suppressPackageStartupMessages({
  library(Matrix)
})

#' Estimators supported by entropy_estimator_matrix() / calculate_entropy()
ENTROPY_ESTIMATORS <- c("plugin", "chao_shen", "chao_shen_cwj", "miller_madow")

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

#' Coerce an expression matrix to a structurally-clean dgCMatrix
#'
#' `Matrix::drop0()` is not cosmetic here. Seurat's `subset()` and the gene
#' universe filter both leave structurally-stored explicit zeros behind, and
#' `diff(@p)` counts stored entries rather than non-zero entries — so K, f1 and
#' f2 are silently wrong without it, in a way that biases exactly the shallow
#' spots the correction is supposed to rescue.
#'
#' @param expr_mat Gene x spot matrix
#' @return A dgCMatrix with no explicitly-stored zeros
as_clean_dgc <- function(expr_mat) {
  if (!inherits(expr_mat, "dgCMatrix")) {
    expr_mat <- methods::as(
      methods::as(methods::as(expr_mat, "dMatrix"), "generalMatrix"),
      "CsparseMatrix"
    )
  }
  Matrix::drop0(expr_mat)
}

#' Per-spot diversity estimators in one column-blocked pass
#'
#' For each spot the kernel needs four statistics — N (total counts), K
#' (observed genes), f1 (singletons), f2 (doubletons) — all of which come off
#' the sparse slots directly, plus the observed proportions p_i = y_i / N. Every
#' supported estimator is a function of those, so they are computed together:
#' the four-estimator bake-off (Stage A3) and the subsampling ladder (Stage A4)
#' would otherwise make four full passes over an ~18,500 x ~11,700 matrix for
#' statistics that are shared.
#'
#' Estimators, all in bits:
#'
#' * `plugin` — H = -sum_i p_i log2(p_i) over observed genes. Downward-biased at
#'   finite depth; carried as the baseline (D2).
#' * `chao_shen` — Chao & Shen (2003). Coverage C = 1 - f1/N (Good-Turing),
#'   adjusted probabilities p~_i = C * p_i, and a Horvitz-Thompson inclusion
#'   correction that up-weights each observed gene by one over its probability
#'   of having been seen at all:
#'       H = -sum_i [ p~_i log2(p~_i) ] / [ 1 - (1 - p~_i)^N ]
#' * `chao_shen_cwj` — the same estimator with the Chao-Wang-Jost (2013)
#'   coverage term, which uses the doubletons:
#'       C = 1 - (f1/N) * [ (N-1) f1 / ((N-1) f1 + 2 f2) ]
#'   Note for the write-up: this is Chao-Shen with an improved coverage
#'   estimate, *not* the full Chao-Wang-Jost entropy estimator, which is a
#'   distinct and more involved formula. The label must not overclaim.
#' * `miller_madow` — H_plugin + (K - 1) / (2 N ln 2). The ln 2 is not optional:
#'   the correction is derived in nats and the metric is in bits, and stating it
#'   without the conversion overstates it by a factor of 1/ln2.
#'
#' Two numerical points. The inclusion probability is computed as
#' `-expm1(N * log1p(-p~))` rather than `1 - (1-p~)^N`: at p~ ~ 1e-4 and
#' N ~ 1e4 the naive form loses most of its significant digits. And f1 == N
#' (every observed gene a singleton) drives C to zero and the estimator to total
#' degeneracy, so the standard f1 := N - 1 substitution is applied — rare at
#' these depths, but it happens in the shallow tail and an unguarded NaN
#' propagating into a spatial plot is expensive to notice.
#'
#' @param expr_mat Gene x spot matrix of integer counts (dgCMatrix or dense)
#' @param estimators Character vector, any of `ENTROPY_ESTIMATORS`
#' @param check_counts Logical, assert integrality and non-negativity (default TRUE)
#' @param return_stats Logical, also return the per-spot N, K, f1, f2 columns
#' @return data.frame with one row per spot, rownames = colnames(expr_mat)
entropy_estimator_matrix <- function(expr_mat,
                                     estimators = ENTROPY_ESTIMATORS,
                                     check_counts = TRUE,
                                     return_stats = FALSE) {
  estimators <- match.arg(estimators, ENTROPY_ESTIMATORS, several.ok = TRUE)
  expr_mat <- as_clean_dgc(expr_mat)

  ncells <- ncol(expr_mat)
  col_sums <- Matrix::colSums(expr_mat)

  # Seed the frame with the right number of rows even when the matrix carries no
  # column names (the single-vector helper builds one such matrix).
  out <- data.frame(row.names = if (!is.null(colnames(expr_mat))) {
    colnames(expr_mat)
  } else {
    seq_len(ncells)
  })
  if (return_stats) {
    out$N <- col_sums
    out$K <- numeric(ncells)
    out$f1 <- numeric(ncells)
    out$f2 <- numeric(ncells)
  }
  for (e in estimators) out[[e]] <- numeric(ncells)

  nnz_budget <- getOption("entropy.block.nnz", 1e6)
  nnz_per_cell <- max(1, length(expr_mat@x) / max(1, ncells))
  block_cells <- max(1L, min(ncells, as.integer(floor(nnz_budget / nnz_per_cell))))
  ln2 <- log(2)

  for (start in seq.int(1L, ncells, by = block_cells)) {
    cols <- seq.int(start, min(start + block_cells - 1L, ncells))
    nc <- length(cols)
    sub <- expr_mat[, cols, drop = FALSE]
    y <- sub@x

    # Assert rather than assume. Checked per block so the temporaries stay
    # inside the nnz budget instead of allocating a copy of the whole matrix.
    if (check_counts && length(y) > 0) {
      if (anyNA(y) || any(y < 0)) {
        stop("entropy_estimator_matrix(): matrix contains negative or missing values; ",
             "these estimators are defined on raw counts only.")
      }
      if (any(abs(y - round(y)) > 1e-8)) {
        stop("entropy_estimator_matrix(): matrix is not integer-valued. The singleton ",
             "and doubleton counts f1/f2 are undefined on a normalized layer; ",
             "pass the raw `counts` layer (D2).")
      }
    }

    K_blk <- diff(sub@p)
    grp <- rep.int(seq_len(nc), K_blk)
    N_blk <- col_sums[cols]
    f1_blk <- tabulate(grp[y == 1], nbins = nc)
    f2_blk <- tabulate(grp[y == 2], nbins = nc)
    nonempty <- N_blk > 0

    p <- y / N_blk[grp]

    # Per-column sums of a term vector aligned with the sparse slots. Going
    # through colSums() on a copy of the block keeps this numerically identical
    # to plugin_entropy_matrix(), so `entropy_plugin` and the pipeline's
    # `entropy_raw_plugin` column agree exactly rather than approximately.
    # No !is.finite() sweep here on purpose: drop0() guarantees y > 0 and both
    # coverage estimates are strictly positive, so every term is finite by
    # construction. Zeroing defensively would convert a future arithmetic bug
    # into a plausible-looking entropy value instead of a visible NaN.
    col_sum_terms <- function(terms) {
      blk <- sub
      blk@x <- terms
      as.numeric(Matrix::colSums(blk))
    }

    H_plugin <- col_sum_terms(-p * log2(p))

    # Coverage. After the f1 := N - 1 substitution both coverage estimates are
    # strictly positive: Good-Turing gives at least 1/N, and the CWJ bracket is
    # bounded by 1 so C_cwj >= 1 - f1/N by the same argument.
    f1_use <- as.numeric(f1_blk)
    degenerate <- nonempty & (f1_use == N_blk)
    f1_use[degenerate] <- N_blk[degenerate] - 1

    chao_shen_from_coverage <- function(C_col) {
      p_tilde <- C_col[grp] * p
      num <- -p_tilde * log2(p_tilde)
      # Horvitz-Thompson inclusion probability 1 - (1 - p~)^N, in stable form.
      lambda <- -expm1(N_blk[grp] * log1p(-p_tilde))
      terms <- num / lambda
      terms[num == 0] <- 0        # p~ == 1 gives 0/1; keep it exact
      col_sum_terms(terms)
    }

    if ("plugin" %in% estimators) out[["plugin"]][cols] <- H_plugin

    if ("chao_shen" %in% estimators) {
      C_gt <- ifelse(nonempty, 1 - f1_use / N_blk, 0)
      out[["chao_shen"]][cols] <- chao_shen_from_coverage(C_gt)
    }

    if ("chao_shen_cwj" %in% estimators) {
      # bracket = (N-1) f1 / ((N-1) f1 + 2 f2); it is 0/0 when f1 and f2 are both
      # zero, which is the fully-saturated spot, and the right limit there is a
      # zero correction (C = 1) rather than NaN.
      a <- (N_blk - 1) * f1_use
      b <- a + 2 * f2_blk
      bracket <- ifelse(b > 0, a / b, 0)
      C_cwj <- ifelse(nonempty, 1 - (f1_use / N_blk) * bracket, 0)
      out[["chao_shen_cwj"]][cols] <- chao_shen_from_coverage(C_cwj)
    }

    if ("miller_madow" %in% estimators) {
      out[["miller_madow"]][cols] <- ifelse(
        nonempty, H_plugin + (K_blk - 1) / (2 * N_blk * ln2), 0
      )
    }

    if (return_stats) {
      out$K[cols] <- K_blk
      out$f1[cols] <- f1_blk
      out$f2[cols] <- f2_blk
    }
  }

  out
}

#' Diversity estimators for a single count vector
#'
#' Thin convenience wrapper that routes a plain count vector through the
#' production kernel, so validation against a reference implementation exercises
#' the same code path the cohort run uses rather than a parallel transcription
#' of the formulas.
#'
#' @param y Numeric vector of counts (zeros allowed)
#' @param estimators Character vector, any of `ENTROPY_ESTIMATORS`
#' @return Named numeric vector of entropy estimates in bits
entropy_estimates <- function(y, estimators = ENTROPY_ESTIMATORS) {
  m <- Matrix::Matrix(matrix(as.numeric(y), ncol = 1L), sparse = TRUE)
  res <- entropy_estimator_matrix(m, estimators = estimators)
  unlist(res[1, , drop = TRUE])
}

#' Ordering invariant across the estimators
#'
#' Both corrections are upward, so `H_plugin <= H_miller_madow` and
#' `H_plugin <= H_chao_shen` hold by construction and a violation is a bug.
#'
#' `H_miller_madow <= H_chao_shen` is *not* a theorem — the two corrections are
#' derived differently and nothing forces them to order — so it is reported as
#' an observation rather than asserted. It is expected to hold in these data
#' because the spots are heavily undersampled relative to an 18,500-gene
#' universe, which is the regime where Chao-Shen corrects hardest; a spot where
#' it fails is informative about coverage, not broken.
#'
#' @param est data.frame from entropy_estimator_matrix()
#' @param tol Absolute tolerance for the strict invariants
#' @param stop_on_violation Logical, error rather than warn (default TRUE)
#' @return Invisible data.frame of violation counts per comparison
check_entropy_ordering <- function(est, tol = 1e-9, stop_on_violation = TRUE) {
  has <- function(...) all(c(...) %in% colnames(est))
  rows <- list()

  strict <- list(
    c("plugin", "miller_madow"),
    c("plugin", "chao_shen"),
    c("plugin", "chao_shen_cwj")
  )
  for (pair in strict) {
    if (!has(pair[1], pair[2])) next
    n_bad <- sum(est[[pair[2]]] < est[[pair[1]]] - tol, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      Comparison = sprintf("%s <= %s", pair[1], pair[2]),
      Invariant = TRUE, N_Violations = n_bad, N_Spots = nrow(est)
    )
    if (n_bad > 0) {
      msg <- sprintf("check_entropy_ordering(): %d/%d spots violate %s <= %s; both corrections are upward, so this is a bug.",
                     n_bad, nrow(est), pair[1], pair[2])
      if (stop_on_violation) stop(msg) else warning(msg)
    }
  }

  if (has("miller_madow", "chao_shen")) {
    n_bad <- sum(est$chao_shen < est$miller_madow - tol, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      Comparison = "miller_madow <= chao_shen",
      Invariant = FALSE, N_Violations = n_bad, N_Spots = nrow(est)
    )
  }

  invisible(do.call(rbind, rows))
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

#' Calculate one or more diversity estimators across spots
#'
#' Seurat wrapper over `entropy_estimator_matrix()`. Passing several estimators
#' computes them in a single blocked pass over the matrix.
#'
#' Column naming: the default is `entropy_<estimator>`, so `estimator = "plugin"`
#' writes `entropy_plugin`, which is numerically identical to the pipeline's
#' `entropy_raw_plugin` when run on the `counts` layer with the same exclusion
#' pattern but is deliberately *not* given that name here. The rename of the
#' pipeline column to whichever estimator is chosen happens once, at Stage C,
#' together with the seven `*_Raw_Plugin_Entropy` QC distribution columns in D1;
#' silently aliasing it now would put two names on one column mid-bake-off.
#'
#' @param seurat_obj Seurat object containing spatial data
#' @param estimator One or more of `ENTROPY_ESTIMATORS` (default: "plugin")
#' @param assay Assay name (default: NULL -> DefaultAssay)
#' @param layer Layer name to extract counts from (default: "counts")
#' @param col.name Metadata column name(s); default `entropy_<estimator>`
#' @param exclude_pattern Regex pattern of gene families to exclude
#' @param verbose Logical, print progress information (default: TRUE)
#' @return Seurat object with one metadata column per requested estimator
calculate_entropy <- function(seurat_obj,
                              estimator = ENTROPY_ESTIMATORS,
                              assay = NULL,
                              layer = "counts",
                              col.name = NULL,
                              exclude_pattern = "^(MT-|RP[SL])",
                              verbose = TRUE) {
  # match.arg(several.ok = TRUE) returns the whole default vector when the
  # argument is untouched, which would silently write four columns; be explicit.
  estimator <- if (missing(estimator)) {
    "plugin"
  } else {
    match.arg(estimator, ENTROPY_ESTIMATORS, several.ok = TRUE)
  }

  if (is.null(col.name)) {
    col.name <- paste0("entropy_", estimator)
  } else if (length(col.name) != length(estimator)) {
    stop("calculate_entropy(): `col.name` must have one entry per estimator.")
  }

  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)

  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  if (is.null(expr_mat) || nrow(expr_mat) == 0) {
    stop(sprintf("No data found in assay '%s', layer '%s'.", assay, layer))
  }

  expr_mat <- exclude_gene_families(
    expr_mat, exclude_pattern, verbose,
    context = sprintf("entropy estimation (%s layer, %s)", layer, paste(estimator, collapse = ", "))
  )

  est <- entropy_estimator_matrix(expr_mat, estimators = estimator)

  if (verbose) {
    for (i in seq_along(estimator)) {
      cat(sprintf("  %-14s -> %-24s mean %.4f, SD %.4f bits\n",
                  estimator[i], col.name[i],
                  mean(est[[estimator[i]]]), sd(est[[estimator[i]]])))
    }
  }

  for (i in seq_along(estimator)) {
    v <- est[[estimator[i]]]
    names(v) <- rownames(est)
    seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = v, col.name = col.name[i])
  }

  seurat_obj
}
