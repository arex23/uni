# Memory-bounded drop-in replacement for SpaNorm's `logpac` adjustment step.
#
# `SpaNorm:::normaliseLogPAC()` builds the whole ngenes x ncells adjustment in
# one shot. For sample1 (15,556 genes x 11,709 spots) a single dense double
# matrix is 1.36 GB, and the expression keeps roughly ten of them alive at once
# (mu, mu.2, two densified copies of the counts, lb, ub, p and the pmin/pmax,
# qnbinom and log2 temporaries), peaking near 13.6 GB. That is what pushes a
# 15.1 GB machine into the OOM killer.
#
# Given the fitted model, every operation in that function is row-wise: row g of
# the output depends only on row g of the counts, gmean[g], alpha[g, ], psi[g]
# and the shared design matrix W. The only quantity that crosses genes is
# `psi.max`, a median/MAD over the full psi vector, and it is computed here
# from the complete vector before any chunking happens. Splitting the work into
# blocks of genes is therefore exact, not an approximation: the same arithmetic
# is applied to the same values in the same order, and the results are written
# into a preallocated output matrix.
#
# Peak memory becomes the output matrix plus a handful of block-sized
# temporaries instead of ten full-sized ones.

# Row-wise reimplementation of SpaNorm's internal calculateMu(). Kept local so
# the pipeline does not reach into the package namespace for it.
.spanorm_mu_block <- function(gmean, alpha, W) {
  mu <- gmean + Matrix::tcrossprod(alpha, W)
  lmu.max <- matrixStats::rowMedians(mu) + 4 * matrixStats::rowMads(mu)
  exp(pmin(mu, lmu.max))
}

# Signature matches SpaNorm:::normaliseLogPAC() so this can stand in for it.
normaliseLogPAC_lowmem <- function(Y, scale.factor, fit.spanorm) {
  W <- fit.spanorm$W
  alpha <- fit.spanorm$alpha
  gmean <- fit.spanorm$gmean

  isbio <- fit.spanorm$wtype %in% "biology"
  W.bio <- W[, isbio, drop = FALSE]
  alpha.bio <- alpha[, isbio, drop = FALSE]

  # Cross-gene step: computed once on the full psi vector, exactly as upstream.
  psi <- fit.spanorm$psi
  psi.max <- exp(stats::median(log(psi)) + 3 * stats::mad(log(psi)))
  psi <- pmin(psi, psi.max)

  ngenes <- nrow(Y)
  ncells <- ncol(Y)

  # Block size is set from a target element budget so it adapts to spot count:
  # 3e6 doubles is ~24 MB per temporary, ~200 MB for a whole block iteration.
  budget <- getOption("spanorm.block.elements", 3e6)
  block.genes <- max(1L, min(ngenes, as.integer(floor(budget / ncells))))

  normmat <- matrix(0, nrow = ngenes, ncol = ncells)

  for (start in seq.int(1L, ngenes, by = block.genes)) {
    idx <- seq.int(start, min(start + block.genes - 1L, ngenes))
    size <- 1 / psi[idx]

    mu <- .spanorm_mu_block(gmean[idx], alpha[idx, , drop = FALSE], W)
    mu.2 <- .spanorm_mu_block(gmean[idx], alpha.bio[idx, , drop = FALSE], W.bio)

    yb <- as.matrix(Y[idx, , drop = FALSE])
    lb <- stats::pnbinom(yb - 1, mu = mu, size = size)
    ub <- stats::dnbinom(yb, mu = mu, size = size) + lb
    p <- (lb + ub) / 2
    p <- pmax(pmin(p, 0.999), 0.001)

    normmat[idx, ] <- log2(stats::qnbinom(p, mu = scale.factor * mu.2, size = size) + 1)
  }

  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)
  normmat
}

# Swap the chunked kernel into the SpaNorm namespace. The public
# SpaNorm::SpaNorm() entry point, the model fit and the returned values are all
# untouched; only the adjustment step's memory profile changes.
enable_spanorm_lowmem <- function(verbose = TRUE) {
  if (!requireNamespace("SpaNorm", quietly = TRUE)) {
    stop("SpaNorm is not installed.")
  }
  current <- utils::getFromNamespace("normaliseLogPAC", "SpaNorm")
  if (identical(current, normaliseLogPAC_lowmem)) {
    return(invisible(FALSE))
  }
  utils::assignInNamespace("normaliseLogPAC", normaliseLogPAC_lowmem, ns = "SpaNorm")
  if (verbose) {
    cat("SpaNorm logpac adjustment patched with gene-blocked low-memory kernel.\n")
  }
  invisible(TRUE)
}
