# R/spatial_neighbors.R
# Hex-lattice spatial neighbourhood construction and Moran's I for Visium data.
#
# Built for Stage B2 (is there spatially-structured variation of the kind
# SpaNorm's spline could act on?) and reused unchanged by Stage D / D7, where
# the same weights carry the spatial inference. No new dependencies: the
# alternative is `spdep`, which pulls `sf` and therefore GDAL/GEOS/PROJ — 300-500
# MB of system libraries for one statistic, on a machine already tight enough
# that SpaNorm needed a blocked kernel to fit.
#
# Provides:
# 1. hex_neighbor_graph()  — sparse adjacency on the Visium hex lattice.
# 2. row_standardise()     — W with rows summing to 1.
# 3. morans_i()            — Moran's I with a permutation p-value.
# 4. morans_i_table()      — Moran's I over several metadata columns.

suppressPackageStartupMessages({
  library(Matrix)
})

#' Six-neighbour adjacency on the Visium hexagonal lattice
#'
#' Visium spots sit on a known hexagonal lattice indexed by `array_row` and
#' `array_col`. Within a row `array_col` increments by 2, and adjacent rows are
#' offset by 1, so the six neighbours of `(r, c)` are `(r, c-2)`, `(r, c+2)`,
#' `(r-1, c-1)`, `(r-1, c+1)`, `(r+1, c-1)` and `(r+1, c+1)`.
#'
#' Building the graph from lattice indices rather than from pixel coordinates
#' matters. A k-nearest-neighbour or distance-threshold graph on the imaged
#' coordinates has to pick a k or a radius, and that choice is a free parameter
#' nobody can defend; here the geometry is known exactly and there is nothing to
#' tune. It also degrades correctly: spots dropped by the D1 QC floors simply
#' have fewer neighbours, rather than silently acquiring distant ones as their
#' nearest, which is what a kNN graph does at the edge of a hole.
#'
#' @param array_row Integer vector of lattice row indices, one per spot
#' @param array_col Integer vector of lattice column indices, one per spot
#' @param spot_ids Optional character vector of spot names for dimnames
#' @return Symmetric sparse `dgCMatrix` of 0/1 adjacency
hex_neighbor_graph <- function(array_row, array_col, spot_ids = NULL) {
  n <- length(array_row)
  if (length(array_col) != n) {
    stop("hex_neighbor_graph(): array_row and array_col must have the same length.")
  }
  if (anyNA(array_row) || anyNA(array_col)) {
    stop("hex_neighbor_graph(): array_row/array_col contain NA. Spots without ",
         "lattice coordinates cannot be placed on the graph; drop them first.")
  }

  r <- as.integer(array_row)
  c_ <- as.integer(array_col)
  key <- paste(r, c_, sep = ":")
  if (anyDuplicated(key)) {
    stop("hex_neighbor_graph(): duplicated (array_row, array_col) pairs; ",
         "two spots cannot occupy one lattice position.")
  }
  lookup <- setNames(seq_len(n), key)

  # The six lattice offsets. Only three are enumerated: each undirected edge
  # would otherwise be found twice, once from each endpoint.
  offsets <- list(c(0L, 2L), c(1L, -1L), c(1L, 1L))

  i_idx <- integer(0)
  j_idx <- integer(0)
  for (off in offsets) {
    nb_key <- paste(r + off[1], c_ + off[2], sep = ":")
    hit <- lookup[nb_key]
    found <- !is.na(hit)
    i_idx <- c(i_idx, which(found))
    j_idx <- c(j_idx, as.integer(hit[found]))
  }

  adj <- Matrix::sparseMatrix(
    i = c(i_idx, j_idx),
    j = c(j_idx, i_idx),
    x = 1,
    dims = c(n, n),
    dimnames = if (is.null(spot_ids)) NULL else list(spot_ids, spot_ids)
  )
  methods::as(adj, "generalMatrix")
}

#' Build the hex adjacency from a Seurat object's metadata
#'
#' @param seurat_obj Seurat object carrying `array_row` / `array_col` (D1)
#' @return Symmetric sparse adjacency over the object's spots, in column order
hex_neighbor_graph_from_seurat <- function(seurat_obj) {
  meta <- seurat_obj@meta.data
  if (!all(c("array_row", "array_col") %in% colnames(meta))) {
    stop("hex_neighbor_graph_from_seurat(): object has no array_row/array_col in ",
         "metadata. These are written by analyze_entropy.R / load_qc_sample() (D1).")
  }
  hex_neighbor_graph(meta$array_row, meta$array_col, rownames(meta))
}

#' Row-standardise a weights matrix
#'
#' Rows sum to 1, except isolated spots (no neighbours), whose rows stay all
#' zero. Isolated spots then contribute nothing to the numerator of Moran's I
#' but still contribute their own variance to the denominator, which is the
#' conservative direction: an isolated spot cannot manufacture autocorrelation.
#'
#' @param adj Sparse adjacency matrix
#' @return Row-standardised sparse matrix
row_standardise <- function(adj) {
  rs <- Matrix::rowSums(adj)
  rs[rs == 0] <- 1
  # Diagonal() %*% drops dimnames; restore them so callers can index by spot id.
  W <- Matrix::Diagonal(x = 1 / rs) %*% adj
  dimnames(W) <- dimnames(adj)
  W
}

#' Moran's I with a permutation p-value
#'
#' The general form is `I = (n / S0) * (z' W z) / (z' z)` with `z` the centred
#' variable and `S0 = sum(W)`. Under row-standardisation `S0` equals the number
#' of *non-isolated* spots, so the `n / S0` factor drops out only when every
#' spot has at least one neighbour — which is the usual case on an intact
#' lattice but not guaranteed once the D1 QC floors have punched holes in it.
#' The factor is therefore carried explicitly rather than assumed away; it
#' cancels in the permutation p-value but not in the reported I.
#'
#' The permutation null is the right one here: the analytic normal
#' approximation assumes the variable is normally distributed, which
#' `percent.mt` and the entropy metrics are not.
#'
#' Note what the permutation tests. Spot values are reshuffled across a *fixed*
#' lattice, so the null is "this variable's values are arranged at random over
#' this tissue's geometry" — the geometry, including any holes left by the D1 QC
#' floors, is held constant and is not part of what is being tested.
#'
#' @param x Numeric vector, one value per spot, in the same order as `adj`
#' @param adj Sparse adjacency (unstandardised); row-standardised internally
#' @param n_perm Number of permutations (default 999)
#' @param seed Integer seed, or NULL to leave the RNG alone
#' @return List with observed I, expected I under the null, permutation mean/SD,
#'   a pseudo z-score, the two-sided permutation p-value and the counts behind it
morans_i <- function(x, adj, n_perm = 999, seed = 23) {
  if (length(x) != nrow(adj)) {
    stop("morans_i(): length(x) does not match nrow(adj).")
  }
  keep <- is.finite(x)
  if (!all(keep)) {
    adj <- adj[keep, keep, drop = FALSE]
    x <- x[keep]
  }
  n <- length(x)

  out <- list(
    N = n, I = NA_real_, S0 = NA_real_,
    Expected_I = if (n > 1) -1 / (n - 1) else NA_real_,
    Perm_Mean = NA_real_, Perm_SD = NA_real_, Z = NA_real_,
    P_Perm = NA_real_, N_Perm = n_perm, N_Extreme = NA_integer_,
    N_Isolated = NA_integer_
  )
  if (n < 3) return(out)

  W <- row_standardise(adj)
  out$N_Isolated <- sum(Matrix::rowSums(adj) == 0)
  out$S0 <- sum(Matrix::rowSums(W))

  z <- x - mean(x)
  denom <- sum(z * z)
  if (denom <= 0) return(out)          # constant variable: I undefined

  S0 <- sum(Matrix::rowSums(W))
  if (S0 <= 0) return(out)
  scale_factor <- n / S0
  moran <- function(zz) {
    scale_factor * as.numeric(crossprod(zz, W %*% zz)) / sum(zz * zz)
  }
  out$I <- moran(z)

  if (n_perm > 0) {
    if (!is.null(seed)) set.seed(seed)
    perm <- vapply(seq_len(n_perm), function(i) moran(z[sample.int(n)]), numeric(1))
    out$Perm_Mean <- mean(perm)
    out$Perm_SD <- sd(perm)
    out$Z <- if (out$Perm_SD > 0) (out$I - out$Perm_Mean) / out$Perm_SD else NA_real_
    # Two-sided, counting the observed value in both numerator and denominator:
    # the standard (n_extreme + 1) / (n_perm + 1) construction, which cannot
    # return 0 and so never claims more evidence than the permutations support.
    n_extreme <- sum(abs(perm - out$Perm_Mean) >= abs(out$I - out$Perm_Mean))
    out$N_Extreme <- n_extreme
    out$P_Perm <- (n_extreme + 1) / (n_perm + 1)
  }

  out
}

#' Moran's I over several metadata columns of a Seurat object
#'
#' @param seurat_obj Seurat object with `array_row` / `array_col` in metadata
#' @param cols Character vector of metadata column names to test
#' @param sample_name Sample label written into the output
#' @param n_perm Number of permutations (default 999)
#' @param seed Integer seed
#' @param adj Optional precomputed adjacency, to avoid rebuilding it per column
#' @return data.frame, one row per column
morans_i_table <- function(seurat_obj, cols, sample_name = "Sample",
                           n_perm = 999, seed = 23, adj = NULL) {
  if (is.null(adj)) adj <- hex_neighbor_graph_from_seurat(seurat_obj)
  meta <- seurat_obj@meta.data
  cols <- cols[cols %in% colnames(meta)]

  rows <- lapply(cols, function(cl) {
    res <- morans_i(as.numeric(meta[[cl]]), adj, n_perm = n_perm, seed = seed)
    data.frame(
      Sample = sample_name, Variable = cl, N = res$N,
      Morans_I = res$I, Expected_I = res$Expected_I, S0 = res$S0,
      Perm_Mean = res$Perm_Mean, Perm_SD = res$Perm_SD, Z_Perm = res$Z,
      P_Perm = res$P_Perm, N_Perm = res$N_Perm, N_Isolated = res$N_Isolated,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Neighbour-count summary for the lattice graph
#'
#' Reported alongside Moran's I as a sanity check on the geometry. On an intact
#' Visium lattice the interior is uniformly 6-connected and only the tissue
#' boundary has fewer; a large isolated-spot count means the lattice
#' coordinates are wrong, not that the tissue is unusual.
#'
#' @param adj Sparse adjacency matrix
#' @param sample_name Sample label written into the output
#' @return One-row data.frame
neighbor_summary <- function(adj, sample_name = "Sample") {
  deg <- Matrix::rowSums(adj)
  data.frame(
    Sample = sample_name, N_Spots = length(deg),
    Mean_Neighbors = mean(deg), Median_Neighbors = median(deg),
    Min_Neighbors = min(deg), Max_Neighbors = max(deg),
    N_Isolated = sum(deg == 0), Pct_Fully_Connected = 100 * mean(deg == 6),
    N_Edges = sum(deg) / 2, stringsAsFactors = FALSE
  )
}
