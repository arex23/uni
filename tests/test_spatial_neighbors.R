# tests/test_spatial_neighbors.R
#
# Checks for the hex-lattice spatial machinery in `R/spatial_neighbors.R`
# (Stage B1 of docs/ESTIMATOR_PIVOT_PLAN.md; reused by Stage D / D7).
#
# Run from the repository root:
#     Rscript tests/test_spatial_neighbors.R
#
# Moran's I is validated against a naive double-loop evaluation of the textbook
# formula rather than against `spdep`. That is deliberate and is the same trade
# made elsewhere in this project: `spdep` pulls `sf` and therefore GDAL/GEOS/PROJ,
# 300-500 MB of system libraries, to check one statistic. The general formula
#
#     I = (n / S0) * sum_ij w_ij (x_i - xbar)(x_j - xbar) / sum_i (x_i - xbar)^2
#
# is short enough to transcribe independently, and a naive O(n^2) loop shares no
# code with the sparse implementation, so agreement between them is real evidence
# rather than a tautology.

source("R/spatial_neighbors.R")

TOL <- 1e-10
n_pass <- 0L
n_fail <- 0L

check <- function(label, ok, detail = "") {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
  } else {
    n_fail <<- n_fail + 1L
    cat(sprintf("FAIL  %s%s\n", label, if (nchar(detail)) paste0("  [", detail, "]") else ""))
  }
}

check_close <- function(label, actual, expected, tol = TOL) {
  d <- max(abs(actual - expected))
  check(label, is.finite(d) && d <= tol, sprintf("max abs diff %.3e", d))
}

# A rectangular patch of the Visium lattice: array_col increments by 2 within a
# row and adjacent rows are offset by 1, so odd rows carry odd column indices.
make_lattice <- function(n_rows, n_cols) {
  g <- do.call(rbind, lapply(seq_len(n_rows) - 1L, function(r) {
    data.frame(array_row = r, array_col = (seq_len(n_cols) - 1L) * 2L + (r %% 2L))
  }))
  g$id <- paste0("spot", seq_len(nrow(g)))
  g
}

# ---------------------------------------------------------------------------
# 1. Lattice geometry
# ---------------------------------------------------------------------------

cat("--- 1. hex lattice geometry ---\n")

g <- make_lattice(9L, 9L)
adj <- hex_neighbor_graph(g$array_row, g$array_col, g$id)

check("adjacency is symmetric", Matrix::isSymmetric(adj))
check("no self-loops", all(Matrix::diag(adj) == 0))
check("adjacency is binary", all(adj@x == 1))
check("dimnames carried through", identical(rownames(adj), g$id))

deg <- Matrix::rowSums(adj)
# Interior spots — those with a full ring of six lattice positions occupied.
interior <- g$array_row > 0 & g$array_row < 8 &
  g$array_col > 1 & g$array_col < max(g$array_col) - 1
check("every interior spot has exactly 6 neighbours", all(deg[interior] == 6),
      sprintf("degrees seen: %s", paste(sort(unique(deg[interior])), collapse = ",")))
check("no spot exceeds 6 neighbours", max(deg) == 6)
check("no isolated spots on an intact lattice", sum(deg == 0) == 0)

# The six offsets, verified explicitly on one interior spot rather than only in
# aggregate: a transposed or sign-flipped offset would still give degree 6.
centre <- which(g$array_row == 4 & g$array_col == 8)
nb <- g[which(adj[centre, ] > 0), c("array_row", "array_col")]
expected_nb <- data.frame(
  array_row = c(4L, 4L, 3L, 3L, 5L, 5L),
  array_col = c(6L, 10L, 7L, 9L, 7L, 9L)
)
key <- function(d) sort(paste(d$array_row, d$array_col, sep = ":"))
check("neighbour offsets are (r,c+-2) and (r+-1,c+-1)", identical(key(nb), key(expected_nb)))

# A hole in the lattice must reduce the degree of its neighbours, not silently
# reassign more distant spots as adjacent the way a kNN graph would.
g_hole <- g[-centre, ]
adj_hole <- hex_neighbor_graph(g_hole$array_row, g_hole$array_col, g_hole$id)
ring <- which(paste(g_hole$array_row, g_hole$array_col, sep = ":") %in%
                paste(expected_nb$array_row, expected_nb$array_col, sep = ":"))
check("removing a spot lowers its neighbours' degree by exactly 1",
      all(Matrix::rowSums(adj_hole)[ring] == deg[which(g$id %in% g_hole$id[ring])] - 1))

# A lone spot far from the patch must be isolated, not attached to the nearest.
g_far <- rbind(g[, c("array_row", "array_col", "id")],
               data.frame(array_row = 100L, array_col = 100L, id = "far"))
adj_far <- hex_neighbor_graph(g_far$array_row, g_far$array_col, g_far$id)
check("a distant spot is isolated, not nearest-neighbour attached",
      Matrix::rowSums(adj_far)["far"] == 0)

check("duplicate lattice positions are rejected",
      inherits(try(hex_neighbor_graph(c(1, 1), c(2, 2)), silent = TRUE), "try-error"))
check("NA lattice coordinates are rejected",
      inherits(try(hex_neighbor_graph(c(1, NA), c(2, 4)), silent = TRUE), "try-error"))

# ---------------------------------------------------------------------------
# 2. Row standardisation
# ---------------------------------------------------------------------------

cat("--- 2. row standardisation ---\n")

W <- row_standardise(adj)
check_close("connected rows sum to 1", Matrix::rowSums(W), rep(1, nrow(W)))

W_far <- row_standardise(adj_far)
rs_far <- Matrix::rowSums(W_far)
check("isolated rows stay all-zero rather than dividing by zero",
      rs_far[["far"]] == 0 && all(is.finite(rs_far)))

# ---------------------------------------------------------------------------
# 3. Moran's I against an independent implementation
# ---------------------------------------------------------------------------

cat("--- 3. Moran's I ---\n")

# Naive O(n^2) transcription of the textbook formula. Shares no code with
# morans_i(): dense, double-looped, and built from the unstandardised adjacency.
morans_i_naive <- function(x, adj) {
  n <- length(x)
  A <- as.matrix(adj)
  Wd <- A
  for (i in seq_len(n)) {
    s <- sum(A[i, ])
    if (s > 0) Wd[i, ] <- A[i, ] / s
  }
  z <- x - mean(x)
  num <- 0
  for (i in seq_len(n)) for (j in seq_len(n)) num <- num + Wd[i, j] * z[i] * z[j]
  (n / sum(Wd)) * num / sum(z^2)
}

set.seed(23)
x_rand <- rnorm(nrow(g))
x_grad <- g$array_row + 0.25 * g$array_col              # smooth spatial gradient
x_check <- ifelse((g$array_row + g$array_col %/% 2) %% 2 == 0, 1, -1)  # alternating

for (nm in c("random", "gradient", "checker")) {
  xv <- switch(nm, random = x_rand, gradient = x_grad, checker = x_check)
  check_close(sprintf("Moran's I matches naive formula (%s)", nm),
              morans_i(xv, adj, n_perm = 0)$I, morans_i_naive(xv, adj))
}

# Same agreement with an isolated spot present, where the n/S0 factor no longer
# cancels. This is the case the plan's "I reduces to z'Wz / z'z" shortcut gets
# wrong, so it is checked explicitly.
x_far <- c(x_rand, 0)
res_far <- morans_i(x_far, adj_far, n_perm = 0)
check_close("Moran's I matches naive formula with an isolated spot",
            res_far$I, morans_i_naive(x_far, adj_far))
check("n/S0 factor is not 1 when a spot is isolated", abs(res_far$S0 - length(x_far)) > 0.5)

# Sign and magnitude behaviour.
I_grad <- morans_i(x_grad, adj, n_perm = 0)$I
I_rand <- morans_i(x_rand, adj, n_perm = 0)$I
I_check <- morans_i(x_check, adj, n_perm = 0)$I
check("smooth gradient gives strong positive I", I_grad > 0.8, sprintf("I = %.4f", I_grad))
check("alternating pattern gives negative I", I_check < 0, sprintf("I = %.4f", I_check))
check("random noise sits near the null expectation", abs(I_rand) < 0.15, sprintf("I = %.4f", I_rand))

# ---------------------------------------------------------------------------
# 4. Permutation inference
# ---------------------------------------------------------------------------

cat("--- 4. permutation p-values ---\n")

res_grad <- morans_i(x_grad, adj, n_perm = 999, seed = 23)
res_rand <- morans_i(x_rand, adj, n_perm = 999, seed = 23)

check("permutation p is bounded below by 1/(n_perm+1)", res_grad$P_Perm >= 1 / 1000 - TOL)
check("permutation p never reaches 0", res_grad$P_Perm > 0)
check("strong structure is significant", res_grad$P_Perm <= 0.01, sprintf("p = %.4f", res_grad$P_Perm))
check("random noise is not significant", res_rand$P_Perm > 0.05, sprintf("p = %.4f", res_rand$P_Perm))
check("permutation mean sits near -1/(n-1)",
      abs(res_grad$Perm_Mean - res_grad$Expected_I) < 0.02,
      sprintf("mean %.4f vs E[I] %.4f", res_grad$Perm_Mean, res_grad$Expected_I))
check("seeding makes the p-value reproducible",
      identical(morans_i(x_grad, adj, n_perm = 199, seed = 7)$P_Perm,
                morans_i(x_grad, adj, n_perm = 199, seed = 7)$P_Perm))
check("permutation does not change the observed I",
      abs(res_grad$I - morans_i(x_grad, adj, n_perm = 0)$I) < TOL)

# Degenerate inputs must return NA rather than a number.
check("constant variable returns NA", is.na(morans_i(rep(1, nrow(g)), adj, n_perm = 0)$I))
res_na <- morans_i(replace(x_rand, 1:5, NA), adj, n_perm = 0)
check("non-finite values are dropped, not propagated",
      is.finite(res_na$I) && res_na$N == nrow(g) - 5)

# ---------------------------------------------------------------------------
# 5. Summaries
# ---------------------------------------------------------------------------

cat("--- 5. summaries ---\n")

ns <- neighbor_summary(adj, "test")
check("neighbour summary reports the intact lattice", ns$N_Isolated == 0 && ns$Max_Neighbors == 6)
check("edge count matches sum(deg)/2", ns$N_Edges == sum(Matrix::rowSums(adj)) / 2)

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
