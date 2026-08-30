#!/usr/bin/env Rscript

# diagnose_entropy_scaling.R
# Stage 2 Diagnostic: Disentangle (1) log-transform category error from
# (2) depth-driven support ceiling in spatial Shannon entropy estimation.

set.seed(23)

library(Seurat)
library(ggplot2)
library(patchwork)

source("R/shannon_entropy.R")
source("R/entropy_correlation.R")

args <- commandArgs(trailingOnly = TRUE)
target_sample <- if (length(args) > 0) args[1] else "sample1"

seurat_file <- file.path("results", "seurat_objects", paste0(target_sample, "_spatial_obj.rds"))
if (!file.exists(seurat_file)) {
  stop(sprintf("Seurat object '%s' not found. Please run analyze_entropy.R on '%s' first.", seurat_file, target_sample))
}

cat("=================================================================\n")
cat("Running Stage 2 Entropy Scaling Diagnostic for:", target_sample, "\n")
cat("=================================================================\n")

spatial_obj <- readRDS(seurat_file)
cat("Loaded Seurat object with", ncol(spatial_obj), "spots and", nrow(spatial_obj), "features.\n")

counts_mat <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
data_log_mat <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "data")

# 1. Per-spot detected-gene support, on the same MT/ribo-excluded gene set the
#    entropy uses. Two supports exist and they are NOT interchangeable:
#      n_detected_counts - genes with a non-zero raw count (what rarefaction acts on)
#      n_detected_norm   - genes with a non-zero SpaNorm logpac value. SpaNorm's
#                          qnbinom step imputes genes with zero raw counts and
#                          zeroes genes with non-zero counts, so this support is
#                          much narrower than the raw one.
#    Entropy computed on the `data` layer sums over n_detected_norm, so that is
#    the support the log-scale degeneracy must be tested against.
keep_genes <- !grepl("^(MT-|RP[SL])", rownames(counts_mat), ignore.case = TRUE)
counts_filtered <- counts_mat[keep_genes, , drop = FALSE]
data_filtered <- data_log_mat[keep_genes, , drop = FALSE]

spatial_obj$n_detected_counts <- Matrix::colSums(counts_filtered > 0)
spatial_obj$n_detected_norm <- Matrix::colSums(data_filtered > 0)
spatial_obj$log2_n_detected_norm <- log2(spatial_obj$n_detected_norm)

# 2. Raw plug-in entropy and 3. log-scale SpaNorm entropy (the existing pipeline
#    category error) are already columns on the saved object, computed by
#    analyze_entropy.R with this same exclude_pattern. Reuse them so the
#    diagnostic reports exactly what the pipeline shipped; recompute only if a
#    future object arrives without them.
if (!"shannon_entropy_raw" %in% colnames(spatial_obj@meta.data)) {
  spatial_obj <- calculate_shannon_entropy(
    spatial_obj,
    assay = "Spatial",
    layer = "counts",
    col.name = "shannon_entropy_raw",
    exclude_pattern = "^(MT-|RP[SL])"
  )
}

if ("shannon_entropy" %in% colnames(spatial_obj@meta.data)) {
  spatial_obj$shannon_entropy_log <- spatial_obj$shannon_entropy
} else {
  spatial_obj <- calculate_shannon_entropy(
    spatial_obj,
    assay = "Spatial",
    layer = "data",
    col.name = "shannon_entropy_log",
    exclude_pattern = "^(MT-|RP[SL])"
  )
}

# 4. Linear back-transformed SpaNorm entropy: 2^x - 1.
#    SpaNorm's logpac adjustment is log2(qnbinom(...) + 1) (SpaNorm:::normaliseLogPAC),
#    so 2^x - 1 is its exact inverse, not an approximation. Mutating @x alone is
#    safe because x = 0 maps to 0, leaving structural zeros untouched.
data_lin_mat <- data_log_mat
data_lin_mat@x <- (2^(data_log_mat@x)) - 1
spatial_obj[["Spatial"]]$linear_data <- data_lin_mat

spatial_obj <- calculate_shannon_entropy(
  spatial_obj,
  assay = "Spatial",
  layer = "linear_data",
  col.name = "shannon_entropy_linear",
  exclude_pattern = "^(MT-|RP[SL])"
)

out_dir <- file.path("results", "statistical_tests")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Targets for diagnostic evaluation. log2 is carried only for the normalized
# support: a log2 transform is monotone, so its Spearman rho is identical to the
# untransformed column by construction and only the Pearson r adds information.
targets_diag <- c(
  nCounts = "nCount_Spatial",
  nFeatures = "nFeature_Spatial",
  n_detected_counts = "n_detected_counts",
  n_detected_norm = "n_detected_norm",
  log2_n_detected_norm = "log2_n_detected_norm"
)

entropy_metrics_diag <- c(
  "shannon_entropy_raw",
  "shannon_entropy_log",
  "shannon_entropy_linear"
)

diag_res <- calculate_entropy_correlations(
  seurat_obj = spatial_obj,
  entropy_cols = entropy_metrics_diag,
  targets = targets_diag,
  sample_name = target_sample,
  output_dir = out_dir,
  file_suffix = "_entropy_scaling_diagnostic",
  save_outputs = TRUE
)

# 5. Degeneracy fit. If the log transform has flattened p_i towards 1/K, entropy
#    collapses onto log2(K) with unit slope. Emitting the fit (rather than
#    asserting it in the docs) makes the H ~ log2(K) claim reproducible.
log2_K <- spatial_obj$log2_n_detected_norm
degeneracy_rows <- lapply(entropy_metrics_diag, function(e_col) {
  H <- spatial_obj@meta.data[[e_col]]
  fit <- lm(H ~ log2_K)
  data.frame(
    Sample = target_sample,
    Entropy_Metric = e_col,
    N = length(H),
    Pearson_r_vs_log2K = cor(H, log2_K),
    Slope_vs_log2K = unname(coef(fit)[2]),
    Intercept_vs_log2K = unname(coef(fit)[1]),
    R_squared = summary(fit)$r.squared,
    Mean_H_over_Mean_log2K = mean(H) / mean(log2_K),
    stringsAsFactors = FALSE
  )
})
degeneracy_df <- do.call(rbind, degeneracy_rows)

# Observed spread of the stored non-zero logpac values, on the same gene set.
# This is the compression that drives the degeneracy; recorded so the range
# quoted in D2 is measured rather than assumed.
nz_vals <- data_filtered@x
value_range_df <- data.frame(
  Sample = target_sample,
  Layer = c("data (logpac, log2 scale)", "linear (2^x - 1)"),
  Min = c(min(nz_vals), min(2^nz_vals - 1)),
  Q01 = c(quantile(nz_vals, 0.01), quantile(2^nz_vals - 1, 0.01)),
  Median = c(median(nz_vals), median(2^nz_vals - 1)),
  Q99 = c(quantile(nz_vals, 0.99), quantile(2^nz_vals - 1, 0.99)),
  Max = c(max(nz_vals), max(2^nz_vals - 1)),
  N_nonzero = length(nz_vals),
  stringsAsFactors = FALSE
)

# Spread of the two supports. SpaNorm's qnbinom step both imputes genes with
# zero raw counts and zeroes genes with non-zero counts, so n_detected_norm is
# far narrower than n_detected_counts. Recorded because it is what separates the
# log-scale degeneracy (driven by n_detected_norm) from the raw depth ceiling
# (driven by n_detected_counts) that Stage 3 rarefaction targets.
support_df <- data.frame(
  Sample = target_sample,
  Support = c("n_detected_counts", "n_detected_norm"),
  Min = c(min(spatial_obj$n_detected_counts), min(spatial_obj$n_detected_norm)),
  Median = c(median(spatial_obj$n_detected_counts), median(spatial_obj$n_detected_norm)),
  Max = c(max(spatial_obj$n_detected_counts), max(spatial_obj$n_detected_norm)),
  Mean_Abs_Diff_Between_Supports = mean(abs(spatial_obj$n_detected_norm - spatial_obj$n_detected_counts)),
  stringsAsFactors = FALSE
)

write.csv(degeneracy_df, file.path(out_dir, paste0(target_sample, "_entropy_scaling_degeneracy.csv")), row.names = FALSE)
write.csv(support_df, file.path(out_dir, paste0(target_sample, "_entropy_scaling_support.csv")), row.names = FALSE)
write.csv(value_range_df, file.path(out_dir, paste0(target_sample, "_entropy_scaling_value_range.csv")), row.names = FALSE)

cat("\n=================================================================\n")
cat("Diagnostic Summary Table:\n")
cat("=================================================================\n")
print(diag_res$summary_table[, c("Entropy_Metric", "Target_Variable", "Pearson_r", "Pearson_p", "Spearman_rho", "Spearman_p")])

cat("\n=================================================================\n")
cat("Degeneracy fit against log2(n_detected_norm):\n")
cat("=================================================================\n")
print(degeneracy_df)

cat("\n=================================================================\n")
cat("Stored non-zero value spread (MT/ribo excluded):\n")
cat("=================================================================\n")
print(value_range_df)

cat("\n=================================================================\n")
cat("Detected-gene support spread (MT/ribo excluded):\n")
cat("=================================================================\n")
print(support_df)

cat(sprintf("\nSaved diagnostic table to: %s\n", file.path(out_dir, paste0(target_sample, "_entropy_scaling_diagnostic.csv"))))
cat(sprintf("Saved diagnostic plot to: %s\n", file.path(out_dir, paste0(target_sample, "_entropy_scaling_diagnostic_plot.png"))))
cat(sprintf("Saved degeneracy fit to: %s\n", file.path(out_dir, paste0(target_sample, "_entropy_scaling_degeneracy.csv"))))
cat(sprintf("Saved value range to: %s\n", file.path(out_dir, paste0(target_sample, "_entropy_scaling_value_range.csv"))))
cat(sprintf("Saved support spread to: %s\n", file.path(out_dir, paste0(target_sample, "_entropy_scaling_support.csv"))))
cat("Diagnostic completed successfully.\n")
