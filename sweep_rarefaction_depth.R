#!/usr/bin/env Rscript

# sweep_rarefaction_depth.R
# Stage 3 Diagnostic: Sweep candidate downsampling depths D in {1000, 2000, 3000, 5000, 8000}
# to evaluate spot retention, depth decoupling, and ranking stability.
#
# Reads raw counts from the sample h5 rather than the Seurat object written by
# analyze_entropy.R. That object has already been filtered at nCount >= D, so
# using it would make the retention column circular (every D <= the pipeline
# threshold would trivially retain ~100% of spots) and would make this script
# depend on the very output it is meant to inform. The spot set here is the
# pre-Stage-3 candidate set: on-tissue, nFeature >= 250, nCount >= 500.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

source("R/shannon_entropy.R")
source("R/gene_universe.R")

args <- commandArgs(trailingOnly = TRUE)
target_sample <- if (length(args) > 0) args[1] else "sample1"

sample_dir <- file.path("data", target_sample)
h5_files <- list.files(sample_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
tp_files <- list.files(sample_dir, pattern = "tissue_positions.*\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(h5_files) == 0) stop(sprintf("No .h5 file found for sample '%s'", target_sample))
if (length(tp_files) == 0) stop(sprintf("No tissue_positions.csv found for sample '%s'", target_sample))

cat("=================================================================\n")
cat(sprintf("Running Rarefaction Depth Sweep for: %s\n", target_sample))
cat("=================================================================\n")

# Candidate spot set: on-tissue, nFeature >= 250, nCount >= 500 (the pre-Stage-3
# floor), so that every candidate D can be evaluated on the same starting point.
counts_all <- Read10X_h5(h5_files[1])
rownames(counts_all) <- sanitize_feature_names(rownames(counts_all))
cat(sprintf("Loaded %d features x %d grid barcodes.\n", nrow(counts_all), ncol(counts_all)))

tp <- read.csv(tp_files[1])
ontissue_bcs <- if ("in_tissue" %in% colnames(tp)) tp$barcode[tp$in_tissue == 1] else tp[[1]][tp[[2]] == 1]
counts_all <- counts_all[, intersect(colnames(counts_all), ontissue_bcs), drop = FALSE]
n_spots_ontissue <- ncol(counts_all)

nCount_all <- Matrix::colSums(counts_all)
nFeature_all <- Matrix::colSums(counts_all > 0)
candidate_spots <- colnames(counts_all)[nCount_all >= 500 & nFeature_all >= 250]
counts_all <- counts_all[, candidate_spots, drop = FALSE]
total_spots <- ncol(counts_all)
cat(sprintf("On-tissue spots: %d; candidate spots (nCount >= 500, nFeature >= 250): %d\n",
            n_spots_ontissue, total_spots))

# Feature pipeline: cohort gene universe as-is (D1/D5), followed by MT/ribosomal exclusion.
# Rarefaction standardizes depth across spots, eliminating the need for an ad-hoc per-sample
# spot detection threshold.
counts_all <- filter_by_gene_universe(counts_all)
counts_entropy <- exclude_gene_families(counts_all, "^(MT-|RP[SL])",
                                        context = "the depth sweep")
cat(sprintf("Entropy matrix: %d features x %d spots.\n", nrow(counts_entropy), ncol(counts_entropy)))

entropy_depth <- Matrix::colSums(counts_entropy)
n_det_counts <- Matrix::colSums(counts_entropy > 0)
nCount_all <- nCount_all[candidate_spots]
nFeature_all <- nFeature_all[candidate_spots]

candidate_depths <- c(1000, 2000, 3000, 5000, 8000)

sweep_rows <- list()
entropy_vectors <- list()

for (D in candidate_depths) {
  # Retention is measured on the entropy matrix's own depth, which is what
  # rarefaction standardises -- not on nCount over all features.
  valid_spots <- names(which(entropy_depth >= D))
  n_valid <- length(valid_spots)
  pct_retained <- round(100 * n_valid / total_spots, 2)

  if (n_valid < 10) {
    warning(sprintf("Depth D = %d retains only %d spots. Skipping.", D, n_valid))
    next
  }

  h_vals <- rarefied_entropy_matrix(
    counts_entropy[, valid_spots, drop = FALSE],
    depth = D,
    n_draws = 5,
    seed = 23,
    allow_shallow = FALSE,
    verbose = FALSE
  )
  entropy_vectors[[as.character(D)]] <- h_vals

  nCount_vals <- nCount_all[valid_spots]
  nFeature_vals <- nFeature_all[valid_spots]
  ndet_vals <- n_det_counts[valid_spots]

  r_nCount <- cor(h_vals, nCount_vals)
  rho_nCount <- cor(h_vals, nCount_vals, method = "spearman")
  r_nFeat <- cor(h_vals, nFeature_vals)
  rho_nFeat <- cor(h_vals, nFeature_vals, method = "spearman")
  r_ndet <- cor(h_vals, ndet_vals)
  rho_ndet <- cor(h_vals, ndet_vals, method = "spearman")

  sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
    Sample = target_sample,
    Depth_D = D,
    Spots_On_Tissue = n_spots_ontissue,
    Candidate_Spots = total_spots,
    Retained_Spots = n_valid,
    Pct_Retained = pct_retained,
    Pearson_r_nCount = r_nCount,
    Spearman_rho_nCount = rho_nCount,
    Pearson_r_nFeature = r_nFeat,
    Spearman_rho_nFeature = rho_nFeat,
    Pearson_r_ndetected = r_ndet,
    Spearman_rho_ndetected = rho_ndet,
    Mean_Rarefied_Entropy = mean(h_vals),
    SD_Rarefied_Entropy = sd(h_vals),
    Median_Rarefied_Entropy = median(h_vals),
    IQR_Rarefied_Entropy = IQR(h_vals),
    Min_Rarefied_Entropy = min(h_vals),
    Max_Rarefied_Entropy = max(h_vals),
    Range_Rarefied_Entropy = max(h_vals) - min(h_vals),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  D = %5d | Retained: %5d (%5.2f%%) | cor(H, nCount): r=%+.4f, rho=%+.4f | cor(H, nFeat): r=%+.4f, rho=%+.4f\n",
              D, n_valid, pct_retained, r_nCount, rho_nCount, r_nFeat, rho_nFeat))
}

sweep_df <- do.call(rbind, sweep_rows)

# Pairwise rank correlation on spots shared by every tested depth
tested_depths <- as.numeric(names(entropy_vectors))
shared_spots <- Reduce(intersect, lapply(entropy_vectors, names))

pairwise_mat <- matrix(NA, nrow = length(tested_depths), ncol = length(tested_depths),
                       dimnames = list(paste0("D_", tested_depths), paste0("D_", tested_depths)))

for (i in seq_along(tested_depths)) {
  v1 <- entropy_vectors[[as.character(tested_depths[i])]][shared_spots]
  for (j in seq_along(tested_depths)) {
    v2 <- entropy_vectors[[as.character(tested_depths[j])]][shared_spots]
    pairwise_mat[i, j] <- cor(v1, v2, method = "spearman")
  }
}
pairwise_df <- as.data.frame(pairwise_mat)
pairwise_df <- cbind(Depth = rownames(pairwise_df), pairwise_df)

out_dir <- file.path("results", "statistical_tests")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

sweep_file <- file.path(out_dir, paste0(target_sample, "_rarefaction_sweep.csv"))
pairwise_file <- file.path(out_dir, paste0(target_sample, "_rarefaction_sweep_pairwise.csv"))
plot_file <- file.path(out_dir, paste0(target_sample, "_rarefaction_sweep_plot.png"))

write.csv(sweep_df, sweep_file, row.names = FALSE)
write.csv(pairwise_df, pairwise_file, row.names = FALSE)

# Generate diagnostic summary plot
p1 <- ggplot(sweep_df, aes(x = factor(Depth_D), y = Pct_Retained)) +
  geom_bar(stat = "identity", fill = "#34495e", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(N=%d)", Pct_Retained, Retained_Spots)), vjust = -0.3, size = 3) +
  labs(title = "Spot Retention vs Rarefaction Depth", x = "Downsampling Depth D", y = "% Candidate Spots Retained") +
  ylim(0, 115) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

plot_corr_df <- data.frame(
  Depth_D = rep(factor(sweep_df$Depth_D), 2),
  Target = rep(c("nCount (Library Size)", "nFeature (Detected Genes)"), each = nrow(sweep_df)),
  Spearman_rho = c(sweep_df$Spearman_rho_nCount, sweep_df$Spearman_rho_nFeature)
)

p2 <- ggplot(plot_corr_df, aes(x = Depth_D, y = Spearman_rho, group = Target, color = Target)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Depth Correlation Decoupling across D", x = "Downsampling Depth D", y = "Spearman Correlation (ρ)") +
  scale_color_manual(values = c("#e74c3c", "#2980b9")) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position = "bottom"
  )

p_combined <- p1 + p2 + plot_layout(ncol = 2) +
  plot_annotation(
    title = paste("Rarefaction Depth Sweep Diagnostic —", target_sample),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(plot_file, plot = p_combined, width = 11, height = 5, dpi = 300)

cat("\n=================================================================\n")
cat("Pairwise Spearman Correlation Matrix on Shared Spots (N =", length(shared_spots), "):\n")
cat("=================================================================\n")
print(round(pairwise_mat, 4))
cat("\nSaved sweep results to:\n")
cat("  -", sweep_file, "\n")
cat("  -", pairwise_file, "\n")
cat("  -", plot_file, "\n")
cat("Depth sweep completed successfully.\n")
