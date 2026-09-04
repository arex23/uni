#!/usr/bin/env Rscript

# build_gene_universe.R
# Stage 0 / Stage 3: Construct the deterministic cohort gene universe across all 22 samples.
#
# Identifies all features with >= 1 count in >= 1 sample on-tissue, filtering out
# structurally-zero probe-set artifacts (~18,067 unrepresented GRCh38 rows).

suppressPackageStartupMessages({
  library(hdf5r)
})

source("R/cohort.R")

cat("=================================================================\n")
cat("Building Cohort Gene Universe across all samples\n")
cat("=================================================================\n")

# Deliberately all 22 samples, not cohort_samples(). The universe is a permissive
# union over the full slide set (D5): the probe panel is identical across samples
# by construction, so including the 6 D6-excluded samples only adds features that
# are rare everywhere and contribute negligibly to entropy. Keeping the union at
# 22 also means the D6 exclusion does not silently redefine the feature space the
# sample1 numbers in D1/D2 were recorded on.
samples <- all_samples()

cat(sprintf("Found %d samples in data/\n", length(samples)))

gene_names <- NULL
gene_ids <- NULL
n_genes <- 0
gene_total_counts <- NULL
sample_detection_counts <- NULL
sample_ontissue_spots <- numeric(length(samples))
names(sample_ontissue_spots) <- samples

for (s in samples) {
  sample_dir <- file.path("data", s)
  h5_files <- list.files(path = sample_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  if (length(h5_files) == 0) {
    stop(sprintf("No .h5 file found for sample '%s'", s))
  }
  tp_files <- list.files(path = sample_dir, pattern = "tissue_positions.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(tp_files) == 0) {
    stop(sprintf("No tissue_positions.csv found for sample '%s'", s))
  }

  tp <- read.csv(tp_files[1])
  ontissue_bcs <- if ("in_tissue" %in% colnames(tp)) {
    tp$barcode[tp$in_tissue == 1]
  } else {
    tp[[1]][tp[[2]] == 1]
  }
  sample_ontissue_spots[[s]] <- length(ontissue_bcs)

  h5 <- H5File$new(h5_files[1], mode = "r")
  barcodes <- h5[["matrix/barcodes"]][]
  features <- h5[["matrix/features/name"]][]
  features_id <- h5[["matrix/features/id"]][]

  if (is.null(gene_names)) {
    gene_names <- features
    gene_ids <- features_id
    n_genes <- length(features)
    gene_total_counts <- numeric(n_genes)
    sample_detection_counts <- numeric(n_genes)
  } else {
    if (!identical(gene_names, features)) {
      stop(sprintf("Feature list mismatch in sample '%s'", s))
    }
  }

  data_vals <- h5[["matrix/data"]][]
  indices <- h5[["matrix/indices"]][] + 1 # 1-indexed gene index
  indptr <- h5[["matrix/indptr"]][]       # 0-indexed column pointers

  on_tissue_col_idx <- which(barcodes %in% ontissue_bcs)

  sample_gene_counts <- numeric(n_genes)
  for (c_idx in on_tissue_col_idx) {
    start_i <- indptr[c_idx] + 1
    end_i <- indptr[c_idx + 1]
    if (end_i >= start_i) {
      g_inds <- indices[start_i:end_i]
      vals <- data_vals[start_i:end_i]
      sample_gene_counts[g_inds] <- sample_gene_counts[g_inds] + vals
    }
  }
  h5$close_all()

  gene_total_counts <- gene_total_counts + sample_gene_counts
  sample_detection_counts <- sample_detection_counts + (sample_gene_counts > 0)
  cat(sprintf("  Processed %s: %d on-tissue spots, %d genes detected\n",
              s, length(on_tissue_col_idx), sum(sample_gene_counts > 0)))
}

# Deprecated probe rows (`DEPRECATED_ENSG*`) are retired probes that 10x still
# emits in the feature list. They carry real counts, so the ">= 1 count" rule
# alone would admit them, but they are not interpretable gene-level measurements
# and are excluded from the universe explicitly. This must be explicit: Seurat
# rewrites underscores to dashes on load, so a universe keyed on the raw h5 names
# would drop these rows by accident rather than by decision (D5).
is_deprecated_probe <- grepl("^DEPRECATED_", gene_names)
is_detected <- sample_detection_counts > 0
is_universe <- is_detected & !is_deprecated_probe

universe_df <- data.frame(
  gene_id = gene_ids,
  gene_name = gene_names,
  samples_detected = sample_detection_counts,
  total_counts = gene_total_counts,
  is_deprecated_probe = is_deprecated_probe,
  is_universe = is_universe,
  stringsAsFactors = FALSE
)

out_dir <- file.path("results", "cohort_qc")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, "gene_universe.csv")
write.csv(universe_df, out_file, row.names = FALSE)

n_total <- nrow(universe_df)
n_univ <- sum(is_universe)
n_zero <- sum(!is_detected)
n_dep <- sum(is_deprecated_probe & is_detected)
n_all22 <- sum(is_universe & sample_detection_counts == length(samples))

cat("\n=================================================================\n")
cat("Cohort Gene Universe Summary:\n")
cat("=================================================================\n")
cat(sprintf("Total features in probe/reference grid: %d\n", n_total))
cat(sprintf("Features in gene universe (>=1 count in >=1 sample, non-deprecated): %d (%.2f%%)\n", n_univ, 100 * n_univ / n_total))
cat(sprintf("Features structurally zero across all %d samples: %d (%.2f%%)\n", length(samples), n_zero, 100 * n_zero / n_total))
cat(sprintf("Deprecated probe rows with counts, excluded from universe: %d (%.2f%% of cohort counts)\n",
            n_dep, 100 * sum(gene_total_counts[is_deprecated_probe]) / sum(gene_total_counts)))
cat(sprintf("Universe features detected in all %d samples: %d (%.2f%% of universe)\n", length(samples), n_all22, 100 * n_all22 / n_univ))
cat(sprintf("\nSaved cohort gene universe to: %s\n", out_file))
cat("=================================================================\n")
