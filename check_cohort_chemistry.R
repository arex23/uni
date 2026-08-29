#!/usr/bin/env Rscript

# check_cohort_chemistry.R
# Stage 0: Check cohort chemistry, feature identity, and quality metrics across all Visium samples.

library(hdf5r)
library(Matrix)

# Locate all sample directories sorted naturally/numerically
sample_dirs <- list.dirs("data", full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("sample[0-9]+$", basename(sample_dirs))]
sample_nums <- as.numeric(gsub(".*sample", "", basename(sample_dirs)))
sample_dirs <- sample_dirs[order(sample_nums)]

if (length(sample_dirs) == 0) {
  stop("No sample directories found under data/")
}

results <- list()
ref_features <- NULL

cat("=================================================================\n")
cat("Checking cohort chemistry for", length(sample_dirs), "samples (on-tissue metrics)...\n")
cat("=================================================================\n\n")

for (i in seq_along(sample_dirs)) {
  sd <- sample_dirs[i]
  sname <- basename(sd)
  
  # Find .h5 file and tissue_positions.csv
  h5_files <- list.files(sd, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  if (length(h5_files) == 0) {
    warning(sprintf("No .h5 file found in %s, skipping.", sd))
    next
  }
  h5_path <- h5_files[1]
  
  tp_files <- list.files(sd, pattern = "tissue_positions.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(tp_files) == 0) {
    warning(sprintf("No tissue_positions.csv found in %s, skipping.", sd))
    next
  }
  tp_path <- tp_files[1]
  
  cat("-----------------------------------------------------------------\n")
  cat("Sample:", sname, "\n")
  cat("H5 file:", h5_path, "\n")
  cat("Tissue positions:", tp_path, "\n")
  cat("H5 root attributes:\n")
  
  # Open H5 file and read attributes
  f <- H5File$new(h5_path, mode = "r")
  attr_names <- h5attr_names(f)
  for (a in attr_names) {
    cat("  ", a, ":", paste(h5attr(f, a), collapse = ", "), "\n")
  }
  
  chem_desc <- if ("chemistry_description" %in% attr_names) paste(h5attr(f, "chemistry_description"), collapse = ", ") else NA_character_
  soft_ver <- if ("software_version" %in% attr_names) paste(h5attr(f, "software_version"), collapse = ", ") else NA_character_
  
  # Read sparse matrix components
  shape <- f[["matrix/shape"]]$read()
  features <- f[["matrix/features/name"]]$read()
  barcodes <- f[["matrix/barcodes"]]$read()
  data_vec <- f[["matrix/data"]]$read()
  indices <- f[["matrix/indices"]]$read()
  indptr <- f[["matrix/indptr"]]$read()
  f$close_all()
  
  # Set reference features from first sample (sample1) and check identity
  if (is.null(ref_features)) {
    ref_features <- features
    feat_identical <- TRUE
  } else {
    feat_identical <- identical(features, ref_features)
  }
  
  mat <- sparseMatrix(
    i = indices + 1,
    p = indptr,
    x = as.numeric(data_vec),
    dims = shape,
    dimnames = list(features, barcodes),
    index1 = TRUE
  )
  
  raw_spots_grid <- ncol(mat)
  n_features <- nrow(mat)
  
  # Read on-tissue barcodes
  tp <- read.csv(tp_path)
  ontissue_barcodes <- tp$barcode[tp$in_tissue == 1]
  ontissue_in_mat <- intersect(barcodes, ontissue_barcodes)
  n_ontissue <- length(ontissue_in_mat)
  
  mat_ontissue <- mat[, ontissue_in_mat, drop = FALSE]
  
  # On-tissue feature and count metrics
  row_sums_ontissue <- Matrix::rowSums(mat_ontissue)
  n_genes_all_zero_ontissue <- sum(row_sums_ontissue == 0)
  
  total_counts_ontissue <- sum(mat_ontissue)
  
  # Mitochondrial metrics on-tissue
  mt_idx <- grep("^MT-", features)
  n_MT_nonzero_ontissue <- sum(row_sums_ontissue[mt_idx] > 0)
  pct_MT_ontissue_total <- round((sum(mat_ontissue[mt_idx, , drop = FALSE]) / total_counts_ontissue) * 100, 4)
  
  col_counts_ontissue <- Matrix::colSums(mat_ontissue)
  col_features_ontissue <- Matrix::colSums(mat_ontissue > 0)
  
  # Spot-level %MT (handle possible 0-count spots safely)
  spot_pct_mt <- ifelse(col_counts_ontissue > 0, (Matrix::colSums(mat_ontissue[mt_idx, , drop = FALSE]) / col_counts_ontissue) * 100, NA_real_)
  mean_spot_pct_MT_ontissue <- round(mean(spot_pct_mt, na.rm = TRUE), 4)
  
  # Ribosomal metrics on-tissue
  rpsl_idx <- grep("^RP[SL]", features)
  n_RPSL_rows <- length(rpsl_idx)
  n_RPSL_nonzero_ontissue <- sum(row_sums_ontissue[rpsl_idx] > 0)
  
  # Depth QC filter: nCount >= 500 & nFeature >= 250 on-tissue spots
  qc_pass <- (col_counts_ontissue >= 500) & (col_features_ontissue >= 250)
  n_ontissue_postqc <- sum(qc_pass)
  pct_ontissue_retained_postqc <- round((n_ontissue_postqc / n_ontissue) * 100, 2)
  
  counts_postqc <- col_counts_ontissue[qc_pass]
  mean_spot_pct_MT_postqc <- round(mean(spot_pct_mt[qc_pass], na.rm = TRUE), 4)
  
  q_postqc <- quantile(counts_postqc, probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75))
  
  results[[sname]] <- data.frame(
    sample = sname,
    chemistry_description = chem_desc,
    software_version = soft_ver,
    features_identical_to_sample1 = feat_identical,
    n_features = n_features,
    raw_spots_grid = raw_spots_grid,
    n_ontissue = n_ontissue,
    n_ontissue_postqc = n_ontissue_postqc,
    pct_ontissue_retained_postqc = pct_ontissue_retained_postqc,
    n_genes_all_zero_ontissue = n_genes_all_zero_ontissue,
    n_MT_nonzero_ontissue = n_MT_nonzero_ontissue,
    pct_MT_ontissue_total = pct_MT_ontissue_total,
    mean_spot_pct_MT_ontissue = mean_spot_pct_MT_ontissue,
    mean_spot_pct_MT_postqc = mean_spot_pct_MT_postqc,
    n_RPSL_rows = n_RPSL_rows,
    n_RPSL_nonzero_ontissue = n_RPSL_nonzero_ontissue,
    total_counts_ontissue = total_counts_ontissue,
    nCount_min_postqc = unname(q_postqc[1]),
    nCount_q05_postqc = unname(q_postqc[2]),
    nCount_q10_postqc = unname(q_postqc[3]),
    nCount_q25_postqc = unname(q_postqc[4]),
    nCount_median_postqc = unname(q_postqc[5]),
    nCount_q75_postqc = unname(q_postqc[6]),
    stringsAsFactors = FALSE
  )
}

df_summary <- do.call(rbind, results)

out_dir <- file.path("results", "cohort_qc")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_csv <- file.path(out_dir, "chemistry_check.csv")
write.csv(df_summary, out_csv, row.names = FALSE)

cat("\n=================================================================\n")
cat("Cohort Chemistry Summary (On-Tissue & Post-QC Quantiles):\n")
cat("=================================================================\n")
print(df_summary, row.names = FALSE)
cat("\nSaved updated summary to:", out_csv, "\n")
