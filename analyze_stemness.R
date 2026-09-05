set.seed(23)

library(Seurat)
library(ggplot2)
library(patchwork)

source("R/entropy_correlation.R")
source("R/quality_confound.R")
source("R/cohort.R")

# Define reference stemness markers panel
stemness_genes <- list(c("PROM1", "SOX2", "POU5F1", "NANOG", "NES", "CD44", "MYC"))

samples <- cohort_samples()

analyze_stemness_sample <- function(sample_name) {
  cat("==========================================\n")
  cat("Processing stemness for sample:", sample_name, "\n")

  seurat_file <- file.path("results", "seurat_objects", paste0(sample_name, "_spatial_obj.rds"))
  if (!file.exists(seurat_file)) {
    stop(sprintf("Seurat object not found: %s. Run analyze_entropy.R on %s first.", seurat_file, sample_name))
  }

  spatial_obj <- readRDS(seurat_file)

  # Check marker detection before module scoring
  counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  panel <- stemness_genes[[1]]

  genes_in_data <- panel[panel %in% rownames(counts)]
  genes_missing <- setdiff(panel, genes_in_data)

  detection_rows <- list()
  if (length(genes_in_data) > 0) {
    sub_counts <- counts[genes_in_data, , drop = FALSE]
    pct_spots <- round(100 * Matrix::rowSums(sub_counts > 0) / ncol(counts), 2)
    mean_counts <- round(Matrix::rowMeans(sub_counts), 3)
    for (g in genes_in_data) {
      pct_val <- pct_spots[[g]]
      # Labels match the detection tiers reported in D4. They describe how often
      # a probe is seen, not how trustworthy it is: POU5F1 and PROM1 land in the
      # moderate tier but are flagged for pseudogene / cross-hybridisation checks.
      det_level <- if (pct_val >= 50.0) {
        "robust (>=50%)"
      } else if (pct_val >= 5.0) {
        "moderate (5-50%)"
      } else if (pct_val >= 1.0) {
        "low (1-5%)"
      } else {
        "rare (<1%)"
      }

      detection_rows[[g]] <- data.frame(
        gene = g,
        pct_spots = pct_val,
        mean_count = mean_counts[[g]],
        detection_level = det_level,
        status = "detected",
        stringsAsFactors = FALSE
      )
    }
  }
  for (g in genes_missing) {
    detection_rows[[g]] <- data.frame(
      gene = g,
      pct_spots = 0,
      mean_count = 0,
      detection_level = "absent (0%)",
      status = "missing/filtered",
      stringsAsFactors = FALSE
    )
  }
  detection_df <- do.call(rbind, detection_rows[panel])

  out_dir <- file.path("results", "stemness_analysis")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Save detailed per-gene detection table
  write.csv(detection_df, file.path(out_dir, paste0(sample_name, "_stemness_marker_detection.csv")), row.names = FALSE)

  # Save sample-level stemness QC summary
  stemness_qc_df <- data.frame(
    Sample = sample_name,
    Target_Genes_Total = length(panel),
    Present_Genes_Count = length(genes_in_data),
    Present_Genes = paste(genes_in_data, collapse = ";"),
    Missing_Genes = paste(genes_missing, collapse = ";"),
    stringsAsFactors = FALSE
  )
  write.csv(stemness_qc_df, file.path(out_dir, paste0(sample_name, "_stemness_qc.csv")), row.names = FALSE)
  cat(sprintf("Logged stemness marker QC to %s\n", out_dir))

  if (length(genes_in_data) == 0) {
    warning("No stemness markers found in this sample. Skipping module score calculation...")
    return(invisible(NULL))
  }

  # Report individual marker detection rates
  cat("Stemness marker detection summary:\n")
  for (g in genes_in_data) {
    cat(sprintf("  - %s: %s of spots (mean count: %s, %s)\n",
                g,
                paste0(detection_df[g, "pct_spots"], "%"),
                detection_df[g, "mean_count"],
                detection_df[g, "detection_level"]))
  }

  # Calculate stemness module score (deterministic with set.seed(23))
  spatial_obj <- AddModuleScore(
    object = spatial_obj,
    features = list(genes_in_data),
    name = "Stemness_Score"
  )

  ent_cols <- intersect(c("entropy_raw_plugin", "entropy_spanorm_plugin"), colnames(spatial_obj@meta.data))
  if (length(ent_cols) > 0) {
    cor_res <- calculate_entropy_correlations(
      seurat_obj = spatial_obj,
      entropy_cols = ent_cols,
      targets = c(Stemness = "Stemness_Score1"),
      sample_name = sample_name,
      output_dir = out_dir,
      file_suffix = "_stemness_correlations",
      save_outputs = TRUE
    )

    # Spatial plots side-by-side (Entropy vs Stemness Score)
    primary_ent <- if ("entropy_raw_plugin" %in% ent_cols) "entropy_raw_plugin" else ent_cols[1]
    p_spatial <- suppressMessages(
      SpatialFeaturePlot(spatial_obj, features = c(primary_ent, "Stemness_Score1")) +
        plot_layout(guides = "collect")
    )

    ggsave(file.path(out_dir, paste0(sample_name, "_spatial_comparison.png")),
           plot = p_spatial, width = 12, height = 5, dpi = 300)
    cat("Spatial comparison plot saved to", out_dir, "\n")

    # Quality confound gate (D3). Spot quality is a live alternative explanation
    # for any entropy-stemness association: a degraded, high-MT spot has flatter,
    # more ambient-like non-MT signal, and low-depth spots are the high-MT ones
    # (cohort rho = -0.41 between depth and MT, Stage 0). The association has to be
    # shown to survive adjustment for percent.mt and log(nCount) before it is
    # treated as biology. The sample1 numbers recorded in D3 were measured on the
    # abandoned rarefied metric and have to be re-run on whatever metric replaces
    # it; the check itself is metric-agnostic and runs on every entropy column.
    run_quality_confound_check(
      seurat_obj = spatial_obj,
      entropy_cols = ent_cols,
      score_col = "Stemness_Score1",
      mt_col = "percent.mt",
      count_col = "nCount_Spatial",
      feature_col = "nFeature_Spatial",
      mt_threshold = 10,
      sample_name = sample_name,
      output_dir = out_dir,
      save_outputs = TRUE
    )
  } else {
    warning("No entropy columns found in metadata. Cannot compute stemness correlation.")
  }

  # Save updated object back to disk for downstream pipeline steps (find_entropy_markers.R)
  saveRDS(spatial_obj, file = seurat_file)
  cat(sprintf("Saved updated spatial object with Stemness_Score1 to %s\n", seurat_file))

  cat("Done stemness analysis for:", sample_name, "\n")
  cat("==========================================\n")
  invisible(spatial_obj)
}

args <- commandArgs(trailingOnly = TRUE)
target_sample <- resolve_target_sample(args)

analyze_stemness_sample(target_sample)
