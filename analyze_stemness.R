set.seed(23)

library(Seurat)
library(ggplot2)
library(patchwork)

source("R/entropy_correlation.R")

# Define reference stemness markers panel
stemness_genes <- list(c("PROM1", "SOX2", "POU5F1", "NANOG", "NES", "CD44", "MYC"))

samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
samples <- samples[samples != ""]

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
      detection_rows[[g]] <- data.frame(
        gene = g,
        pct_spots = pct_spots[[g]],
        mean_count = mean_counts[[g]],
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
      status = "missing/filtered",
      stringsAsFactors = FALSE
    )
  }
  detection_df <- do.call(rbind, detection_rows[panel])

  cat("Stemness marker detection summary:\n")
  print(detection_df)

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

  # Note if low-abundance markers (POU5F1, NANOG, PROM1) are poorly detected
  low_abund <- intersect(c("POU5F1", "NANOG", "PROM1"), panel)
  low_detected <- low_abund[low_abund %in% genes_in_data]
  if (length(low_detected) > 0) {
    low_pcts <- detection_df[low_detected, "pct_spots"]
    if (all(low_pcts < 2.0)) {
      cat("Note: Low-abundance pluripotency markers (POU5F1/NANOG/PROM1) detected in <2% of spots. Module score is predominantly driven by robust markers (CD44/NES/MYC).\n")
    }
  }

  # Calculate stemness module score (deterministic with set.seed(23))
  spatial_obj <- AddModuleScore(
    object = spatial_obj,
    features = list(genes_in_data),
    name = "Stemness_Score"
  )

  if ("shannon_entropy" %in% colnames(spatial_obj@meta.data)) {
    cor_res <- calculate_entropy_correlations(
      seurat_obj = spatial_obj,
      entropy_cols = "shannon_entropy",
      targets = c(Stemness = "Stemness_Score1"),
      sample_name = sample_name,
      output_dir = out_dir,
      file_suffix = "_stemness_correlations",
      save_outputs = TRUE
    )

    # Spatial plots side-by-side
    p_spatial <- suppressMessages(
      SpatialFeaturePlot(spatial_obj, features = c("shannon_entropy", "Stemness_Score1")) +
        plot_layout(guides = "collect")
    )

    ggsave(file.path(out_dir, paste0(sample_name, "_spatial_comparison.png")),
           plot = p_spatial, width = 12, height = 5, dpi = 300)
    cat("Spatial comparison plot saved to", out_dir, "\n")
  } else {
    warning("Column 'shannon_entropy' not found in metadata. Cannot compare.")
  }

  # Save updated object back to disk for downstream pipeline steps (find_entropy_markers.R)
  saveRDS(spatial_obj, file = seurat_file)
  cat(sprintf("Saved updated spatial object with Stemness_Score1 to %s\n", seurat_file))

  cat("Done stemness analysis for:", sample_name, "\n")
  cat("==========================================\n")
  invisible(spatial_obj)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  target_sample <- args[1]
  if (!target_sample %in% samples) {
    stop(sprintf("Sample '%s' not found in data/. Available samples: %s", target_sample, paste(samples, collapse = ", ")))
  }
} else {
  target_sample <- samples[1]
}

analyze_stemness_sample(target_sample)
