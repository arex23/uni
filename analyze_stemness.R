library(Seurat)
library(ggplot2)
library(patchwork)

# Define stemness markers
stemness_genes <- list(c("PROM1", "SOX2", "POU5F1", "NANOG", "NES", "CD44", "MYC"))

# Get all seurat objects in results/seurat_objects
seurat_files <- list.files("results/seurat_objects", pattern = "\\.rds$", full.names = TRUE)

if (length(seurat_files) == 0) {
  stop("No Seurat object found in results/seurat_objects.")
}

out_dir <- "results/stemness_analysis"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (seurat_file in seurat_files) {
  sample_name <- gsub("_spatial_obj\\.rds$", "", basename(seurat_file))
  cat("==========================================\n")
  cat("Processing sample:", sample_name, "\n")
  
  spatial_obj <- readRDS(seurat_file)
  
  # Check which genes are present in the dataset
  all_genes <- rownames(spatial_obj)
  present_genes <- stemness_genes[[1]][stemness_genes[[1]] %in% all_genes]
  
  cat("Present stemness genes in dataset:", paste(present_genes, collapse=", "), "\n")
  missing_genes <- setdiff(stemness_genes[[1]], present_genes)
  if (length(missing_genes) > 0) {
    cat("Missing stemness genes:", paste(missing_genes, collapse=", "), "\n")
  }
  
  if (length(present_genes) == 0) {
    cat("No stemness markers found in this sample. Skipping...\n")
    next
  }
  
  # Calculate stemness module score
  spatial_obj <- AddModuleScore(
    object = spatial_obj,
    features = list(present_genes),
    name = "Stemness_Score"
  )
  
  if ("shannon_entropy" %in% colnames(spatial_obj@meta.data)) {
    entropy_vals <- spatial_obj$shannon_entropy
    stemness_vals <- spatial_obj$Stemness_Score1
    
    # Correlation between entropy and stemness score
    cor_res <- cor.test(entropy_vals, stemness_vals, method = "spearman")
    cat("Spearman correlation between entropy and stemness score:\n")
    cat("  Rho: ", cor_res$estimate, "\n")
    cat("  p-value: ", cor_res$p.value, "\n")
    
    # Save correlation result to text
    sink(file.path(out_dir, paste0(sample_name, "_correlation.txt")))
    print(cor_res)
    sink()
    
    # Scatter plot
    p_scatter <- ggplot(spatial_obj@meta.data, aes(x = shannon_entropy, y = Stemness_Score1)) +
      geom_point(alpha = 0.5, color = "blue") +
      geom_smooth(method = "lm", color = "red") +
      theme_minimal() +
      labs(
        title = paste("Entropy vs Stemness Score -", sample_name),
        x = "Shannon Entropy",
        y = "Stemness Module Score"
      ) +
      annotate("text", x = min(entropy_vals, na.rm=TRUE), 
               y = max(stemness_vals, na.rm=TRUE), 
               label = paste("Rho =", round(cor_res$estimate, 3), "\np =", signif(cor_res$p.value, 3)),
               hjust = 0, vjust = 1, color = "darkred")
    
    ggsave(file.path(out_dir, paste0(sample_name, "_scatter.png")), p_scatter, width = 6, height = 5)
    
    # Spatial plots side-by-side
    p_spatial <- SpatialFeaturePlot(spatial_obj, features = c("shannon_entropy", "Stemness_Score1")) +
      plot_layout(guides = "collect")
      
    ggsave(file.path(out_dir, paste0(sample_name, "_spatial_comparison.png")), p_spatial, width = 12, height = 5)
    cat("Plots saved to", out_dir, "\n")
  } else {
    cat("Column 'shannon_entropy' not found in metadata. Cannot compare.\n")
  }
}

cat("==========================================\n")
cat("Stemness analysis complete.\n")
