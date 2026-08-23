library(ggplot2)
library(patchwork)

calculate_entropy_correlations <- function(seurat_obj,
                                            entropy_cols = "shannon_entropy",
                                            count_col = NULL,
                                            feature_col = NULL,
                                            sample_name = "Sample",
                                            output_dir = file.path("results", "statistical_tests"),
                                            save_outputs = TRUE) {

  meta <- seurat_obj@meta.data

  if (is.null(count_col)) count_col <- grep("^nCount", colnames(meta), value = TRUE)[1]
  if (is.null(feature_col)) feature_col <- grep("^nFeature", colnames(meta), value = TRUE)[1]

  sample_prefix <- if (!is.null(sample_name) && nchar(sample_name) > 0) sample_name else "Sample"

  targets <- c(nCounts = count_col, nFeatures = feature_col)
  summary_rows <- list()
  plot_list <- list()

  for (e_col in entropy_cols) {
    e_label <- if (e_col == "shannon_entropy") "Shannon Entropy" else e_col

    for (t_name in names(targets)) {
      t_col <- targets[[t_name]]

      x <- meta[[t_col]]
      y <- meta[[e_col]]

      valid <- !is.na(x) & !is.na(y)
      x_clean <- x[valid]
      y_clean <- y[valid]

      p_test <- cor.test(x_clean, y_clean, method = "pearson")
      s_test <- suppressWarnings(cor.test(x_clean, y_clean, method = "spearman", exact = FALSE))

      p_r <- unname(p_test$estimate)
      p_pval <- p_test$p.value
      s_rho <- unname(s_test$estimate)
      s_pval <- s_test$p.value

      p_str <- if (p_pval < 0.0001) sprintf("%.2e", p_pval) else sprintf("%.4f", p_pval)
      s_str <- if (s_pval < 0.0001) sprintf("%.2e", s_pval) else sprintf("%.4f", s_pval)

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Sample = sample_prefix,
        Entropy_Metric = e_label,
        Target_Variable = t_name,
        Pearson_r = round(p_r, 4),
        Pearson_p = p_str,
        Spearman_rho = round(s_rho, 4),
        Spearman_p = s_str,
        stringsAsFactors = FALSE
      )

      anno_text <- sprintf("Pearson r: %.4f (p = %s)\nSpearman \u03c1: %.4f (p = %s)", p_r, p_str, s_rho, s_str)

      p_scatter <- ggplot(meta, aes(x = .data[[t_col]], y = .data[[e_col]])) +
        geom_point(alpha = 0.4, color = "#2c3e50", size = 1.2) +
        geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
        labs(
          title = paste(e_label, "vs", t_name),
          subtitle = paste("Sample:", sample_prefix),
          x = paste(t_name, paste0("(", t_col, ")")),
          y = e_label
        ) +
        annotate(
          "label", x = Inf, y = -Inf, label = anno_text,
          hjust = 1.05, vjust = -0.2, size = 3.5,
          fill = alpha("white", 0.85),
          label.padding = grid::unit(0.4, "lines")
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
          axis.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank()
        )

      plot_list[[paste(e_col, t_name, sep = "_")]] <- p_scatter
    }
  }

  summary_table <- do.call(rbind, summary_rows)
  combined_plot <- wrap_plots(plot_list, ncol = 2)

  if (save_outputs) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    write.csv(summary_table, file.path(output_dir, paste0(sample_prefix, "_entropy_correlations.csv")), row.names = FALSE)
    ggsave(file.path(output_dir, paste0(sample_prefix, "_entropy_correlations_plot.png")),
           plot = combined_plot, width = 10, height = 6, dpi = 300)
  }

  return(list(
    summary_table = summary_table,
    plots = plot_list,
    combined_plot = combined_plot
  ))
}
