library(ggplot2)
library(patchwork)

calculate_entropy_correlations <- function(seurat_obj,
                                            entropy_cols = "shannon_entropy",
                                            targets = NULL,
                                            count_col = NULL,
                                            feature_col = NULL,
                                            sample_name = "Sample",
                                            output_dir = file.path("results", "statistical_tests"),
                                            file_suffix = "_entropy_correlations",
                                            save_outputs = TRUE) {

  meta <- seurat_obj@meta.data

  if (is.null(targets)) {
    if (is.null(count_col)) count_col <- grep("^nCount", colnames(meta), value = TRUE)[1]
    if (is.null(feature_col)) feature_col <- grep("^nFeature", colnames(meta), value = TRUE)[1]
    mt_col <- grep("^(percent\\.mt|percent_mt|percent_mito)", colnames(meta), value = TRUE, ignore.case = TRUE)[1]
    ribo_col <- grep("^(percent\\.ribo|percent_ribo|percent\\.rp)", colnames(meta), value = TRUE, ignore.case = TRUE)[1]

    targets <- c(nCounts = count_col, nFeatures = feature_col)
    if (!is.na(mt_col) && !is.null(mt_col)) targets["Percent_MT"] <- mt_col
    if (!is.na(ribo_col) && !is.null(ribo_col)) targets["Percent_Ribo"] <- ribo_col
  }

  sample_prefix <- if (!is.null(sample_name) && nchar(sample_name) > 0) sample_name else "Sample"

  summary_rows <- list()
  plot_list <- list()

  for (e_col in entropy_cols) {
    e_label <- if (e_col == "shannon_entropy") {
      "Shannon Entropy (Normalized)"
    } else if (e_col == "shannon_entropy_raw") {
      "Shannon Entropy (Raw)"
    } else {
      e_col
    }

    for (t_name in names(targets)) {
      t_col <- targets[[t_name]]

      if (is.na(t_col) || is.null(t_col) || !t_col %in% colnames(meta)) {
        warning(sprintf("Target column '%s' not found in metadata. Skipping.", t_col))
        next
      }
      if (!e_col %in% colnames(meta)) {
        warning(sprintf("Entropy column '%s' not found in metadata. Skipping.", e_col))
        next
      }

      x <- meta[[t_col]]
      y <- meta[[e_col]]

      valid <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
      x_clean <- x[valid]
      y_clean <- y[valid]

      if (length(x_clean) < 3) {
        warning(sprintf("Not enough valid data points for '%s' vs '%s'. Skipping.", e_col, t_name))
        next
      }

      has_variance <- (sd(x_clean) > 0) && (sd(y_clean) > 0)

      if (has_variance) {
        p_test <- suppressWarnings(cor.test(x_clean, y_clean, method = "pearson"))
        s_test <- suppressWarnings(cor.test(x_clean, y_clean, method = "spearman", exact = FALSE))

        p_r <- if (!is.null(p_test$estimate) && !is.na(p_test$estimate)) unname(p_test$estimate) else NA_real_
        p_pval <- if (!is.null(p_test$p.value) && !is.na(p_test$p.value)) p_test$p.value else NA_real_
        s_rho <- if (!is.null(s_test$estimate) && !is.na(s_test$estimate)) unname(s_test$estimate) else NA_real_
        s_pval <- if (!is.null(s_test$p.value) && !is.na(s_test$p.value)) s_test$p.value else NA_real_
      } else {
        p_r <- NA_real_
        p_pval <- NA_real_
        s_rho <- NA_real_
        s_pval <- NA_real_
      }

      p_r_str <- if (is.na(p_r) || is.nan(p_r)) "NA" else sprintf("%.4f", round(p_r, 4))
      s_rho_str <- if (is.na(s_rho) || is.nan(s_rho)) "NA" else sprintf("%.4f", round(s_rho, 4))

      p_str <- if (is.na(p_pval) || is.nan(p_pval)) {
        "NA"
      } else if (p_pval < 0.0001) {
        sprintf("%.2e", p_pval)
      } else {
        sprintf("%.4f", p_pval)
      }

      s_str <- if (is.na(s_pval) || is.nan(s_pval)) {
        "NA"
      } else if (s_pval < 0.0001) {
        sprintf("%.2e", s_pval)
      } else {
        sprintf("%.4f", s_pval)
      }

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Sample = sample_prefix,
        Entropy_Metric = e_label,
        Target_Variable = t_name,
        Pearson_r = if (is.na(p_r) || is.nan(p_r)) NA_real_ else round(p_r, 4),
        Pearson_p = p_str,
        Spearman_rho = if (is.na(s_rho) || is.nan(s_rho)) NA_real_ else round(s_rho, 4),
        Spearman_p = s_str,
        stringsAsFactors = FALSE
      )

      anno_text <- if (has_variance) {
        sprintf("Pearson r: %s (p = %s)\nSpearman \u03c1: %s (p = %s)", p_r_str, p_str, s_rho_str, s_str)
      } else {
        sprintf("Pearson r: NA\nSpearman \u03c1: NA\n(Zero variance in %s)", if (sd(x_clean) == 0) t_name else e_label)
      }

      p_scatter <- ggplot(meta, aes(x = .data[[t_col]], y = .data[[e_col]])) +
        geom_point(alpha = 0.4, color = "#2c3e50", size = 1.2)

      if (has_variance) {
        p_scatter <- p_scatter +
          geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2)
      }

      p_scatter <- p_scatter +
        labs(
          title = paste(e_label, "vs", t_name),
          subtitle = paste("Sample:", sample_prefix),
          x = paste(t_name, paste0("(", t_col, ")")),
          y = e_label
        ) +
        annotate(
          "label", x = Inf, y = -Inf, label = anno_text,
          hjust = 1.05, vjust = -0.2, size = 3.2,
          fill = alpha("white", 0.85),
          label.padding = grid::unit(0.35, "lines")
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray30"),
          axis.title = element_text(face = "bold", size = 9.5),
          panel.grid.minor = element_blank()
        )

      plot_list[[paste(e_col, t_name, sep = "_")]] <- p_scatter
    }
  }

  summary_table <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()
  combined_plot <- if (length(plot_list) > 1) {
    wrap_plots(plot_list, ncol = 2)
  } else if (length(plot_list) == 1) {
    plot_list[[1]]
  } else {
    NULL
  }

  if (save_outputs && nrow(summary_table) > 0) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    write.csv(summary_table, file.path(output_dir, paste0(sample_prefix, file_suffix, ".csv")), row.names = FALSE)
    if (!is.null(combined_plot)) {
      n_plots <- length(plot_list)
      plot_ncol <- if (n_plots > 1) 2 else 1
      plot_nrow <- ceiling(n_plots / plot_ncol)
      plot_width <- if (plot_ncol == 2) 11 else 6
      plot_height <- if (plot_nrow == 1) 5 else if (plot_nrow == 2) 8 else 3.5 * plot_nrow

      ggsave(file.path(output_dir, paste0(sample_prefix, file_suffix, "_plot.png")),
             plot = combined_plot, width = plot_width, height = plot_height, dpi = 300)
    }
  }

  return(list(
    summary_table = summary_table,
    plots = plot_list,
    combined_plot = combined_plot
  ))
}
