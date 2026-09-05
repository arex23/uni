library(ggplot2)
library(patchwork)

#' Pearson and Spearman statistics for one pair of vectors
#'
#' Single source of truth for the correlation numbers reported across the
#' project: effect size, exact double-precision p-value, log10 p-value computed
#' in log space (so it does not saturate at 0 for the ~11.5k-spot samples), and
#' a 95% confidence interval. Spearman CIs use the Fisher-z transform with the
#' Bonett-Wright standard error sqrt(1.06 / (N - 3)) (D3).
#'
#' @param x,y Numeric vectors of equal length; non-finite pairs are dropped
#' @return List with N, sd_x, sd_y, has_variance and the Pearson / Spearman
#'   estimate, p-value, log10 p-value and CI bounds. All statistics are NA when
#'   either vector is constant or fewer than three complete pairs remain.
pair_correlation_stats <- function(x, y) {
  valid <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  x_clean <- x[valid]
  y_clean <- y[valid]
  N <- length(x_clean)

  out <- list(
    N = N,
    sd_x = if (N > 1) sd(x_clean) else NA_real_,
    sd_y = if (N > 1) sd(y_clean) else NA_real_,
    has_variance = FALSE,
    pearson_r = NA_real_,
    pearson_p = NA_real_,
    pearson_log10_p = NA_real_,
    pearson_ci_low = NA_real_,
    pearson_ci_high = NA_real_,
    spearman_rho = NA_real_,
    spearman_p = NA_real_,
    spearman_log10_p = NA_real_,
    spearman_ci_low = NA_real_,
    spearman_ci_high = NA_real_
  )

  if (N < 3) return(out)

  out$has_variance <- (out$sd_x > 0) && (out$sd_y > 0)
  if (!out$has_variance) return(out)

  p_test <- suppressWarnings(cor.test(x_clean, y_clean, method = "pearson"))
  s_test <- suppressWarnings(cor.test(x_clean, y_clean, method = "spearman", exact = FALSE))

  out$pearson_r <- if (!is.null(p_test$estimate) && !is.na(p_test$estimate)) unname(p_test$estimate) else NA_real_
  out$pearson_p <- if (!is.null(p_test$p.value) && !is.na(p_test$p.value)) unname(p_test$p.value) else NA_real_

  # Pearson log10 p-value via log-t distribution to avoid zero-saturation
  t_stat <- if (!is.null(p_test$statistic) && !is.na(p_test$statistic)) unname(p_test$statistic) else NA_real_
  df_val <- if (!is.null(p_test$parameter) && !is.na(p_test$parameter)) unname(p_test$parameter) else NA_real_
  out$pearson_log10_p <- if (!is.na(t_stat) && !is.na(df_val) && df_val > 0) {
    (log(2) + pt(-abs(t_stat), df = df_val, log.p = TRUE)) / log(10)
  } else {
    NA_real_
  }

  # Pearson 95% Confidence Interval
  if (!is.null(p_test$conf.int)) {
    out$pearson_ci_low <- unname(p_test$conf.int[1])
    out$pearson_ci_high <- unname(p_test$conf.int[2])
  }

  s_rho <- if (!is.null(s_test$estimate) && !is.na(s_test$estimate)) unname(s_test$estimate) else NA_real_
  out$spearman_rho <- s_rho
  out$spearman_p <- if (!is.null(s_test$p.value) && !is.na(s_test$p.value)) unname(s_test$p.value) else NA_real_

  # Spearman log10 p-value via the same t approximation cor.test uses (exact = FALSE)
  out$spearman_log10_p <- if (!is.na(s_rho) && N > 2 && abs(s_rho) < 1) {
    t_s <- s_rho * sqrt((N - 2) / (1 - s_rho^2))
    (log(2) + pt(-abs(t_s), df = N - 2, log.p = TRUE)) / log(10)
  } else {
    NA_real_
  }

  # Spearman 95% Confidence Interval via Fisher-z transformation.
  # SE uses the Bonett-Wright correction sqrt(1.06 / (N - 3)); the plain
  # Pearson SE 1 / sqrt(N - 3) understates the width for rank correlations.
  if (!is.na(s_rho) && N > 3 && abs(s_rho) < 1) {
    z_val <- atanh(s_rho)
    se_z <- sqrt(1.06 / (N - 3))
    crit_z <- qnorm(0.975)
    out$spearman_ci_low <- tanh(z_val - crit_z * se_z)
    out$spearman_ci_high <- tanh(z_val + crit_z * se_z)
  } else if (!is.na(s_rho) && abs(s_rho) >= 1) {
    out$spearman_ci_low <- s_rho
    out$spearman_ci_high <- s_rho
  }

  out
}


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

    # percent.ribo is deliberately not a default target: the probe panel carries
    # no cytoplasmic ribosomal probes, so the covariate is identically zero in
    # every sample (D2, results/cohort_qc/chemistry_check.csv) and correlating
    # against it only emits all-NA rows and empty plot panels. Pass it via
    # `targets` explicitly if a future chemistry ever populates those rows.
    targets <- c(nCounts = count_col, nFeatures = feature_col)
    if (!is.na(mt_col) && !is.null(mt_col)) targets["Percent_MT"] <- mt_col
  }

  sample_prefix <- if (!is.null(sample_name) && nchar(sample_name) > 0) sample_name else "Sample"

  summary_rows <- list()
  plot_list <- list()

  for (e_col in entropy_cols) {
    e_label <- if (e_col == "entropy_raw_plugin") {
      "Shannon Entropy (Raw Plug-in)"
    } else if (e_col == "entropy_spanorm_plugin") {
      "Shannon Entropy (SpaNorm Plug-in)"
    } else if (e_col == "shannon_entropy") {
      "Shannon Entropy (Normalized)"
    } else if (e_col == "shannon_entropy_raw") {
      "Shannon Entropy (Raw)"
    } else if (e_col == "shannon_entropy_log") {
      "Shannon Entropy (SpaNorm Log-Scale)"
    } else if (e_col == "shannon_entropy_linear") {
      "Shannon Entropy (SpaNorm Linear)"
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

      stats <- pair_correlation_stats(meta[[t_col]], meta[[e_col]])
      N <- stats$N

      if (N < 3) {
        warning(sprintf("Not enough valid data points for '%s' vs '%s'. Skipping.", e_col, t_name))
        next
      }

      has_variance <- stats$has_variance
      p_r <- stats$pearson_r
      p_pval <- stats$pearson_p
      p_log10_pval <- stats$pearson_log10_p
      p_ci_low <- stats$pearson_ci_low
      p_ci_high <- stats$pearson_ci_high
      s_rho <- stats$spearman_rho
      s_pval <- stats$spearman_p
      s_log10_pval <- stats$spearman_log10_p
      s_ci_low <- stats$spearman_ci_low
      s_ci_high <- stats$spearman_ci_high

      # Format display strings for plot annotations ONLY
      p_r_str <- if (is.na(p_r) || is.nan(p_r)) "NA" else sprintf("%.4f", p_r)
      s_rho_str <- if (is.na(s_rho) || is.nan(s_rho)) "NA" else sprintf("%.4f", s_rho)

      p_str <- if (is.na(p_pval) || is.nan(p_pval)) {
        "NA"
      } else if (!is.na(p_log10_pval) && p_log10_pval < -300) {
        sprintf("10^%.1f", p_log10_pval)
      } else if (p_pval < 0.0001) {
        sprintf("%.2e", p_pval)
      } else {
        sprintf("%.4f", p_pval)
      }

      s_str <- if (is.na(s_pval) || is.nan(s_pval)) {
        "NA"
      } else if (!is.na(s_log10_pval) && s_log10_pval < -300) {
        sprintf("10^%.1f", s_log10_pval)
      } else if (s_pval < 0.0001) {
        sprintf("%.2e", s_pval)
      } else {
        sprintf("%.4f", s_pval)
      }

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Sample = sample_prefix,
        Entropy_Metric = e_label,
        Target_Variable = t_name,
        N = N,
        Pearson_r = p_r,
        Pearson_p = p_pval,
        Pearson_log10_p = p_log10_pval,
        Pearson_CI_low = p_ci_low,
        Pearson_CI_high = p_ci_high,
        Spearman_rho = s_rho,
        Spearman_p = s_pval,
        Spearman_log10_p = s_log10_pval,
        Spearman_CI_low = s_ci_low,
        Spearman_CI_high = s_ci_high,
        stringsAsFactors = FALSE
      )

      anno_text <- if (has_variance) {
        sprintf("N: %d\nPearson r: %s (p = %s)\nSpearman \u03c1: %s (p = %s)", N, p_r_str, p_str, s_rho_str, s_str)
      } else {
        sprintf("N: %d\nPearson r: NA\nSpearman \u03c1: NA\n(Zero variance in %s)", N, if (!is.na(stats$sd_x) && stats$sd_x == 0) t_name else e_label)
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
