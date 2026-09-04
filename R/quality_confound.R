# R/quality_confound.R
# Gate: is the entropy-stemness association driven by spot quality?
#
# The mechanism to rule out: a degraded, high-MT spot has flatter, more
# ambient-like non-MT signal, which reads as higher transcriptional diversity, and
# shallow spots are the high-MT ones (cohort rho = -0.41 between depth and MT,
# Stage 0). If `Stemness_Score1` also tracks degradation, then any
# entropy-stemness correlation is a quality artefact rather than biology.
#
# The check is metric-agnostic: it runs on whichever entropy columns it is handed
# and makes no assumption about how they were estimated.
#
# Provides three checks (D3):
# 1. covariate_correlation_table() - zero-order correlations of the stemness score
#    and each entropy metric against percent.mt, nCount and nFeature. cor(score,
#    percent.mt) and cor(score, nCount) were not in any output file before this.
# 2. partial_correlation_table() - entropy vs stemness controlling for percent.mt
#    and log(nCount), by residualisation on the covariates (base `lm.fit`, no new
#    dependency; ppcor::pcor.test computes the same quantity).
# 3. subset_correlation_table() - the same zero-order correlations restricted to
#    clean spots (percent.mt < threshold) and, as the mirror image, to the high-MT
#    spots. If the association survives in the clean subset it is not degradation.
#
# Requires pair_correlation_stats() from R/entropy_correlation.R.

#' Label a metadata column for reporting
confound_var_label <- function(col) {
  labels <- c(
    entropy_raw_plugin = "Shannon Entropy (Raw Plug-in)",
    entropy_spanorm_plugin = "Shannon Entropy (SpaNorm Plug-in)",
    Stemness_Score1 = "Stemness Score",
    percent.mt = "Percent_MT",
    percent.ribo = "Percent_Ribo",
    nCount_Spatial = "nCounts",
    nFeature_Spatial = "nFeatures",
    log_nCount_Spatial = "log(nCounts)"
  )
  if (col %in% names(labels)) unname(labels[[col]]) else col
}

#' Zero-order correlations for a list of variable pairs
#'
#' @param data Data frame holding every column named in `pairs`
#' @param pairs List of length-2 character vectors, c(x_col, y_col)
#' @param sample_name Sample identifier written into the table
#' @param subset_label Value for the `Subset` column
#' @return Tidy data frame, one row per pair, with the same statistics
#'   (estimate, exact p, log10 p, 95% CI) reported everywhere else in the project
correlation_pair_table <- function(data, pairs, sample_name = "Sample", subset_label = "all_spots") {
  rows <- list()
  for (pr in pairs) {
    x_col <- pr[[1]]
    y_col <- pr[[2]]
    if (!all(c(x_col, y_col) %in% colnames(data))) {
      warning(sprintf("Skipping pair '%s' vs '%s': column not present.", x_col, y_col))
      next
    }
    st <- pair_correlation_stats(data[[x_col]], data[[y_col]])
    rows[[length(rows) + 1]] <- data.frame(
      Sample = sample_name,
      Subset = subset_label,
      X_Variable = confound_var_label(x_col),
      Y_Variable = confound_var_label(y_col),
      X_Column = x_col,
      Y_Column = y_col,
      N = st$N,
      Pearson_r = st$pearson_r,
      Pearson_p = st$pearson_p,
      Pearson_log10_p = st$pearson_log10_p,
      Pearson_CI_low = st$pearson_ci_low,
      Pearson_CI_high = st$pearson_ci_high,
      Spearman_rho = st$spearman_rho,
      Spearman_p = st$spearman_p,
      Spearman_log10_p = st$spearman_log10_p,
      Spearman_CI_low = st$spearman_ci_low,
      Spearman_CI_high = st$spearman_ci_high,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

#' Partial correlation of two variables given covariates
#'
#' Both variables are regressed on the covariates (plus intercept) and the
#' residuals correlated; for `method = "spearman"` every variable is rank
#' transformed first, which is the standard Spearman partial correlation.
#' Significance uses the partial-correlation t statistic
#' \eqn{t = r_{xy \cdot z} \sqrt{(N - 2 - k) / (1 - r^2)}} on \eqn{N - 2 - k}
#' degrees of freedom, and the CI the Fisher-z transform with
#' \eqn{SE = 1 / \sqrt{N - 3 - k}} (Bonett-Wright \eqn{\sqrt{1.06/(N-3-k)}} for
#' the rank version, matching the zero-order CIs in D3). \eqn{k} is the rank of
#' the covariate design, so collinear covariates do not inflate the df.
#'
#' @return List with N, k, zero-order estimate, partial estimate, p, log10 p, CI
#'   and the two residual vectors (for plotting)
partial_correlation_stats <- function(data, x_col, y_col, covariate_cols, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  vars <- c(x_col, y_col, covariate_cols)

  d <- data[, vars, drop = FALSE]
  keep <- Reduce(`&`, lapply(d, function(v) !is.na(v) & is.finite(v)))
  d <- d[keep, , drop = FALSE]
  N <- nrow(d)

  out <- list(
    N = N, K = length(covariate_cols), Estimate_zero_order = NA_real_,
    Estimate = NA_real_, P = NA_real_, Log10_p = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_, DF = NA_real_,
    residual_x = NULL, residual_y = NULL
  )

  if (N < length(covariate_cols) + 5) {
    warning(sprintf("Too few complete observations (%d) for partial correlation of '%s' vs '%s'.", N, x_col, y_col))
    return(out)
  }

  if (method == "spearman") {
    d[] <- lapply(d, rank)
  }

  x <- d[[x_col]]
  y <- d[[y_col]]
  if (sd(x) == 0 || sd(y) == 0) return(out)

  out$Estimate_zero_order <- cor(x, y)

  design <- cbind(`(Intercept)` = rep(1, N), as.matrix(d[, covariate_cols, drop = FALSE]))
  fit_x <- stats::lm.fit(design, x)
  fit_y <- stats::lm.fit(design, y)
  res_x <- fit_x$residuals
  res_y <- fit_y$residuals

  # Rank of the design minus the intercept: the number of covariate dimensions
  # actually removed, which is what the residual degrees of freedom depend on.
  k_used <- fit_x$rank - 1L
  df <- N - 2 - k_used
  if (df < 1 || sd(res_x) == 0 || sd(res_y) == 0) return(out)

  r <- cor(res_x, res_y)
  out$Estimate <- r
  out$K <- k_used
  out$DF <- df
  out$residual_x <- res_x
  out$residual_y <- res_y

  if (abs(r) < 1) {
    t_stat <- r * sqrt(df / (1 - r^2))
    out$P <- 2 * pt(-abs(t_stat), df = df)
    out$Log10_p <- (log(2) + pt(-abs(t_stat), df = df, log.p = TRUE)) / log(10)
  } else {
    out$P <- 0
    out$Log10_p <- -Inf
  }

  n_eff <- N - 3 - k_used
  if (n_eff > 0 && abs(r) < 1) {
    se_z <- if (method == "spearman") sqrt(1.06 / n_eff) else 1 / sqrt(n_eff)
    crit_z <- qnorm(0.975)
    out$CI_low <- tanh(atanh(r) - crit_z * se_z)
    out$CI_high <- tanh(atanh(r) + crit_z * se_z)
  }

  out
}

#' Partial correlation table over several covariate sets and both methods
#'
#' @param covariate_sets Named list of character vectors of covariate columns
partial_correlation_table <- function(data, x_col, y_col, covariate_sets, sample_name = "Sample") {
  rows <- list()
  for (set_name in names(covariate_sets)) {
    cov_cols <- covariate_sets[[set_name]]
    missing_cov <- setdiff(cov_cols, colnames(data))
    if (length(missing_cov) > 0) {
      warning(sprintf("Skipping covariate set '%s': missing %s.", set_name, paste(missing_cov, collapse = ", ")))
      next
    }
    for (m in c("pearson", "spearman")) {
      st <- partial_correlation_stats(data, x_col, y_col, cov_cols, method = m)
      rows[[length(rows) + 1]] <- data.frame(
        Sample = sample_name,
        X_Variable = confound_var_label(x_col),
        Y_Variable = confound_var_label(y_col),
        Covariate_Set = set_name,
        Covariates = paste(cov_cols, collapse = " + "),
        Method = m,
        N = st$N,
        K_Covariates = st$K,
        DF = st$DF,
        Zero_Order_Estimate = st$Estimate_zero_order,
        Partial_Estimate = st$Estimate,
        Partial_p = st$P,
        Partial_log10_p = st$Log10_p,
        Partial_CI_low = st$CI_low,
        Partial_CI_high = st$CI_high,
        # Attenuation is a ratio, so it is only reported when the zero-order
        # estimate is far enough from zero to make the ratio meaningful; the raw
        # plug-in metric sits at r = +0.008 against stemness and would otherwise
        # emit four-digit percentages that mean nothing. It is also undefined
        # when conditioning flips the sign -- the raw plug-in goes from r = -0.036
        # zero-order to +0.152 given log(nCount), which is a reversal, not an
        # attenuation, and the ratio renders it as "523.79% attenuated".
        Pct_Attenuation = if (!is.na(st$Estimate_zero_order) &&
                              !is.na(st$Estimate) &&
                              abs(st$Estimate_zero_order) >= 0.02 &&
                              sign(st$Estimate) == sign(st$Estimate_zero_order)) {
          round(100 * (1 - st$Estimate / st$Estimate_zero_order), 2)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

#' Run the full quality-confound check for one sample
#'
#' Writes three tidy CSVs and a four-panel diagnostic figure to `output_dir`:
#'   <sample>_quality_confound_covariates.csv - check 1
#'   <sample>_quality_confound_partial.csv    - check 2
#'   <sample>_quality_confound_subsets.csv    - check 3
#'   <sample>_quality_confound_plot.png
#'
#' @param seurat_obj Seurat object carrying the entropy columns, the stemness
#'   score and the QC covariates in `@meta.data`
#' @param mt_threshold percent.mt cut defining the "clean" subset
#' @return Named list of the three tables (invisibly NULL if prerequisites are missing)
run_quality_confound_check <- function(seurat_obj,
                                       entropy_cols = c("entropy_raw_plugin", "entropy_spanorm_plugin"),
                                       score_col = "Stemness_Score1",
                                       mt_col = "percent.mt",
                                       count_col = "nCount_Spatial",
                                       feature_col = "nFeature_Spatial",
                                       mt_threshold = 10,
                                       sample_name = "Sample",
                                       output_dir = file.path("results", "stemness_analysis"),
                                       save_outputs = TRUE) {

  meta <- seurat_obj@meta.data
  entropy_cols <- entropy_cols[entropy_cols %in% colnames(meta)]

  if (length(entropy_cols) == 0 || !score_col %in% colnames(meta)) {
    warning("Quality confound check skipped: no entropy column or no stemness score in metadata.")
    return(invisible(NULL))
  }
  if (!mt_col %in% colnames(meta)) {
    warning(sprintf("Quality confound check skipped: '%s' not in metadata.", mt_col))
    return(invisible(NULL))
  }

  # log(nCount) rather than nCount: depth is right-skewed and enters the support
  # ceiling logarithmically, so the linear covariate would under-adjust.
  log_count_col <- paste0("log_", count_col)
  meta[[log_count_col]] <- log(meta[[count_col]])

  primary_ent <- entropy_cols[1]

  # --- Check 1: zero-order covariate correlations -----------------------------
  qc_targets <- intersect(c(mt_col, count_col, feature_col), colnames(meta))
  cov_pairs <- list()
  for (v in c(score_col, entropy_cols)) {
    for (t in qc_targets) cov_pairs[[length(cov_pairs) + 1]] <- c(v, t)
  }
  for (e in entropy_cols) cov_pairs[[length(cov_pairs) + 1]] <- c(e, score_col)

  covariate_table <- correlation_pair_table(meta, cov_pairs, sample_name = sample_name, subset_label = "all_spots")

  # --- Check 2: partial correlations ------------------------------------------
  covariate_sets <- list(
    percent_mt = mt_col,
    log_nCount = log_count_col,
    percent_mt_and_log_nCount = c(mt_col, log_count_col)
  )
  partial_table <- do.call(rbind, lapply(entropy_cols, function(e) {
    partial_correlation_table(meta, e, score_col, covariate_sets, sample_name = sample_name)
  }))

  # --- Check 3: MT subset sensitivity -----------------------------------------
  n_all <- nrow(meta)
  subsets <- list(
    all_spots = rep(TRUE, n_all)
  )
  clean_label <- sprintf("percent_mt_lt_%g", mt_threshold)
  dirty_label <- sprintf("percent_mt_ge_%g", mt_threshold)
  subsets[[clean_label]] <- !is.na(meta[[mt_col]]) & meta[[mt_col]] < mt_threshold
  subsets[[dirty_label]] <- !is.na(meta[[mt_col]]) & meta[[mt_col]] >= mt_threshold

  # The absolute threshold is a weak test on its own: sample1 averages 3.7% MT,
  # so `percent.mt < 10` retains 98.6% of spots and the "clean subset" is very
  # nearly the full set. Quartile strata of percent.mt are added so the check can
  # actually discriminate — if the entropy-stemness association is degradation,
  # it should be concentrated in Q4 and absent in Q1.
  mt_vals <- meta[[mt_col]]
  mt_breaks <- unique(quantile(mt_vals, probs = seq(0, 1, 0.25), na.rm = TRUE))
  if (length(mt_breaks) >= 3) {
    mt_bin <- cut(mt_vals, breaks = mt_breaks, include.lowest = TRUE, labels = FALSE)
    for (q in sort(unique(mt_bin[!is.na(mt_bin)]))) {
      subsets[[sprintf("percent_mt_Q%d", q)]] <- !is.na(mt_bin) & mt_bin == q
    }
  }

  subset_pairs <- list()
  for (e in entropy_cols) {
    subset_pairs[[length(subset_pairs) + 1]] <- c(e, score_col)
    subset_pairs[[length(subset_pairs) + 1]] <- c(e, mt_col)
  }
  subset_pairs[[length(subset_pairs) + 1]] <- c(score_col, mt_col)

  subset_table <- do.call(rbind, lapply(names(subsets), function(s_name) {
    idx <- subsets[[s_name]]
    if (sum(idx) < 3) {
      warning(sprintf("Subset '%s' holds %d spots; skipping.", s_name, sum(idx)))
      return(NULL)
    }
    tab <- correlation_pair_table(meta[idx, , drop = FALSE], subset_pairs,
                                  sample_name = sample_name, subset_label = s_name)
    if (nrow(tab) == 0) return(NULL)
    mt_in_subset <- meta[[mt_col]][idx]
    cbind(tab[, 1:2, drop = FALSE],
          Spots_In_Subset = sum(idx),
          Pct_Spots_In_Subset = round(100 * sum(idx) / n_all, 2),
          MT_Min = round(min(mt_in_subset, na.rm = TRUE), 3),
          MT_Median = round(median(mt_in_subset, na.rm = TRUE), 3),
          MT_Max = round(max(mt_in_subset, na.rm = TRUE), 3),
          tab[, -(1:2), drop = FALSE])
  }))

  # --- Console summary ---------------------------------------------------------
  cat("Quality confound check (", sample_name, "):\n", sep = "")
  fmt <- function(v) if (length(v) == 0 || is.na(v)) "NA" else sprintf("%+.4f", v)
  for (t in qc_targets) {
    row <- covariate_table[covariate_table$X_Column == score_col & covariate_table$Y_Column == t, ]
    cat(sprintf("  %s vs %-16s r = %s, rho = %s\n", "Stemness", t, fmt(row$Pearson_r), fmt(row$Spearman_rho)))
  }
  if (!is.null(partial_table) && nrow(partial_table) > 0) {
    for (m in c("pearson", "spearman")) {
      row <- partial_table[partial_table$X_Variable == confound_var_label(primary_ent) &
                             partial_table$Covariate_Set == "percent_mt_and_log_nCount" &
                             partial_table$Method == m, ]
      cat(sprintf("  %s partial (entropy ~ stemness | %s + log %s): %s -> %s\n",
                  m, mt_col, count_col, fmt(row$Zero_Order_Estimate), fmt(row$Partial_Estimate)))
    }
  }
  if (!is.null(subset_table) && nrow(subset_table) > 0) {
    for (s_name in names(subsets)) {
      row <- subset_table[subset_table$Subset == s_name &
                            subset_table$X_Column == primary_ent &
                            subset_table$Y_Column == score_col, ]
      if (nrow(row) == 0) next
      cat(sprintf("  %-22s N = %5d, entropy vs stemness r = %s, rho = %s\n",
                  s_name, row$N, fmt(row$Pearson_r), fmt(row$Spearman_rho)))
    }
  }

  if (!save_outputs) {
    return(invisible(list(covariates = covariate_table, partial = partial_table, subsets = subset_table)))
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(covariate_table, file.path(output_dir, paste0(sample_name, "_quality_confound_covariates.csv")), row.names = FALSE)
  write.csv(partial_table, file.path(output_dir, paste0(sample_name, "_quality_confound_partial.csv")), row.names = FALSE)
  write.csv(subset_table, file.path(output_dir, paste0(sample_name, "_quality_confound_subsets.csv")), row.names = FALSE)

  # --- Figure ------------------------------------------------------------------
  anno <- function(tab, x_col, y_col, subset_label = "all_spots") {
    row <- tab[tab$X_Column == x_col & tab$Y_Column == y_col & tab$Subset == subset_label, ]
    if (nrow(row) == 0) return("")
    sprintf("N: %d\nPearson r: %.4f\nSpearman ρ: %.4f", row$N[1], row$Pearson_r[1], row$Spearman_rho[1])
  }
  panel <- function(df, x_col, y_col, title, label, color_col = NULL) {
    p <- if (is.null(color_col)) {
      ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
        geom_point(alpha = 0.4, color = "#2c3e50", size = 1.2)
    } else {
      ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_col]])) +
        geom_point(alpha = 0.6, size = 1.2) +
        scale_color_viridis_c(option = "magma", name = "percent.mt")
    }
    p + geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
      labs(title = title, subtitle = paste("Sample:", sample_name), x = x_col, y = y_col) +
      annotate("label", x = Inf, y = -Inf, label = label, hjust = 1.05, vjust = -0.2, size = 3.2,
               fill = alpha("white", 0.85), label.padding = grid::unit(0.35, "lines")) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
            plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray30"),
            axis.title = element_text(face = "bold", size = 9.5),
            panel.grid.minor = element_blank())
  }

  p1 <- panel(meta, mt_col, score_col, "Stemness Score vs Percent MT", anno(covariate_table, score_col, mt_col))
  ent_label <- confound_var_label(primary_ent)
  p2 <- panel(meta, mt_col, primary_ent, paste(ent_label, "vs Percent MT"), anno(covariate_table, primary_ent, mt_col))
  p3 <- panel(meta, score_col, primary_ent, paste(ent_label, "vs Stemness Score"),
              anno(covariate_table, primary_ent, score_col), color_col = mt_col)

  res_stats <- partial_correlation_stats(meta, primary_ent, score_col, c(mt_col, log_count_col), method = "pearson")
  p4 <- if (!is.null(res_stats$residual_x)) {
    res_df <- data.frame(residual_stemness = res_stats$residual_y, residual_entropy = res_stats$residual_x)
    panel(res_df, "residual_stemness", "residual_entropy",
          sprintf("Residuals after %s + log(%s)", mt_col, count_col),
          sprintf("N: %d\nPartial r: %.4f\nlog10(p): %.1f", res_stats$N, res_stats$Estimate, res_stats$Log10_p))
  } else {
    NULL
  }

  panels <- Filter(Negate(is.null), list(p1, p2, p3, p4))
  combined <- wrap_plots(panels, ncol = 2)
  ggsave(file.path(output_dir, paste0(sample_name, "_quality_confound_plot.png")),
         plot = combined, width = 11, height = 8, dpi = 300)
  cat(sprintf("Quality confound tables and figure written to %s\n", output_dir))

  invisible(list(covariates = covariate_table, partial = partial_table, subsets = subset_table))
}
