#!/usr/bin/env Rscript

# diagnose_correction_routes.R
# Stage A5a: do the two routes to a depth-corrected entropy agree?
#
# A3 and A4 left Stage A5 between two branches of the plan rather than on one of
# them. Chao-Shen is the best of the three corrections on every measured axis but
# does not clear the full A0 bar, so the choice is between correcting in the
# METRIC (adopt Chao-Shen) and correcting in the INFERENCE (keep the plug-in and
# make the log(nCount)-adjusted partial the headline number, D2's third branch).
#
# THE QUESTION. Route M's depth-corrected object is `chao_shen` itself. Route I's
# is the plug-in residualised on log(nCount). If those two objects rank spots the
# same way, the routes are doing the same work, the choice is low-stakes, and the
# agreement is itself a robustness result. If they diverge, they are not the same
# correction and one of them has to be argued for.
#
# This is the metric-level analogue of the hypothesis D3 records from the
# rarefaction era -- that partialling log(nCount) out of the plug-in reproduced
# what the corrected metric reported, i.e. that the depth confound rather than
# the estimator was doing the work. It needs no stemness score, so it can be
# answered now instead of waiting for Stage C.
#
# THE SECOND MEASUREMENT, and the reason this is not just one correlation.
# Entropy saturates logarithmically in N, so a LINEAR log(nCount) term can leave
# structured curvature behind. A high route agreement would not settle anything
# if Route I's residual still carried a depth trend. So each estimator's residual
# after partialling is also regressed on a spline of log(nCount): the R2 of that
# fit is residual depth structure the linear adjustment did not reach. `splines`
# is base R -- no new dependency (D2, D3).
#
# WHAT IS NOT RECOMPUTED. The reliability of a depth-partialled metric is already
# measured: `Split_Half_Rho_Given_Depth` in <sample>_subsampling_reliability.csv,
# from two independent 50% thinnings with both partialled on the spot's full-depth
# total. It is joined in here, not recomputed. Note that the rank partial is
# invariant to monotone transforms of the covariate, so partialling on N and on
# log(N) give the same number -- that is not a bug to be fixed.
#
# Usage:
#   Rscript diagnose_correction_routes.R [sample] [--all]
#
#   --all   run the frozen 16-sample cohort (D6) and write cohort_correction_*.csv

set.seed(23)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(splines)
})

source("R/shannon_entropy.R")
source("R/entropy_correlation.R")
source("R/quality_confound.R")   # partial_correlation_stats()
source("R/gene_universe.R")
source("R/load_sample.R")
source("R/cohort.R")
source("R/diagnostic_sweep.R")

ENTROPY_EXCLUDE <- "^(MT-|RP[SL])"

# A5a decision rule, fixed in docs/DECISIONS.md D2 BEFORE this script was first
# run, for the same reason the A0 criteria were. Read from here, never edited to
# match a result. If a number lands near a boundary it is reported as
# near-boundary; it is not re-thresholded.
#
#   Routes agree      rho(Route M, Route I) >= 0.80 on >= 12 of 16 samples.
#                     -> the routes correct the same thing; adopt Chao-Shen as
#                        the primary metric with the plug-in co-reported, and
#                        record the agreement as a robustness result.
#   Routes diverge    rho < 0.80 on more than 4 samples.
#                     -> not the same correction; the residual depth structure
#                        below decides which route to believe.
#   Independently     if the plug-in's linearly-partialled residual keeps spline
#                     R2 > 0.05 on a majority of samples while Chao-Shen's stays
#                     below it, linear partialling under-corrects and the
#                     metric-side correction is required regardless.
A5A_MIN_ROUTE_AGREEMENT   <- 0.80
A5A_MIN_AGREEING_SAMPLES  <- 12
A5A_MAX_RESIDUAL_SPLINE_R2 <- 0.05
A5A_SPLINE_DF             <- 4

# The depth covariate is nCount_Spatial, not the estimator's own MT/ribo-excluded
# total. The two are nearly identical, but nCount_Spatial is what the
# inference-side route would actually regress on downstream
# (run_quality_confound_check(), D3), and the point of Route I is to reproduce
# that route faithfully rather than an idealised version of it.
DEPTH_COL <- "nCount_Spatial"

#' Run the Stage A5a correction-route comparison on one sample
#'
#' @param target_sample Sample name
#' @return Named list of the tables written, for the cohort rbind
run_correction_routes <- function(target_sample) {

cat("=================================================================\n")
cat("Stage A5a correction routes:", target_sample, "\n")
cat("=================================================================\n")

spatial_obj <- load_qc_sample(target_sample)$obj
meta <- spatial_obj@meta.data
counts_mat <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
counts_filtered <- exclude_gene_families(counts_mat, ENTROPY_EXCLUDE,
                                         context = "correction routes")
rm(spatial_obj, counts_mat)
gc(verbose = FALSE)

est <- entropy_estimator_matrix(counts_filtered, return_stats = TRUE)
rm(counts_filtered)
gc(verbose = FALSE)

d <- data.frame(
  log_nCount = log(meta[[DEPTH_COL]]),
  nCount     = meta[[DEPTH_COL]],
  percent_mt = meta$percent.mt,
  stringsAsFactors = FALSE
)
for (e in ENTROPY_ESTIMATORS) d[[e]] <- est[[e]]

# Drop incomplete rows once, up front, so every call below residualises the same
# spots. partial_correlation_stats() does its own complete-case filtering, so
# without this the one-covariate and two-covariate residuals could come back on
# different row sets and be silently incomparable.
keep <- Reduce(`&`, lapply(d, function(v) !is.na(v) & is.finite(v)))
if (any(!keep)) {
  cat(sprintf("Dropping %d spot(s) with non-finite depth or percent.mt.\n", sum(!keep)))
  d <- d[keep, , drop = FALSE]
}

cat(sprintf("Spots: %d   median depth: %.0f\n", nrow(d), median(d$nCount)))

# Spline basis for log(nCount), as plain columns so partial_correlation_stats()
# can take them as covariates unchanged rather than growing a formula interface.
spline_basis <- splines::ns(d$log_nCount, df = A5A_SPLINE_DF)
SPLINE_COLS <- paste0("ns_logN_", seq_len(ncol(spline_basis)))
for (j in seq_along(SPLINE_COLS)) d[[SPLINE_COLS[j]]] <- spline_basis[, j]

# --- 1. Residualise on depth (Route I) --------------------------------------
#
# partial_correlation_stats() already does the residualisation, the rank
# transform and the degrees-of-freedom bookkeeping (R/quality_confound.R). It is
# called with the estimator as x and the plug-in as y so that a single call
# yields BOTH residuals plus their partial correlation, rather than
# re-implementing lm.fit here.

# The function residualises BOTH of its variables, so pairing each candidate
# estimator with the plug-in yields the candidate's residual, the plug-in's
# residual, and their partial correlation from one call. Pairing an estimator
# with `log_nCount` itself would not work: the covariate regressed on itself has
# zero residual variance and the function correctly bails out.
resid_depth <- list()        # residual on log(nCount)
resid_depth_mt <- list()     # residual on log(nCount) + percent.mt
resid_depth_rank <- list()
both_adjusted <- list()      # partial correlation vs the plug-in, both adjusted

candidates <- setdiff(ENTROPY_ESTIMATORS, "plugin")
for (e in candidates) {
  st      <- partial_correlation_stats(d, e, "plugin", "log_nCount", method = "pearson")
  st_mt   <- partial_correlation_stats(d, e, "plugin", c("log_nCount", "percent_mt"), method = "pearson")
  st_rank <- partial_correlation_stats(d, e, "plugin", "log_nCount", method = "spearman")

  resid_depth[[e]]      <- st$residual_x
  resid_depth_mt[[e]]   <- st_mt$residual_x
  resid_depth_rank[[e]] <- st_rank$residual_x

  # Identical on every iteration; assigned from the last pairing.
  resid_depth[["plugin"]]      <- st$residual_y
  resid_depth_mt[["plugin"]]   <- st_mt$residual_y
  resid_depth_rank[["plugin"]] <- st_rank$residual_y

  both_adjusted[[e]] <- list(pearson = st$Estimate, spearman = st_rank$Estimate)
}

stopifnot(!is.null(resid_depth[["plugin"]]),
          length(resid_depth[["plugin"]]) == nrow(d))

# --- 2. Route agreement ------------------------------------------------------
#
# Primary comparison: Route M's corrected object (chao_shen, unadjusted) against
# Route I's corrected object (the plug-in residual). Secondary: both adjusted,
# which isolates whether Chao-Shen carries structure the partial does not reach.

route_i        <- resid_depth[["plugin"]]
route_i_mt     <- resid_depth_mt[["plugin"]]
route_i_rank   <- resid_depth_rank[["plugin"]]

agreement <- do.call(rbind, lapply(candidates, function(e) {
  route_m <- d[[e]]
  data.frame(
    Sample = target_sample, Estimator = e, N = nrow(d),
    Pearson_r_RouteM_vs_RouteI  = cor(route_m, route_i),
    Spearman_rho_RouteM_vs_RouteI = cor(rank(route_m), route_i_rank),
    Pearson_r_RouteM_vs_RouteI_mt = cor(route_m, route_i_mt),
    # Both adjusted: ~1 means the metric-side correction adds nothing once the
    # partial has been taken; < 1 means it carries structure the partial misses.
    # Taken from partial_correlation_stats() rather than recomputed.
    Pearson_r_BothAdjusted  = both_adjusted[[e]]$pearson,
    Spearman_rho_BothAdjusted = both_adjusted[[e]]$spearman,
    stringsAsFactors = FALSE
  )
}))

# --- 3. Residual depth structure --------------------------------------------
#
# Spline R2 of the residual on log(nCount) AFTER the linear term has already been
# removed: whatever it picks up is curvature the linear adjustment could not
# reach. The raw column is the same fit on the unadjusted metric, as context for
# how much depth there was to remove in the first place.

spline_r2 <- function(y) {
  fit <- stats::lm(y ~ splines::ns(d$log_nCount, df = A5A_SPLINE_DF))
  summary(fit)$r.squared
}

structure_df <- do.call(rbind, lapply(ENTROPY_ESTIMATORS, function(e) {
  r <- resid_depth[[e]]
  # The percent.mt association at zero order AND with depth removed. A0
  # criterion 5(a) checks only the zero-order form, and the cohort run showed
  # that is not enough: on samples where the plug-in's depth bias runs opposite
  # to a real H-MT relation the two cancel, and the estimator is credited with a
  # clean MT association it does not have. The partial is what the association
  # actually is once depth is not doing the masking, so both are recorded.
  mt_partial <- partial_correlation_stats(d, e, "percent_mt", "log_nCount",
                                          method = "pearson")
  # And again against a SPLINE of depth. The linear partial is not a fair basis
  # for comparing estimators on this axis: it leaves most of the plug-in's depth
  # structure in place (Spline_R2_Residual, median 0.60) while leaving almost
  # none of Chao-Shen's (median 0.04), and because percent.mt is itself coupled
  # to depth (cohort rho = -0.41, Stage 0) that leftover depth leaks into the
  # plug-in's MT partial and suppresses it. Adjusting both on the same spline
  # basis removes the asymmetry, so this is the column to compare across
  # estimators; the linear one is kept because it is what the inference-side
  # route would actually apply.
  mt_partial_spline <- partial_correlation_stats(d, e, "percent_mt", SPLINE_COLS,
                                                 method = "pearson")
  data.frame(
    Sample = target_sample, Estimator = e, N = nrow(d),
    Spearman_rho_raw_vs_depth      = cor(d[[e]], d$nCount, method = "spearman"),
    Spearman_rho_residual_vs_depth = cor(r, d$nCount, method = "spearman"),
    Spline_R2_Raw       = spline_r2(d[[e]]),
    Spline_R2_Residual  = spline_r2(r),
    Pearson_r_mt            = cor(d[[e]], d$percent_mt),
    Pearson_r_mt_given_depth = mt_partial$Estimate,
    Pearson_r_mt_given_depth_spline = mt_partial_spline$Estimate,
    SD_Raw      = sd(d[[e]]),
    SD_Residual = sd(r),
    stringsAsFactors = FALSE
  )
}))

# --- 4. Residual shape by depth decile (for the figure) ---------------------

decile <- cut(d$log_nCount, breaks = quantile(d$log_nCount, probs = seq(0, 1, 0.1)),
              include.lowest = TRUE, labels = FALSE)
deciles_df <- do.call(rbind, lapply(ENTROPY_ESTIMATORS, function(e) {
  r <- resid_depth[[e]]
  data.frame(
    Sample = target_sample, Estimator = e, Decile = seq_len(10),
    Median_nCount = as.numeric(tapply(d$nCount, decile, median)),
    Mean_Residual = as.numeric(tapply(r, decile, mean)),
    SE_Residual = as.numeric(tapply(r, decile, function(v) sd(v) / sqrt(length(v)))),
    stringsAsFactors = FALSE
  )
}))

# --- 5. Join the reliability already measured at A4 --------------------------

out_dir <- file.path("results", "statistical_tests")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

rel_file <- file.path(out_dir, paste0(target_sample, "_subsampling_reliability.csv"))
if (file.exists(rel_file)) {
  rel <- read.csv(rel_file, stringsAsFactors = FALSE)
  structure_df <- merge(structure_df,
                        rel[, c("Metric", "Split_Half_Rho_Given_Depth")],
                        by.x = "Estimator", by.y = "Metric", all.x = TRUE, sort = FALSE)
} else {
  cat(sprintf("NOTE: %s absent; run diagnose_subsampling_stability.R for the reliability column.\n",
              basename(rel_file)))
  structure_df$Split_Half_Rho_Given_Depth <- NA_real_
}

routes_df <- merge(structure_df, agreement, by = c("Sample", "Estimator", "N"),
                   all.x = TRUE, sort = FALSE)
routes_df$Routes_Agree <- routes_df$Spearman_rho_RouteM_vs_RouteI >= A5A_MIN_ROUTE_AGREEMENT
routes_df$Residual_Depth_Structure <- routes_df$Spline_R2_Residual > A5A_MAX_RESIDUAL_SPLINE_R2

write.csv(routes_df, file.path(out_dir, paste0(target_sample, "_correction_routes.csv")), row.names = FALSE)
write.csv(deciles_df, file.path(out_dir, paste0(target_sample, "_correction_routes_deciles.csv")), row.names = FALSE)

# --- 6. Figure ---------------------------------------------------------------

scatter_df <- do.call(rbind, lapply(candidates, function(e) {
  data.frame(Estimator = e, Route_M = d[[e]], Route_I = route_i, stringsAsFactors = FALSE)
}))

p_scatter <- ggplot(scatter_df, aes(x = Route_I, y = Route_M)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "#2c3e50") +
  facet_wrap(~Estimator, nrow = 1, scales = "free_y") +
  labs(title = sprintf("Do the two correction routes agree? -- %s", target_sample),
       subtitle = "x: plug-in residualised on log(nCount) (inference-side).  y: the corrected estimator itself (metric-side).",
       x = "Route I: plug-in | log(nCount)", y = "Route M: corrected estimator (bits)") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12),
                     strip.text = element_text(face = "bold"))

deciles_df$Estimator <- factor(deciles_df$Estimator, levels = ENTROPY_ESTIMATORS)
p_dec <- ggplot(deciles_df, aes(x = Median_nCount, y = Mean_Residual, colour = Estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = Mean_Residual - 2 * SE_Residual,
                    ymax = Mean_Residual + 2 * SE_Residual), width = 0, alpha = 0.6) +
  geom_line() + geom_point(size = 1.6) +
  scale_x_log10() +
  labs(title = "Residual depth structure the linear adjustment did not reach",
       subtitle = "Mean residual per depth decile after regressing on log(nCount). A flat line at zero means fully adjusted.",
       x = "Median nCount in decile (log scale)", y = "Mean residual (bits)") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(out_dir, paste0(target_sample, "_correction_routes_plot.png")),
       plot = p_scatter / p_dec, width = 13, height = 9, dpi = 200)

# --- 7. Report ---------------------------------------------------------------

cat("\n=================================================================\n")
cat("Route agreement (Route M = corrected estimator, Route I = plug-in | log nCount):\n")
cat("=================================================================\n")
print(agreement[, c("Estimator", "Pearson_r_RouteM_vs_RouteI",
                    "Spearman_rho_RouteM_vs_RouteI", "Pearson_r_RouteM_vs_RouteI_mt",
                    "Spearman_rho_BothAdjusted")], row.names = FALSE, digits = 4)

cat("\n=================================================================\n")
cat(sprintf("Residual depth structure (spline df %d; flagged above R2 %.2f):\n",
            A5A_SPLINE_DF, A5A_MAX_RESIDUAL_SPLINE_R2))
cat("=================================================================\n")
print(routes_df[, c("Estimator", "Spearman_rho_raw_vs_depth",
                    "Spearman_rho_residual_vs_depth", "Spline_R2_Raw",
                    "Spline_R2_Residual", "Residual_Depth_Structure",
                    "Split_Half_Rho_Given_Depth")], row.names = FALSE, digits = 4)

cat("\npercent.mt association: zero-order, linear depth partial, spline depth partial:\n")
print(routes_df[, c("Estimator", "Pearson_r_mt", "Pearson_r_mt_given_depth",
                    "Pearson_r_mt_given_depth_spline")],
      row.names = FALSE, digits = 4)

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n",
            file.path(out_dir, paste0(target_sample, "_correction_routes.csv")),
            file.path(out_dir, paste0(target_sample, "_correction_routes_deciles.csv")),
            file.path(out_dir, paste0(target_sample, "_correction_routes_plot.png"))))

invisible(list(correction_routes = routes_df,
               correction_routes_deciles = deciles_df))
}

# --- Dispatch ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
sweep_spec <- resolve_sweep_samples(args)

sweep <- run_sweep(sweep_spec$samples, run_correction_routes, label = "A5a routes")

if (sweep_spec$all) {
  out_dir <- file.path("results", "statistical_tests")
  cat("\n=================================================================\n")
  cat("Cohort tables:\n")
  cat("=================================================================\n")
  cohort <- write_cohort_tables(
    sweep, c("correction_routes", "correction_routes_deciles"), out_dir)

  routes <- cohort$correction_routes
  if (!is.null(routes)) {
    cat("\n=================================================================\n")
    cat("A5a decision rule (fixed before this run; see D2)\n")
    cat("=================================================================\n")

    verdict <- do.call(rbind, lapply(unique(routes$Estimator), function(e) {
      sub <- routes[routes$Estimator == e, ]
      n_agree <- sum(sub$Routes_Agree, na.rm = TRUE)
      n_struct <- sum(sub$Residual_Depth_Structure, na.rm = TRUE)
      data.frame(
        Estimator = e, N_Samples = nrow(sub),
        Median_Rho_RouteM_vs_RouteI = median(sub$Spearman_rho_RouteM_vs_RouteI, na.rm = TRUE),
        Min_Rho_RouteM_vs_RouteI = min(sub$Spearman_rho_RouteM_vs_RouteI, na.rm = TRUE),
        N_Samples_Routes_Agree = n_agree,
        Median_Spline_R2_Residual = median(sub$Spline_R2_Residual, na.rm = TRUE),
        N_Samples_Residual_Depth_Structure = n_struct,
        Median_Split_Half_Rho_Given_Depth = median(sub$Split_Half_Rho_Given_Depth, na.rm = TRUE),
        Median_Abs_r_mt = median(abs(sub$Pearson_r_mt), na.rm = TRUE),
        Median_Abs_r_mt_Given_Depth = median(abs(sub$Pearson_r_mt_given_depth), na.rm = TRUE),
        Max_Abs_r_mt_Given_Depth = max(abs(sub$Pearson_r_mt_given_depth), na.rm = TRUE),
        Median_Abs_r_mt_Given_Depth_Spline = median(abs(sub$Pearson_r_mt_given_depth_spline), na.rm = TRUE),
        Max_Abs_r_mt_Given_Depth_Spline = max(abs(sub$Pearson_r_mt_given_depth_spline), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))

    plugin_struct <- routes[routes$Estimator == "plugin", ]
    verdict_plugin_n <- sum(plugin_struct$Residual_Depth_Structure, na.rm = TRUE)

    write.csv(verdict, file.path(out_dir, "cohort_correction_routes_verdict.csv"), row.names = FALSE)
    cat(sprintf("  %s\n\n", file.path(out_dir, "cohort_correction_routes_verdict.csv")))
    print(verdict, row.names = FALSE, digits = 4)

    cs <- verdict[verdict$Estimator == "chao_shen", ]
    if (nrow(cs) == 1) {
      cat("\n-----------------------------------------------------------------\n")
      cat(sprintf("Chao-Shen vs Route I: agree on %d/%d samples (rule needs >= %d).  -> %s\n",
                  cs$N_Samples_Routes_Agree, cs$N_Samples, A5A_MIN_AGREEING_SAMPLES,
                  if (cs$N_Samples_Routes_Agree >= A5A_MIN_AGREEING_SAMPLES) "ROUTES AGREE" else "ROUTES DIVERGE"))
      cat(sprintf("Plug-in residual keeps depth structure on %d/%d samples; Chao-Shen's on %d/%d.\n",
                  verdict_plugin_n, nrow(plugin_struct),
                  cs$N_Samples_Residual_Depth_Structure, cs$N_Samples))
      if (verdict_plugin_n > nrow(plugin_struct) / 2 &&
          cs$N_Samples_Residual_Depth_Structure <= cs$N_Samples / 2) {
        cat("  -> Linear partialling UNDER-CORRECTS: metric-side correction required.\n")
      }
      cat("-----------------------------------------------------------------\n")
    }
  }
}

report_sweep_failures(sweep, length(sweep_spec$samples))
cat("Stage A5a diagnostic completed.\n")
