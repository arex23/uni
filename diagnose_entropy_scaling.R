#!/usr/bin/env Rscript

# diagnose_entropy_scaling.R
# Stage A3: estimator bake-off. Compares the four diversity estimators of D2
# against sequencing depth and the quality covariates, and evaluates the
# pre-registered A0 acceptance criteria that can be checked without downstream
# input (criteria 1, 2, 4 and 5a).
#
# This script previously answered the Stage 2 question -- raw vs SpaNorm logpac
# vs linear back-transform -- which D2 has settled: entropy on a log-scale layer
# is a relabelled gene counter. Those outputs are kept on disk under
# results/statistical_tests/<sample>_entropy_scaling_*.csv and are what D2's
# diagnostic table cites; they are not regenerated here.
#
# Deliberately does NOT run SpaNorm. Every estimator in the bake-off reads raw
# counts, so the script loads and QCs the sample directly through
# load_qc_sample() (D1 steps 1-3, R/load_sample.R). That is what makes the
# three-sample requirement affordable: ~15 seconds and ~1 GB per sample instead
# of a SpaNorm fit at ~7.6 GB peak whose output would go unused.
#
# Usage:
#   Rscript diagnose_entropy_scaling.R [sample] [--all] [--from-rds]
#
#   --all        run the frozen 16-sample cohort (D6) instead of one sample, and
#                additionally write cohort_estimator_*.csv. A0 criteria are then
#                summarised as "passes on k of 16" -- the k/16 counts are Stage
#                A5 evidence to read, not a threshold to tune.
#   --from-rds   load results/seurat_objects/<sample>_spatial_obj.rds instead of
#                re-running the QC. Only useful for cross-checking a sample that
#                has already been through analyze_entropy.R.

set.seed(23)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

source("R/shannon_entropy.R")
source("R/entropy_correlation.R")
source("R/gene_universe.R")
source("R/load_sample.R")
source("R/cohort.R")
source("R/diagnostic_sweep.R")

ENTROPY_EXCLUDE <- "^(MT-|RP[SL])"

# A0 criteria, frozen in docs/DECISIONS.md D2 before this script was written.
# Read from here, never edited to match a result.
A0_MAX_ABS_RHO_DEPTH   <- 0.30   # criterion 1
A0_MIN_RHO_NFEATURE    <- -0.15  # criterion 2 (no sign flip)
A0_MIN_IQR             <- 0.10   # criterion 4
A0_MAX_MT_EXCESS       <- 0.10   # criterion 5a, vs the plug-in baseline

# RETIRED: A0_MIN_CV <- 0.02, the other half of criterion 4.
#
# Retired on a validity argument, not because it was inconvenient -- see D2.
# CV = SD/mean is a ratio-scale statistic. Entropy in bits has no meaningful
# zero: the origin is set by the size of the gene universe (log2 K), which is a
# property of the reference, not of the tissue. Rescaling the universe shifts
# every H by a constant and changes the CV while changing nothing biological.
# The mean therefore has no interpretation as a unit, and SD/mean has none
# either. IQR in bits is on the scale the metric is actually read on and is the
# defensible half.
#
# The CV half is still COMPUTED and reported below so the record is intact and
# the original verdict is reproducible from this script. It just no longer
# gates. Both verdicts appear in the output.
A0_MIN_CV_RETIRED      <- 0.02

#' Run the Stage A3 estimator bake-off on one sample
#'
#' @param target_sample Sample name
#' @param from_rds Load the saved Seurat object instead of re-running D1 QC
#' @return Named list of the tables written, for the cohort rbind
run_estimator_bakeoff <- function(target_sample, from_rds = FALSE) {

cat("=================================================================\n")
cat("Stage A3 estimator bake-off:", target_sample, "\n")
cat("=================================================================\n")

# --- 1. Load ---------------------------------------------------------------

if (from_rds) {
  seurat_file <- file.path("results", "seurat_objects", paste0(target_sample, "_spatial_obj.rds"))
  if (!file.exists(seurat_file)) {
    stop(sprintf("Seurat object '%s' not found. Drop --from-rds to load from data/ instead.", seurat_file))
  }
  cat("Loading saved object:", seurat_file, "\n")
  spatial_obj <- readRDS(seurat_file)
} else {
  cat("Loading and QC'ing from data/ (D1 steps 1-3, no SpaNorm)\n")
  spatial_obj <- load_qc_sample(target_sample)$obj
}
cat(sprintf("Spots: %d   Features: %d\n", ncol(spatial_obj), nrow(spatial_obj)))

# --- 2. Estimators, one blocked pass ---------------------------------------

counts_mat <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
counts_filtered <- exclude_gene_families(counts_mat, ENTROPY_EXCLUDE,
                                         context = "estimator bake-off")
rm(counts_mat)
gc(verbose = FALSE)

cat("\nComputing all four estimators in one pass...\n")
t0 <- Sys.time()
est <- entropy_estimator_matrix(counts_filtered, return_stats = TRUE)
cat(sprintf("  done in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# n_detected is the support the estimators actually see: genes with a non-zero
# raw count on the MT/ribo-excluded set. nFeature_Spatial counts the full gene
# set including MT, so the two are close but not identical, and both are carried
# as targets because the plug-in ceiling is driven by the former.
spatial_obj$n_detected <- est$K
spatial_obj$spot_total_counts <- est$N
spatial_obj$singleton_fraction <- ifelse(est$N > 0, est$f1 / est$N, NA_real_)
spatial_obj$coverage_gt <- 1 - spatial_obj$singleton_fraction

est_cols <- character(0)
for (e in ENTROPY_ESTIMATORS) {
  cl <- paste0("entropy_", e)
  spatial_obj[[cl]] <- est[[e]]
  est_cols <- c(est_cols, cl)
}

out_dir <- file.path("results", "statistical_tests")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 3. Ordering invariants on real data (A2's real-data half) -------------

cat("\nOrdering invariants:\n")
ordering <- check_entropy_ordering(est, stop_on_violation = TRUE)
ordering$Sample <- target_sample
print(ordering)
write.csv(ordering[, c("Sample", "Comparison", "Invariant", "N_Violations", "N_Spots")],
          file.path(out_dir, paste0(target_sample, "_estimator_ordering.csv")), row.names = FALSE)

# --- 4. Correlations against depth and quality covariates ------------------

targets_a3 <- c(
  nCounts = "nCount_Spatial",
  nFeatures = "nFeature_Spatial",
  n_detected = "n_detected",
  Percent_MT = "percent.mt"
)

comparison <- calculate_entropy_correlations(
  seurat_obj = spatial_obj,
  entropy_cols = est_cols,
  targets = targets_a3,
  sample_name = target_sample,
  output_dir = out_dir,
  file_suffix = "_estimator_comparison",
  save_outputs = TRUE
)
summary_table <- comparison$summary_table

# --- 5. Distribution block (criterion 4) -----------------------------------

dist_df <- do.call(rbind, lapply(ENTROPY_ESTIMATORS, function(e) {
  H <- est[[e]]
  m <- mean(H)
  data.frame(
    Sample = target_sample, Estimator = e, N = length(H),
    Mean = m, SD = sd(H), CV = if (m != 0) sd(H) / m else NA_real_,
    Median = median(H), IQR = IQR(H),
    Min = min(H), Max = max(H), Range = diff(range(H)),
    stringsAsFactors = FALSE
  )
}))
write.csv(dist_df, file.path(out_dir, paste0(target_sample, "_estimator_distribution.csv")), row.names = FALSE)

# Per-spot count statistics. The singleton fraction f1/N is what drives the
# Chao-Shen correction magnitude, so its own coupling to depth is the mechanism
# behind the overcorrection risk this plan names: shallow spots have a higher
# f1/N and are therefore corrected harder. Recorded so that a sign flip in
# criterion 2 can be attributed rather than merely observed.
cov_stats <- pair_correlation_stats(spatial_obj$singleton_fraction, spatial_obj$nCount_Spatial)
spot_stats_df <- data.frame(
  Sample = target_sample, N = nrow(est),
  Median_N = median(est$N), Median_K = median(est$K),
  Median_f1 = median(est$f1), Median_f2 = median(est$f2),
  Median_Singleton_Fraction = median(spatial_obj$singleton_fraction, na.rm = TRUE),
  Median_Coverage_GT = median(spatial_obj$coverage_gt, na.rm = TRUE),
  Min_Coverage_GT = min(spatial_obj$coverage_gt, na.rm = TRUE),
  N_Spots_f1_eq_N = sum(est$f1 == est$N & est$N > 0),
  Spearman_SingletonFrac_vs_nCount = cov_stats$spearman_rho,
  stringsAsFactors = FALSE
)
write.csv(spot_stats_df, file.path(out_dir, paste0(target_sample, "_estimator_spot_stats.csv")), row.names = FALSE)

# --- 6. A0 criteria ---------------------------------------------------------

get_stat <- function(estimator, target, column) {
  lab <- entropy_col_label(paste0("entropy_", estimator))
  row <- summary_table[summary_table$Entropy_Metric == lab &
                         summary_table$Target_Variable == target, ]
  if (nrow(row) == 0) NA_real_ else row[[column]][1]
}

mt_plugin <- abs(get_stat("plugin", "Percent_MT", "Pearson_r"))

criteria_df <- do.call(rbind, lapply(ENTROPY_ESTIMATORS, function(e) {
  rho_ncount <- get_stat(e, "nCounts", "Spearman_rho")
  rho_nfeat  <- get_stat(e, "nFeatures", "Spearman_rho")
  r_mt       <- get_stat(e, "Percent_MT", "Pearson_r")
  d <- dist_df[dist_df$Estimator == e, ]

  c1 <- abs(rho_ncount) < A0_MAX_ABS_RHO_DEPTH && abs(rho_nfeat) < A0_MAX_ABS_RHO_DEPTH
  c2 <- rho_nfeat > A0_MIN_RHO_NFEATURE
  c4 <- d$IQR >= A0_MIN_IQR
  c4_orig <- d$CV >= A0_MIN_CV_RETIRED && d$IQR >= A0_MIN_IQR
  c5a <- abs(r_mt) <= mt_plugin + A0_MAX_MT_EXCESS

  data.frame(
    Sample = target_sample, Estimator = e,
    Spearman_rho_nCount = rho_ncount, Spearman_rho_nFeature = rho_nfeat,
    CV = d$CV, IQR = d$IQR, Pearson_r_percent_mt = r_mt,
    Pearson_r_percent_mt_plugin = mt_plugin,
    C1_Depth_Decoupling = c1, C2_No_Sign_Flip = c2,
    C4_Dynamic_Range = c4, C5a_MT_Not_Worse = c5a,
    Passes_A3_Criteria = c1 && c2 && c4 && c5a,
    C4_Dynamic_Range_CV_Retired = c4_orig,
    Passes_A3_Criteria_As_Preregistered = c1 && c2 && c4_orig && c5a,
    stringsAsFactors = FALSE
  )
}))
write.csv(criteria_df, file.path(out_dir, paste0(target_sample, "_estimator_criteria.csv")), row.names = FALSE)

# --- 7. Report --------------------------------------------------------------

cat("\n=================================================================\n")
cat("Correlations (Spearman rho):\n")
cat("=================================================================\n")
wide <- reshape(summary_table[, c("Entropy_Metric", "Target_Variable", "Spearman_rho")],
                idvar = "Entropy_Metric", timevar = "Target_Variable", direction = "wide")
colnames(wide) <- sub("^Spearman_rho\\.", "", colnames(wide))
print(wide, row.names = FALSE, digits = 4)

cat("\n=================================================================\n")
cat("Distributions (criterion 4: IQR >= 0.10 bits; CV reported, retired as a gate):\n")
cat("=================================================================\n")
print(dist_df[, c("Estimator", "Mean", "SD", "CV", "Median", "IQR", "Range")], row.names = FALSE, digits = 4)

cat("\n=================================================================\n")
cat("Per-spot count statistics:\n")
cat("=================================================================\n")
print(spot_stats_df, row.names = FALSE, digits = 4)

cat("\n=================================================================\n")
cat("A0 criteria checkable at A3 (criterion 3 needs A4, 5b needs Stage C):\n")
cat("=================================================================\n")
print(criteria_df[, c("Estimator", "Spearman_rho_nCount", "Spearman_rho_nFeature", "IQR",
                      "Pearson_r_percent_mt", "C1_Depth_Decoupling", "C2_No_Sign_Flip",
                      "C4_Dynamic_Range", "C5a_MT_Not_Worse", "Passes_A3_Criteria")],
      row.names = FALSE, digits = 4)

if (!identical(criteria_df$Passes_A3_Criteria,
               criteria_df$Passes_A3_Criteria_As_Preregistered)) {
  cat("\nRetiring the CV half of criterion 4 changed the verdict for:\n")
  ch <- criteria_df[criteria_df$Passes_A3_Criteria !=
                      criteria_df$Passes_A3_Criteria_As_Preregistered, ]
  print(ch[, c("Estimator", "CV", "IQR", "C4_Dynamic_Range_CV_Retired",
               "C4_Dynamic_Range", "Passes_A3_Criteria_As_Preregistered",
               "Passes_A3_Criteria")], row.names = FALSE, digits = 4)
}

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n  %s\n  %s\n",
            file.path(out_dir, paste0(target_sample, "_estimator_comparison.csv")),
            file.path(out_dir, paste0(target_sample, "_estimator_distribution.csv")),
            file.path(out_dir, paste0(target_sample, "_estimator_spot_stats.csv")),
            file.path(out_dir, paste0(target_sample, "_estimator_ordering.csv")),
            file.path(out_dir, paste0(target_sample, "_estimator_criteria.csv"))))

invisible(list(
  estimator_criteria = criteria_df,
  estimator_distribution = dist_df,
  estimator_spot_stats = spot_stats_df,
  estimator_ordering = ordering[, c("Sample", "Comparison", "Invariant",
                                    "N_Violations", "N_Spots")],
  estimator_comparison = summary_table
))
}

# --- Dispatch ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
from_rds <- "--from-rds" %in% args
args <- args[args != "--from-rds"]
sweep_spec <- resolve_sweep_samples(args)

sweep <- run_sweep(sweep_spec$samples,
                   function(s) run_estimator_bakeoff(s, from_rds = from_rds),
                   label = "A3 bake-off")

if (sweep_spec$all) {
  out_dir <- file.path("results", "statistical_tests")
  cat("\n=================================================================\n")
  cat("Cohort tables:\n")
  cat("=================================================================\n")
  cohort <- write_cohort_tables(
    sweep,
    c("estimator_criteria", "estimator_distribution", "estimator_spot_stats",
      "estimator_ordering", "estimator_comparison"),
    out_dir)

  if (!is.null(cohort$estimator_criteria)) {
    crit_summary <- summarise_criteria(
      cohort$estimator_criteria,
      c("C1_Depth_Decoupling", "C2_No_Sign_Flip", "C4_Dynamic_Range",
        "C5a_MT_Not_Worse", "Passes_A3_Criteria",
        "C4_Dynamic_Range_CV_Retired", "Passes_A3_Criteria_As_Preregistered"))
    write.csv(crit_summary,
              file.path(out_dir, "cohort_estimator_criteria_summary.csv"),
              row.names = FALSE)
    cat(sprintf("  %s\n", file.path(out_dir, "cohort_estimator_criteria_summary.csv")))

    cat("\n=================================================================\n")
    cat("A0 criteria across the cohort: passes on k of n\n")
    cat("=================================================================\n")
    print(crit_summary[, c("Estimator", "Criterion", "N_Pass", "N_Samples",
                           "Pct_Pass", "Failing_Samples")],
          row.names = FALSE)
  }
}

report_sweep_failures(sweep, length(sweep_spec$samples))
cat("Stage A3 diagnostic completed.\n")
