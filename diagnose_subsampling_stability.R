#!/usr/bin/env Rscript

# diagnose_subsampling_stability.R
# Stage A4: the decisive depth-invariance test, and the evidence for A0
# criterion 3.
#
# Every spot's counts are downsampled by a COMMON PROPORTION and all four
# estimators of D2 are recomputed at each level. A depth-unbiased estimator
# returns approximately the same value at 50% depth as at 100%; a downward-
# biased one drops systematically. This is the test that would have caught
# rarefaction's problem in advance.
#
# TWO FRAMING POINTS, which belong in D2 verbatim:
#
# 1. This uses downsampling as a MEASURING INSTRUMENT, not as the estimator.
#    The project is not rarefying; it is perturbing depth in a controlled way to
#    see which estimator is invariant to it. That is the answer to "didn't you
#    say rarefaction was inappropriate?" -- and it is a better answer than a
#    citation, because it is measured on this project's own data.
#
# 2. Common-proportion downsampling preserves the RELATIVE depth structure
#    across spots. Rarefaction to a common depth destroys it: every spot is
#    forced onto the same N, which is the assumption D2 rejected. Halving every
#    spot leaves the deep spots deep and the shallow spots shallow, so this
#    diagnostic does not smuggle the rejected assumption back in.
#
# HOW THE DOWNSAMPLING IS DONE. Independent binomial thinning of each stored
# count, y' ~ Binomial(y, prop). `scuttle::downsampleMatrix()`, which the removed
# rarefaction code used, instead samples without replacement per column so that
# each spot lands on exactly round(prop * N) counts; `scuttle` is no longer a
# dependency (D2) and thinning is the better fit regardless. The realized total
# is then random with mean prop * N -- at N ~ 1e4 that is about 1% of variation,
# immaterial -- and the point here is to perturb depth, not to hit a target
# depth. Hitting a target depth is precisely the rarefaction assumption.
#
# Thinning is applied to the MT/ribo-excluded matrix the estimators actually
# read. Because binomial thinning is independent per entry, that is
# distributionally identical to thinning the full matrix and excluding
# afterwards, and it avoids a second copy of the largest object in the script.
#
# ONE DRAW PER LEVEL, seeded. Averaging over draws would suppress sampling
# noise, but a single draw is what an experiment at that depth would actually
# have produced. The signed median dH separates systematic bias from noise
# anyway -- noise is symmetric about zero, bias is not -- which is why both the
# absolute and the signed medians are reported.
#
# Usage:
#   Rscript diagnose_subsampling_stability.R [sample] [--all] [--props=0.75,0.5,0.25]
#
#   --all   run the frozen 16-sample cohort (D6) instead of one sample, and
#           additionally write cohort_subsampling_*.csv. The RNG is reseeded per
#           sample by run_sweep(), so a sample's thinning draw does not depend on
#           its position in the loop and --all reproduces the single-sample run.

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
source("R/quality_confound.R")   # partial_correlation_stats()
source("R/cohort.R")
source("R/diagnostic_sweep.R")

ENTROPY_EXCLUDE <- "^(MT-|RP[SL])"

# A0 criterion 3, frozen in docs/DECISIONS.md D2. Read from here, never edited
# to match a result.
A0_MAX_DELTA_OVER_SD <- 0.20   # median |dH| at 50% depth, as a fraction of full-depth SD
A0_REFERENCE_PROP    <- 0.50

# RETIRED: A0_MIN_RANK_RHO <- 0.90, the rank half of criterion 3.
#
# Retired on a validity argument, not because it was inconvenient -- see D2.
# Common-proportion thinning preserves the relative depth ordering of the spots
# almost exactly (binomial noise on a spot's total is ~0.5*sqrt(N), negligible
# against the spread of N across spots). So ANY metric that is a monotone
# function of depth reproduces its own ranks under thinning, and scores ~1.
# The criterion is therefore maximised by exactly the failure mode criterion 1
# exists to punish, and this script measures that directly: the row labelled
# DEPTH in the ceiling table below uses each spot's own total counts as a
# pseudo-metric and scores ~0.9998 on all three samples. Plug-in's 0.994-0.996
# is not a sign it is a good estimator; it is a restatement of its 0.96-0.98
# correlation with depth.
#
# Still COMPUTED and reported so the record is intact and the original verdict
# is reproducible from this script. It no longer gates.
A0_MIN_RANK_RHO_RETIRED <- 0.90

# What replaced it is DESCRIPTIVE, with no threshold. A gate invented after
# seeing results would violate the D6 pre-registration discipline; removing a
# structurally invalid one does not, because the defect is demonstrable without
# reference to which estimator wins. See `split_half` below.

#' Run the Stage A4 subsampling ladder on one sample
#'
#' @param target_sample Sample name
#' @param props Retained depth fractions to thin to
#' @return Named list of the tables written, for the cohort rbind
run_subsampling_stability <- function(target_sample, props = c(0.75, 0.50, 0.25)) {

cat("=================================================================\n")
cat("Stage A4 subsampling stability:", target_sample, "\n")
cat("Proportions:", paste(props, collapse = ", "), "\n")
cat("=================================================================\n")

spatial_obj <- load_qc_sample(target_sample)$obj
meta <- spatial_obj@meta.data
counts_mat <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
counts_filtered <- exclude_gene_families(counts_mat, ENTROPY_EXCLUDE,
                                         context = "subsampling stability")
counts_filtered <- Matrix::drop0(counts_filtered)
rm(spatial_obj, counts_mat)
gc(verbose = FALSE)

cat(sprintf("Spots: %d   Features: %d   Stored counts: %s\n",
            ncol(counts_filtered), nrow(counts_filtered),
            format(length(counts_filtered@x), big.mark = ",")))

#' Binomial thinning of a sparse count matrix
#'
#' Zeros produced by thinning are dropped so that K, f1 and f2 are read off a
#' structurally clean matrix -- the same drop0 requirement the estimator kernel
#' has, and the reason a thinned matrix must never be passed on unfiltered.
thin_counts <- function(m, prop) {
  m@x <- as.numeric(rbinom(length(m@x), size = as.integer(m@x), prob = prop))
  Matrix::drop0(m)
}

cat("\nFull-depth estimators...\n")
full_est <- entropy_estimator_matrix(counts_filtered, return_stats = TRUE)
full_sd <- vapply(ENTROPY_ESTIMATORS, function(e) sd(full_est[[e]]), numeric(1))
full_N <- full_est$N

rows <- list()
delta_long <- list()

for (prop in props) {
  cat(sprintf("\nThinning to %.0f%% ...\n", 100 * prop))
  t0 <- Sys.time()
  thinned <- thin_counts(counts_filtered, prop)
  sub_est <- entropy_estimator_matrix(thinned, return_stats = TRUE)
  realized <- sum(sub_est$N) / sum(full_N)
  cat(sprintf("  realized depth fraction %.4f, %.1f s\n",
              realized, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  rm(thinned)
  gc(verbose = FALSE)

  # Ceiling for the rank statistic: a pure depth readout. Thinning every spot by
  # the same proportion barely disturbs the depth ordering, so this is very
  # close to 1 -- which is the whole reason the rank half of criterion 3 is not
  # a test of estimator quality.
  depth_ceiling <- suppressWarnings(cor(sub_est$N, full_N, method = "spearman"))
  cat(sprintf("  rank rho ceiling (nCount itself, a pure depth metric): %.5f\n",
              depth_ceiling))

  for (e in ENTROPY_ESTIMATORS) {
    d <- sub_est[[e]] - full_est[[e]]
    rho <- suppressWarnings(cor(sub_est[[e]], full_est[[e]], method = "spearman"))
    rows[[length(rows) + 1]] <- data.frame(
      Sample = target_sample, Estimator = e, Prop = prop,
      Realized_Depth_Fraction = realized,
      N = length(d),
      Median_Abs_dH = median(abs(d)),
      Median_Signed_dH = median(d),
      Mean_Signed_dH = mean(d),
      SD_dH = sd(d),
      Full_Depth_SD = full_sd[[e]],
      Median_Abs_dH_over_SD = median(abs(d)) / full_sd[[e]],
      Spearman_rho_vs_full = rho,
      Depth_Ceiling_Rank_Rho = depth_ceiling,
      Rank_Rho_Gap_To_Ceiling = depth_ceiling - rho,
      Pearson_r_vs_full = suppressWarnings(cor(sub_est[[e]], full_est[[e]])),
      stringsAsFactors = FALSE
    )
    delta_long[[length(delta_long) + 1]] <- data.frame(
      Estimator = e, Prop = prop, dH = d,
      H_full = full_est[[e]], H_sub = sub_est[[e]],
      stringsAsFactors = FALSE
    )
  }
  rm(sub_est)
  gc(verbose = FALSE)
}

stability <- do.call(rbind, rows)
delta_df <- do.call(rbind, delta_long)

# --- Depth-free reliability (replaces the retired rank half of criterion 3) --
#
# The question the rank half was TRYING to ask is "does this metric rank spots
# reproducibly?". Correlating a thinned replicate against the full-depth values
# cannot answer it, because both carry the same depth signal and depth is what
# thinning preserves.
#
# Two INDEPENDENT thinnings at the same proportion give a split-half replicate
# pair at matched depth. The zero-order correlation between them is still
# inflated by depth for the same reason. Partialling both on the spot's
# FULL-DEPTH total counts removes it, and what survives is the reproducible
# part of the metric that is not a restatement of how deeply the spot was
# sequenced -- which is precisely the quantity the project needs, since D3's
# open question is whether the depth confound was doing all the work.
#
# Two reference rows make the table self-interpreting:
#   DEPTH  -- the spot's own total counts as a pseudo-metric. Split-half ~1,
#             depth-free ~0. All of its reproducibility is depth.
#   NOISE  -- independent normal draws. Both ~0. No reproducibility at all.
# An estimator carrying real, depth-free spatial signal must separate from
# both: high split-half AND clearly non-zero depth-free.
#
# This is reported WITHOUT a pass/fail threshold. Inventing a gate after seeing
# results is what D6 forbids; retiring an invalid one is not.

cat("\n=================================================================\n")
cat(sprintf("Depth-free reliability: two independent %.0f%% thinnings\n",
            100 * A0_REFERENCE_PROP))
cat("=================================================================\n")

rep_a <- entropy_estimator_matrix(thin_counts(counts_filtered, A0_REFERENCE_PROP))
gc(verbose = FALSE)
rep_b <- entropy_estimator_matrix(thin_counts(counts_filtered, A0_REFERENCE_PROP))
gc(verbose = FALSE)

split_rows <- function(label, a, b) {
  df <- data.frame(a = a, b = b, depth = full_N)
  st <- partial_correlation_stats(df, "a", "b", "depth", method = "spearman")
  data.frame(Sample = target_sample, Metric = label, Prop = A0_REFERENCE_PROP,
             Split_Half_Rho = st$Estimate_zero_order,
             Split_Half_Rho_Given_Depth = st$Estimate,
             stringsAsFactors = FALSE)
}

split_half <- do.call(rbind, c(
  list(split_rows("DEPTH (nCount itself)",
                  Matrix::colSums(thin_counts(counts_filtered, A0_REFERENCE_PROP)),
                  Matrix::colSums(thin_counts(counts_filtered, A0_REFERENCE_PROP))),
       split_rows("NOISE (independent normals)", rnorm(nrow(rep_a)), rnorm(nrow(rep_a)))),
  lapply(ENTROPY_ESTIMATORS, function(e) split_rows(e, rep_a[[e]], rep_b[[e]]))))

rm(rep_a, rep_b)
gc(verbose = FALSE)
print(split_half[, c("Metric", "Split_Half_Rho", "Split_Half_Rho_Given_Depth")],
      row.names = FALSE, digits = 4)



out_dir <- file.path("results", "statistical_tests")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(stability, file.path(out_dir, paste0(target_sample, "_subsampling_stability.csv")), row.names = FALSE)
write.csv(split_half, file.path(out_dir, paste0(target_sample, "_subsampling_reliability.csv")), row.names = FALSE)

# --- A0 criterion 3 ---------------------------------------------------------

ref <- stability[abs(stability$Prop - A0_REFERENCE_PROP) < 1e-9, ]
criteria3 <- data.frame(
  Sample = target_sample, Estimator = ref$Estimator, Prop = ref$Prop,
  Median_Abs_dH = ref$Median_Abs_dH, Full_Depth_SD = ref$Full_Depth_SD,
  Median_Abs_dH_over_SD = ref$Median_Abs_dH_over_SD,
  Spearman_rho_vs_full = ref$Spearman_rho_vs_full,
  Depth_Ceiling_Rank_Rho = ref$Depth_Ceiling_Rank_Rho,
  Median_Signed_dH = ref$Median_Signed_dH,
  C3_Magnitude = ref$Median_Abs_dH_over_SD < A0_MAX_DELTA_OVER_SD,
  C3_Rank_Preserved_Retired = ref$Spearman_rho_vs_full > A0_MIN_RANK_RHO_RETIRED,
  stringsAsFactors = FALSE
)
# Criterion 3 is now its magnitude half alone. The rank half is carried for the
# record only; see A0_MIN_RANK_RHO_RETIRED above.
criteria3$C3_Subsampling_Stability <- criteria3$C3_Magnitude
criteria3$C3_Subsampling_Stability_As_Preregistered <-
  criteria3$C3_Magnitude & criteria3$C3_Rank_Preserved_Retired
criteria3 <- merge(criteria3,
                   split_half[, c("Metric", "Split_Half_Rho", "Split_Half_Rho_Given_Depth")],
                   by.x = "Estimator", by.y = "Metric", all.x = TRUE, sort = FALSE)
write.csv(criteria3, file.path(out_dir, paste0(target_sample, "_subsampling_criteria.csv")), row.names = FALSE)

# --- Figure -----------------------------------------------------------------

delta_df$Estimator <- factor(delta_df$Estimator, levels = ENTROPY_ESTIMATORS)
delta_df$Level <- factor(sprintf("%.0f%%", 100 * delta_df$Prop),
                         levels = sprintf("%.0f%%", 100 * sort(props, decreasing = TRUE)))

p_delta <- ggplot(delta_df, aes(x = Level, y = dH, fill = Estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  facet_wrap(~Estimator, nrow = 1) +
  labs(title = sprintf("Depth invariance: H(thinned) - H(full depth) -- %s", target_sample),
       subtitle = "A depth-unbiased estimator stays on the dashed line as depth is reduced",
       x = "Retained depth", y = expression(Delta*H~"(bits)")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(face = "bold", size = 12),
                     strip.text = element_text(face = "bold"))

sub50 <- delta_df[abs(delta_df$Prop - A0_REFERENCE_PROP) < 1e-9, ]
p_scatter <- ggplot(sub50, aes(x = H_full, y = H_sub)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.15, size = 0.5, colour = "#2c3e50") +
  facet_wrap(~Estimator, nrow = 1, scales = "free") +
  labs(title = "Spot-wise agreement at 50% depth",
       subtitle = "Points on the dashed line are unchanged by halving every spot's depth",
       x = "H at full depth (bits)", y = "H at 50% depth (bits)") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12),
                     strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, paste0(target_sample, "_subsampling_stability_plot.png")),
       plot = p_delta / p_scatter, width = 14, height = 9, dpi = 200)

# --- Report -----------------------------------------------------------------

cat("\n=================================================================\n")
cat("Subsampling ladder:\n")
cat("=================================================================\n")
print(stability[, c("Estimator", "Prop", "Median_Abs_dH", "Median_Signed_dH",
                    "Median_Abs_dH_over_SD", "Spearman_rho_vs_full",
                    "Depth_Ceiling_Rank_Rho")],
      row.names = FALSE, digits = 4)

cat("\n=================================================================\n")
cat(sprintf("A0 criterion 3 at %.0f%% depth (|dH| < %.2f x SD):\n",
            100 * A0_REFERENCE_PROP, A0_MAX_DELTA_OVER_SD))
cat("=================================================================\n")
print(criteria3[, c("Estimator", "Median_Abs_dH", "Full_Depth_SD", "Median_Abs_dH_over_SD",
                    "C3_Magnitude", "C3_Subsampling_Stability",
                    "Split_Half_Rho_Given_Depth")], row.names = FALSE, digits = 4)

cat(sprintf("\nRetired rank half, for the record only (rho > %.2f; pure-depth ceiling %.5f):\n",
            A0_MIN_RANK_RHO_RETIRED, criteria3$Depth_Ceiling_Rank_Rho[1]))
print(criteria3[, c("Estimator", "Spearman_rho_vs_full", "C3_Rank_Preserved_Retired",
                    "C3_Subsampling_Stability_As_Preregistered")],
      row.names = FALSE, digits = 4)

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n  %s\n",
            file.path(out_dir, paste0(target_sample, "_subsampling_stability.csv")),
            file.path(out_dir, paste0(target_sample, "_subsampling_reliability.csv")),
            file.path(out_dir, paste0(target_sample, "_subsampling_criteria.csv")),
            file.path(out_dir, paste0(target_sample, "_subsampling_stability_plot.png"))))

invisible(list(
  subsampling_criteria = criteria3,
  subsampling_stability = stability,
  subsampling_reliability = split_half
))
}

# --- Dispatch ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
prop_arg <- grep("^--props=", args, value = TRUE)
props <- if (length(prop_arg) > 0) {
  as.numeric(strsplit(sub("^--props=", "", prop_arg[1]), ",")[[1]])
} else {
  c(0.75, 0.50, 0.25)
}
args <- args[!grepl("^--props=", args)]
sweep_spec <- resolve_sweep_samples(args)

sweep <- run_sweep(sweep_spec$samples,
                   function(s) run_subsampling_stability(s, props = props),
                   label = "A4 subsampling")

if (sweep_spec$all) {
  out_dir <- file.path("results", "statistical_tests")
  cat("\n=================================================================\n")
  cat("Cohort tables:\n")
  cat("=================================================================\n")
  cohort <- write_cohort_tables(
    sweep,
    c("subsampling_criteria", "subsampling_stability", "subsampling_reliability"),
    out_dir)

  if (!is.null(cohort$subsampling_criteria)) {
    crit_summary <- summarise_criteria(
      cohort$subsampling_criteria,
      c("C3_Magnitude", "C3_Subsampling_Stability",
        "C3_Rank_Preserved_Retired", "C3_Subsampling_Stability_As_Preregistered"))
    write.csv(crit_summary,
              file.path(out_dir, "cohort_subsampling_criteria_summary.csv"),
              row.names = FALSE)
    cat(sprintf("  %s\n", file.path(out_dir, "cohort_subsampling_criteria_summary.csv")))

    cat("\n=================================================================\n")
    cat("Criterion 3 across the cohort: passes on k of n\n")
    cat("=================================================================\n")
    print(crit_summary[, c("Estimator", "Criterion", "N_Pass", "N_Samples",
                           "Pct_Pass", "Failing_Samples")],
          row.names = FALSE)

    cat("\nDepth-free split-half reliability (descriptive, no threshold):\n")
    rel <- cohort$subsampling_reliability
    agg <- aggregate(Split_Half_Rho_Given_Depth ~ Metric, data = rel,
                     FUN = function(v) c(median = median(v), min = min(v), max = max(v)))
    agg <- cbind(Metric = agg$Metric, as.data.frame(agg$Split_Half_Rho_Given_Depth))
    print(agg, row.names = FALSE, digits = 4)
  }
}

report_sweep_failures(sweep, length(sweep_spec$samples))
cat("Stage A4 diagnostic completed.\n")
