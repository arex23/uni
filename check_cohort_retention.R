#!/usr/bin/env Rscript
# Cohort retention gate (D6), run after run_cohort.sh.
#
# The spot QC gates drop each sample's shallow tail, and that tail is not the same
# size everywhere. Dropped spots are spatially clustered rather than scattered, so
# each sample enters downstream spatial inference with a differently-holed hex
# lattice, and a sample that lost a third of its tissue is not comparable to one
# that lost none.
#
# The threshold (COHORT_MIN_RETENTION_PCT in R/cohort.R) is frozen before the
# cohort run for the same reason the %MT threshold is: a retention floor chosen
# after seeing entropy-stemness correlations is a selection-bias generator, which
# is exactly what D6 rejected alternative 1 forbids.
#
# This script only reports. It does not edit R/cohort.R -- promoting a sample out
# of the cohort is a decision to record in D6, by hand, with the number attached.

source("R/cohort.R")

qc_dir <- file.path("results", "analyze_entropy")
out_file <- file.path("results", "cohort_qc", "retention_gate.csv")

cohort <- cohort_samples()

rows <- list()
missing <- character(0)

for (s in cohort) {
  f <- file.path(qc_dir, paste0(s, "_qc_metrics.csv"))
  if (!file.exists(f)) {
    missing <- c(missing, s)
    next
  }
  qc <- read.csv(f, stringsAsFactors = FALSE)
  rows[[s]] <- data.frame(
    Sample = s,
    Spots_On_Tissue = qc$Spots_On_Tissue[1],
    Spots_Post_Depth_QC = qc$Spots_Post_Depth_QC[1],
    Spots_Final = qc$Spots_Final[1],
    Pct_OnTissue_Retained = qc$Pct_OnTissue_Retained[1],
    Mean_Percent_MT = qc$Mean_Percent_MT[1],
    Mean_Raw_Plugin_Entropy = qc$Mean_Raw_Plugin_Entropy[1],
    SD_Raw_Plugin_Entropy = qc$SD_Raw_Plugin_Entropy[1],
    stringsAsFactors = FALSE
  )
}

if (length(missing) > 0) {
  cat(sprintf("WARNING: no QC table for %d cohort sample(s): %s\n",
              length(missing), paste(missing, collapse = ", ")))
  cat("Run ./run_cohort.sh entropy first; the gate below covers only what exists.\n\n")
}

if (length(rows) == 0) {
  stop("No per-sample QC tables found in results/analyze_entropy/. Nothing to gate.")
}

df <- do.call(rbind, rows)
df$Passes_Retention_Gate <- df$Pct_OnTissue_Retained >= COHORT_MIN_RETENTION_PCT
df <- df[order(df$Pct_OnTissue_Retained), ]
rownames(df) <- NULL

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.csv(df, out_file, row.names = FALSE)

cat("=================================================================\n")
cat(sprintf("Cohort retention gate: Pct_OnTissue_Retained >= %g%%\n", COHORT_MIN_RETENTION_PCT))
cat("=================================================================\n")
print(df[, c("Sample", "Spots_On_Tissue", "Spots_Final",
             "Pct_OnTissue_Retained", "Passes_Retention_Gate")], row.names = FALSE)

failing <- df$Sample[!df$Passes_Retention_Gate]
cat("\n")
if (length(failing) > 0) {
  cat(sprintf("BELOW THRESHOLD (%d): %s\n", length(failing), paste(failing, collapse = ", ")))
  cat("Record the decision in D6 before running any cross-sample inference.\n")
} else if (length(missing) > 0) {
  cat(sprintf("All %d sample(s) with a QC table pass, but %d of %d are still missing.\n",
              nrow(df), length(missing), length(cohort)))
  cat("The gate is not satisfied until the whole cohort has run.\n")
} else {
  cat(sprintf("All %d cohort samples pass. Comparable spot coverage; Stage 4 can proceed.\n", nrow(df)))
}

# Cross-sample spread of the entropy baseline. A large spread across samples is a
# batch effect to model in Stage 4, not a reason to re-gate anything here. Note
# that the plug-in metric is depth-coupled by construction, so part of this spread
# just tracks the samples' differing sequencing depths.
cat(sprintf("\nMean plug-in entropy across cohort: %.4f - %.4f (spread %.4f bits)\n",
            min(df$Mean_Raw_Plugin_Entropy, na.rm = TRUE),
            max(df$Mean_Raw_Plugin_Entropy, na.rm = TRUE),
            diff(range(df$Mean_Raw_Plugin_Entropy, na.rm = TRUE))))
cat(sprintf("Written: %s\n", out_file))
