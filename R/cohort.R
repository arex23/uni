#' Frozen cohort definition (D6)
#'
#' The 22 sequenced samples are enumerated from `data/`, but only 16 enter
#' downstream biological analysis. The 6 exclusions are frozen here, upfront,
#' on an orthogonal technical criterion (mean spot %MT >= 10% measured in
#' `check_cohort_chemistry.R`, Stage 0) and must never be revisited after
#' looking at entropy or stemness results -- see D6, rejected alternative 1.
#'
#' Scripts that build cohort-wide reference data (`check_cohort_chemistry.R`,
#' `build_gene_universe.R`) deliberately run over all 22: the chemistry check is
#' what produces the exclusion criterion, and the gene universe is defined as a
#' permissive union over the full slide set (D5). Every script that produces a
#' *biological* result runs over `cohort_samples()` only.

# Excluded for tissue degradation / mitochondrial contamination (D6).
# Values are mean spot %MT on post-QC on-tissue spots, from
# results/cohort_qc/chemistry_check.csv.
COHORT_EXCLUDED <- c(
  sample2  = 11.50,
  sample3  = 11.52,
  sample12 = 10.54,
  sample16 = 14.74,
  sample19 = 10.80,
  sample20 = 16.84
)

# Frozen quality threshold behind COHORT_EXCLUDED, in percent.
COHORT_MT_THRESHOLD <- 10

# Minimum percentage of on-tissue spots a sample must retain through the QC
# gates to stay in the cohort. Frozen before the cohort run so the decision
# cannot be made after seeing correlations. Applied by
# `check_cohort_retention.R` against the per-sample QC tables.
#
# Under the nCount >= 500 / nFeature >= 250 floors this is a guard rather than a
# live discriminator: every cohort sample is expected to retain close to 100% of
# its on-tissue spots. It earns its keep if a future estimator reintroduces a
# depth-dependent gate, or if a sample turns out to be far shallower than Stage 0
# predicted -- a sample that has lost a third of its tissue area is not comparable
# to one that lost none.
COHORT_MIN_RETENTION_PCT <- 60

#' All samples present in data/, numerically ordered
#'
#' @param data_dir Directory holding one subdirectory per sample (default "data")
#' @return Character vector of sample names
all_samples <- function(data_dir = "data") {
  s <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
  s <- s[grepl("^sample[0-9]+$", s)]
  if (length(s) == 0) {
    stop(sprintf("No sample directories found in '%s/'.", data_dir))
  }
  s[order(as.numeric(gsub("sample", "", s)))]
}

#' The analysis cohort: all samples minus the frozen D6 exclusions
#'
#' @param data_dir Directory holding one subdirectory per sample (default "data")
#' @return Character vector of sample names in the analysis cohort
cohort_samples <- function(data_dir = "data") {
  s <- all_samples(data_dir)
  missing_excl <- setdiff(names(COHORT_EXCLUDED), s)
  if (length(missing_excl) > 0) {
    stop(sprintf("D6 excludes samples not present in '%s/': %s. The frozen cohort no longer matches the data.",
                 data_dir, paste(missing_excl, collapse = ", ")))
  }
  setdiff(s, names(COHORT_EXCLUDED))
}

#' Resolve the sample a script should run on
#'
#' Takes the first command-line argument if given, otherwise the first cohort
#' sample. Refuses samples excluded by D6 unless `allow_excluded = TRUE`.
#'
#' @param args Character vector of command-line arguments
#' @param allow_excluded Logical, permit running on a D6-excluded sample
#' @param data_dir Directory holding one subdirectory per sample (default "data")
#' @return A single sample name
resolve_target_sample <- function(args, allow_excluded = FALSE, data_dir = "data") {
  cohort <- cohort_samples(data_dir)
  if (length(args) == 0) return(cohort[1])

  target <- args[1]
  if (target %in% names(COHORT_EXCLUDED)) {
    if (!allow_excluded) {
      stop(sprintf(paste0("Sample '%s' is excluded from the analysis cohort (D6): mean spot %%MT = %.2f%% ",
                          ">= %g%%. Re-run with allow_excluded = TRUE only for diagnostics, never for results."),
                   target, COHORT_EXCLUDED[[target]], COHORT_MT_THRESHOLD))
    }
    cat(sprintf("WARNING: '%s' is a D6-excluded sample (mean spot %%MT = %.2f%%). Diagnostic run only.\n",
                target, COHORT_EXCLUDED[[target]]))
    return(target)
  }
  if (!target %in% cohort) {
    stop(sprintf("Sample '%s' not found in %s/. Cohort samples: %s",
                 target, data_dir, paste(cohort, collapse = ", ")))
  }
  target
}
