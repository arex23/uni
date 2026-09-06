#' Cohort sweep driver for the off-pipeline diagnostics
#'
#' The Stage A diagnostics (`diagnose_entropy_scaling.R`,
#' `diagnose_subsampling_stability.R`, `diagnose_correction_routes.R`) all have
#' the same shape: a per-sample function returning a named list of data frames,
#' run either on one sample or on the frozen cohort. That driver lives here once
#' rather than three times.
#'
#' Requires `R/cohort.R` to be sourced first (`cohort_samples()`,
#' `resolve_target_sample()`).

#' Resolve the sample set a diagnostic should run on
#'
#' `--all` runs the frozen 16-sample analysis cohort (D6) via `cohort_samples()`
#' -- never `list.dirs("data")`, which would silently pull in the six excluded
#' samples. Without it, the usual single-sample resolution applies.
#'
#' The unconsumed arguments are returned so callers can keep parsing their own
#' flags (`--from-rds`, `--props=`) without having to know about `--all`.
#'
#' @param args Character vector of command-line arguments
#' @param data_dir Directory holding one subdirectory per sample
#' @return List with `samples`, the remaining `args`, and the `all` flag
resolve_sweep_samples <- function(args, data_dir = "data") {
  all_flag <- "--all" %in% args
  args <- args[args != "--all"]
  samples <- if (all_flag) {
    cohort_samples(data_dir)
  } else {
    resolve_target_sample(args, data_dir = data_dir)
  }
  list(samples = samples, args = args, all = all_flag)
}

#' Run a per-sample diagnostic across a sample set
#'
#' A sample that errors is recorded and the sweep continues. An unattended
#' 16-sample run is 30-45 minutes; losing all of it because one sample tripped
#' an assertion is worse than losing one sample and being told which. The
#' failures are returned and printed, never swallowed -- `check_entropy_ordering()`
#' aborting would be a real bug worth surfacing loudly.
#'
#' The RNG is reseeded before each sample so a sample's draws depend on the
#' sample rather than on its position in the loop. Without this, running
#' `--all` and running one sample alone would give different numbers for the
#' same sample, which would make the subsampling ladder irreproducible.
#'
#' @param samples Character vector of sample names
#' @param fn Function of one sample name, returning a named list of data frames
#' @param label Short name of the diagnostic, used in progress and error lines
#' @param seed Integer reseeded before each sample
#' @return List with `results` (named by sample) and `failures` (named by sample)
run_sweep <- function(samples, fn, label, seed = 23) {
  results <- list()
  failures <- character(0)

  for (i in seq_along(samples)) {
    s <- samples[i]
    cat(sprintf("\n\n#################################################################\n"))
    cat(sprintf("# %s  [%d/%d]  %s\n", label, i, length(samples), s))
    cat(sprintf("#################################################################\n"))
    t0 <- Sys.time()
    set.seed(seed)

    res <- tryCatch(fn(s), error = function(e) {
      cat(sprintf("\n*** %s FAILED on %s: %s\n", label, s, conditionMessage(e)))
      structure(conditionMessage(e), class = "sweep_failure")
    })

    if (inherits(res, "sweep_failure")) {
      failures[[s]] <- as.character(res)
    } else {
      results[[s]] <- res
      cat(sprintf("\n[%s done in %.1f s]\n", s,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    gc(verbose = FALSE)
  }

  list(results = results, failures = failures)
}

#' Bind one named element across a sweep's per-sample results
#'
#' @param results The `results` list from `run_sweep()`
#' @param element Name of the data frame to extract from each sample
#' @return A single data frame, or NULL when no sample produced that element
bind_cohort <- function(results, element) {
  parts <- lapply(results, function(r) r[[element]])
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) return(NULL)
  do.call(rbind, parts)
}

#' Write the cohort-level rbind of every element a sweep produced
#'
#' Per-sample CSVs are written by the diagnostics themselves and keep their
#' names; this adds one `cohort_<element>.csv` per element, written only when
#' the sweep actually ran on more than one sample.
#'
#' @param sweep The list returned by `run_sweep()`
#' @param elements Character vector of element names to bind
#' @param out_dir Output directory
#' @param prefix File name prefix (default "cohort")
#' @return Named list of the written data frames
write_cohort_tables <- function(sweep, elements, out_dir, prefix = "cohort") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  written <- list()
  for (el in elements) {
    df <- bind_cohort(sweep$results, el)
    if (is.null(df)) next
    path <- file.path(out_dir, sprintf("%s_%s.csv", prefix, el))
    write.csv(df, path, row.names = FALSE)
    cat(sprintf("  %s  (%d rows)\n", path, nrow(df)))
    written[[el]] <- df
  }
  written
}

#' Report which samples failed, and stop pretending the sweep was complete
#'
#' @param sweep The list returned by `run_sweep()`
#' @param n_requested Number of samples the sweep was asked to run
report_sweep_failures <- function(sweep, n_requested) {
  n_ok <- length(sweep$results)
  cat(sprintf("\nSweep completed on %d/%d samples.\n", n_ok, n_requested))
  if (length(sweep$failures) > 0) {
    cat("FAILED samples (their rows are absent from every cohort table):\n")
    for (s in names(sweep$failures)) {
      cat(sprintf("  %-10s %s\n", s, sweep$failures[[s]]))
    }
  }
  invisible(NULL)
}

#' Count criterion passes across the cohort
#'
#' A0 criteria are per-sample booleans; the cohort-level evidence is "passes on
#' k of n", not a single collapsed verdict. Deliberately returns the counts and
#' the sample lists rather than a pass/fail -- k/n is evidence to read, not a
#' new threshold to tune.
#'
#' @param criteria_df Cohort criteria table with a `Sample`, `Estimator` and
#'   boolean criterion columns
#' @param criterion_cols Names of the boolean columns to summarise
#' @return Data frame with one row per (estimator x criterion)
summarise_criteria <- function(criteria_df, criterion_cols) {
  ests <- unique(criteria_df$Estimator)
  rows <- list()
  for (e in ests) {
    sub <- criteria_df[criteria_df$Estimator == e, ]
    for (cc in criterion_cols) {
      v <- sub[[cc]]
      rows[[length(rows) + 1]] <- data.frame(
        Estimator = e, Criterion = cc,
        N_Samples = sum(!is.na(v)),
        N_Pass = sum(v, na.rm = TRUE),
        Pct_Pass = round(100 * mean(v, na.rm = TRUE), 1),
        Failing_Samples = paste(sub$Sample[!is.na(v) & !v], collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
