# R/gene_universe.R
# Helper functions for loading and applying the cohort gene universe.

#' Apply Seurat's feature-name mangling
#'
#' `CreateSeuratObject()`/`CreateAssayObject()` reject underscores in feature names
#' and silently rewrite them to dashes ("Feature names cannot have underscores
#' ('_'), replacing with dashes ('-')"). The gene universe CSV stores the raw
#' `matrix/features/name` strings from the h5, so any match against Seurat rownames
#' has to mangle the universe the same way or those features are silently dropped.
#'
#' @param x Character vector of feature names
#' @return The names with underscores replaced by dashes
sanitize_feature_names <- function(x) {
  gsub("_", "-", x)
}

#' Load the cohort gene universe from CSV
#'
#' @param path Path to gene universe CSV (default: results/cohort_qc/gene_universe.csv)
#' @return Character vector of gene symbols in the universe
load_gene_universe <- function(path = file.path("results", "cohort_qc", "gene_universe.csv")) {
  if (!file.exists(path)) {
    stop(sprintf("Gene universe file '%s' not found. Run build_gene_universe.R first.", path))
  }
  gu_df <- read.csv(path, stringsAsFactors = FALSE)
  if ("is_universe" %in% colnames(gu_df)) {
    gu_df <- gu_df[gu_df$is_universe == TRUE, ]
  }
  gene_col <- if ("gene_name" %in% colnames(gu_df)) "gene_name" else colnames(gu_df)[1]
  return(as.character(gu_df[[gene_col]]))
}

#' Filter a Seurat object or count matrix to the cohort gene universe
#'
#' Matching is done on Seurat-sanitized names. Universe entries that do not appear
#' in the object are reported: with the deprecated probe rows already excluded at
#' build time, the only expected shortfall is duplicated gene symbols, which
#' `make.unique()` renames to `SYMBOL.1` in the object but not in the universe.
#'
#' @param obj Seurat object or dgCMatrix/matrix
#' @param universe Character vector of gene symbols. If NULL, loaded from results/cohort_qc/gene_universe.csv
#' @param assay Assay name if obj is a Seurat object (default: "Spatial")
#' @param verbose Logical, report the feature accounting (default: TRUE)
#' @return Filtered object with features restricted to the universe
filter_by_gene_universe <- function(obj, universe = NULL, assay = "Spatial", verbose = TRUE) {
  if (is.null(universe)) {
    universe <- load_gene_universe()
  }
  universe <- unique(sanitize_feature_names(universe))

  if (inherits(obj, "Seurat")) {
    current_features <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  } else if (inherits(obj, c("Matrix", "matrix", "dgCMatrix"))) {
    current_features <- rownames(obj)
  } else {
    stop("Unsupported object type for filter_by_gene_universe. Must be Seurat or matrix/dgCMatrix.")
  }

  valid_features <- intersect(current_features, universe)
  unmatched <- setdiff(universe, current_features)

  if (verbose) {
    cat(sprintf("Filtering from %d to %d features in gene universe (%d dropped as outside the universe).\n",
                length(current_features), length(valid_features), length(current_features) - length(valid_features)))
    if (length(unmatched) > 0) {
      cat(sprintf("  %d universe features absent from this object (duplicated symbols renamed by make.unique): %s\n",
                  length(unmatched), paste(utils::head(unmatched, 5), collapse = ", ")))
    }
  }

  if (inherits(obj, "Seurat")) {
    return(subset(obj, features = valid_features))
  }
  obj[valid_features, , drop = FALSE]
}
