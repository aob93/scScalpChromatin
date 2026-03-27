#!/usr/bin/env Rscript

trim_or_empty <- function(x) {
  val <- Sys.getenv(x, unset = "")
  trimws(val)
}

split_csv <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) {
    return(character())
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

require_env <- function(name) {
  value <- trim_or_empty(name)
  if (!nzchar(value)) {
    stop(sprintf("Required environment variable %s is not set.", name), call. = FALSE)
  }
  value
}

message(sprintf("SCMORE_PROJECT_ROOT: %s", trim_or_empty("SCMORE_PROJECT_ROOT")))
message(sprintf("SCMORE_RESULTS_ROOT: %s", trim_or_empty("SCMORE_RESULTS_ROOT")))
message(sprintf("SCMORE_TRAIT_NAME: %s", trim_or_empty("SCMORE_TRAIT_NAME")))

raw_sumstats <- require_env("SCMORE_RAW_SUMSTATS")
if (!file.exists(raw_sumstats)) {
  stop(sprintf("Raw GWAS summary-statistics file not found: %s", raw_sumstats), call. = FALSE)
}
message(sprintf("Found raw GWAS summary statistics: %s", raw_sumstats))

raw_gene_results <- require_env("SCMORE_RAW_GENE_RESULTS")
if (!file.exists(raw_gene_results)) {
  stop(sprintf("Raw gene-level results file not found: %s", raw_gene_results), call. = FALSE)
}
message(sprintf("Found raw gene-level results: %s", raw_gene_results))

archr_dir <- trim_or_empty("SCMORE_ARCHR_PROJECT_DIR")
if (nzchar(archr_dir)) {
  if (!dir.exists(archr_dir)) {
    stop(sprintf("SCMORE_ARCHR_PROJECT_DIR does not exist: %s", archr_dir), call. = FALSE)
  }
  message(sprintf("Found ArchR project directory: %s", archr_dir))
} else {
  message("SCMORE_ARCHR_PROJECT_DIR is not set yet. That is fine if you are only preparing GWAS / MAGMA inputs now.")
}

keep_status <- split_csv(trim_or_empty("SCMORE_KEEP_DISEASE_STATUS"))
if (length(keep_status)) {
  message(sprintf("Disease-status filter: %s", paste(keep_status, collapse = ", ")))
}

message("Validation complete.")
