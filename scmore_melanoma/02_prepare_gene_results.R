#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

trim_or_empty <- function(x) {
  trimws(Sys.getenv(x, unset = ""))
}

require_env <- function(name) {
  value <- trim_or_empty(name)
  if (!nzchar(value)) {
    stop(sprintf("Required environment variable %s is not set.", name), call. = FALSE)
  }
  value
}

get_optional_col <- function(dt, env_name) {
  col_name <- trim_or_empty(env_name)
  if (!nzchar(col_name)) {
    return(NULL)
  }
  if (!col_name %in% names(dt)) {
    stop(sprintf("Column %s from %s was not found in the gene-results file.", col_name, env_name), call. = FALSE)
  }
  col_name
}

get_required_col <- function(dt, env_name) {
  col_name <- require_env(env_name)
  if (!col_name %in% names(dt)) {
    stop(sprintf("Column %s from %s was not found in the gene-results file.", col_name, env_name), call. = FALSE)
  }
  col_name
}

results_root <- trim_or_empty("SCMORE_RESULTS_ROOT")
if (!nzchar(results_root)) {
  results_root <- file.path(trim_or_empty("SCMORE_PROJECT_ROOT"), "scmore_melanoma", "results")
}
trait_name <- require_env("SCMORE_TRAIT_NAME")
input_file <- require_env("SCMORE_RAW_GENE_RESULTS")

dt <- fread(input_file)

gene_col <- get_required_col(dt, "SCMORE_GENE_ID_COL")
chr_col <- get_required_col(dt, "SCMORE_GENE_CHR_COL")
start_col <- get_required_col(dt, "SCMORE_GENE_START_COL")
stop_col <- get_required_col(dt, "SCMORE_GENE_STOP_COL")
nsnps_col <- get_required_col(dt, "SCMORE_GENE_NSNPS_COL")
nparam_col <- get_required_col(dt, "SCMORE_GENE_NPARAM_COL")
n_col <- get_required_col(dt, "SCMORE_GENE_N_COL")
zstat_col <- get_required_col(dt, "SCMORE_GENE_ZSTAT_COL")
p_col <- get_required_col(dt, "SCMORE_GENE_P_COL")
symbol_col <- get_optional_col(dt, "SCMORE_GENE_SYMBOL_COL")

out <- data.table(
  GENE = as.character(dt[[gene_col]]),
  CHR = as.character(dt[[chr_col]]),
  START = as.integer(dt[[start_col]]),
  STOP = as.integer(dt[[stop_col]]),
  NSNPS = as.numeric(dt[[nsnps_col]]),
  NPARAM = as.numeric(dt[[nparam_col]]),
  N = as.numeric(dt[[n_col]]),
  ZSTAT = as.numeric(dt[[zstat_col]]),
  P = as.numeric(dt[[p_col]])
)

if (!is.null(symbol_col)) {
  out[, SYMBOL := as.character(dt[[symbol_col]])]
}

out <- out[!is.na(GENE) & !is.na(CHR) & !is.na(START) & !is.na(STOP) & !is.na(P)]
out <- unique(out, by = "GENE")

dir.create(file.path(results_root, "inputs"), recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(results_root, "inputs", sprintf("%s_scmore_gene_results.tsv.gz", trait_name))
fwrite(out, output_file, sep = "\t")

message(sprintf("Wrote normalized scMORE gene results to %s", output_file))
message(sprintf("Rows retained: %s", nrow(out)))
