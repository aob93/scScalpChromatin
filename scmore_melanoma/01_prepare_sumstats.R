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
    stop(sprintf("Column %s from %s was not found in the summary-statistics file.", col_name, env_name), call. = FALSE)
  }
  col_name
}

get_required_col <- function(dt, env_name) {
  col_name <- require_env(env_name)
  if (!col_name %in% names(dt)) {
    stop(sprintf("Column %s from %s was not found in the summary-statistics file.", col_name, env_name), call. = FALSE)
  }
  col_name
}

results_root <- trim_or_empty("SCMORE_RESULTS_ROOT")
if (!nzchar(results_root)) {
  results_root <- file.path(trim_or_empty("SCMORE_PROJECT_ROOT"), "scmore_melanoma", "results")
}
trait_name <- require_env("SCMORE_TRAIT_NAME")
input_file <- require_env("SCMORE_RAW_SUMSTATS")

dt <- fread(input_file)

chr_col <- get_required_col(dt, "SCMORE_SUMSTAT_CHR_COL")
pos_col <- get_required_col(dt, "SCMORE_SUMSTAT_POS_COL")
beta_col <- get_optional_col(dt, "SCMORE_SUMSTAT_BETA_COL")
or_col <- get_optional_col(dt, "SCMORE_SUMSTAT_OR_COL")
se_col <- get_required_col(dt, "SCMORE_SUMSTAT_SE_COL")
p_col <- get_optional_col(dt, "SCMORE_SUMSTAT_P_COL")
lp_col <- get_optional_col(dt, "SCMORE_SUMSTAT_LP_COL")
af_col <- get_required_col(dt, "SCMORE_SUMSTAT_AF_COL")
n_col <- get_required_col(dt, "SCMORE_SUMSTAT_N_COL")
snp_col <- get_required_col(dt, "SCMORE_SUMSTAT_SNP_COL")

if (is.null(beta_col) && is.null(or_col)) {
  stop("One of SCMORE_SUMSTAT_BETA_COL or SCMORE_SUMSTAT_OR_COL must be set.", call. = FALSE)
}
if (is.null(p_col) && is.null(lp_col)) {
  stop("One of SCMORE_SUMSTAT_P_COL or SCMORE_SUMSTAT_LP_COL must be set.", call. = FALSE)
}

es <- if (!is.null(beta_col)) {
  as.numeric(dt[[beta_col]])
} else {
  log(as.numeric(dt[[or_col]]))
}

lp <- if (!is.null(lp_col)) {
  as.numeric(dt[[lp_col]])
} else {
  pvals <- as.numeric(dt[[p_col]])
  pvals[is.na(pvals) | pvals <= 0] <- .Machine$double.xmin
  -log10(pvals)
}

out <- data.table(
  CHR = as.character(dt[[chr_col]]),
  POS = as.integer(dt[[pos_col]]),
  ES = es,
  SE = as.numeric(dt[[se_col]]),
  LP = lp,
  AF = as.numeric(dt[[af_col]]),
  SZ = as.numeric(dt[[n_col]]),
  SNP = as.character(dt[[snp_col]])
)

out <- out[!is.na(CHR) & !is.na(POS) & !is.na(ES) & !is.na(SE) & !is.na(LP) & !is.na(AF) & !is.na(SZ) & !is.na(SNP)]
out <- unique(out, by = "SNP")

dir.create(file.path(results_root, "inputs"), recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(results_root, "inputs", sprintf("%s_scmore_sumstats.tsv.gz", trait_name))
fwrite(out, output_file, sep = "\t")

message(sprintf("Wrote normalized scMORE summary statistics to %s", output_file))
message(sprintf("Rows retained: %s", nrow(out)))
