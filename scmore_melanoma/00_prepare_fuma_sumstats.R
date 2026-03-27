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

normalize_chr <- function(x) {
  chr <- toupper(as.character(x))
  chr <- sub("^CHR", "", chr)
  chr
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
snp_col <- get_optional_col(dt, "SCMORE_SUMSTAT_SNP_COL")
a1_col <- get_optional_col(dt, "SCMORE_SUMSTAT_A1_COL")
a2_col <- get_optional_col(dt, "SCMORE_SUMSTAT_A2_COL")
beta_col <- get_optional_col(dt, "SCMORE_SUMSTAT_BETA_COL")
or_col <- get_optional_col(dt, "SCMORE_SUMSTAT_OR_COL")
se_col <- get_optional_col(dt, "SCMORE_SUMSTAT_SE_COL")
p_col <- get_optional_col(dt, "SCMORE_SUMSTAT_P_COL")
lp_col <- get_optional_col(dt, "SCMORE_SUMSTAT_LP_COL")
n_col <- get_optional_col(dt, "SCMORE_SUMSTAT_N_COL")

if (is.null(p_col) && is.null(lp_col)) {
  stop("One of SCMORE_SUMSTAT_P_COL or SCMORE_SUMSTAT_LP_COL must be set for FUMA formatting.", call. = FALSE)
}

pvals <- if (!is.null(p_col)) {
  as.numeric(dt[[p_col]])
} else {
  10^(-as.numeric(dt[[lp_col]]))
}
pvals[is.na(pvals) | pvals <= 0] <- .Machine$double.xmin

out <- data.table(
  CHR = normalize_chr(dt[[chr_col]]),
  BP = as.integer(dt[[pos_col]]),
  P = pvals
)

if (!is.null(snp_col)) {
  out[, SNP := as.character(dt[[snp_col]])]
}
if (!is.null(a1_col)) {
  out[, A1 := as.character(dt[[a1_col]])]
}
if (!is.null(a2_col)) {
  out[, A2 := as.character(dt[[a2_col]])]
}
if (!is.null(beta_col)) {
  out[, BETA := as.numeric(dt[[beta_col]])]
}
if (!is.null(or_col)) {
  out[, OR := as.numeric(dt[[or_col]])]
}
if (!is.null(se_col)) {
  out[, SE := as.numeric(dt[[se_col]])]
}
if (!is.null(n_col)) {
  out[, N := as.numeric(dt[[n_col]])]
}

preferred_order <- c("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "OR", "SE", "N")
existing_order <- preferred_order[preferred_order %in% names(out)]
setcolorder(out, existing_order)

required_core <- c("CHR", "BP", "P")
out <- out[complete.cases(out[, ..required_core])]

if ("SNP" %in% names(out)) {
  out <- unique(out, by = "SNP")
} else {
  out <- unique(out, by = c("CHR", "BP"))
}

dir.create(file.path(results_root, "inputs"), recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(results_root, "inputs", sprintf("%s_fuma_sumstats.tsv.gz", trait_name))
fwrite(out, output_file, sep = "\t", quote = FALSE)

message(sprintf("Wrote FUMA-formatted summary statistics to %s", output_file))
message(sprintf("Rows retained: %s", nrow(out)))
message(sprintf("Genome build to choose on FUMA upload: %s", trim_or_empty("SCMORE_SUMSTAT_BUILD")))
