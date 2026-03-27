#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scMORE)
  library(Seurat)
  library(Signac)
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

results_root <- trim_or_empty("SCMORE_RESULTS_ROOT")
if (!nzchar(results_root)) {
  results_root <- file.path(trim_or_empty("SCMORE_PROJECT_ROOT"), "scmore_melanoma", "results")
}
trait_name <- require_env("SCMORE_TRAIT_NAME")
multiome_rds <- trim_or_empty("SCMORE_MULTIOME_RDS")
if (!nzchar(multiome_rds)) {
  multiome_rds <- file.path(results_root, "inputs", sprintf("%s_controls_multimodal.rds", trait_name))
}
sumstats_file <- file.path(results_root, "inputs", sprintf("%s_scmore_sumstats.tsv.gz", trait_name))
gene_file <- file.path(results_root, "inputs", sprintf("%s_scmore_gene_results.tsv.gz", trait_name))

if (!file.exists(multiome_rds)) {
  stop(sprintf("Multimodal Seurat object not found: %s", multiome_rds), call. = FALSE)
}
if (!file.exists(sumstats_file)) {
  stop(sprintf("Normalized scMORE summary statistics not found: %s", sumstats_file), call. = FALSE)
}
if (!file.exists(gene_file)) {
  stop(sprintf("Normalized scMORE gene results not found: %s", gene_file), call. = FALSE)
}

single_cell <- readRDS(multiome_rds)
if (!"cell_type" %in% colnames(single_cell[[]])) {
  stop("The multimodal Seurat object is missing a cell_type column.", call. = FALSE)
}
Idents(single_cell) <- single_cell$cell_type

snp_info <- fread(sumstats_file)
gene_info <- fread(gene_file)

dir.create(file.path(results_root, "outputs"), recursive = TRUE, showWarnings = FALSE)
output_rds <- file.path(results_root, "outputs", sprintf("%s_scmore_result.rds", trait_name))

result <- scMore(
  single_cell = single_cell,
  snp_info = snp_info,
  gene_info = gene_info,
  n_targets = as.integer(require_env("SCMORE_N_TARGETS")),
  perm_n = as.integer(require_env("SCMORE_PERM_N")),
  theta = as.numeric(require_env("SCMORE_THETA")),
  alpha = as.numeric(require_env("SCMORE_ALPHA")),
  buffer = as.numeric(require_env("SCMORE_BUFFER")),
  top_n = as.integer(require_env("SCMORE_TOP_N")),
  p1 = as.numeric(require_env("SCMORE_P1")),
  p2 = as.numeric(require_env("SCMORE_P2")),
  p3 = as.numeric(require_env("SCMORE_P3")),
  peak2gene_method = require_env("SCMORE_PEAK2GENE_METHOD"),
  infer_method = require_env("SCMORE_INFER_METHOD"),
  method = require_env("SCMORE_SPECIFICITY_METHOD"),
  nSeed = as.integer(require_env("SCMORE_SEED"))
)

saveRDS(result, output_rds)
message(sprintf("Wrote scMORE result to %s", output_rds))
