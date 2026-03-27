#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(GenomicRanges)
  library(SummarizedExperiment)
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

split_csv <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) {
    return(character())
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

get_assay_matrix <- function(se) {
  assay_names <- SummarizedExperiment::assayNames(se)
  if (!length(assay_names)) {
    stop("No assays were found in the ArchR matrix object.", call. = FALSE)
  }
  SummarizedExperiment::assay(se, assay_names[[1]])
}

archr_dir <- require_env("SCMORE_ARCHR_PROJECT_DIR")
results_root <- trim_or_empty("SCMORE_RESULTS_ROOT")
if (!nzchar(results_root)) {
  results_root <- file.path(trim_or_empty("SCMORE_PROJECT_ROOT"), "scmore_melanoma", "results")
}
multiome_rds <- trim_or_empty("SCMORE_MULTIOME_RDS")
if (!nzchar(multiome_rds)) {
  multiome_rds <- file.path(results_root, "inputs", sprintf("%s_controls_multimodal.rds", require_env("SCMORE_TRAIT_NAME")))
}

celltype_col <- require_env("SCMORE_CELLTYPE_COLUMN")
keep_status <- split_csv(trim_or_empty("SCMORE_KEEP_DISEASE_STATUS"))
keep_samples <- split_csv(trim_or_empty("SCMORE_KEEP_SAMPLES"))

addArchRThreads(threads = 8)
proj <- loadArchRProject(archr_dir, force = TRUE)

meta <- as.data.frame(proj@cellColData)
meta$cell_id <- rownames(meta)

if (length(keep_status)) {
  if (!"diseaseStatus" %in% colnames(meta)) {
    stop("ArchR project is missing diseaseStatus metadata required for SCMORE_KEEP_DISEASE_STATUS filtering.", call. = FALSE)
  }
  meta <- meta[meta$diseaseStatus %in% keep_status, , drop = FALSE]
}

if (length(keep_samples)) {
  sample_col <- if ("Sample2" %in% colnames(meta)) "Sample2" else if ("Sample" %in% colnames(meta)) "Sample" else NULL
  if (is.null(sample_col)) {
    stop("ArchR project is missing Sample/Sample2 metadata required for SCMORE_KEEP_SAMPLES filtering.", call. = FALSE)
  }
  meta <- meta[meta[[sample_col]] %in% keep_samples, , drop = FALSE]
}

if (!nrow(meta)) {
  stop("No cells remain after applying the requested scMORE filters.", call. = FALSE)
}

if (!celltype_col %in% colnames(meta)) {
  stop(sprintf("Cell-type column %s is not present in the ArchR project metadata.", celltype_col), call. = FALSE)
}

keep_cells <- meta$cell_id

peak_se <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
peak_counts <- get_assay_matrix(peak_se)[, keep_cells, drop = FALSE]
peak_set <- getPeakSet(proj)
peak_names <- rownames(peak_counts)
if (is.null(names(peak_set)) || !length(names(peak_set))) {
  names(peak_set) <- peak_names[seq_len(length(peak_set))]
}
peak_idx <- match(peak_names, names(peak_set))
if (any(is.na(peak_idx))) {
  stop("Could not align PeakMatrix rows to the ArchR peak set.", call. = FALSE)
}
peak_ranges <- peak_set[peak_idx]
names(peak_ranges) <- peak_names

gene_se <- getMatrixFromProject(proj, useMatrix = "GeneIntegrationMatrix")
gene_counts <- get_assay_matrix(gene_se)[, keep_cells, drop = FALSE]
gene_names <- SummarizedExperiment::rowData(gene_se)$name
if (is.null(gene_names) || !length(gene_names)) {
  gene_names <- rownames(gene_counts)
}
rownames(gene_counts) <- make.unique(as.character(gene_names))

meta <- meta[match(keep_cells, meta$cell_id), , drop = FALSE]
rownames(meta) <- meta$cell_id

seu <- CreateSeuratObject(
  counts = gene_counts,
  assay = "RNA",
  meta.data = meta
)

seu[["peaks"]] <- CreateChromatinAssay(
  counts = peak_counts,
  ranges = peak_ranges
)

DefaultAssay(seu) <- "RNA"
Idents(seu) <- seu[[celltype_col, drop = TRUE]]
seu$cell_type <- Idents(seu)

dir.create(dirname(multiome_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(seu, multiome_rds)

message(sprintf("Wrote scMORE multimodal Seurat object to %s", multiome_rds))
message(sprintf("Cells retained: %s", ncol(seu)))
message(sprintf("Genes retained: %s", nrow(seu[["RNA"]])))
message(sprintf("Peaks retained: %s", nrow(seu[["peaks"]])))
