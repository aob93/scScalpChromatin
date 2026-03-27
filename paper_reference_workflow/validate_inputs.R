#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)), "..", "pipeline_config.R"))

stage <- commandArgs(trailingOnly = TRUE)
stage <- if (length(stage)) stage[[1]] else "all"

needs_full_reference <- stage %in% c("integrate", "all")
needs_subcluster_references <- stage %in% c("subcluster", "peaks", "advanced", "all")

message(sprintf("Project root: %s", scscalp_cfg$project_root))
message(sprintf("Results root: %s", scscalp_cfg$results$root))

fragment_files <- scscalp_find_atac_fragment_files(scscalp_cfg$inputs$raw_atac_fragments)
message(sprintf("Found %s ATAC fragment files under %s", length(fragment_files), scscalp_cfg$inputs$raw_atac_fragments))
message(sprintf("ATAC samples: %s", paste(names(fragment_files), collapse = ", ")))

if (needs_full_reference) {
  full_path <- scscalp_resolve_rna_reference_path(level = "full")
  full_ref <- scscalp_load_rna_reference(level = "full")
  message(sprintf("Validated full RNA reference: %s", full_path))
  message(sprintf("Full RNA cells: %s", ncol(full_ref)))
}

if (needs_subcluster_references) {
  subgroups <- c("Endothelial", "Fibroblasts", "Keratinocytes", "Lymphoid", "Myeloid")
  for (subgroup in subgroups) {
    subgroup_path <- scscalp_resolve_rna_reference_path(level = "subcluster", subgroup = subgroup)
    subgroup_ref <- scscalp_load_rna_reference(level = "subcluster", subgroup = subgroup)
    message(sprintf("Validated %s RNA reference: %s", subgroup, subgroup_path))
    message(sprintf("%s RNA cells: %s", subgroup, ncol(subgroup_ref)))
  }
}

message(sprintf("Input validation complete for stage: %s", stage))
