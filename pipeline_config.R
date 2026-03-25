`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x)) || !nzchar(x[1])) {
    return(y)
  }
  x
}

scscalp_locate_project_root <- function() {
  env_root <- Sys.getenv("SCSCALP_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = FALSE))
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    start_dir <- dirname(sub("^--file=", "", file_arg[[1]]))
  } else {
    start_dir <- getwd()
  }
  start_dir <- normalizePath(start_dir, mustWork = FALSE)

  current <- start_dir
  repeat {
    sentinel_files <- c("README.md", "plotting_config.R", "sample_metadata.R")
    if (all(file.exists(file.path(current, sentinel_files)))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not determine project root. Set SCSCALP_PROJECT_ROOT.", call. = FALSE)
    }
    current <- parent
  }
}

scscalp_cfg <- local({
  project_root <- scscalp_locate_project_root()
  results_root <- normalizePath(
    Sys.getenv("SCSCALP_RESULTS_ROOT", unset = file.path(project_root, "results")),
    mustWork = FALSE
  )

  list(
    project_root = project_root,
    results = list(
      root = results_root,
      rna = normalizePath(
        Sys.getenv("SCSCALP_RNA_RESULTS_DIR", unset = file.path(results_root, "scRNA_preprocessing")),
        mustWork = FALSE
      ),
      atac = normalizePath(
        Sys.getenv("SCSCALP_ATAC_RESULTS_DIR", unset = file.path(results_root, "scATAC_preprocessing")),
        mustWork = FALSE
      )
    ),
    inputs = list(
      raw_rna = normalizePath(
        Sys.getenv(
          "SCSCALP_RAW_RNA_DIR",
          unset = file.path(project_root, "data", "scRNA", "filtered_feature_bc_matrices")
        ),
        mustWork = FALSE
      ),
      raw_atac_fragments = normalizePath(
        Sys.getenv(
          "SCSCALP_RAW_ATAC_FRAGMENT_DIR",
          unset = file.path(project_root, "data", "scATAC", "fragment_files")
        ),
        mustWork = FALSE
      ),
      demuxlet_best = Filter(
        nzchar,
        trimws(strsplit(Sys.getenv("SCSCALP_DEMUXLET_BEST", unset = ""), ",", fixed = TRUE)[[1]])
      )
    )
  )
})

scscalp_asset_path <- function(...) {
  file.path(scscalp_cfg$project_root, ...)
}

scscalp_rna_preprocess_dir <- function(...) {
  file.path(scscalp_cfg$results$rna, "preprocessing_output", ...)
}

scscalp_rna_subcluster_dir <- function(subgroup = NULL, ...) {
  pieces <- c(scscalp_cfg$results$rna, "harmonized_subclustering", subgroup, list(...))
  do.call(file.path, pieces[!vapply(pieces, is.null, logical(1))])
}

scscalp_atac_baseline_dir <- function(...) {
  file.path(scscalp_cfg$results$atac, "baseline_preprocessing", ...)
}

scscalp_atac_subcluster_dir <- function(subgroup = NULL, ...) {
  pieces <- c(scscalp_cfg$results$atac, if (!is.null(subgroup)) sprintf("subclustered_%s", subgroup), list(...))
  do.call(file.path, pieces)
}

scscalp_atac_fine_dir <- function(...) {
  file.path(scscalp_cfg$results$atac, "fine_clustered", ...)
}

scscalp_parse_version_spec <- function(spec) {
  trimmed <- trimws(spec)
  match <- regexec("^([<>]=?|==?)?\\s*([0-9][0-9A-Za-z.\\-]*)$", trimmed)
  groups <- regmatches(trimmed, match)[[1]]
  if (!length(groups)) {
    stop(sprintf("Invalid version requirement: %s", spec), call. = FALSE)
  }
  list(op = groups[[2]] %||% "==", version = groups[[3]])
}

scscalp_compare_versions <- function(current, requirement) {
  parsed <- scscalp_parse_version_spec(requirement)
  cmp <- utils::compareVersion(current, parsed$version)
  switch(
    parsed$op,
    ">" = cmp > 0,
    ">=" = cmp >= 0,
    "<" = cmp < 0,
    "<=" = cmp <= 0,
    "=" = cmp == 0,
    "==" = cmp == 0,
    stop(sprintf("Unsupported version operator: %s", parsed$op), call. = FALSE)
  )
}

scscalp_check_requested_package_versions <- function(pkgs = loadedNamespaces()) {
  pkg_names <- sort(unique(pkgs))
  for (pkg in pkg_names) {
    env_name <- paste0("SCSCALP_PKG_", gsub("[^A-Z0-9]+", "_", toupper(pkg)))
    requested <- Sys.getenv(env_name, unset = "")
    if (!nzchar(requested)) {
      next
    }
    current <- as.character(utils::packageVersion(pkg))
    if (!scscalp_compare_versions(current, requested)) {
      stop(
        sprintf("Package %s version %s does not satisfy %s=%s", pkg, current, env_name, requested),
        call. = FALSE
      )
    }
  }
}
