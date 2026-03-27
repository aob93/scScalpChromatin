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

scscalp_find_10x_sample_dirs <- function(base_dir) {
  sample_dirs <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)
  sample_names <- basename(sample_dirs)

  nested_matrix_dirs <- file.path(sample_dirs, "outs", "filtered_feature_bc_matrix")
  nested_ok <- dir.exists(nested_matrix_dirs)

  direct_ok <- file.exists(file.path(sample_dirs, "matrix.mtx")) |
    file.exists(file.path(sample_dirs, "matrix.mtx.gz"))

  resolved_dirs <- c(
    stats::setNames(nested_matrix_dirs[nested_ok], sample_names[nested_ok]),
    stats::setNames(sample_dirs[direct_ok & !nested_ok], sample_names[direct_ok & !nested_ok])
  )

  if (!length(resolved_dirs)) {
    stop(
      sprintf(
        "No 10x matrix directories found under %s. Expected sample/outs/filtered_feature_bc_matrix or direct matrix directories.",
        base_dir
      ),
      call. = FALSE
    )
  }

  resolved_dirs
}

scscalp_find_atac_fragment_files <- function(base_dir) {
  sample_dirs <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)
  sample_names <- basename(sample_dirs)

  direct_fragment_files <- c("fragments.tsv.gz", "atac_fragments.tsv.gz")
  nested_fragment_paths <- c(
    file.path("outs", "fragments.tsv.gz"),
    file.path("outs", "atac_fragments.tsv.gz")
  )

  resolved <- character()

  if (length(sample_dirs)) {
    for (i in seq_along(sample_dirs)) {
      sample_dir <- sample_dirs[[i]]
      sample_name <- sample_names[[i]]

      direct_match <- file.path(sample_dir, direct_fragment_files)
      direct_match <- direct_match[file.exists(direct_match)]
      nested_match <- file.path(sample_dir, nested_fragment_paths)
      nested_match <- nested_match[file.exists(nested_match)]

      chosen <- c(direct_match, nested_match)
      if (length(chosen)) {
        resolved[[sample_name]] <- chosen[[1]]
      }
    }
  }

  root_level_matches <- file.path(base_dir, direct_fragment_files)
  root_level_matches <- root_level_matches[file.exists(root_level_matches)]
  if (length(root_level_matches)) {
    resolved[[basename(base_dir)]] <- root_level_matches[[1]]
  }

  if (!length(resolved)) {
    stop(
      sprintf(
        "No ATAC fragment files found under %s. Expected fragments.tsv.gz or atac_fragments.tsv.gz directly under each sample directory or under sample/outs/.",
        base_dir
      ),
      call. = FALSE
    )
  }

  resolved
}

scscalp_rna_preprocess_dir <- function(...) {
  file.path(scscalp_cfg$results$rna, "preprocessing_output", ...)
}

scscalp_rna_subcluster_dir <- function(subgroup = NULL, ...) {
  pieces <- c(scscalp_cfg$results$rna, "harmonized_subclustering", subgroup, list(...))
  do.call(file.path, pieces[!vapply(pieces, is.null, logical(1))])
}

scscalp_env_key <- function(x) {
  key <- gsub("[^A-Z0-9]+", "_", toupper(x))
  gsub("^_+|_+$", "", key)
}

scscalp_resolve_rna_reference_path <- function(level = c("full", "subcluster"), subgroup = NULL) {
  level <- match.arg(level)

  if (identical(level, "full")) {
    override <- Sys.getenv("SCSCALP_RNA_REFERENCE_RDS", unset = "")
    path <- if (nzchar(override)) override else scscalp_rna_preprocess_dir("scalp.rds")
    return(normalizePath(path, mustWork = FALSE))
  }

  if (is.null(subgroup) || !nzchar(subgroup)) {
    stop("subgroup must be provided when resolving a subcluster RNA reference.", call. = FALSE)
  }

  subgroup_env <- paste0("SCSCALP_RNA_REFERENCE_SUBGROUP_", scscalp_env_key(subgroup), "_RDS")
  override <- Sys.getenv(subgroup_env, unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, mustWork = FALSE))
  }

  ref_dir <- Sys.getenv("SCSCALP_RNA_REFERENCE_SUBCLUSTER_DIR", unset = "")
  if (nzchar(ref_dir)) {
    candidates <- c(
      file.path(ref_dir, subgroup, sprintf("%s.rds", subgroup)),
      file.path(ref_dir, sprintf("%s.rds", subgroup)),
      file.path(ref_dir, subgroup, "scalp.rds")
    )
    existing <- candidates[file.exists(candidates)]
    chosen <- if (length(existing)) existing[[1]] else candidates[[1]]
    return(normalizePath(chosen, mustWork = FALSE))
  }

  normalizePath(scscalp_rna_subcluster_dir(subgroup, sprintf("%s.rds", subgroup)), mustWork = FALSE)
}

scscalp_derive_disease_status <- function(sample_ids) {
  status <- rep(NA_character_, length(sample_ids))
  status[grepl("C_SD", sample_ids)] <- "C_SD"
  status[grepl("C_PB", sample_ids)] <- "C_PB"
  status[grepl("AA", sample_ids)] <- "AA"
  status
}

scscalp_require_umap_reduction <- function(object, path, context) {
  reduction_names <- names(object@reductions)
  if ("umap" %in% reduction_names) {
    return(object)
  }

  matching_umap <- reduction_names[tolower(reduction_names) == "umap"]
  if (length(matching_umap)) {
    object[["umap"]] <- object[[matching_umap[[1]]]]
    return(object)
  }

  stop(
    sprintf(
      "RNA reference at %s is missing a 'umap' reduction required for %s plotting.",
      path,
      context
    ),
    call. = FALSE
  )
}

scscalp_prepare_rna_reference <- function(object, path, level = c("full", "subcluster"), subgroup = NULL) {
  level <- match.arg(level)
  context <- if (identical(level, "full")) {
    "full RNA integration"
  } else {
    sprintf("%s RNA integration", subgroup)
  }

  if (!inherits(object, "Seurat")) {
    stop(sprintf("RNA reference at %s is not a Seurat object.", path), call. = FALSE)
  }

  meta <- object[[]]

  if (!"Sample" %in% colnames(meta) && "orig.ident" %in% colnames(meta)) {
    object$Sample <- object$orig.ident
    meta <- object[[]]
  }

  if (!"diseaseStatus" %in% colnames(meta) && "Sample" %in% colnames(meta)) {
    object$diseaseStatus <- scscalp_derive_disease_status(object$Sample)
    meta <- object[[]]
  }

  if (identical(level, "full") && !"BroadClust" %in% colnames(meta) && "NamedClust" %in% colnames(meta)) {
    broad_clust <- gsub("[0-9]+", "", object$NamedClust)
    object$BroadClust <- sub(".", "", broad_clust)
    meta <- object[[]]
  }

  required_cols <- if (identical(level, "full")) {
    c("NamedClust", "BroadClust", "Sample", "diseaseStatus")
  } else {
    c("FineClust", "Sample", "diseaseStatus")
  }

  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols)) {
    stop(
      sprintf(
        "RNA reference at %s is missing required metadata columns for %s: %s",
        path,
        context,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  scscalp_require_umap_reduction(object, path, context)
}

scscalp_load_rna_reference <- function(level = c("full", "subcluster"), subgroup = NULL) {
  level <- match.arg(level)
  path <- scscalp_resolve_rna_reference_path(level = level, subgroup = subgroup)

  if (!file.exists(path)) {
    stop(
      sprintf(
        "RNA reference RDS not found at %s. Configure SCSCALP_RNA_REFERENCE_RDS or SCSCALP_RNA_REFERENCE_SUBCLUSTER_DIR / subgroup-specific overrides.",
        path
      ),
      call. = FALSE
    )
  }

  object <- readRDS(path)
  scscalp_prepare_rna_reference(object = object, path = path, level = level, subgroup = subgroup)
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

scscalp_normalize_sample_id <- function(sample_ids) {
  normalized <- as.character(sample_ids)
  normalized <- sub("^AA_([0-9]+)$", "AA\\1", normalized)
  normalized <- sub("^C_PB_([0-9]+)$", "C_PB\\1", normalized)
  normalized <- sub("^C_SD_([0-9]+)$", "C_SD\\1", normalized)
  normalized
}

scscalp_lookup_sample_metadata <- function(
  sample_ids,
  metadata,
  field_name = "sample metadata",
  allow_missing = FALSE,
  missing_value = NA
) {
  normalized_ids <- scscalp_normalize_sample_id(sample_ids)
  matched <- metadata[normalized_ids]
  missing <- unique(normalized_ids[vapply(matched, is.null, logical(1))])
  if (length(missing)) {
    if (!allow_missing) {
      stop(
        sprintf(
          "Missing %s for sample IDs: %s",
          field_name,
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    matched[vapply(matched, is.null, logical(1))] <- list(missing_value)
  }
  unname(unlist(matched, use.names = FALSE))
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

scscalp_require_hg38_bsgenome <- function() {
  pkg <- "BSgenome.Hsapiens.UCSC.hg38"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("Required package %s is not installed in this R environment.", pkg),
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
  invisible(TRUE)
}

scscalp_start_logging <- function(logfile, append = FALSE) {
  tee_flag <- if (append) "-a" else ""
  output_cmd <- trimws(sprintf("tee %s %s", tee_flag, shQuote(logfile)))
  message_cmd <- sprintf("%s >&2", output_cmd)

  output_con <- pipe(output_cmd, open = "wt")
  message_con <- pipe(message_cmd, open = "wt")

  sink(output_con, type = "output")
  sink(message_con, type = "message")

  list(
    logfile = logfile,
    output_con = output_con,
    message_con = message_con
  )
}

scscalp_stop_logging <- function(log_state) {
  if (!is.null(log_state$message_con)) {
    sink(type = "message")
  }
  if (!is.null(log_state$output_con)) {
    sink(type = "output")
  }
  if (!is.null(log_state$message_con) && isOpen(log_state$message_con)) {
    close(log_state$message_con)
  }
  if (!is.null(log_state$output_con) && isOpen(log_state$output_con)) {
    close(log_state$output_con)
  }
}

scscalp_enable_seurat_slot_compat <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  if (utils::compareVersion(as.character(utils::packageVersion("SeuratObject")), "5.0.0") < 0) {
    return(invisible(FALSE))
  }

  ns <- asNamespace("SeuratObject")

  get_assay_data_compat <- function(object, assay = NULL, layer = NULL, slot = deprecated(), ...) {
    CheckDots(...)
    if (!missing(slot) && is.null(layer)) {
      layer <- slot
    }
    object <- UpdateSlots(object = object)
    assay <- assay %||% DefaultAssay(object = object)
    assay <- arg_match(arg = assay, values = Assays(object = object))
    GetAssayData(object = object[[assay]], layer = layer)
  }

  set_assay_data_compat <- function(object, layer = "data", new.data, slot = deprecated(), assay = NULL, ...) {
    CheckDots(...)
    if (!missing(slot) && is.null(layer)) {
      layer <- slot
    }
    object <- UpdateSlots(object = object)
    assay <- assay %||% DefaultAssay(object = object)
    object[[assay]] <- SetAssayData(object = object[[assay]], layer = layer, new.data = new.data, ...)
    object
  }

  environment(get_assay_data_compat) <- ns
  environment(set_assay_data_compat) <- ns

  base::registerS3method("GetAssayData", "Seurat", get_assay_data_compat, envir = ns)
  base::registerS3method("SetAssayData", "Seurat", set_assay_data_compat, envir = ns)
  invisible(TRUE)
}

scscalp_parse_size_bytes <- function(x) {
  spec <- trimws(toupper(x))
  match <- regexec("^([0-9]*\\.?[0-9]+)\\s*([KMGTP]?)B?$", spec)
  groups <- regmatches(spec, match)[[1]]
  if (!length(groups)) {
    stop(sprintf("Invalid size specification: %s", x), call. = FALSE)
  }

  value <- as.numeric(groups[[2]])
  unit <- groups[[3]]
  multipliers <- stats::setNames(
    c(1, 1024^1, 1024^2, 1024^3, 1024^4, 1024^5),
    c("", "K", "M", "G", "T", "P")
  )
  if (!unit %in% names(multipliers)) {
    stop(sprintf("Unsupported size unit in: %s", x), call. = FALSE)
  }
  multiplier <- unname(multipliers[[unit]])
  value * multiplier
}

scscalp_future_maxsize_bytes <- function() {
  spec <- Sys.getenv("SCSCALP_FUTURE_GLOBALS_MAXSIZE", unset = "")
  if (!nzchar(spec)) {
    legacy_gb <- Sys.getenv("SCSCALP_FUTURE_GLOBALS_MAXSIZE_GB", unset = "")
    spec <- if (nzchar(legacy_gb)) paste0(legacy_gb, "G") else "20G"
  }
  scscalp_parse_size_bytes(spec)
}

scscalp_configure_future <- function() {
  options(future.globals.maxSize = scscalp_future_maxsize_bytes())
}

scscalp_set_future_plan <- function(strategy, workers, ...) {
  if (is.character(strategy)) {
    future::plan(
      strategy,
      workers = workers,
      maxSizeOfObjects = scscalp_future_maxsize_bytes(),
      ...
    )
  } else {
    future::plan(
      future::tweak(
        strategy,
        workers = workers,
        maxSizeOfObjects = scscalp_future_maxsize_bytes(),
        ...
      )
    )
  }
}

scscalp_configure_future()
scscalp_enable_seurat_slot_compat()

scscalp_get_assay_data <- function(object, layer, assay = NULL, ...) {
  args <- list(object = object, ...)
  if (!is.null(assay)) {
    args$assay <- assay
  }

  get_assay_data <- get("GetAssayData", envir = asNamespace("SeuratObject"), inherits = FALSE)
  object_class <- class(object)[[1]]
  method_formals <- tryCatch(
    names(formals(getS3method("GetAssayData", object_class))),
    error = function(e) names(formals(get_assay_data))
  )
  if ("layer" %in% method_formals) {
    args$layer <- layer
  } else {
    args$slot <- layer
  }

  tryCatch(
    do.call(get_assay_data, args),
    error = function(e) {
      if (!inherits(object, "Seurat") || !"layer" %in% names(args) || !grepl("multiple layers", conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }

      joined_object <- SeuratObject::JoinLayers(
        object = object,
        assay = assay %||% SeuratObject::DefaultAssay(object),
        layers = layer,
        new = layer
      )
      args$object <- joined_object
      do.call(get_assay_data, args)
    }
  )
}

scscalp_set_assay_data <- function(object, layer, new.data, assay = NULL, ...) {
  args <- list(object = object, new.data = new.data, ...)
  if (!is.null(assay)) {
    args$assay <- assay
  }

  set_assay_data <- get("SetAssayData", envir = asNamespace("SeuratObject"), inherits = FALSE)
  object_class <- class(object)[[1]]
  method_formals <- tryCatch(
    names(formals(getS3method("SetAssayData", object_class))),
    error = function(e) names(formals(set_assay_data))
  )
  if ("layer" %in% method_formals) {
    args$layer <- layer
  } else {
    args$slot <- layer
  }

  do.call(set_assay_data, args)
}

scscalp_diet_seurat <- function(object, assays = NULL, layers = NULL, features = NULL, dimreducs = NULL, ...) {
  diet_seurat <- get("DietSeurat", envir = asNamespace("Seurat"), inherits = FALSE)
  args <- list(
    object = object,
    assays = assays,
    features = features,
    dimreducs = dimreducs,
    ...
  )

  if ("layers" %in% names(formals(diet_seurat))) {
    args$layers <- layers
  } else {
    if (!is.null(layers)) {
      args$counts <- "counts" %in% layers
      args$data <- "data" %in% layers
      args$scale.data <- "scale.data" %in% layers
    }
  }

  do.call(diet_seurat, args)
}

scscalp_should_use_knn_only <- function() {
  flag <- tolower(trimws(Sys.getenv("SCSCALP_CLUSTER_USE_KNN_ONLY", unset = "")))
  flag %in% c("1", "true", "t", "yes", "y")
}

scscalp_is_neighbor_memory_error <- function(err) {
  msg <- conditionMessage(err)
  patterns <- c(
    "std::bad_alloc",
    "cannot allocate",
    "vector memory exhausted",
    "sparse->dense coercion",
    "exceeding 2\\^31"
  )
  any(vapply(patterns, grepl, logical(1), x = msg, ignore.case = TRUE))
}

scscalp_find_neighbors_and_clusters <- function(
  object,
  reduction,
  dims,
  resolution,
  k.param = 20,
  prune.SNN = 1 / 15,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "euclidean",
  nn.eps = 0,
  l2.norm = FALSE,
  algorithm = 1,
  random.seed = 1,
  group.singletons = TRUE,
  verbose = TRUE,
  use_knn_only = scscalp_should_use_knn_only()
) {
  assay <- SeuratObject::DefaultAssay(object[[reduction]])
  data.use <- Seurat::Embeddings(object[[reduction]])
  data.use <- data.use[, dims, drop = FALSE]
  graph_names <- c(
    nn = paste0(reduction, "_nn"),
    snn = paste0(reduction, "_snn")
  )
  build_graphs <- function(compute.SNN) {
    Seurat::FindNeighbors(
      object = data.use,
      k.param = k.param,
      compute.SNN = compute.SNN,
      prune.SNN = prune.SNN,
      nn.method = nn.method,
      n.trees = n.trees,
      annoy.metric = annoy.metric,
      nn.eps = nn.eps,
      verbose = verbose,
      l2.norm = l2.norm
    )
  }
  store_graphs <- function(target, graphs) {
    for (nm in names(graphs)) {
      if (inherits(graphs[[nm]], "Graph")) {
        SeuratObject::DefaultAssay(graphs[[nm]]) <- assay
      }
      target[[graph_names[[nm]]]] <- graphs[[nm]]
    }
    target
  }

  cluster_graph <- unname(graph_names[["snn"]])
  if (use_knn_only) {
    message(sprintf("Using kNN graph only for clustering on %s.", reduction))
    object <- store_graphs(object, build_graphs(compute.SNN = FALSE))
    cluster_graph <- unname(graph_names[["nn"]])
  } else {
    object <- tryCatch(
      store_graphs(object, build_graphs(compute.SNN = TRUE)),
      error = function(err) {
        if (!scscalp_is_neighbor_memory_error(err)) {
          stop(err)
        }
        message(
          sprintf(
            "SNN graph construction failed for %s (%s). Falling back to kNN graph only.",
            reduction,
            conditionMessage(err)
          )
        )
        store_graphs(object, build_graphs(compute.SNN = FALSE))
      }
    )

    if (!cluster_graph %in% names(object)) {
      cluster_graph <- unname(graph_names[["nn"]])
    }
  }

  Seurat::FindClusters(
    object = object,
    graph.name = cluster_graph,
    resolution = resolution,
    algorithm = algorithm,
    random.seed = random.seed,
    group.singletons = group.singletons,
    verbose = verbose
  )
}
