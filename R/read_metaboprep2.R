#' @title Read Exported Metaboprep2 Data
#'
#' @description
#' Reads data exported by `export_data()` from a specified directory and returns the necessary data to
#' populate a `Metabolites` object.
#' This includes loading the primary data matrix, feature and sample metadata, feature/sample summaries,
#' feature tree (if present), and additional metadata from YAML configuration files.
#'
#' @param dirpath A character string giving the path to the directory containing the exported data.
#' This should be the top-level directory created by `export_data()`, typically named like
#' `"metaboprep_release_<date>"`.
#'
#' @returns A reconstructed `Metabolites` object populated with all available layers and metadata.
#'
#' @export
#'
#' @importFrom yaml read_yaml
#'
#' @examples
#' \dontrun{
#' # Import from a previously exported directory
#' metabolites <- read_metaboprep2("path/to/metaboprep_release_2024_12_31")
#' }
read_metaboprep2 <- function(dirpath) {

  # testing
  if (FALSE) {
    dirpath <- "/Users/xx20081/git/metaboprep2/vignettes/metaboprep_release_2025_04_21"
  }

  # read optional matrix from a TSV
  read_optional_matrix <- function(path, row_id) {
    if (!file.exists(path)) return(NULL)
    df <- data.table::fread(path)
    row_ids <- df[[row_id]]
    df[[row_id]] <- NULL
    mat <- as.matrix(df)
    rownames(mat) <- row_ids
    mat
  }

  # read transposed matrix (feature_summary)
  read_transposed_matrix <- function(path, row_id) {
    if (!file.exists(path)) return(NULL)
    df <- data.table::fread(path)
    ids <- df[[row_id]]
    df[[row_id]] <- NULL
    mat <- t(as.matrix(df))
    colnames(mat) <- ids
    mat
  }

  # stack 3D arrays with union of all rows/cols
  stack_3d_array <- function(data_list, all_rows=NULL, all_cols=NULL) {
    if (!length(data_list)) return(integer(0L))

    layers   <- names(data_list)

    order_by_number <- function(x) {
      nums <- as.numeric(gsub("[^0-9]", "", x))
      order(nums)
    }

    if (is.null(all_rows)) {
      all_rows <- unique(unlist(lapply(data_list, function(m) rownames(m))))
      all_rows <- all_rows[order_by_number(all_rows)]
    }

    if (is.null(all_cols)) {
      all_cols <- unique(unlist(lapply(data_list, function(m) colnames(m))))
      all_cols <- all_cols[order_by_number(all_cols)]
    }

    arr <- array(NA_real_, dim = c(length(all_rows), length(all_cols), length(layers)),
                 dimnames = list(all_rows, all_cols, layers))

    for (layer in layers) {
      m <- data_list[[layer]]
      arr[rownames(m), colnames(m), layer] <- as.matrix(m)
    }

    return(arr)
  }

  # layers from the subdirectories
  layer_dirs <- list.dirs(dirpath, recursive = FALSE, full.names = TRUE)
  layers     <- basename(layer_dirs)
  layer_idx  <- sapply(layer_dirs, function(x) yaml::read_yaml(file.path(x, "config.yml"))$layer)
  layers     <- layers[layer_idx]
  layer_dirs <- layer_dirs[layer_idx]

  # initialize containers
  lists <- list(
    data_list            = list(),
    feature_summary_list = list(),
    sample_summary_list  = list(),
    pc_list              = list(),
    prob_pc_list         = list(),
    var_exp_list         = list(),
    feature_tree_list    = list(),
    config_props         = list()
  )
  features <- NULL
  samples  <- NULL

  for (i in seq_along(layers)) {
    layer <- layers[i]
    path  <- layer_dirs[i]

    # config
    config <- yaml::read_yaml(file.path(path, "config.yml"))
    config[["layer"]] <- NULL
    config[["feature_tree"]] <- NULL
    for (name in names(config)) {
      lists$config_props[[name]][[layer]] <- config[[name]]
    }

    # feature tree
    tree_path <- file.path(path, "feature_tree.RDS")
    if (file.exists(tree_path)) {
      lists$feature_tree_list[[layer]] <- readRDS(tree_path)
    }

    # shared features and samples (same across layers)
    features_tmp <- data.table::fread(file.path(path, "features.tsv"))
    samples_tmp  <- data.table::fread(file.path(path, "samples.tsv"))
    if (is.null(features)) {
      features <- features_tmp
    } else {
      features <- rbind(features, features_tmp[!feature_id %in% features$feature_id])
    }
    if (is.null(samples)) {
      samples <- samples_tmp
    } else {
      samples <- rbind(samples, samples_tmp[!sample_id %in% samples$sample_id])
    }

    # core data matrices
    dat <- data.table::fread(file.path(path, "data.tsv"), header = TRUE)
    row_ids <- dat$sample_id
    dat$sample_id <- NULL
    mat <- as.matrix(dat)
    rownames(mat) <- row_ids
    lists$data_list[[layer]] <- mat

    # optional matrices
    lists$feature_summary_list[[layer]] <- read_transposed_matrix(file.path(path, "feature_summary.tsv"), "feature_id")
    lists$sample_summary_list[[layer]]  <- read_optional_matrix(file.path(path, "sample_summary.tsv"), "sample_id")
    lists$pc_list[[layer]]              <- read_optional_matrix(file.path(path, "pcs.tsv"), "sample_id")
    lists$prob_pc_list[[layer]]         <- read_optional_matrix(file.path(path, "prob_pcs.tsv"), "sample_id")
    lists$var_exp_list[[layer]]         <- read_optional_matrix(file.path(path, "var_exp.tsv"), "pc")
  }

  feature_ids <- features$feature_id
  sample_ids  <- samples$sample_id

  out <- list(
    data            = stack_3d_array(lists$data_list, sample_ids, feature_ids),
    features        = features,
    samples         = samples,
    feature_summary = stack_3d_array(lists$feature_summary_list, all_cols = feature_ids),
    sample_summary  = stack_3d_array(lists$sample_summary_list, all_rows = sample_ids),
    pcs             = stack_3d_array(lists$pc_list, all_rows = sample_ids),
    prob_pcs        = stack_3d_array(lists$prob_pc_list, all_rows = sample_ids),
    var_exp         = stack_3d_array(lists$var_exp_list),
    feature_tree    = lists$feature_tree_list,
    config          = lists$config_props
  )

  return(out)
}
