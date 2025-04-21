# Silence R CMD check
globalVariables(c("derived_feature", "feature_id", "pathway", "sample_id",
                  "i.derived_feature", "exclude", "clean_id",
                  "i.chemical_name", "i.super_pathway", "platform", "hmdb_id",
                  "synonym", "row_id__", "key_value", "anno_fuzzy_dist", "join_key",
                  "idx", "raw_feature_ids"), package = "metaboprep2")

#' @title Annotate Features with Common IDs
#'
#' @description
#' This function annotates a `data.table` by mapping feature identifiers (`id_col`)
#' to standardized annotations. Optionally, pathway and platform information can be included.
#'
#' @param features `data.table`
#'   The dataset containing features to be annotated.
#'
#' @param fixed_match_cols `list`
#'   The name of the column containing feature IDs - can be a vector of column names if multiple columns are to be tried.
#'   This can be a column with identifiers such as `comp_ids`, `chem_ids`,
#'   or identifiers from Nightingale Excel or flat text files. If multiple columns are given then the resulting match precedence will
#'   be in order of the columns given.
#'
#' @param fuzzy_match_cols `list`,
#'   The name of the column containing pathway information.
#'   If provided, source pathway data from `features` will be included in
#'   the annotation with standardised column name `pathway`.
#' @param fuzzy_method method for fuzzy matching see fuzzyjoin::stringdist_inner_join
#' @param fuzzy_dist distance for fuzzy matching see fuzzyjoin::stringdist_inner_join
#'
#' @return `data.table`
#'   A modified version of `features` with annotations added based on the provided IDs.
#'
#' @importFrom fuzzyjoin stringdist_inner_join
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' # Example dataset
#' features <- data.table::data.table(comp_ids = c("63042", "62919", "62921"))
#'
#' # Annotate features using 'comp_ids'
#' annotated_features <- annotate_features(features, fixed_match_cols = list(comp_id = 'comp_ids'))
#' annotated_features
#'
annotate_features <- function(features,
                              fixed_match_cols = list(pubchem_id = NULL, hmdb_id = NULL, kegg_id = NULL, chebi_id = NULL, metlin_id = NULL, name = NULL, smiles = NULL, inchi_key = NULL, chem_id = NULL, comp_id = NULL),
                              fuzzy_match_cols = list(name = NULL),
                              fuzzy_method = "jw",
                              fuzzy_dist = 0.1) {

  # testing
  if (FALSE) {
    file_paths <- c(
      rdata_path    = file.path(Sys.getenv("BBS_METABOLON_DIR"),  "ReportData.Rdata"),
      log_path      = file.path(Sys.getenv("BBS_METABOLON_DIR"),  "bbs_only_2024_06_25_logfile.txt"),
      clinical_path = file.path(Sys.getenv("BBS_CLINICAL_DIR"),   "03_sample_clinical_data_all.csv")
    )
    metabolites <- Metabolites(project_name = "bbse",
                               format       = "metaboprep1",
                               filepath     = file_paths)
    metabolites <- import_data(metabolites)
    features <- get_features(metabolites, layer = "raw")
    features[, comp_id := sub(".*?([0-9]+).*", "\\1", feature_id)]
    fixed_match_cols <- list(pubchem_id = NULL, hmdb_id = "hmdb", kegg_id = "kegg", chebi_id = NULL, metlin_id = NULL, name = "chemical_name", smiles = "smiles", inchi_key = "inchikey", chem_id = "chem_id", comp_id = "comp_id")
    fuzzy_match_cols <- list(name = NULL) #"chemical_name")
    fuzzy_method <- "jw"
    fuzzy_dist <- 0.1
  }

  # checks
  fixed_match_cols <- fixed_match_cols[!sapply(fixed_match_cols, is.null)]
  fuzzy_match_cols <- fuzzy_match_cols[!sapply(fuzzy_match_cols, is.null)]
  stopifnot("features must contain all the column(s) in `fixed_match_cols` and `fuzzy_match_cols`" = all(c(fixed_match_cols, fuzzy_match_cols) %in% names(features)))
  anno_cols <- c("pubchem_id", "hmdb_id", "kegg_id", "chebi_id", "metlin_id", "name", "smiles", "inchi_key", "chem_id", "comp_id")
  for (ann_col in names(c(fixed_match_cols, fuzzy_match_cols))) {
    msg <- paste0("Annotation column `", ann_col, "` must be one of ", paste0(anno_cols, collapse=", "))
    if (!ann_col %in% anno_cols) {
      stop(msg)
    }
  }
  fuzzy_anno_cols <- c("name")
  for (ann_col in names(fuzzy_match_cols)) {
    msg <- paste0("Fuzzy annotation column `", ann_col, "` must be ", paste0(fuzzy_anno_cols, collapse=", "))
    if (!ann_col %in% fuzzy_anno_cols) {
      stop(msg)
    }
  }

  # get the Metaboanalyst database
  ma_anno <- get_metaboanalyser_db()
  ma_anno[, `:=`(chem_id = NA_character_, comp_id = NA_character_, derived_feature = FALSE)]

  # get the additional annotations file
  anno_file <- system.file("extdata", "feature_annotations.tsv", package = "metaboprep2")
  local_anno <- data.table::fread(anno_file, na.strings=c("NA", ""))
  local_anno <- local_anno[, list(hmdb_id = unlist(tstrsplit(hmdb_id, ",", fixed=TRUE))), by=setdiff(names(local_anno), "hmdb_id")]

  # combine annotation tables
  anno <- rbind(ma_anno, local_anno, fill=TRUE)

  # expand the synonyms if searching by name
  if ("name" %in% names(c(fixed_match_cols, fuzzy_match_cols))) {
    anno <- anno[, list(name = unlist(tstrsplit(synonym, "; ?"))), by=setdiff(names(anno), "name")]
  }

  # unique
  anno <- unique(anno)

  # key pairs
  features[, row_id__ := .I]
  features[, join_key := NA_character_]


  # cycle possible fixed joins
  for (k in seq_along(fixed_match_cols)) {

    key         <- fixed_match_cols[k]
    feature_key <- key[[1]]
    anno_key    <- names(key)[1]
    cat("Attempting join", paste0("[features$", feature_key, "]"), "~", paste0("[annotation$", anno_key, "]"), "\n")

    # unmatched and has a local key
    unmatched <- features[is.na(join_key) & !is.na(get(feature_key)), list(row_id__ = row_id__, key_value = get(feature_key))]
    if (anno_key != "name") { # dont split names which might actually have commas
      unmatched <- unmatched[,
                             list(row_id__ = row_id__,
                               key_value = trimws(unlist(tstrsplit(key_value, "[,;]", fixed = FALSE)))),
                             by = .SD
      ][!is.na(key_value)]
    }
    if (nrow(unmatched) == 0) break  # all done

    # has an anno key
    filtered_anno <- anno[!is.na(get(anno_key))]

    # join
    anno_key_prefixed <- paste0("anno_", anno_key)
    data.table::setnames(filtered_anno, names(filtered_anno), paste0("anno_", names(filtered_anno)))
    unmatched[, key_value := as.character(key_value)]
    filtered_anno[, (anno_key_prefixed) := as.character(get(anno_key_prefixed))]
    joined <- filtered_anno[unmatched, on=stats::setNames("key_value", anno_key_prefixed), nomatch=NULL, mult="first"][, anno_fuzzy_dist := NA_real_]
    if (nrow(joined) == 0) next

    # add to features
    cols_to_assign <- grep("^anno_", names(joined), value = TRUE)

    for (col in cols_to_assign) {
      features[joined$row_id__, (col) := joined[[col]]]
    }
    features[joined$row_id__, join_key := paste0(feature_key, " ~ ", anno_key)]
  }


  # cycle possible fuzzy joins
  for (k in seq_along(fuzzy_match_cols)) {

    key         <- fuzzy_match_cols[k]
    feature_key <- key[[1]]
    anno_key    <- names(key)[1]
    cat("Attempting fuzzy join", paste0("[features$", feature_key, "]"), "~", paste0("[annotation$", anno_key, "]"), "\n")

    # unmatched and has a local key
    unmatched <- features[is.na(join_key) & !is.na(get(feature_key)), list(row_id__ = row_id__, key_value = get(feature_key))]
    if (nrow(unmatched) == 0) break  # all done

    # has an anno key
    filtered_anno <- anno[!is.na(get(anno_key))]

    # join
    anno_key_prefixed <- paste0("anno_", anno_key)
    data.table::setnames(filtered_anno, names(filtered_anno), paste0("anno_", names(filtered_anno)))
    joined <- fuzzyjoin::stringdist_inner_join(
      unmatched,
      filtered_anno,
      by = stats::setNames(anno_key_prefixed, "key_value"),
      method = fuzzy_method,
      max_dist = fuzzy_dist,
      ignore_case = TRUE,
      distance_col = "anno_fuzzy_dist"
    )
    joined <- as.data.table(joined)[, .SD[which.min(anno_fuzzy_dist)], by = row_id__]
    if (nrow(joined) == 0) next

    # add to features
    cols_to_assign <- grep("^anno_", names(joined), value = TRUE)

    for (col in cols_to_assign) {
      features[joined$row_id__, (col) := joined[[col]]]
    }
    features[joined$row_id__, join_key := paste0(feature_key, " ~fuzzy~ ", anno_key)]
  }

  # clean up
  features[, row_id__ := NULL]

  # ensure all anno columns exist in features, even if none matched
  expected_anno_cols <- paste0("anno_", names(anno))
  missing_cols <- setdiff(expected_anno_cols, names(features))
  if (length(missing_cols) > 0) {
    features[, (missing_cols) := NA_character_]
  }

  # warn number of unmatched
  if (any(is.na(features$join_key))) {
    warning(paste0(sum(is.na(features$join_key)), " unannotated metabolites remain after merging with metabolomics database\n"))
  }

  return(features)
}



#' @title Download MetaboAnalyseR Database
#' @returns data.table
#' @importFrom qs qread
#' @importFrom data.table as.data.table
#' @importFrom utils download.file
#' @export
#'
get_metaboanalyser_db <- function() {

  # caching
  cache_dir <- tools::R_user_dir("metaboprep2", which = "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_dir, "master_compound_db.qs")

  # check if already downloaded
  if (!file.exists(cache_file)) {
    message("Downloading MetaboAnalyst DB...")
    tryCatch(
      expr = {
        url <- "https://www.metaboanalyst.ca/resources/libs/master_compound_db.qs"
        utils::download.file(url, cache_file, method = "curl")
      },
      error = function(e) {
        stop("Failed to download MetaboAnalyst DB: ", e$message)
      }
    )
  } else {
    message("Using cached MetaboAnalyst DB.")
  }

  # read
  db <- tryCatch(
    expr = {
      data.table::as.data.table(qs::qread(cache_file))
    },
    error = function(e) {
      stop("Failed to read cached MetaboAnalyst DB: ", e$message)
    }
  )

  return(db)
}
