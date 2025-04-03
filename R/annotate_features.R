# Silence R CMD check
globalVariables(c("derived_feature", "feature_id", "pathway", "sample_id",
                  "i.derived_feature", "exclude", "clean_id",
                  "i.chemical_name", "i.super_pathway", "platform",
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
#' @param id_col `character`
#'   The name of the column containing feature IDs.
#'   This can be a column with identifiers such as `comp_ids`, `chem_ids`,
#'   or identifiers from Nightingale Excel or flat text files.
#'
#' @param pathway_col `character`, *optional*
#'   The name of the column containing pathway information.
#'   If provided, pathway data will be included in the annotation.
#'
#' @param platform_col `character`, *optional*
#'   The name of the column specifying the platform type.
#'   If supplied, platform-specific annotations will be added.
#'
#' @return `data.table`
#'   A modified version of `features` with annotations added based on the provided IDs.
#'
#' @export
#'
#' @examples
#' # Example dataset
#' features <- data.table::data.table(comp_ids = c("63042", "62919", "62921"))
#'
#' # Annotate features using 'comp_ids'
#' annotated_features <- annotate_features(features, id_col = "comp_ids")
#'
#' # Annotate with additional pathway and platform information
#' annotated_features <- annotate_features(features,
#'                                         id_col = "comp_ids",
#'                                         pathway_col = "pathway",
#'                                         platform_col = "platform")
#'
annotate_features <- function(features,
                              id_col,
                              pathway_col=NULL,
                              platform_col=NULL) {

  # check
  stopifnot("features must contain the `id_col`" = id_col %in% names(features))

  # get the annotations file
  anno_file <- system.file("extdata", "feature_annotations.tsv", package = "metaboprep2")
  anno <- data.table::fread(anno_file)[exclude==FALSE, ]
  keys <- names(anno)[grepl("_id", names(anno)) & !grepl("old", names(anno))]

  # initialize feature_id column as "unmapped" by default
  features[, `:=`(feature_id      = paste0("unmapped_", get(id_col)),
                  platform        = NA_character_,
                  derived_feature = FALSE,
                  idx             = .I)]

  # cycle the possible join keys to find a match
  features[, clean_id := clean_names(get(id_col))]
  data.table::setkey(features, clean_id)
  for (key in rev(keys)) { # rev ensures that order of precedence is chem_id>comp_id>ng1_id>ng2_id
    anno[, (key) := clean_names(get(key))]
    data.table::setkeyv(anno, key)
    features[anno, `:=`(feature_id      = paste0(key, "_", get(paste0("i.",key))),
                        derived_feature = i.derived_feature,
                        chemical_name   = i.chemical_name)]
    if (is.null(pathway_col)) {
      features[anno, pathway := i.super_pathway]
    }
    if (is.null(platform_col)) {
      features[anno, platform := NA_character_]
    }
  }

  # ensure same order as entry
  data.table::setcolorder(features, "feature_id")
  features <- features[order(idx)]
  features[, idx := NULL]
  features[, clean_id := NULL]

  # warn if unmapped
  unmapped <- features[grepl("unmapped", feature_id)]
  if (nrow(unmapped) > 0) {
    warning(
      paste0(nrow(unmapped), " unmapped features:\n - ",
             paste0(unmapped[, get(id_col)], collapse=", "))
    )
  }

  return(features)
}
