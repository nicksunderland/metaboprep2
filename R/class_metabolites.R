# Silence R CMD check
globalVariables(c(""), package = "metaboprep2")

#' @title Metabolites object
#' @description
#' A Metabolites object is a container for matrices of metabolites data and associated metadata.
#'
#' @param project_name character, name for this project
#' @param filepath character, path to the original raw data file
#' @param source character, a valid metabolomic data source
#'
#' @slot project_name character, name for this project
#' @slot filepath character, path to the original raw data file
#' @slot source character, a valid metabolomic data source
#'
#' @return an S7 class metaboprep2::Metabolites object
#'
#' @import S7
#' @import data.table
#'
#' @export
Metabolites <- new_class(
  name    = "Metabolites",
  package = "metaboprep2",
  properties = list(
    project_name = class_character,
    filepath     = class_character,
    source       = class_character,
    samples      = class_data.frame,
    features     = class_data.frame,
    data         = new_S3_class('array')
  ),
  constructor = function(filepath, project_name, source="metabolon_v1") {

    data_list <- switch(source,
                        metabolon_v1 = read_metabolon_v1(filepath))
                        #metabolon_v2 = read_metabolon_v2(filepath))

    object <- new_object(S7::S7_object(),
                         project_name = project_name,
                         filepath     = filepath,
                         source       = source,
                         samples      = data_list[["samples"]],
                         features     = data_list[["features"]],
                         data         = data_list[["data"]])

    return(object)
  },
  validator = function(self) {
    NULL
  }
)

#' @title get_data
#' @param metabolites an object of class Metabolites
#' @param type character, type of data to extract
#' @param metabolite_ids optional, a vector of row names to extract (default is NULL, meaning all rows)
#' @param sample_ids optional, a vector of column names to extract (default is NULL, meaning all columns)
#' @param as_df logical, whether to return the result as a data.frame (default is FALSE, meaning it returns a matrix)
#' @export
get_data <- new_generic("get_data", c("metabolites", "type"), function(metabolites, type, metabolite_ids = NULL, sample_ids = NULL, as_df = FALSE) { S7_dispatch() })
#' @name get_data
method(get_data, list(Metabolites, class_character)) <- function(metabolites, type, metabolite_ids = NULL, sample_ids = NULL, as_df = FALSE) {

  data_names <- dimnames(metabolites@data)[[3]]
  if (!(type %in% data_names)) {
    error_msg <- paste("Error: '", type, "' is not a valid type. Valid options are: ", paste(data_names, collapse = ", "), ".", sep = "")
    stop(error_msg)
  }

  data <- metabolites@data[, , type]

  if (!is.null(metabolite_ids)) {
    if (!all(metabolite_ids %in% colnames(data))) {
      stop("Error: Some of the specified col names are not found in the data.")
    }
    data <- data[, metabolite_ids, drop = FALSE]
  }

  if (!is.null(sample_ids)) {
    if (!all(sample_ids %in% rownames(data))) {
      stop("Error: Some of the specified row names are not found in the data.")
    }
    data <- data[sample_ids, , drop = FALSE]
  }

  if (as_df) {
    data <- as.data.frame(data)
  }

  return(data)
}


#' @title write_data
#' @param metabolites an object of class Metabolites
#' @param directory character, path to the write directory
#' @importFrom data.table fwrite
#' @export
write_data <- new_generic("write_data", c("metabolites", "directory"), function(metabolites, directory) { S7_dispatch() })
#' @name write_data
method(write_data, list(Metabolites, class_character)) <- function(metabolites, directory) {

  today <- gsub("-", "_", Sys.Date())
  dir   <- file.path(sub("/$", "", directory), paste0("metaboprep_release_", today), "raw_data")
  dir.create(dir, showWarnings = FALSE)

  samples_path <- file.path(dir, paste(metabolites@project_name, today, "Metabolon_sampledata.txt", sep="_"))
  data.table::fwrite(metabolites@samples, samples_path, sep="\t", quote=FALSE)

  features_path      <- file.path(dir, paste(metabolites@project_name, today, "Metabolon_featuredata.txt", sep="_"))
  features           <- as.data.frame(metabolites@features)
  rownames(features) <- features[, "feature_names"]
  data.table::fwrite(features, features_path, row.names=TRUE, sep="\t", quote=TRUE)

  sheet_names <- dimnames(metabolites@data)[[3]]
  for (sheet in sheet_names) {
    data_path <- file.path(dir, paste(metabolites@project_name, today, sheet, "Metabolon_metabolitedata.txt", sep="_"))
    tbl       <- data.table::as.data.table(metabolites@data[, , sheet])
    data.table::fwrite(tbl, data_path, row.names=TRUE, col.names = TRUE, sep="\t", quote=FALSE)
  }

  invisible(metabolites)
}
