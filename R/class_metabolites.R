# Silence R CMD check
globalVariables(c(""), package = "metaboprep2")

#' @title Metabolites object
#' @description
#' A Metabolites object is a container for matrices of metabolites data and associated metadata.
#'
#' @param project_name character, name for this project
#' @param filepath character, path to the original raw data file
#' @param format character, a valid metabolomic data format
#' @param ... optional other slots to fill, see below.
#'
#' @slot project_name character, name for this project
#' @slot filepath character, path to the original raw data file
#' @slot source character, a valid metabolomic data source
#' @slot samples data.table,
#' @slot features data.table,
#' @slot data multi-dimentional matrix,
#' @slot feature_missingness numeric 0-1 default 0.2, proportion of missinginess on features used as a QC threshold
#' @slot sample_missingness  numeric 0-1 default 0.2, proportion of missinginess on samples used as a QC threshold
#' @slot total_peak_area_sd numeric default 5, number of standard deviation (SD) units from the mean to perform a sample QC based on total abundance. If you would like to exclude total sum abundance as a QC parameter set the value to NA
#' @slot outlier_udist numeric default 5, outlier values for features. The number of interquartile range (IQR) unit distances from the median to call a value an outlier. "Outliers" here are extreme values that are the product of error.
#' @slot outlier_treatment character default 'leave_be', what to do with outliers for the purposes of the PCA, and for the PCA only. Options: "leave_be" if you would like no action on outliers; "winsorize" if you would like outliers to be winsorized to the 100th quantile of all remaining (non-outlying) values, at a feature; "turn_NA" if you would like outliers converted to NA. This means they will be imputed to the median for the purposes of the PCA.
#' @slot tree_cut_height numeric default 0.5, feature|metabolite independence. To identify "independent" features in your data set, how similar should clustered|grouped features be? tree_cut_height = 1-absolute(spearman's rho), such that a tree_cut_height of 0.2 would indicate a intra-cluster similarity of >0.8, and a tree_cut_height of 0.8 would indicate a intra-cluster similarity of >0.2. Larger tree_cut_height values yield more clustering and thus fewer "representative" or "independent" features.
#' @slot pc_outlier_sd numeric default 5, the number of standard deviation (SD) units from the mean to perform a sample QC based on principal components. If you would like to exclude PC exclusions as a QC parameter set the value to NA. Setting to NA may be advisable if you expect significant structure among individuals for example your data is derived from different tissues, different geographies/ecologies/environments.
#' @slot derived_var_exclusion logical default TRUE, Nightingale derived variable exclusion. Derived variables in Nightingale are all variables derived from two or more variables already present in the data set. In this instance it represents all ratios that Nightingale supply in their data releases.
#' @slot xenobiotics_var_exclusion logical default TRUE, whether to exclude Xenobiotics from the missingness estimate
#' @slot feature_tree object
#' @slot pcs matrix
#' @slot prob_pcs matrix
#' @slot var_exp matrix
#' @slot acceleration_factor numeric
#' @slot n_parallel numeric
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
    project_name          = new_property(class_character, default="", validator = function(value) if (value!="") NULL else "empty - please specify a project name"),
    created               = new_property(class_POSIXct, default=quote(Sys.time())),
    filepath              = new_property(class_character, validator = function(value) if (file.exists(value)) NULL else "does not appear to exist"),
    format                = new_property(class_character, default="metabolon_v1", validator = function(value) if (value %in% available_data_formats()) NULL else paste("must be one of:", paste(available_data_formats(), collapse=", "))),
    samples               = class_data.frame,
    features              = class_data.frame,
    data                  = class_numeric,
    feature_missingness   = new_property(class_numeric, default=0.2),
    sample_missingness    = new_property(class_numeric, default=0.2),
    total_peak_area_sd    = new_property(class_numeric, default=5.0),
    outlier_udist         = new_property(class_numeric, default=5.0),
    outlier_treatment     = new_property(class_character, default="leave_be", validator = function(value) if (value %in% c("leave_be", "winsorize", "turn_NA")) NULL else "@outlier_treatment should be one of leave_be|winsorize|turn_NA"),
    tree_cut_height       = new_property(class_numeric, default=0.5),
    pc_outlier_sd         = new_property(class_numeric, default=5.0),
    derived_var_exclusion     = new_property(class_logical, default=TRUE),
    xenobiotics_var_exclusion = new_property(class_logical, default=TRUE),
    feature_tree          = class_list,
    pcs                   = class_numeric,
    prob_pcs              = class_numeric,
    var_exp               = class_numeric,
    acceleration_factor   = class_numeric,
    n_parallel            = class_numeric
  ),
  validator = function(self) {
    if ((nrow(self@features)>0 & length(self@data)>0) && (nrow(self@features) != ncol(self@data))) {
      sprintf("Number of @features (%i) must equal the number metabolites in @data (%i)", nrow(self@features), ncol(self@data))
    }
    if ((nrow(self@samples)>0 & length(self@data)>0) && (nrow(self@samples) != nrow(self@data))) {
      sprintf("Number of @samples (%i) must equal the number metabolite samples in @data (%i)", nrow(self@samples), nrow(self@data))
    }
    if ((nrow(self@features)>0 & length(self@data)>0) && !check_features_table(self@features)) {
      "DEVELOPER NOTE: @features does not contain the required columns, adjust your read function"
    }

  }
)



#' @title Import Data
#' @param metabolites an object of class Metabolites
#' @export
import_data <- new_generic("import_data", c("metabolites"), function(metabolites) { S7_dispatch() })
#' @name import_data
method(import_data, Metabolites) <- function(metabolites) {

  format <- match.arg(metabolites@format, choices=available_data_formats())

  data_list <- switch(format,
                      metabolon_v1 = read_metabolon_v1(filepath))
  #metabolon_v2 = read_metabolon_v2(filepath))

  metabolites@samples  <- data_list[["samples"]]
  metabolites@features <- data_list[["features"]]
  metabolites@data     <- data_list[["data"]]

  return(metabolites)
}




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
