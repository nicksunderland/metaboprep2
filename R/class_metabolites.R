# Silence R CMD check
globalVariables(c(""), package = "metaboprep2")

#' @title Metabolites object
#' @description
#' A Metabolites object is a container for matrices of metabolites data and associated metadata.
#'
#' @param project_name character, name for this project
#' @param filepath character, path to the original raw data file
#' @param format character, a valid metabolomic data format
#' @param feature_missingness numeric 0-1 default 0.2, proportion of missinginess on features used as a QC threshold
#' @param sample_missingness  numeric 0-1 default 0.2, proportion of missinginess on samples used as a QC threshold
#' @param total_peak_area_sd numeric default 5, number of standard deviation (SD) units from the mean to perform a sample QC based on total abundance. If you would like to exclude total sum abundance as a QC parameter set the value to NA
#' @param outlier_udist numeric default 5, outlier values for features. The number of interquartile range (IQR) unit distances from the median to call a value an outlier. "Outliers" here are extreme values that are the product of error.
#' @param outlier_treatment character default 'leave_be', what to do with outliers for the purposes of the PCA, and for the PCA only. Options: "leave_be" if you would like no action on outliers; "winsorize" if you would like outliers to be winsorized to the 100th quantile of all remaining (non-outlying) values, at a feature; "turn_NA" if you would like outliers converted to NA. This means they will be imputed to the median for the purposes of the PCA.
#' @param winsorize_quantile numeric defailt 1, winsorize to what quantile of remaining (not outlier) values
#' @param tree_cut_height numeric default 0.5, feature|metabolite independence. To identify "independent" features in your data set, how similar should clustered|grouped features be? tree_cut_height = 1-absolute(spearman's rho), such that a tree_cut_height of 0.2 would indicate a intra-cluster similarity of >0.8, and a tree_cut_height of 0.8 would indicate a intra-cluster similarity of >0.2. Larger tree_cut_height values yield more clustering and thus fewer "representative" or "independent" features.
#' @param pc_outlier_sd numeric default 5, the number of standard deviation (SD) units from the mean to perform a sample QC based on principal components. If you would like to exclude PC exclusions as a QC parameter set the value to NA. Setting to NA may be advisable if you expect significant structure among individuals for example your data is derived from different tissues, different geographies/ecologies/environments.
#' @param derived_var_exclusion logical default TRUE, Nightingale derived variable exclusion. Derived variables in Nightingale are all variables derived from two or more variables already present in the data set. In this instance it represents all ratios that Nightingale supply in their data releases.
#' @param xenobiotics_var_exclusion logical default TRUE, whether to exclude Xenobiotics from the missingness estimate
#' @param created POSIXct, when this object was created
#' @param source named character, source for this data layer e.g. 'raw' may be the source for 'post_qc'
#' @param samples data.table, not to be set directly
#' @param features data.table, not to be set directly
#' @param data numeric matrix, not to be set directly
#' @param exclusions multidimensional sparse character matrix, exclusion codes (mask for data), not to be set directly
#' @param feature_summary numeric matrix, not to be set directly
#' @param feature_tree named list of tree objects matrix, not to be set directly
#' @param sample_summary numeric matrix, not to be set directly
#' @param pcs numeric matrix, not to be set directly
#' @param prob_pcs numeric matrix, not to be set directly
#' @param var_exp named numeric, not to be set directly
#' @param acceleration_factor named numeric, not to be set directly
#' @param n_parallel named numeric, not to be set directly
#'
#' @slot project_name character, name for this project
#' @slot filepath character, path to the original raw data file
#' @slot source character, a valid metabolomic data source
#' @slot samples data.frame / data.table
#' @slot features data.frame / data.table
#' @slot data multi-dimentional matrix,
#' @slot exclusions matrix???
#' @slot feature_summary matrix
#' @slot sample_summary matrix
#' @slot feature_missingness numeric 0-1 default 0.2, proportion of missinginess on features used as a QC threshold
#' @slot sample_missingness  numeric 0-1 default 0.2, proportion of missinginess on samples used as a QC threshold
#' @slot total_peak_area_sd numeric default 5, number of standard deviation (SD) units from the mean to perform a sample QC based on total abundance. If you would like to exclude total sum abundance as a QC parameter set the value to NA
#' @slot outlier_udist numeric default 5, outlier values for features. The number of interquartile range (IQR) unit distances from the median to call a value an outlier. "Outliers" here are extreme values that are the product of error.
#' @slot outlier_treatment character default 'leave_be', what to do with outliers for the purposes of the PCA, and for the PCA only. Options: "leave_be" if you would like no action on outliers; "winsorize" if you would like outliers to be winsorized to the 100th quantile of all remaining (non-outlying) values, at a feature; "turn_NA" if you would like outliers converted to NA. This means they will be imputed to the median for the purposes of the PCA.
#' @slot winsorize_quantile numeric defailt 1, winsorize to what quantile of remaining (not outlier) values
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
    # parameters
    project_name              = new_property(class_character, default="", validator = function(value) if (value!="") NULL else "empty - please specify a project name"),
    created                   = new_property(class_POSIXct, default=quote(Sys.time())),
    filepath                  = new_property(class_character, validator = function(value) if (file.exists(value)) NULL else "does not appear to exist"),
    format                    = new_property(class_character, default="metabolon_v1", validator = function(value) if (value %in% available_data_formats()) NULL else paste("must be one of:", paste(available_data_formats(), collapse=", "))),
    source                    = new_property(class_character),
    feature_missingness       = new_property(class_numeric, default=0.2, validator = function(value) if (data.table::between(value, 0, 1)) NULL else "should be a value between 0 and 1, inclusive"),
    sample_missingness        = new_property(class_numeric, default=0.2, validator = function(value) if (data.table::between(value, 0, 1)) NULL else "should be a value between 0 and 1, inclusive"),
    total_peak_area_sd        = new_property(class_numeric, default=5.0, validator = function(value) if (value>0) NULL else "should be a value >0"),
    outlier_udist             = new_property(class_numeric, default=5.0, validator = function(value) if (value>0) NULL else "should be a value >0"),
    outlier_treatment         = new_property(class_character, default="leave_be", validator = function(value) if (value %in% c("leave_be", "winsorize", "turn_NA")) NULL else "should be one of leave_be|winsorize|turn_NA"),
    winsorize_quantile        = new_property(class_numeric, default=1.0),
    tree_cut_height           = new_property(class_numeric, default=0.5, validator = function(value) if (value>=0) NULL else "should be a value >=0"),
    pc_outlier_sd             = new_property(class_numeric, default=5.0, validator = function(value) if (is.na(value) || value>0) NULL else "should be a value >0 or NA"),
    derived_var_exclusion     = new_property(class_logical, default=TRUE),
    xenobiotics_var_exclusion = new_property(class_logical, default=TRUE),
    # data objects
    samples               = class_data.frame,
    features              = class_data.frame,
    data                  = class_numeric,
    exclusions            = new_property(class_any, default=SparseArray::SparseArray(matrix(character())), validator = function(value) if (SparseArray::is_sparse(value)) NULL else "should be a SparseArray"),
    feature_summary       = class_numeric,
    feature_tree          = class_list,
    sample_summary        = class_numeric,
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
                      metabolon_v1 = read_metabolon_v1(metabolites@filepath))
  #metabolon_v2 = read_metabolon_v2(metabolites@filepath))

  metabolites@samples    <- data_list[["samples"]]
  metabolites@features   <- data_list[["features"]]
  metabolites@data       <- data_list[["data"]]

  return(metabolites)
}



#' Get Exclusion Codes
#' @param verbose logical, printing
#' @returns named list, list("code_description" = code)
#' @importFrom utils str
#' @export
#' @examples
#' get_exclusion_codes()
get_exclusion_codes <- function(verbose=TRUE) {

  l <- list("derived_feature" = 1,
            "xenobiotic_feature" = 2,
            "extreme_sample_missingness" = 3,
            "extreme_feature_missingness" = 4,
             "user_defined_sample_missingness" = 5,
             "user_defined_feature_missingness" = 6,
             "user_defined_sample_totalpeakarea" = 7,
             "user_defined_sample_pca_outlier" = 8,
             "outlier_udist_turned_na" = 9,
             "outlier_udist_winsorized" = 10)

  if (verbose) {
    utils::str(l)
  }

  invisible(l)
}

#' @title update_exclusions
#' @param metabolites an object of class Metabolites
#' @param type character, type of data
#' @param code integer, exclusions code
#' @param feature_ids optional, a vector of row names to extract (default is NULL, meaning all rows)
#' @param sample_ids optional, a vector of column names to extract (default is NULL, meaning all columns)
#' @param arr_ids logical, whether ids are array indices
#' @importFrom SparseArray SparseArray
#' @export
update_exclusions <- new_generic("update_exclusions", c("metabolites", "type"), function(metabolites, type, code, sample_ids = NULL, feature_ids = NULL, arr_ids = FALSE) { S7_dispatch() })
#' @name update_exclusions
method(update_exclusions, list(Metabolites, class_character)) <- function(metabolites, type, code, sample_ids = NULL, feature_ids = NULL, arr_ids = FALSE) {

  # create exclusions matrix if doesn't exist
  if (length(metabolites@exclusions) == 0 || !type %in% dimnames(metabolites@exclusions)[[3]]) {
    excl_mat <- SparseArray::SparseArray(array("", dim = c(dim(metabolites@data)[1:2], 1),
                                                   dimnames = c(dimnames(metabolites@data)[1:2], type)),
                                         type="character")
    metabolites@exclusions <- add_layer(current    = metabolites@exclusions,
                                        layer      = excl_mat,
                                        layer_name = type)
  }

  if (arr_ids) {
    stopifnot("Indicated sample_ids and feature_ids are array IDs, but they are of different lengths" = !is.null(feature_ids) && !is.null(sample_ids) && length(sample_ids)==length(feature_ids))

    # exclude some samples for some features and vice versa
    for (fi in feature_ids) {
      for (si in sample_ids) {
        metabolites@exclusions[si, fi, type] <- sub("^,","", paste(metabolites@exclusions[si, fi, type], code, sep=","))
      }
    }

  } else {

    # exclude samples for all features
    if (!is.null(sample_ids) && length(sample_ids) > 0) {
      s <- which( rownames(metabolites@exclusions) %in% sample_ids )
      f <- seq_len(ncol(metabolites@exclusions))
      for (fi in f) {
        for (si in s) {
          metabolites@exclusions[si, fi, type] <- sub("^,","", paste(metabolites@exclusions[si, fi, type], code, sep=","))
        }
      }
    }

    # exclude features for all samples
    if (!is.null(feature_ids) && length(feature_ids) > 0) {
      f <- which( colnames(metabolites@exclusions) %in% feature_ids )
      s <- seq_len(nrow(metabolites@exclusions))
      for (fi in f) {
        for (si in s) {
          metabolites@exclusions[si, fi, type] <- sub("^,","", paste(metabolites@exclusions[si, fi, type], code, sep=","))
        }
      }
    }

  }
  return(metabolites)
}


#' @title get_data
#' @param metabolites an object of class Metabolites
#' @param type character, type of data to extract
#' @param apply_exclusions logical, whether to exclude samples/features based on the current `@exclusions` matrix
#' @param feature_ids optional, a vector of row names to extract (default is NULL, meaning all rows)
#' @param sample_ids optional, a vector of column names to extract (default is NULL, meaning all columns)
#' @param as_df logical, whether to return the result as a data.frame (default is FALSE, meaning it returns a matrix)
#' @export
get_data <- new_generic("get_data", c("metabolites", "type"), function(metabolites, type, apply_exclusions = TRUE, feature_ids = NULL, sample_ids = NULL, as_df = FALSE) { S7_dispatch() })
#' @name get_data
method(get_data, list(Metabolites, class_character)) <- function(metabolites, type, apply_exclusions = TRUE, feature_ids = NULL, sample_ids = NULL, as_df = FALSE) {

  data_names <- dimnames(metabolites@data)[[3]]
  if (!(type %in% data_names)) {
    error_msg <- paste("Error: '", type, "' is not a valid type. Valid options are: ", paste(data_names, collapse = ", "), ".", sep = "")
    stop(error_msg)
  }

  data <- metabolites@data[, , type]

  # apply exclusions if requested and present
  if (apply_exclusions==TRUE & type %in% dimnames(metabolites@exclusions)[[3]]) {
    excl <- metabolites@exclusions[, , type]
    ex_f <- which(apply(excl, 2, function(x) all(x!=""))) # all empty, i.e. no exclusions
    ex_s <- which(apply(excl, 1, function(x) all(x!=""))) # all empty, i.e. no exclusions
    if (length(ex_s) > 0) data <- data[-ex_s, , drop = FALSE]
    if (length(ex_f) > 0) data <- data[, -ex_f, drop = FALSE]
  }

  if (!is.null(feature_ids)) {
    if (!all(feature_ids %in% colnames(data))) {
      stop("Error: Some of the specified col names are not found in the data. If exclusons==TRUE, check they are not excluded for some reason.")
    }
    data <- data[, feature_ids, drop = FALSE]
  }

  if (!is.null(sample_ids)) {
    if (!all(sample_ids %in% rownames(data))) {
      stop("Error: Some of the specified row names are not found in the data. If exclusons==TRUE, check they are not excluded for some reason.")
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
