# Silence R CMD check
globalVariables(c(""), package = "metaboprep2")

#' @title Metabolites Object
#' @description
#' A `Metabolites` object is a container for matrices of metabolite data, along with associated metadata.
#' It allows for efficient storage and manipulation of data, supporting quality control, transformations,
#' and various analyses. This object facilitates easy access to data layers, sample and feature summaries,
#' outlier treatment, and more.
#'
#' @param project_name character, the name of the project associated with the data.
#' @param created POSIXct, the timestamp when the `Metabolites` object was created.
#' @param filepath character, the file path to the original raw data file.
#' @param format character, the format of the metabolomic data (e.g., "metabolon_v1").
#' @param feature_missingness numeric, the proportion (0-1) of missingness in features, used as a QC threshold. Default is 0.2.
#' @param sample_missingness numeric, the proportion (0-1) of missingness in samples, used as a QC threshold. Default is 0.2.
#' @param total_peak_area_sd numeric, the number of standard deviations from the mean to use for sample QC based on total abundance. Default is 5. If set to `NA`, total sum abundance will not be used for QC.
#' @param outlier_udist numeric, the number of interquartile range (IQR) unit distances from the median used to identify outliers for features. Default is 5.
#' @param outlier_treatment character, specifies how to treat outliers for PCA analysis. Options:
#'     "leave_be" (no action on outliers), "winsorize" (winsorize to the 100th quantile),
#'     "turn_NA" (convert outliers to NA, which will be imputed to the median).
#' @param winsorize_quantile numeric, the quantile of remaining values to use for winsorizing outliers. Default is 1.
#' @param tree_cut_height numeric, the threshold for feature independence in hierarchical clustering. Default is 0.5.
#' @param pc_outlier_sd numeric, the number of standard deviations from the mean used for PCA-based sample QC. Default is 5. Set to `NA` to exclude PCA-based QC.
#' @param derived_var_exclusion logical, whether to exclude derived variables (from multiple existing variables) in Nightingale data. Default is `TRUE`.
#' @param xenobiotics_var_exclusion logical, whether to exclude xenobiotics from missingness estimates. Default is `TRUE`.
#' @param source named character, the source for this data layer (e.g., 'raw' for 'post_qc').
#' @param samples data.table, a data table containing sample-related information (not to be set directly).
#' @param features data.table, a data table containing feature-related information (not to be set directly).
#' @param data numeric matrix, the data matrix containing metabolite values (not to be set directly).
#' @param exclusions list, holds exclusion codes for data masking (not to be set directly).
#' @param feature_summary numeric matrix, summary statistics for features (not to be set directly).
#' @param feature_tree named list of tree objects, hierarchical clustering results for features (not to be set directly).
#' @param sample_summary numeric matrix, summary statistics for samples (not to be set directly).
#' @param pcs numeric matrix, principal component scores for the data (not to be set directly).
#' @param prob_pcs numeric matrix, probabilities associated with principal components (not to be set directly).
#' @param var_exp named numeric, variance explained by each principal component (not to be set directly).
#' @param acceleration_factor named numeric, the acceleration factor for parallel analysis (not to be set directly).
#' @param n_parallel named numeric, the number of parallel components in the analysis (not to be set directly).
#'
#' @slot project_name character, name for this project.
#' @slot created POSIXct, the timestamp when the `Metabolites` object was created.
#' @slot filepath character, path to the original raw data file.
#' @slot format character, a valid metabolomic data source.
#' @slot samples data.table, the samples data table.
#' @slot features data.table, the features data table.
#' @slot data numeric matrix, the metabolite data.
#' @slot source named character, the source for this data layer (e.g., 'raw' for 'post_qc').
#' @slot exclusions list, exclusion codes (mask for data).
#' @slot feature_summary numeric matrix, feature summary statistics.
#' @slot sample_summary numeric matrix, sample summary statistics.
#' @slot feature_missingness numeric, proportion of missingness in features (QC threshold).
#' @slot sample_missingness numeric, proportion of missingness in samples (QC threshold).
#' @slot total_peak_area_sd numeric, standard deviation for sample QC based on total abundance.
#' @slot outlier_udist numeric, outlier threshold in terms of IQR for features.
#' @slot outlier_treatment character, outlier treatment for PCA analysis.
#' @slot winsorize_quantile numeric, quantile used for winsorizing outliers.
#' @slot tree_cut_height numeric, threshold for feature independence in clustering.
#' @slot pc_outlier_sd numeric, standard deviation for PCA-based sample QC.
#' @slot derived_var_exclusion logical, whether to exclude derived variables from Nightingale data.
#' @slot xenobiotics_var_exclusion logical, whether to exclude xenobiotics from missingness estimate.
#' @slot feature_tree object, hierarchical clustering results for features.
#' @slot pcs matrix, principal component scores.
#' @slot prob_pcs matrix, probabilities for principal components.
#' @slot var_exp matrix, variance explained by each principal component.
#' @slot acceleration_factor numeric, acceleration factor for parallel analysis.
#' @slot n_parallel numeric, number of parallel components.
#'
#' @return An object of class `Metabolites` (S7 class).
#'
#' @import S7
#' @import data.table
#'
#' @examples
#' # Example of importing data into a Metabolites object and running QC
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' metabolites
#'
#' @export
Metabolites <- new_class(
  name    = "Metabolites",
  package = "metaboprep2",
  properties = list(
    # parameters
    project_name              = new_property(class_character, default="", validator = function(value) if (value!="") NULL else "empty - please specify a project name"),
    created                   = new_property(class_POSIXct, default=quote(Sys.time())),
    filepath                  = new_property(class_character, validator = function(value) if (file.exists(value) || dir.exists(value)) NULL else "does not appear to exist"),
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
    samples               = new_property(class_data.frame, default=quote(data.table::data.table(sample_id=character())), validator = function(value) if (length(setdiff(c("sample_id"), names(value))) > 0) "should have a column sample_id" else NULL),
    features              = new_property(class_data.frame, default=quote(data.table::data.table(feature_id=character(), platform=character(), pathway=character(), derived_feature=logical())), validator = function(value) if (length(setdiff(c("feature_id", "platform", "pathway", "derived_feature"), names(value))) > 0) "should have columns feature_id, platform, pathway, derived_feature" else NULL),
    data                  = class_numeric,
    exclusions            = new_property(class_list,
                                         default=list(raw = list(samples  = list(extreme_sample_missingness        = character(),
                                                                                 user_defined_sample_missingness   = character(),
                                                                                 user_defined_sample_totalpeakarea = character(),
                                                                                 user_defined_sample_pca_outlier   = character()),
                                                                 features = list(extreme_feature_missingness       = character(),
                                                                                 user_defined_feature_missingness  = character()))),
                                         validator = function(value) {
                                           if (length(value)==0) return(NULL)
                                           check_names <- all(sapply(value, function(layer) identical(names(layer), c("samples","features"))))
                                           check_samp  <- all(sapply(value, function(layer) identical(names(layer[["samples"]]), c("extreme_sample_missingness","user_defined_sample_missingness", "user_defined_sample_totalpeakarea", "user_defined_sample_pca_outlier"))))
                                           check_feat  <- all(sapply(value, function(layer) identical(names(layer[["features"]]), c("extreme_feature_missingness","user_defined_feature_missingness"))))
                                           if (check_names && check_samp && check_feat) {
                                              return(NULL)
                                           } else {
                                              return("should be a list of lists with names 'samples' and 'features'. In addition, 'samples' should be a list of character vectors named 'extreme_sample_missingness','user_defined_sample_missingness', 'user_defined_sample_totalpeakarea', 'user_defined_sample_pca_outlier'; and 'features' a list of character vectors named 'extreme_feature_missingness', 'user_defined_feature_missingness'")
                                           }
                                         }),
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
      sprintf("Number of @features (%i) must equal the number of features in @data (%i)", nrow(self@features), ncol(self@data))
    }
    if ((nrow(self@samples)>0 & length(self@data)>0) && (nrow(self@samples) != nrow(self@data))) {
      sprintf("Number of @samples (%i) must equal the number of samples in @data (%i)", nrow(self@samples), nrow(self@data))
    }
    if ((nrow(self@samples)>0 & length(self@data)>0) && !identical(self@samples[, sample_id], rownames(self@data))) {
      "Column `sample_id` in @samples must be identical to the rownames of @data"
    }
    if ((nrow(self@features)>0 & length(self@data)>0) && !identical(self@features[, feature_id], colnames(self@data))) {
      "Column `feature_id` in @features must be identical to the colnames of @data"
    }

  }
)



#' @title Import Data
#'
#' @description
#' This function is used to import data into a `Metabolites` object. It processes and loads the data
#' in a format that is compatible with subsequent analysis steps in the `metaboprep2` package.
#'
#' @param metabolites An object of class `Metabolites`. This object contains the necessary metadata
#'                    and structure for the data being imported.
#'
#' @return The function returns the updated `Metabolites` object with the imported data integrated.
#'
#' @examples
#' # Example of importing data into a Metabolites object
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#'
#' @export
import_data <- new_generic("import_data", c("metabolites"), function(metabolites) { S7_dispatch() })
#' @name import_data
method(import_data, Metabolites) <- function(metabolites) {

  format <- match.arg(metabolites@format, choices=available_data_formats())

  if (format != "metaboprep2") {

    data_list <- switch(format,
                        metabolon_v1   = read_metabolon_v1(metabolites@filepath),
                        metabolon_v2   = read_metabolon_v2(metabolites@filepath),
                        nightingale_v1 = read_nightingale_v1(metabolites@filepath),
                        nightingale_v2 = read_nightingale_v2(metabolites@filepath))

    metabolites@samples    <- data_list[["samples"]]
    metabolites@features   <- data_list[["features"]]
    metabolites@data       <- data_list[["data"]]

  } else {

    data_list <- read_metaboprep2(dirpath = metabolites@filepath)

    metabolites@samples         <- data_list[["samples"]]
    metabolites@features        <- data_list[["features"]]
    metabolites@data            <- data_list[["data"]]
    metabolites@sample_summary  <- data_list[["sample_summary"]]
    metabolites@feature_summary <- data_list[["feature_summary"]]
    metabolites@pcs             <- data_list[["pcs"]]
    metabolites@prob_pcs        <- data_list[["prob_pcs"]]
    metabolites@var_exp         <- data_list[["var_exp"]]
    metabolites@feature_tree    <- data_list[["feature_tree"]]

    for (i in seq_along(data_list$config)) {
      slot_name <- names(data_list$config)[i]
      slot_data <- data_list$config[[i]]
      S7::prop(metabolites, slot_name) <- switch(slot_name,
                                                 "created"             = as.POSIXct(slot_data[[1]]),
                                                 "exclusions"          = as.list(slot_data),
                                                 "source"              = stats::na.omit(unlist(slot_data)),
                                                 "acceleration_factor" = stats::na.omit(unlist(slot_data)),
                                                 "n_parallel"          = stats::na.omit(unlist(slot_data)),
                                                 slot_data[[1]])
    }

  }

  return(metabolites)
}


#' @title Add a Layer of Data (internal use)
#' @description
#' This function adds an additional layer of data along the third dimension to an existing 3D array (or 2D matrix/vector) by stacking a new layer of data.
#' It ensures that the dimensions of the new layer match the first two dimensions of the existing array or matrix.
#' If there is a mismatch in row or column names and the `force` parameter is set to `TRUE`, the function will align the data by filling missing values with `NA`.
#' It is used internally and not intended for routine user use.
#'
#' @param current A vector, matrix, or 3D array representing the current stack of data.
#' @param layer A matrix or array that represents the new layer of data to be added. It should match the dimensions of the first two dimensions of `current`.
#' @param layer_name A character string specifying the name of the new dimension for the 3rd axis. This can be used to annotate the new data layer.
#' @param force A logical value indicating whether to force the join and create `NA` values where row or column names do not match between `current` and `layer`. Default is `FALSE`.
#'
#' @returns A 3D array with the added layer in the third dimension.
#'
#' @examples
#' # Example: Adding a new data layer to an existing 3D array
#' current_stack <- array(1:12,
#'                        dim = c(3, 4, 1),
#'                        dimnames = list(c(1:3), c(1:4), c(1)))
#' new_layer <- matrix(13:16, nrow = 3, ncol = 4)
#' rownames(new_layer) <- 1:3
#' colnames(new_layer) <- 1:4
#' updated_stack <- add_layer(current    = current_stack,
#'                            layer      = new_layer,
#'                            layer_name = "Layer2",
#'                            force      = FALSE)
#' dim(updated_stack)
#'
#' # Example with mismatched row/column names (force = TRUE)
#' rownames(current_stack) <- c("A", "B", "C")
#' rownames(new_layer) <- c("B", "C", "D")
#' updated_stack_with_force <- add_layer(current    = current_stack,
#'                                       layer      = new_layer,
#'                                       layer_name = "Layer2",
#'                                       force      = TRUE)
#' dim(updated_stack_with_force)
#'
#' @export
add_layer <- function(current, layer, layer_name, force=FALSE) {

  if (is.vector(layer)) {
    layer <- array(layer,
                   dim = c(length(layer), 1, 1),
                   dimnames = list(names(layer), "value", layer_name))
  } else if (is.matrix(layer)) {
    layer <- array(layer,
                   dim = c(nrow(layer), ncol(layer), 1),
                   dimnames = list(rownames(layer), colnames(layer), layer_name))
  }

  # If current is empty, initialize it directly with the first layer (no NAs)
  if (length(current) == 0) {
    return(layer)
  }

  # Ensure row and column names match before appending
  if (!identical(rownames(current), rownames(layer)) ||
      !identical(colnames(current), colnames(layer))) {

    # force together anyway (ok for things like PCs or varexp)
    if (force) {
      all_rows <- union(rownames(current), rownames(layer))
      all_cols <- union(colnames(current), colnames(layer))
      current_expanded <- array(NA_real_,
                                dim = c(length(all_rows), length(all_cols), dim(current)[3]),
                                dimnames = list(all_rows, all_cols, dimnames(current)[[3]]))
      current_expanded[rownames(current), colnames(current), ] <- current
      layer_expanded <- array(NA_real_,
                              dim = c(length(all_rows), length(all_cols), 1),
                              dimnames = list(all_rows, all_cols, layer_name))
      layer_expanded[rownames(layer), colnames(layer), ] <- layer
      current <- current_expanded
      layer <- layer_expanded

    } else {
      stop("Error: Row names and column names must match before joining layers.")
    }
  }

  # Check if the layer_name already exists in the current array depth
  existing_layer_index <- which(dimnames(current)[[3]] == layer_name)

  # If the layer exists, overwrite it
  if (length(existing_layer_index) > 0) {

    current[, , existing_layer_index] <- layer

    # If it doesnt append layer
  } else {

    if (is.array(current)) {

      current <- array(c(current, layer),
                       dim = c(dim(layer)[1], dim(layer)[2], dim(current)[3] + 1),
                       dimnames = list(rownames(layer),
                                       colnames(layer),
                                       c(dimnames(current)[[3]], layer_name)))

    } else {
      stop("why is current not an array or sparse array")
    }

  }

  return(current)
}



#' @title Internal Data Extraction Helper
#'
#' @description
#' A helper function for internal use that facilitates data extraction from a `Metabolites` object.
#' This function is utilized by several exposed data getter functions, including:
#' \link{get_data}, \link{get_features}, \link{get_samples}, \link{get_feature_summary},
#' \link{get_sample_summary}, \link{get_feature_tree}, \link{get_pcs}, \link{get_var_exp},
#' \link{get_prob_pcs}, \link{get_acceleration_factor}, and \link{get_n_parallel}.
#'
#' @param metabolites An object of class `Metabolites`.
#' @param slot the name of the metabolites object slot to extract
#' @param layer Character. The layer of data to extract.
#' @param feature_ids Optional. A vector of row names to extract. Defaults to `NULL`, meaning all rows are extracted.
#' @param sample_ids Optional. A vector of column names to extract. Defaults to `NULL`, meaning all columns are extracted.
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE` (returns a matrix).
#' @param drop Logical. Whether to drop all NA rows or columns (excluded samples or features)
#'
#' @return A matrix or a `data.frame` containing the extracted data.
#'
#' @importFrom glue glue
#'
get_ <- new_generic("get_", c("metabolites", "slot", "layer"), function(metabolites, slot, layer, feature_ids = NULL, sample_ids = NULL, drop=FALSE, as_df = FALSE) { S7_dispatch() })
#' @name get_
method(get_, list(Metabolites, class_character, class_character)) <- function(metabolites, slot, layer, feature_ids = NULL, sample_ids = NULL, drop=FALSE, as_df = FALSE) {

  # check
  stopifnot("`as_df` should be a logical" = is.logical(as_df))
  slot <- match.arg(slot, choices = c("data", "samples", "features", "feature_summary", "sample_summary",
                                      "feature_tree", "pcs", "prob_pcs", "var_exp", "acceleration_factor",
                                      "n_parallel"))

  # get the slot data
  data <- slot(metabolites, slot)

  # check if empty
  if (length(data) == 0) {
    error_msg <- glue::glue("Error: slot @{slot} is length 0.")
    stop(error_msg)
  }

  # what kind of slot is this
  slot_type <- class(data)[1]

  # ensure data layer exists - names in the maximal dimension
  data_layers <- switch(slot_type,
                        "data.table" = dimnames(metabolites@data)[[3]],         # samples / features
                        "data.frame" = dimnames(metabolites@data)[[3]],         # samples / features
                        "numeric"    = names(data),          # var_exp
                        "integer"    = names(data),          # n_parallel / acceleration_factor
                        "list"       = names(data),          # feature tree list
                        "array"      = dimnames(data)[[3]])  # 3D data arrays (data, summaries, pcs, var_exp etc)

  # ensure requested type is a layer of data
  layer <- match.arg(layer, choices = data_layers)

  # extract the data layer
  if (slot_type == "array") {
    data <- data[, , layer, drop=TRUE]
  } else if (slot_type != "data.table") {
    data <- data[[layer]]
  }

  # subset matrix by features and samples (cant subset var_exp)
  if (slot_type == "array" && slot != "var_exp") {

    # feature exclusions (if pcs matrix cols are PCs not features, so ignore)
    cols_are_not_features <- c("pcs", "prob_pcs", "samples", "sample_summary")

    # requested feature filtering
    if (!is.null(feature_ids) && !slot %in% cols_are_not_features) {
      if (!all(feature_ids %in% colnames(data))) {
        stop("Error: Some of the specified feature names are not found in colnames(data).")
      }
      data <- data[, feature_ids, drop = FALSE]
    }

    # requested sample filtering
    if (!is.null(sample_ids)) {
      if (!all(sample_ids %in% rownames(data))) {
        stop("Error: Some of the specified sample names are not found in rownames(data).")
      }
      data <- data[sample_ids, , drop = FALSE]
    }

  # subset tables by features or samples
  } else if (grepl("data\\.table|data\\.frame", slot_type)) {

    # requested feature filtering
    if (!is.null(feature_ids) && "feature_id" %in% names(data)) {
      if (!all(feature_ids %in% data[, feature_id])) {
        stop("Error: Some of the specified feature names are not found in data[, feature_id].")
      }
      data <- data[feature_id %in% feature_ids, ]
    }

    # requested sample filtering
    if (!is.null(sample_ids) && "sample_id" %in% names(data)) {
      if (!all(sample_ids %in% data[, sample_id])) {
        stop("Error: Some of the specified sample names are not found in data[, sample_id].")
      }
      data <- data[sample_id %in% sample_ids, ]
    }

  }

  # drop all NAs
  if (drop) {
    if (is.vector(data)) {
      data <- data[!is.na(data)]
    } else if (is.matrix(data)) {
      data <- data[rowSums(is.na(data)) < ncol(data), , drop = FALSE]
      data <- data[, colSums(is.na(data)) < nrow(data), drop = FALSE]
    }
  }

  # return as data.frame
  if (as_df) {
    data <- as.data.frame(data)
  }

  # return
  return(data)
}



#' @title Retrieve Parallel Analysis No. PCs
#' @description
#' Extracts the number of principal components estimated to be informative
#' from the parallel analysis on the data in specified data `layer` from a `Metabolites` object.
#' @param metabolites A `Metabolites` object containing the acceleration factor data.
#' @param layer A character string specifying the data layer from which to extract the acceleration factor.
#' @return An `integer` representing the result of the parallel analysis for the specified `layer`.
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_n_parallel(metabolites, layer = "post_qc")
#' @export
get_n_parallel <- new_generic("get_n_parallel", c("metabolites", "layer"), function(metabolites, layer) { S7_dispatch() })
#' @name get_n_parallel
method(get_n_parallel, list(Metabolites, class_character)) <- function(metabolites, layer) {

  return(get_(metabolites, slot="n_parallel", layer=layer, feature_ids=NULL, sample_ids=NULL, drop=FALSE, as_df=FALSE))

}

#' @title Retrieve the Acceleration Factor
#' @description
#' Extracts the acceleration factor for a specified data `layer` from a `Metabolites` object.
#' @param metabolites A `Metabolites` object containing the acceleration factor data.
#' @param layer A character string specifying the data layer from which to extract the acceleration factor.
#' @return An `integer` representing the acceleration factor for the specified `layer`.
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_acceleration_factor(metabolites, layer = "post_qc")
#' @export
get_acceleration_factor <- new_generic("get_acceleration_factor", c("metabolites", "layer"), function(metabolites, layer) { S7_dispatch() })
#' @name get_acceleration_factor
method(get_acceleration_factor, list(Metabolites, class_character)) <- function(metabolites, layer) {

  return(get_(metabolites, slot="acceleration_factor", layer=layer, feature_ids=NULL, sample_ids=NULL, drop=FALSE, as_df=FALSE))

}


#' @title Retrieve Probabilistic Principal Components
#'
#' @description
#' Extracts the probabilistic principal components (PCs) from a `Metabolites` object
#' for a specified data `layer`, with optional exclusions applied.
#'
#' @param metabolites A `Metabolites` object containing the probabilistic PC data.
#' @param layer A character string specifying the data layer from which to extract the probabilistic PCs.
#' @param sample_ids Optional. A vector of column names to extract. Defaults to `NULL`, meaning all columns are included.
#' @param drop Logical. Whether to drop completely NA rows (i.e. excluded)
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE`, returning a matrix.
#'
#' @return A matrix or `data.frame` (if `as_df = TRUE`) containing the probabilistic principal components.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_prob_pcs(metabolites, layer = "post_qc")
#'
#' @export
get_prob_pcs <- new_generic("get_prob_pcs", c("metabolites", "layer"), function(metabolites, layer, sample_ids = NULL, drop=FALSE, as_df = FALSE) { S7_dispatch() })
#' @name get_prob_pcs
method(get_prob_pcs, list(Metabolites, class_character)) <- function(metabolites, layer, sample_ids = NULL, drop=FALSE, as_df = FALSE) {

  return(get_(metabolites, slot="prob_pcs", layer=layer, feature_ids=NULL, sample_ids=sample_ids, drop=drop, as_df=as_df))

}


#' @title Retrieve Principal Components
#'
#' @description
#' Extracts the principal components (PCs) from a `Metabolites` object
#' for a specified data `layer`, with optional exclusions and subsetting.
#'
#' @param metabolites A `Metabolites` object containing the principal component data.
#' @param layer A character string specifying the data layer from which to extract the PCs.
#' @param sample_ids Optional. A vector of column names to extract. Defaults to `NULL`, meaning all columns are included.
#' @param drop Logical. Whether to drop completely NA rows (i.e. excluded)
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE`, returning a matrix.
#'
#' @return A matrix or `data.frame` (if `as_df = TRUE`) containing the principal components.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_pcs(metabolites, layer = "post_qc")
#'
#' @export
get_pcs <- new_generic("get_pcs", c("metabolites", "layer"), function(metabolites, layer, sample_ids = NULL, drop=FALSE, as_df = FALSE) { S7_dispatch() })
#' @name get_pcs
method(get_pcs, list(Metabolites, class_character)) <- function(metabolites, layer, sample_ids = NULL, drop=FALSE, as_df = FALSE) {

  return(get_(metabolites, slot="pcs", layer=layer, feature_ids=NULL, sample_ids=sample_ids, drop=drop, as_df=as_df))

}


#' @title Retrieve Variance Explained by PCs
#'
#' @description
#' Extracts the variance explained by each principal components (PCs) from a `Metabolites` object
#' for a specified data `layer`, with optional exclusions and subsetting.
#'
#' @param metabolites A `Metabolites` object containing the principal component data.
#' @param layer A character string specifying the data layer from which to extract the PCs.
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE`, returning a matrix.
#'
#' @return A matrix or `data.frame` (if `as_df = TRUE`) containing the principal components.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_var_exp(metabolites, layer = "post_qc")
#'
#' @export
get_var_exp <- new_generic("get_var_exp", c("metabolites", "layer"), function(metabolites, layer, as_df = FALSE) { S7_dispatch() })
#' @name get_var_exp
method(get_var_exp, list(Metabolites, class_character)) <- function(metabolites, layer, as_df = FALSE) {

  return(get_(metabolites, slot="var_exp", layer=layer, feature_ids=NULL, sample_ids=NULL, drop=FALSE, as_df=as_df))

}


#' @title Retrieve Feature Tree
#'
#' @description
#' Extracts the hierarchical feature clustering (feature tree) from a `Metabolites` object
#' for a specified data `layer`, with optional exclusions and subsetting.
#'
#' @param metabolites A `Metabolites` object containing the feature tree data.
#' @param layer A character string specifying the data layer from which to extract the feature tree.
#'
#' @return A hierarchical clustering object (`hclust`) or a `data.frame` (if `as_df = TRUE`) containing feature tree data.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_feature_tree(metabolites, layer = "post_qc")
#'
#' @export
get_feature_tree <- new_generic("get_feature_tree", c("metabolites", "layer"), function(metabolites, layer) { S7_dispatch() })
#' @name get_feature_tree
method(get_feature_tree, list(Metabolites, class_character)) <- function(metabolites, layer) {

  return(get_(metabolites, slot="feature_tree", layer=layer, feature_ids=NULL, sample_ids=NULL, drop=FALSE, as_df=FALSE))

}


#' @title Retrieve Feature Data
#'
#' @description
#' Extracts feature-level data from a `Metabolites` object for a specified data `layer`,
#' with optional exclusions and subsetting.
#'
#' @param metabolites A `Metabolites` object containing feature data.
#' @param layer A character string specifying the data layer from which to extract features.
#' @param feature_ids Optional. A vector of row names to extract. Defaults to `NULL`, meaning all rows are included.
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE`, returning a matrix.
#'
#' @return A data.table or `data.frame` (if `as_df = TRUE`) containing feature data for the specified `layer`.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_features(metabolites, layer = "post_qc")
#'
#' @export

get_features <- new_generic("get_features", c("metabolites", "layer"), function(metabolites, layer, feature_ids = NULL, as_df = FALSE) { S7_dispatch() })
#' @name get_features
method(get_features, list(Metabolites, class_character)) <- function(metabolites, layer, feature_ids = NULL, as_df = FALSE) {

  return(get_(metabolites, slot="features", layer=layer, feature_ids=feature_ids, sample_ids=NULL, drop=FALSE, as_df=as_df))

}


#' @title Retrieve Sample Data
#'
#' @description
#' Extracts sample-level data from a `Metabolites` object for a specified data `layer`,
#' with optional exclusions and subsetting.
#'
#' @param metabolites A `Metabolites` object containing sample data.
#' @param layer A character string specifying the data layer from which to extract samples.
#' @param sample_ids Optional. A vector of column names to extract. Defaults to `NULL`, meaning all columns are included.
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE`, returning a matrix.
#'
#' @return A matrix or `data.frame` (if `as_df = TRUE`) containing sample data for the specified `layer`.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_samples(metabolites, layer = "post_qc")
#'
#' @export
get_samples <- new_generic("get_samples", c("metabolites", "layer"), function(metabolites, layer, sample_ids = NULL, as_df = FALSE) { S7_dispatch() })
#' @name get_samples
method(get_samples, list(Metabolites, class_character)) <- function(metabolites, layer, sample_ids = NULL, as_df = FALSE) {

  return(get_(metabolites, slot="samples", layer=layer, feature_ids=NULL, sample_ids=sample_ids, drop=FALSE, as_df=as_df))

}


#' @title Retrieve Feature Summary Data
#'
#' @description
#' Extracts feature-level summary statistics from a `Metabolites` object for a specified data `layer`,
#' with optional exclusions and subsetting.
#'
#' @param metabolites A `Metabolites` object containing feature summary data.
#' @param layer A character string specifying the data layer from which to extract feature summary data.
#' @param feature_ids Optional. A vector of row names (features) to extract. Defaults to `NULL`, meaning all features are included.
#' @param drop Logical. Whether to drop all NA (excluded features)
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE`, meaning the result is returned as a matrix.
#'
#' @return A matrix or `data.frame` (if `as_df = TRUE`) containing the feature summary data for the specified `layer`.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_feature_summary(metabolites, layer = "post_qc")
#'
#' @export
get_feature_summary <- new_generic("get_feature_summary", c("metabolites", "layer"), function(metabolites, layer, feature_ids = NULL, drop=FALSE, as_df = FALSE) { S7_dispatch() })
#' @name get_feature_summary
method(get_feature_summary, list(Metabolites, class_character)) <- function(metabolites, layer, feature_ids = NULL, drop=FALSE, as_df = FALSE) {

  return(get_(metabolites, slot="feature_summary", layer=layer, feature_ids=feature_ids, sample_ids=NULL, drop=drop, as_df=as_df))

}


#' @title Retrieve Sample Summary Data
#'
#' @description
#' Extracts sample-level summary statistics from a `Metabolites` object for a specified data `layer`,
#' with optional exclusions and subsetting.
#'
#' @param metabolites A `Metabolites` object containing sample summary data.
#' @param layer A character string specifying the data layer from which to extract sample summary data.
#' @param sample_ids Optional. A vector of column names (samples) to extract. Defaults to `NULL`, meaning all samples are included.
#' @param drop Logical. Whether to drop all NA (excluded samples)
#' @param as_df Logical. Whether to return the result as a `data.frame`. Defaults to `FALSE`, meaning the result is returned as a matrix.
#'
#' @return A matrix or `data.frame` (if `as_df = TRUE`) containing the sample summary data for the specified `layer`.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_sample_summary(metabolites, layer = "post_qc")
#'
#' @export
get_sample_summary <- new_generic("get_sample_summary", c("metabolites", "layer"), function(metabolites, layer, sample_ids = NULL, drop=FALSE, as_df = FALSE) { S7_dispatch() })
#' @name get_sample_summary
method(get_sample_summary, list(Metabolites, class_character)) <- function(metabolites, layer, sample_ids = NULL, drop=FALSE, as_df = FALSE) {

  return(get_(metabolites, slot="sample_summary", layer=layer, feature_ids=NULL, sample_ids=sample_ids, drop=drop, as_df=as_df))

}


#' @title Retrieve Data from a Metabolites Object
#'
#' @description
#' Extracts specified data from a `Metabolites` object for a given data layer.
#' This function allows for optional exclusions of data based on the current `@exclusions` matrix
#' and supports subsetting by feature and sample IDs.
#'
#' @param metabolites A `Metabolites` object containing the data to be extracted.
#' @param layer A character string specifying the data layer to extract (e.g., `"raw"`, `"post_qc"`, etc.).
#' @param feature_ids Optional. A vector of row names (features) to extract. Defaults to `NULL`, meaning all features are included.
#' @param sample_ids Optional. A vector of column names (samples) to extract. Defaults to `NULL`, meaning all samples are included.
#' @param drop Logical. Whether to drop all NA rows or columns (excluded samples or features)
#' @param as_df Logical. If `TRUE`, the result is returned as a `data.frame`; otherwise, the result is returned as a matrix. Default is `FALSE`.
#'
#' @return A matrix (or `data.frame`, if `as_df = TRUE`) containing the requested data for the specified `layer` and subsetting conditions.
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' metabolites <- metabolite_qc(metabolites, source="raw")
#' get_data(metabolites, layer = "post_qc")
#'
#' @export
get_data <- new_generic("get_data", c("metabolites", "layer"), function(metabolites, layer, feature_ids = NULL, sample_ids = NULL, drop=FALSE, as_df = FALSE) { S7_dispatch() })
#' @name get_data
method(get_data, list(Metabolites, class_character)) <- function(metabolites, layer, feature_ids = NULL, sample_ids = NULL, drop=FALSE, as_df = FALSE) {

  return(get_(metabolites, slot="data", layer=layer, feature_ids=feature_ids, sample_ids=sample_ids, drop=drop, as_df=as_df))

}




#' @title Export Data from a Metabolites Object
#'
#' @description
#' Exports all data layers from a `Metabolites` object to a structured directory format.
#' For each data layer, the function creates a subdirectory containing:
#' - the primary data matrix (`data.tsv`),
#' - associated feature and sample metadata (`features.tsv`, `samples.tsv`),
#' - feature and sample summaries (if present, `feature_summary.tsv`, `sample_summary.tsv`),
#' - a serialized feature tree (if present),
#' - and a `config.yml` file with additional metadata and processing parameters.
#'
#'
#' @param metabolites A `Metabolites` object containing the data to be exported.
#' @param directory A character string specifying the path to the directory where the data should be written.
#'
#' @importFrom data.table fwrite
#'
#' @return This function does not return a value. It writes the data to the specified directory as a file (e.g., CSV).
#'
#' @examples
#' file <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' metabolites <- Metabolites(project_name = "MYPROJECT",
#'                            format       = "metabolon_v1",
#'                            filepath     = file)
#' metabolites <- import_data(metabolites)
#' \dontrun{
#' # Specify the directory to export the data
#' export_data(metabolites, directory = "path/to/directory")
#' }
#'
#' @importFrom yaml write_yaml
#'
#' @export
export_data <- new_generic("export_data", c("metabolites", "directory"), function(metabolites, directory) { S7_dispatch() })
#' @name export_data
method(export_data, list(Metabolites, class_character)) <- function(metabolites, directory) {

  # make the directories
  today    <- gsub("-", "_", Sys.Date())
  dir      <- file.path(sub("/$", "", directory), paste0("metaboprep_release_", today))
  layers   <- stats::setNames(dimnames(metabolites@data)[[3]], dimnames(metabolites@data)[[3]])
  sub_dirs <- lapply(layers, function(x) file.path(dir, x))
  invisible(lapply(sub_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

  # combine and write the data for each layer
  for (j in seq_along(layers)) {

    layer <- layers[[j]]

    properties <- S7::props(metabolites)

    # deal with feature tree separately
    if (layer %in% names(properties$feature_tree)) {
      tree_path <- file.path(sub_dirs[[layer]], "feature_tree.RDS")
      saveRDS(properties$feature_tree[[layer]], tree_path)
      properties$feature_tree[[layer]] <- tree_path
    }

    # write the yaml config
    config <- list(layer = j)
    for (i in seq_along(properties)) {
      if (!inherits(properties[[i]], c("data.frame", "array", "matrix"))) {
        if (layer %in% names(properties[[i]])) {
          config[[names(properties)[i]]] <- properties[[i]][[layer]]
        } else if (!is.null(names(properties[[i]]))) {
          config[[names(properties)[i]]] <- NA_character_
        } else if (inherits(properties[[i]], c("POSIXct"))) {
          config[[names(properties)[i]]] <- as.character(properties[[i]])
        } else {
          config[[names(properties)[i]]] <- properties[[i]]
        }
      }
    }
    yaml::write_yaml(config, file.path(sub_dirs[[layer]], "config.yml"))

    # exclusions
    excl_feats <- unlist(metabolites@exclusions[[layer]][["features"]])
    excl_samps <- unlist(metabolites@exclusions[[layer]][["samples"]])

    # get the features
    features   <- get_features(metabolites, layer=layer)
    incl_feats <- features[!feature_id %in% excl_feats, feature_id]

    # get the samples
    samples    <- get_samples(metabolites, layer=layer)
    incl_samps <- samples[!sample_id %in% excl_samps, sample_id]

    # write the samples
    data.table::fwrite(samples[sample_id %in% incl_samps], file.path(sub_dirs[[layer]], "samples.tsv"), sep="\t")

    # write the features
    data.table::fwrite(features[feature_id %in% incl_feats], file.path(sub_dirs[[layer]], "features.tsv"), sep="\t")

    # write the feature summary if present
    if (layer %in% dimnames(metabolites@feature_summary)[[3]]) {
      feat_sum <- t(get_feature_summary(metabolites, layer=layer, feature_ids = incl_feats))
      feat_sum <- cbind(
        data.table::data.table(feature_id = rownames(feat_sum)),
        data.table::as.data.table(feat_sum)
      )
      data.table::fwrite(feat_sum, file.path(sub_dirs[[layer]], "feature_summary.tsv"), sep="\t")
    }

    # write the sample summary if present
    if (layer %in% dimnames(metabolites@sample_summary)[[3]]) {
      samp_sum <- get_sample_summary(metabolites, layer=layer, sample_ids = incl_samps)
      samp_sum <- cbind(
        data.table::data.table(sample_id = rownames(samp_sum)),
        data.table::as.data.table(samp_sum)
      )
      data.table::fwrite(samp_sum, file.path(sub_dirs[[layer]], "sample_summary.tsv"), sep="\t")
    }

    # write the data
    if (layer %in% dimnames(metabolites@data)[[3]]) {
      data <- get_data(metabolites, layer=layer, sample_ids = incl_samps, feature_ids = incl_feats)
      data <- cbind(
        data.table::data.table(sample_id = rownames(data)),
        data.table::as.data.table(data)
      )
      data.table::fwrite(data, file.path(sub_dirs[[layer]], "data.tsv"), sep="\t")
    }

    # write the pcs
    if (layer %in% dimnames(metabolites@data)[[3]]) {
      pcs <- get_pcs(metabolites, layer=layer, sample_ids = incl_samps)
      pcs <- cbind(
        data.table::data.table(sample_id = rownames(pcs)),
        data.table::as.data.table(pcs)
      )
      data.table::fwrite(pcs, file.path(sub_dirs[[layer]], "pcs.tsv"), sep="\t")
    }

    # write the prob pcs
    if (layer %in% dimnames(metabolites@data)[[3]]) {
      prob_pcs <- get_prob_pcs(metabolites, layer=layer, sample_ids = incl_samps)
      prob_pcs <- cbind(
        data.table::data.table(sample_id = rownames(prob_pcs)),
        data.table::as.data.table(prob_pcs)
      )
      data.table::fwrite(prob_pcs, file.path(sub_dirs[[layer]], "prob_pcs.tsv"), sep="\t")
    }

    # write the var explained
    if (layer %in% dimnames(metabolites@data)[[3]]) {
      var_exp <- get_var_exp(metabolites, layer=layer)
      var_exp <- data.table::data.table(pc    = names(var_exp),
                                        value = var_exp)
      data.table::fwrite(var_exp, file.path(sub_dirs[[layer]], "var_exp.tsv"), sep="\t")
    }


    cli::cli_alert_success("Exported `{layer}` data layer to directory {.file {sub_dirs[[layer]]}}")
  }

  invisible(NULL)
}
