#' @title Sample Summary Statistics
#' @param metabolites an object of class Metabolites
#' @include class_metabolites.R
#' @export
sample_summary <- new_generic("sample_summary", c("metabolites"), function(metabolites, type="raw") { S7_dispatch() })
#' @name sample_summary
method(sample_summary, Metabolites) <- function(metabolites, type="raw") {

  # features to exclude
  exclude_features <- character()
  if (metabolites@derived_var_exclusion) {
    exclude_features <- c(exclude_features, metabolites@features[derived_feature==TRUE, feature_names]) # derived
  }
  if (metabolites@derived_var_exclusion) {
    exclude_features <- c(exclude_features, metabolites@features[grepl("(?i)xenobiotic", pathway), feature_names]) # xenobiotics
  }

  # missingness
  missing <- missingness(metabolites@data[, , type], by="row", exclude_features = exclude_features)

  # total peak area
  tpa <- total_peak_area(metabolites@data[, , type])

  # count sample outliers
  omat     <- outlier_detection(metabolites@data[, , type], nsd = metabolites@outlier_udist, meansd = FALSE, by = "column")
  sumomat  <- apply(omat, 1, sum)
  outliers <- data.table::data.table(sample_id = names(sumomat),
                                     outlier_count = sumomat)

  # combine
  dt_list <- list(missing, tpa, outliers)
  out <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), dt_list)

  return(out)
}


