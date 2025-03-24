#' @title Sample Summary Statistics
#' @param metabolites an object of class Metabolites
#' @include class_metabolites.R
#' @export
sample_summary <- new_generic("sample_summary", c("metabolites"), function(metabolites) { S7_dispatch() })
#' @name sample_summary
method(sample_summary, Metabolites) <- function(metabolites) {

  # features to exclude
  exclude_features <- character()
  if (metabolites@derived_var_exclusion) {
    exclude_features <- c(exclude_features, metabolites@features[derived_feature==TRUE, feature_names]) # derived
  }
  if (metabolites@derived_var_exclusion) {
    exclude_features <- c(exclude_features, metabolites@features[grepl("(?i)xenobiotic", pathway), feature_names]) # xenobiotics
  }

  # missingness
  missing <- missingness(metabolites@data[, , "raw"], by="row", exclude_features = exclude_features)

  # total peak area
  tpa <- total_peak_area(metabolites@data[, , "raw"])

  # count sample outliers
  omat     <- outlier_detection(metabolites@data[, , "raw"], nsd = metabolites@outlier_udist, meansd = FALSE, by = "column")
  outliers <- data.frame(outlier_count = apply(omat, 1, sum))

  # combine
  output <- cbind(missing, tpa, outliers)
  sample_id <- rownames(output)
  output = cbind(sample_id, output)

  return(data.table::as.data.table(output))
}


