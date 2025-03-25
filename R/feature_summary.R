#' @title Feature Summary Statistics
#' @description
#' This function estimates feature statistics for samples in a matrix of metabolite features.
#' @param metabolites an object of class Metabolites
#' @include class_metabolites.R
#' @export
feature_summary <- new_generic("feature_summary", c("metabolites"), function(metabolites, type="raw") { S7_dispatch() })
#' @name feature_summary
method(feature_summary, Metabolites) <- function(metabolites, type="raw"){

  ## feature missingness
  featuremis = missingness(metabolites@data[, , type], by="column", exclude_features = NA)

  ### distribution discritions
  description = data.table::as.data.table(feature_describe(metabolites@data[, , type]), keep.rownames = TRUE)
  data.table::setnames(description, "rn", "feature_id")

  ### count of sample outliers per feature
  omat     <- outlier_detection(metabolites@data[, , type], nsd = metabolites@outlier_udist, meansd = FALSE, by = "column")
  sumomat  <- apply(omat, 2, sum)
  outliers <- data.table::data.table(feature_id    = names(sumomat),
                                     outlier_count = sumomat)

  ### identify independent features
  # if(50 > nrow(wdata)*0.8 ){
  #   MSS = nrow(wdata) * 0.8  ## this allows 20% missingness on data sets with less than 50 individuals
  # } else {
  #   MSS = 50
  # }
  ## ** ALL FEATURES THAT GO INTO THE DENDROGRAM AND CAN BE REPRESENTITIVE
  ##    FEATURES MUST HAVE > 80% Presence or <= 20% missing
  min_sample_size = floor( nrow(metabolites@data[, , type]) * 0.8 )

  # features to exclude
  exclude_features <- character()
  if (metabolites@derived_var_exclusion) {
    exclude_features <- metabolites@features[derived_feature==TRUE, feature_names] # derived
  }

  ##
  indf = tree_and_independent_features(data,
                                       minimum_samplesize = min_sample_size,
                                       tree_cut_height    = metabolites@tree_cut_height,
                                       exclude_features   = exclude_features )

  dt_list <- list(featuremis, outliers, description, indf$dat)
  out <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), dt_list)

  return(  list( table = out, tree = indf$speartree ) )
}

