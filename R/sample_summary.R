#' @title Sample Summary Statistics
#' @param metabolites an object of class Metabolites
#' @param layer character
#' @param output character
#' @include class_metabolites.R
#' @export
sample_summary <- new_generic("sample_summary", c("metabolites"), function(metabolites, layer="raw", output="object") { S7_dispatch() })
#' @name sample_summary
method(sample_summary, Metabolites) <- function(metabolites, layer="raw", output="object") {

  # the data to work with
  dat <- get_data(metabolites, layer=layer, apply_exclusions=TRUE)

  # features to exclude
  exclude_features <- character()
  if (metabolites@derived_var_exclusion) {
    exclude_features <- c(exclude_features, metabolites@features[derived_feature==TRUE, feature_id]) # derived
  }
  if (metabolites@derived_var_exclusion) {
    exclude_features <- c(exclude_features, metabolites@features[grepl("(?i)xenobiotic", pathway), feature_id]) # xenobiotics
  }

  # missingness
  missing <- missingness(dat, by="row", exclude_features = exclude_features)

  # total peak area
  tpa <- total_peak_area(dat)

  # count sample outliers
  omat     <- outlier_detection(dat, nsd = metabolites@outlier_udist, meansd = FALSE, by = "column")
  sumomat  <- apply(omat, 1, sum)
  outliers <- data.table::data.table(sample_id = names(sumomat),
                                     outlier_count = sumomat)

  # combine
  dt_list <- list(missing, tpa, outliers)
  out <- Reduce(function(x, y) data.table::merge.data.table(x, y, by = "sample_id", all = TRUE), dt_list)

  # ensure correct order
  out <- out[order(match(sample_id, rownames(metabolites@data[, , layer])))]

  # ensure correct order (use the unfiltered data to get the rowname names, inject NAs if absent from filtered data)
  ordered_base_ids <- data.table::data.table(sample_id = rownames(metabolites@data[, , layer]))
  out <- out[ordered_base_ids, on="sample_id", nomatch = NA]

  # as matrix
  mat <- as.matrix(out[, !("sample_id")])
  rownames(mat) <- out$sample_id

  # add to sample_summary matrix
  metabolites@sample_summary <- add_layer(current    = metabolites@sample_summary,
                                          layer      = mat,
                                          layer_name = layer)


  # return desired output
  return(
    switch(output,
           "object"     = metabolites,
           "data.table" = out,
           "matrix"     = mat
    )
  )
}


