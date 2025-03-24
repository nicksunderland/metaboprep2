#' @title Batch Normalisation
#' @param metabolites an object of class Metabolites
#' @include class_metabolites.R
#' @importFrom stats median
#' @export
batch_normalisation <- new_generic("batch_normalisation", c("metabolites"), function(metabolites) { S7_dispatch() })
#' @name batch_normalisation
method(batch_normalisation, Metabolites) <- function(metabolites) {

  # add a copy of the raw data to the back of the matrix stack
  metabolites@data <- array(
    c(metabolites@data, metabolites@data[, , "raw"]),
    dim = c(dim(metabolites@data)[1], dim(metabolites@data)[2], dim(metabolites@data)[3] + 1),
    dimnames = list(
      dimnames(metabolites@data)[[1]],
      dimnames(metabolites@data)[[2]],
      c(dimnames(metabolites@data)[[3]], "batch_normalised")
    )
  )

  platform_ids <- sort(unique(metabolites@features[["platform"]]))

  for (plat in platform_ids) {

    metabolite_index = which(metabolites@features[, platform] == plat)

    batch_ids = unique(metabolites@samples[, get(plat)])

    for (bid in batch_ids) {

      sample_index = which(metabolites@samples[, get(plat)] == bid)

      m = apply(metabolites@data[sample_index, metabolite_index, "raw"],
                2, function(x) {
                  stats::median(x, na.rm = TRUE)
                })

      for (j in 1:length(m)) {

        metabolites@data[sample_index,
                         metabolite_index[j],
                         "batch_normalised"] = metabolites@data[sample_index,
                                                                metabolite_index[j],
                                                                "batch_normalised"]/m[j]

      }
    }
  }

  return(metabolites)
}
