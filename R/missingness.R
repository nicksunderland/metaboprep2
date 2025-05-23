# Silence R CMD check
globalVariables(c("missingness_w_exclusions", "i.missingness_w_exclusions"), package = "metaboprep2")

#' @title Estimate Missingness
#' @description
#' This function estimates missingness in a matrix of data and provides an option to exclude certain columns or features from the analysis, such as xenobiotics (with high missingness rates) in metabolomics data sets.
#' @param data matrix, a numeric matrix with samples in rows and features in columns
#' @param exclude_features character, a vector of feature names (i.e. column names) to exclude from sample (i.e. by="rows") missingness estimates
#' @param by character, whether to calculate missingness by rows (samples) or column (features)
#' @return data.frame, a data frame of missingness estimates for each sample/feature. If a vector of feature names was also passed to the function a second column of missingness estimates will also be returned providing sample missingness estimates to the exclusion of those features provided.
#'
#' @export
#'
#' @examples
#' ## simulate some data
#' set.seed(1110)
#' ex_data = sapply(1:5, function(x){ rnorm(10, 40, 5) })
#' rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
#' colnames(ex_data) = paste0("var", 1:ncol(ex_data))
#' ## add some missingness to the data
#' ex_data[ sample(1:50, 10) ] = NA
#' ## estimate missingness
#' mis_est = missingness(ex_data, by="row")
#' mis_est_v2 = missingness(ex_data, exclude_features = "var5", by="row")
#'
missingness <- function(data, exclude_features = NA, by = "row") {

  by <- match.arg(by, choices = c("row", "column"))

  if (by == "row") {
    margin <- 1
    id_col <- "sample_id"
  } else if (by == "column") {
    margin <- 2
    id_col <- "feature_id"
    stopifnot("`exclude_features` is only relevated when calculating missingness for samples, set NA for feature missingness" =  is.na(exclude_features))
  }

  miss <- apply(data, margin, function(x) { sum(is.na(x)) / length(x) })
  out <- data.table::data.table(ids         = names(miss), #ifelse(is.null(names(miss)), character(), names(miss)),
                                missingness = miss,
                                missingness_w_exclusions = NA_real_)

  if( !all(is.na(exclude_features)) || length(exclude_features)>0 ) {

    r = which(colnames(data) %in% exclude_features)

    if (length(r) > 0) {

      miss2 <- apply(data[, -r, drop = FALSE], margin, function(x) { sum(is.na(x)) / length(x) })

      out2  <- data.table::data.table(ids = names(miss2),
                                      missingness_w_exclusions = miss2)

      out[out2, missingness_w_exclusions := i.missingness_w_exclusions, on="ids"]

    }

  }

  data.table::setnames(out, "ids", id_col)

  return(out)
}
