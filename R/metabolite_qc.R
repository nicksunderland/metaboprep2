#' @title Metabolite Quality Control
#' @description
#' This function is a wrapper function that performs the key quality controls steps on a metabolomics data set.
#' @param metabolites an object of class Metabolites
#' @param source character
#' @param destination character
#' @include class_metabolites.R
#' @importFrom stats quantile
#' @export
metabolite_qc <- new_generic("metabolite_qc", c("metabolites"), function(metabolites, source="raw", destination="post_qc") { S7_dispatch() })
#' @name metabolite_qc
method(metabolite_qc, Metabolites) <- function(metabolites, source="raw", destination="post_qc"){

  source <- match.arg(source, choices = dimnames(metabolites@data)[[3]])
  stopifnot("`destination` parameter is fixed to `post_qc` and only present to help document where the QCd data will land" = destination=="post_qc")

  # principles:
  # 1. keep the underlying data as it is
  # 2. build an exclusion matrix mask, accumulating codes for exclusion reasons
  # -> analysis -> ids of QC fails -> add to exclusion matrix -> apply exclusions to copy of data -> more analysis

  # data to work from and add another layer - we now work with the 'destination data', plus any exclusions, from now on
  dat <- get_data(metabolites, type = source, apply_exclusions = FALSE)
  metabolites@data <- add_layer(current    = metabolites@data,
                                layer      = dat,
                                layer_name = destination)

  # add the source of this new layer
  metabolites@source <- c(metabolites@source, stats::setNames(source, destination))

  # add an exclusions sparse matrix, copy the source exclusions if it exists
  if (source %in% dimnames(metabolites@exclusions)[[3]]) {
    excl_mat <- metabolites@exclusions[, , source, drop=FALSE]
  } else {
    excl_mat <- SparseArray::SparseArray(array("",
                                               dim = dim(metabolites@data[, , destination]),
                                               dimnames = dimnames(metabolites@data[, , destination])[1:2]),
                                         type="character")
  }
  metabolites@exclusions <- add_layer(current    = metabolites@exclusions,
                                      layer      = excl_mat,
                                      layer_name = destination)


  # derived features to exclude
  derived_feature_excl <- character()
  if (metabolites@derived_var_exclusion) {
    derived_feature_excl <- metabolites@features[derived_feature==TRUE, feature_names]
    metabolites <- update_exclusions(metabolites,
                                     type        = destination,
                                     code        = get_exclusion_codes(verbose=FALSE)[["derived_feature"]],
                                     feature_ids = derived_feature_excl)
  }

  # xenobiotic features to exclude
  if (metabolites@derived_var_exclusion) {
    xenobiotic_feature_excl <- metabolites@features[grepl("(?i)xenobiotic", pathway), feature_names]
    metabolites <- update_exclusions(metabolites,
                                     type        = destination,
                                     code        = get_exclusion_codes(verbose=FALSE)[["xenobiotic_feature"]],
                                     feature_ids = xenobiotic_feature_excl)
  }

  # very bad sample missingness
  dat <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  samplemis = missingness(dat, by="row", exclude_features = derived_feature_excl)
  if (!all(is.na(samplemis$missingness_w_exclusions))) {
    sample_ids <- samplemis[missingness_w_exclusions >= 0.8, sample_id]
  } else {
    sample_ids <- samplemis[missingness >= 0.8, sample_id]
  }
  metabolites <- update_exclusions(metabolites,
                                   type       = destination,
                                   code       = get_exclusion_codes(verbose=FALSE)[["extreme_sample_missingness"]],
                                   sample_ids = sample_ids)

  # very bad feature missingness
  dat <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  featuremis  <- missingness(dat, by="column")
  feature_ids <- featuremis[missingness >= 0.8, feature_id]
  metabolites <- update_exclusions(metabolites,
                                   type        = destination,
                                   code        = get_exclusion_codes(verbose=FALSE)[["extreme_feature_missingness"]],
                                   feature_ids = feature_ids)

  # re-estimate sample missingness
  dat <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  samplemis <- missingness(dat, by="row", exclude_features = derived_feature_excl)
  if (!all(is.na(samplemis$missingness_w_exclusions))) {
    sample_ids <- samplemis[missingness_w_exclusions >= metabolites@sample_missingness, sample_id]
  } else {
    sample_ids <- samplemis[missingness >= metabolites@sample_missingness, sample_id]
  }
  metabolites <- update_exclusions(metabolites,
                                   type       = destination,
                                   code       = get_exclusion_codes(verbose=FALSE)[["user_defined_sample_missingness"]],
                                   sample_ids = sample_ids)

  # re-estimate feature missingness
  dat         <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  featuremis  <- missingness(dat, by="column")
  feature_ids <- featuremis[missingness >= metabolites@feature_missingness, feature_id]
  metabolites <- update_exclusions(metabolites,
                                   type        = destination,
                                   code        = get_exclusion_codes(verbose=FALSE)[["user_defined_feature_missingness"]],
                                   feature_ids = feature_ids)

  # total peak area
  dat <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  tpa <- total_peak_area(dat, features_exclude = derived_feature_excl)
  tpa[, `:=`(sdev = sd(tpa_total),
             mean = mean(tpa_total))]
  tpa[, `:=`(UL   = mean + sdev * metabolites@total_peak_area_sd,
             LL   = mean - sdev * metabolites@total_peak_area_sd)]
  sample_ids <- tpa[!data.table::between(tpa_total, LL, UL), sample_id]
  metabolites   <- update_exclusions(metabolites,
                                     type       = destination,
                                     code       = get_exclusion_codes(verbose=FALSE)[["user_defined_sample_totalpeakarea"]],
                                     sample_ids = sample_ids)

  # PCA data
  dat <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  if (metabolites@outlier_treatment != "leave_be") {

    omat <- outlier_detection(dat, nsd = metabolites@outlier_udist, meansd = FALSE, by="column")

    # which samples and features to set NA
    indices <- which(omat == 1, arr.ind = TRUE)
    adjust_row_names <- rownames(omat)[indices[,1]]
    adjust_col_names <- colnames(omat)[indices[,2]]

    if(metabolites@outlier_treatment == "turn_NA") {

      # turn NA in actual destination data
      metabolites@data[adjust_row_names, adjust_col_names, destination] <- NA_real_

      # log in exclusion matrix
      metabolites   <- update_exclusions(metabolites,
                                         type        = destination,
                                         code        = get_exclusion_codes(verbose=FALSE)[["outlier_udist_turned_na"]],
                                         sample_ids  = adjust_row_names,
                                         feature_ids = adjust_col_names,
                                         arr_ids     = TRUE)

    } else if (outlier_treatment == "winsorize") {

      for(feature_id in adjust_col_names) {

        # quantile value of the feature, minus the outliers
        quantile_value = quantile(dat[!rownames(dat) %in% adjust_row_names, feature_id], probs = c(metabolites@winsorize_quantile), na.rm = TRUE)

        # set in the destination data
        metabolites@data[adjust_row_names, feature_id, destination] <- quantile_value

        # log in the exclusion matrix
        metabolites   <- update_exclusions(metabolites,
                                           type        = destination,
                                           code        = get_exclusion_codes(verbose=FALSE)[["outlier_udist_winsorized"]],
                                           sample_ids  = adjust_row_names,
                                           feature_ids = adjust_col_names,
                                           arr_ids     = TRUE)
      }
    }
  }#end dealing with PCA based adjustments

  # re-identify feature independence and PC outliers (feature_summary will use the current exclusions internally)
  metabolites <- feature_summary(metabolites, type = destination)

  # identify PC outliers using the newly populated @feature_summary data
  metabolites <- pc_and_outliers(metabolites, type = destination)

  # extract PCs 1-2 or 1-number of Acceleration factor PCs
  af <- metabolites@acceleration_factor[[destination]]
  if(af < 2) {
    pcs = metabolites@pcs[, 1:2, destination]
  } else {
    pcs = metabolites@pcs[, 1:af, destination]
  }

  # perform exclusion on top PCs to ID outliers
  if (!is.na(metabolites@pc_outlier_sd)) {

    outliers <- outlier_detection(pcs, nsd = metabolites@pc_outlier_sd, meansd = TRUE)
    indices <- which(outliers == 1, arr.ind = TRUE)
    excl_col_names <- unique(colnames(outliers)[indices[,2]])

    # log in exclusion matrix
    metabolites <- update_exclusions(metabolites,
                                     type        = destination,
                                     code        = get_exclusion_codes(verbose=FALSE)[["user_defined_sample_pca_outlier"]],
                                     sample_ids  = excl_col_names)
  }

  # return the metabolites with underlying data (+/- adjustments) and exclusion matrix
  return(metabolites)
}
