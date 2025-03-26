#' @title Metabolite Quality Control
#' @description
#' This function is a wrapper function that performs the key quality controls steps on a metabolomics data set.
#' Key principles:
#'  1. keep the source underlying data as it is
#'  2. copy the source data to a new data layer called post_qc (the destination)
#'  3. build an exclusion destination matrix mask, accumulating codes for exclusion reasons
#'  4. make any adjustments needed in the destination copy of the data, flag these in the exclusion matrix
#'  5. return the Metabolites object with the newly populated data layers
#'
#' @param metabolites an object of class Metabolites
#' @param source character
#' @param destination character
#' @include class_metabolites.R
#' @importFrom stats quantile
#' @importFrom glue glue
#' @import cli
#' @export
metabolite_qc <- new_generic("metabolite_qc", c("metabolites"), function(metabolites, source="raw", destination="post_qc") { S7_dispatch() })
#' @name metabolite_qc
method(metabolite_qc, Metabolites) <- function(metabolites, source="raw", destination="post_qc"){

  cli::cli_h1("Starting Metabolite QC Process")



  # input validation
  cli::cli_alert_info("Validating input parameters...")
  source <- match.arg(source, choices = dimnames(metabolites@data)[[3]])
  stopifnot("`destination` parameter is fixed to `post_qc` and only present to help document where the QCd data will land" = destination=="post_qc")
  cli::cli_alert_success("Input validation complete.")



  # data to work from and add another layer - we now work with the 'destination data', plus any exclusions, from now on
  cli::cli_alert_info(glue::glue("Copying {source} data to new {destination} data layer..."))
  dat <- get_data(metabolites, type = source, apply_exclusions = FALSE)
  metabolites@data <- add_layer(current    = metabolites@data,
                                layer      = dat,
                                layer_name = destination)
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
  cli::cli_alert_success(glue::glue("{destination} data ready for QC."))



  # derived features to exclude
  cli::cli_alert_info("Analysing exclusions...")
  derived_feature_excl <- character()
  if (metabolites@derived_var_exclusion) {
    cli::cli_alert_info(glue::glue("Applying exclusion of derived features..."))
    derived_feature_excl <- metabolites@features[derived_feature==TRUE, feature_id]
    metabolites <- update_exclusions(metabolites,
                                     type        = destination,
                                     code        = get_exclusion_codes(verbose=FALSE)[["derived_feature"]],
                                     feature_ids = derived_feature_excl)
  } else {
    cli::cli_alert_info(glue::glue("Skipping exclusion of derived features..."))
  }
  cli::cli_alert_success(glue::glue("Derived features exclusions assessment complete. ({length(derived_feature_excl)} features excluded)"))



  # xenobiotic features to exclude
  xenobiotic_feature_excl <- character()
  if (metabolites@xenobiotics_var_exclusion) {
    cli::cli_alert_info(glue::glue("Applying exclusion of xenobiotic features..."))
    xenobiotic_feature_excl <- metabolites@features[grepl("(?i)xenobiotic", pathway), feature_id]
    metabolites <- update_exclusions(metabolites,
                                     type        = destination,
                                     code        = get_exclusion_codes(verbose=FALSE)[["xenobiotic_feature"]],
                                     feature_ids = xenobiotic_feature_excl)
  } else {
    cli::cli_alert_info(glue::glue("Skipping exclusion of xenobiotic features..."))
  }
  cli::cli_alert_success(glue::glue("Xenobiotic features exclusions assessment complete. ({length(xenobiotic_feature_excl)} features excluded)"))



  # very bad sample missingness
  cli::cli_alert_info(glue::glue("Assessing for extreme sample missingness >=80%..."))
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
  cli::cli_alert_success(glue::glue("Extreme sample missingness exclusions assessment complete. ({length(sample_ids)} samples excluded)"))



  # very bad feature missingness
  cli::cli_alert_info(glue::glue("Assessing for extreme feature missingness >=80%..."))
  dat <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  featuremis  <- missingness(dat, by="column")
  feature_ids <- featuremis[missingness >= 0.8, feature_id]
  metabolites <- update_exclusions(metabolites,
                                   type        = destination,
                                   code        = get_exclusion_codes(verbose=FALSE)[["extreme_feature_missingness"]],
                                   feature_ids = feature_ids)
  cli::cli_alert_success(glue::glue("Extreme feature missingness exclusions assessment complete. ({length(feature_ids)} features excluded})"))



  # re-estimate sample missingness
  cli::cli_alert_info(glue::glue("Assessing for sample missingness at specified level of >={round(metabolites@sample_missingness*100)}%..."))
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
  cli::cli_alert_success(glue::glue("Sample missingness exclusions assessment complete. ({length(sample_ids)} samples excluded)"))



  # re-estimate feature missingness
  cli::cli_alert_info(glue::glue("Assessing for feature missingness at specified level of >={round(metabolites@feature_missingness*100)}%..."))
  dat         <- get_data(metabolites, type = destination, apply_exclusions = TRUE)
  featuremis  <- missingness(dat, by="column")
  feature_ids <- featuremis[missingness >= metabolites@feature_missingness, feature_id]
  metabolites <- update_exclusions(metabolites,
                                   type        = destination,
                                   code        = get_exclusion_codes(verbose=FALSE)[["user_defined_feature_missingness"]],
                                   feature_ids = feature_ids)
  cli::cli_alert_success(glue::glue("Feature missingness exclusions assessment complete. ({length(feature_ids)} features excluded)"))



  # total peak area
  cli::cli_alert_info(glue::glue("Calculating total peak abundance outliers at +/- {metabolites@total_peak_area_sd} Sdev..."))
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
  cli::cli_alert_success(glue::glue("Total peak abundance outlier exclusions assessment complete. ({length(sample_ids)} samples excluded)"))



  # PCA data
  cli::cli_alert_info(glue::glue("Running principal component outlier analysis at +/- {metabolites@outlier_udist} Sdev..."))
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
  cli::cli_alert_info(glue::glue("PCA outlier analysis complete."))



  # re-identify feature independence and PC outliers (feature_summary will use the current exclusions internally)
  cli::cli_alert_info(glue::glue("Re-identify feature independence and PC outliers..."))
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
  cli::cli_alert_info(glue::glue("PCA outlier analysis 2 complete."))



  cli::cli_h1("Metabolite QC Process Completed")

  # return the metabolites with underlying data (+/- adjustments) and exclusion matrix
  return(metabolites)
}
