#' @title Metabolite Quality Control
#' @description
#' This function is a wrapper function that performs the key quality controls steps on a metabolomics data set.
#' Key principles:
#'  1. keep the source underlying data as it is
#'  2. copy the source data to a new data layer called qcing for processing
#'  3. build an exclusion list, accumulating codes for exclusion reasons
#'  4. make any adjustments needed in the destination copy of the data, flag these in the exclusion list
#'  5. copy the final result to a data layer called post_qc
#'  6. return the Metabolites object with the newly populated data layers
#'
#' @param metabolites an object of class Metabolites
#' @param source character
#' @include class_metabolites.R
#' @importFrom stats quantile
#' @importFrom glue glue
#' @import cli
#' @export
metabolite_qc <- new_generic("metabolite_qc", c("metabolites"), function(metabolites, source="raw") { S7_dispatch() })
#' @name metabolite_qc
method(metabolite_qc, Metabolites) <- function(metabolites, source="raw"){

  cli::cli_h1("Starting Metabolite QC Process")


  # input validation
  cli::cli_alert_info("Validating input parameters...")
  source <- match.arg(source, choices = dimnames(metabolites@data)[[3]])
  cli::cli_alert_success("Input validation complete.")


  # features and samples included in this layer
  incl_samps <- setdiff(get_samples(metabolites, layer=source) [, sample_id], unlist(metabolites@exclusions[[source]][["samples"]]))
  incl_feats <- setdiff(get_features(metabolites, layer=source)[, feature_id], unlist(metabolites@exclusions[[source]][["features"]]))

  # run normalisation
  if (metabolites@batch_normalise==TRUE &&
      !all(is.na(get_features(metabolites, layer=source)[, platform])) && # must have platform specified
      all(unique(get_features(metabolites, layer=source)[, platform]) %in% names(get_samples(metabolites, layer=source)))) { # must have the batches specified
    cli::cli_alert_info(glue::glue("Batch normalising raw data..."))
    metabolites <- batch_normalisation(metabolites)
    source <- "batch_normalised"
    cli::cli_alert_success(glue::glue("Normalised raw data."))
  }

  # xenobiotic features to exclude
  xenobiotic_feats <- character()
  if (metabolites@xenobiotics_var_exclusion) {
    cli::cli_alert_info(glue::glue("Applying exclusion of xenobiotic features..."))
    xenobiotic_feats <- get_features(metabolites, "raw")[grepl("(?i)xenobiotic", pathway), feature_id]
  } else {
    cli::cli_alert_info(glue::glue("Skipping exclusion of xenobiotic features..."))
  }
  cli::cli_alert_success(glue::glue("Xenobiotic features exclusions assessment complete. ({length(xenobiotic_feats)} features excluded)"))

  # features to include
  cli::cli_alert_info("Analysing exclusions...")
  derived_feats <- character()
  incl_feats    <- setdiff(incl_feats, unlist(metabolites@exclusions[[source]][["features"]]))
  if (metabolites@derived_var_exclusion) {
    cli::cli_alert_info(glue::glue("Applying exclusion of derived features..."))
    derived_feats <- get_features(metabolites, "raw")[derived_feature==TRUE, feature_id]
  } else {
    cli::cli_alert_info(glue::glue("Skipping exclusion of derived features..."))
  }
  cli::cli_alert_success(glue::glue("Derived features exclusions assessment complete. ({length(derived_feats)} features excluded)"))


  # run sample summary and feature summary on the raw data
  cli::cli_alert_info(glue::glue("Sample & Feature Summary Statistics for raw data..."))
  metabolites <- sample_summary( metabolites, layer=source, sample_ids=incl_samps, feature_ids=incl_feats, output="object")
  metabolites <- feature_summary(metabolites, layer=source, sample_ids=incl_samps, feature_ids=incl_feats, output="object")
  metabolites <- pc_and_outliers(metabolites, layer=source, sample_ids=incl_samps, feature_ids=incl_feats)
  cli::cli_alert_success(glue::glue("Sample & Feature Summary Statistics for raw data..."))


  # data to work from and add another layer - we now work with the 'destination data', plus any exclusions, from now on
  cli::cli_alert_info(glue::glue("Copying {source} data to new 'qcing' data layer..."))
  dat <- get_data(metabolites, layer=source, sample_ids=incl_samps, feature_ids=incl_feats)
  metabolites@data <- add_layer(current    = metabolites@data,
                                layer      = dat,
                                layer_name = "qcing")
  metabolites@source <- c(metabolites@source, stats::setNames(source, "qcing"))


  # add an exclusions list, copy the source exclusions if it exists
  if (source %in% names(metabolites@exclusions)) {
    metabolites@exclusions[["qcing"]] <- metabolites@exclusions[[source]]
  } else {
    metabolites@exclusions[["qcing"]] <- metabolites@exclusions[["raw"]]
  }
  cli::cli_alert_success(glue::glue("qcing data ready for QC."))




  # very bad sample missingness
  cli::cli_alert_info(glue::glue("Assessing for extreme sample missingness >=80%..."))
  dat        <- get_data(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, c(xenobiotic_feats, derived_feats)))
  samplemis  <- missingness(dat, by="row")
  excl_samps <- samplemis[missingness >= 0.8, sample_id]
  metabolites@exclusions[["qcing"]][["samples"]][["extreme_sample_missingness"]] <- unique(c(metabolites@exclusions[["qcing"]][["samples"]][["extreme_sample_missingness"]], excl_samps))
  incl_samps <- setdiff(incl_samps, excl_samps)
  cli::cli_alert_success(glue::glue("Extreme sample missingness exclusions assessment complete. ({length(excl_samps)} samples excluded)"))


  # very bad feature missingness
  cli::cli_alert_info(glue::glue("Assessing for extreme feature missingness >=80%..."))
  dat        <- get_data(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, xenobiotic_feats))
  featuremis <- missingness(dat, by="column")
  excl_feats <- featuremis[missingness >= 0.8, feature_id]
  metabolites@exclusions[["qcing"]][["features"]][["extreme_feature_missingness"]] <- unique(c(metabolites@exclusions[["qcing"]][["features"]][["extreme_feature_missingness"]], excl_feats))
  incl_feats <- setdiff(incl_feats, excl_feats)
  cli::cli_alert_success(glue::glue("Extreme feature missingness exclusions assessment complete. ({length(excl_feats)} features excluded)"))



  # re-estimate sample missingness and exclude based on user-defined
  cli::cli_alert_info(glue::glue("Assessing for sample missingness at specified level of >={round(metabolites@sample_missingness*100)}%..."))
  dat        <- get_data(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, c(xenobiotic_feats, derived_feats)))
  samplemis  <- missingness(dat, by="row")
  excl_samps <- samplemis[missingness >= metabolites@sample_missingness, sample_id]
  metabolites@exclusions[["qcing"]][["samples"]][["user_defined_sample_missingness"]] <- unique(c(metabolites@exclusions[["qcing"]][["samples"]][["user_defined_sample_missingness"]], excl_samps))
  incl_samps <- setdiff(incl_samps, excl_samps)
  cli::cli_alert_success(glue::glue("Sample missingness exclusions assessment complete. ({length(excl_samps)} samples excluded)"))



  # re-estimate feature missingness
  cli::cli_alert_info(glue::glue("Assessing for feature missingness at specified level of >={round(metabolites@feature_missingness*100)}%..."))
  dat        <- get_data(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, xenobiotic_feats))
  featuremis <- missingness(dat, by="column")
  excl_feats <- featuremis[missingness >= metabolites@feature_missingness, feature_id]
  metabolites@exclusions[["qcing"]][["features"]][["user_defined_feature_missingness"]] <- unique(c(metabolites@exclusions[["qcing"]][["features"]][["user_defined_feature_missingness"]], excl_feats))
  incl_feats <- setdiff(incl_feats, excl_feats)
  cli::cli_alert_success(glue::glue("Feature missingness exclusions assessment complete. ({length(excl_feats)} features excluded)"))



  # total peak area
  cli::cli_alert_info(glue::glue("Calculating total peak abundance outliers at +/- {metabolites@total_peak_area_sd} Sdev..."))
  dat <- get_data(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, c(xenobiotic_feats, derived_feats)))
  tpa <- total_peak_area(dat)
  tpa[, `:=`(sdev = sd(tpa_total),
             mean = mean(tpa_total))]
  tpa[, `:=`(UL   = mean + sdev * metabolites@total_peak_area_sd,
             LL   = mean - sdev * metabolites@total_peak_area_sd)]
  excl_samps <- tpa[!data.table::between(tpa_total, LL, UL), sample_id]
  metabolites@exclusions[["qcing"]][["samples"]][["user_defined_sample_totalpeakarea"]] <- unique(c(metabolites@exclusions[["qcing"]][["samples"]][["user_defined_sample_totalpeakarea"]], excl_samps))
  incl_samps <- setdiff(incl_samps, excl_samps)
  cli::cli_alert_success(glue::glue("Total peak abundance outlier exclusions assessment complete. ({length(excl_samps)} samples excluded)"))



  # PCA data
  # first deal with outlier data depending on option selected
  cli::cli_alert_info(glue::glue("Running data outlier analysis at +/- {metabolites@outlier_udist} Sdev..."))
  dat <- get_data(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, xenobiotic_feats))
  if (metabolites@outlier_treatment != "leave_be") {

    omat <- outlier_detection(dat, nsd = metabolites@outlier_udist, meansd = FALSE, by="column")

    # which samples and features to set NA
    indices <- which(omat == 1, arr.ind = TRUE)
    adjust_samps <- rownames(omat)[indices[,1]]
    adjust_feats <- colnames(omat)[indices[,2]]

    if(metabolites@outlier_treatment == "turn_NA") {

      # turn NA in actual destination data
      metabolites@data[adjust_row_names, adjust_col_names, "qcing"] <- NA_real_
      cli::cli_alert_info(glue::glue("All identified outliers were turned into NA..."))

    } else if (outlier_treatment == "winsorize") {

      for(feature_id in adjust_col_names) {

        # quantile value of the feature, minus the outliers
        quantile_value = quantile(dat[!rownames(dat) %in% adjust_row_names, feature_id], probs = c(metabolites@winsorize_quantile), na.rm = TRUE)

        # set in the destination data
        metabolites@data[adjust_row_names, feature_id, "qcing"] <- quantile_value
        cli::cli_alert_info(glue::glue("Outliers were winsorized to the {winsorize_quantile * 100} quantile of remaining (non outlying) values."))
      }
    }
  }#end dealing with adjustments pre PCA run
  cli::cli_alert_info(glue::glue("Data outlier analysis complete."))

  # re-identify feature independence and PC outliers (feature_summary will use the current exclusions internally)
  cli::cli_alert_info(glue::glue("Re-identify feature independence and PC outliers..."))
  metabolites <- feature_summary(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, c(xenobiotic_feats, derived_feats)))

  # identify PC outliers using the newly populated @feature_summary data
  metabolites <- pc_and_outliers(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=setdiff(incl_feats, xenobiotic_feats))

  # run sample summary with the current exclusions internally
  cli::cli_alert_info(glue::glue("Re-running sample summary using current exclusions..."))
  metabolites <- sample_summary(metabolites, layer = "qcing", sample_ids=incl_samps, feature_ids=incl_feats)
  cli::cli_alert_info(glue::glue("QC sample summary complete."))

  # extract PCs 1-2 or 1-number of Acceleration factor PCs
  af <- metabolites@acceleration_factor[["qcing"]]
  if(af < 2) {
    pcs <- metabolites@pcs[, 1:2, "qcing"]
  } else {
    pcs <- metabolites@pcs[, 1:af, "qcing"]
  }

  # perform exclusion on top PCs to ID outliers
  if (!is.na(metabolites@pc_outlier_sd)) {

    outliers   <- outlier_detection(pcs, nsd = metabolites@pc_outlier_sd, meansd = TRUE)
    excl_samps <- names(which(apply(outliers, 1, function(x) sum(x) > 0)))
    metabolites@exclusions[["qcing"]][["samples"]][["user_defined_sample_pca_outlier"]] <- unique(c(metabolites@exclusions[["qcing"]][["samples"]][["user_defined_sample_pca_outlier"]], excl_samps))
    incl_samps <- setdiff(incl_samps, excl_samps)
  }
  cli::cli_alert_info(glue::glue("PCA outlier analysis 2 complete."))


  # Make final QC dataset
  cli::cli_alert_info(glue::glue("Creating final QC dataset..."))
  incl_feats  <- unique(c(incl_feats, derived_feats, xenobiotic_feats)) # add back
  dat <- get_data(metabolites, layer="qcing", sample_ids=incl_samps, feature_ids=incl_feats)
  metabolites@data <- add_layer(current    = metabolites@data,
                                layer      = dat,
                                layer_name = "post_qc",
                                force      = TRUE)
  metabolites@source <- c(metabolites@source, stats::setNames("qcing", "post_qc"))
  metabolites@exclusions[["post_qc"]] <- metabolites@exclusions[["qcing"]]

  cli::cli_alert_info(glue::glue("Rerunning Sample & Feature Summary Statistics for QC data..."))
  metabolites <- sample_summary( metabolites, layer="post_qc",
                                 sample_ids= incl_samps,
                                 feature_ids=incl_feats, output="object")
  metabolites <- feature_summary(metabolites, layer="post_qc", sample_ids=incl_samps, feature_ids=incl_feats, output="object")
  metabolites <- pc_and_outliers(metabolites, layer="post_qc", sample_ids=incl_samps, feature_ids=incl_feats)
  cli::cli_alert_success(glue::glue("Sample & Feature Summary Statistics for QC data complete."))


  cli::cli_h1("Metabolite QC Process Completed")

  # return the metabolites with underlying data (+/- adjustments) and exclusion matrix
  return(metabolites)
}
