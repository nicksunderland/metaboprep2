#' @title Metabolite Quality Control
#' @description
#' This function is a wrapper function that performs the key quality controls steps on a metabolomics data set.
#' @param metabolites an object of class Metabolites
#' @include class_metabolites.R
#' @importFrom stats quantile
#' @export
metabolite_qc <- new_generic("metabolite_qc", c("metabolites"), function(metabolites, type="raw") { S7_dispatch() })
#' @name metabolite_qc
method(metabolite_qc, Metabolites) <- function(metabolites, type="raw"){

  # principles:
  # 1. keep the underlying data as it is
  # 2. build an exclusion matrix mask, accumulating codes for exclusion reasons
  # -> analysis -> ids of QC fails -> add to exclusion matrix -> apply exclusions to copy of data -> more analysis

  # exclusions codes
  excl_codes <- list("derived_feature" = 1,
                     "xenobiotic_feature" = 2,
                     "extreme_sample_missingness" = 3,
                     "extreme_feature_missingness" = 4,
                     "user_defined_sample_missingness" = 5,
                     "user_defined_feature_missingness" = 6,
                     "user_defined_sample_totalpeakarea" = 7,
                     "user_defined_sample_pca_outlier" = 8)

  # exclusions sparse matrix
  excl_matrix <- matrix(NA_character_, nrow = nrow(metabolites@data), ncol = ncol(metabolites@data),
                        dimnames = dimnames(metabolites@data)[1:2])



  # derived features to exclude
  derived_feature_excl <- character()
  if (metabolites@derived_var_exclusion) {
    derived_feature_excl <- metabolites@features[derived_feature==TRUE, feature_names]
    metabolites <- update_exclusions(metabolites, type=type, code = excl_codes[["derived_feature"]], feature_ids = derived_feature_excl)
  }

  # xenobiotic features to exclude
  if (metabolites@derived_var_exclusion) {
    xenobiotic_feature_excl <- metabolites@features[grepl("(?i)xenobiotic", pathway), feature_names]
    metabolites <- update_exclusions(metabolites, type=type, code = excl_codes[["xenobiotic_feature"]], feature_ids = xenobiotic_feature_excl)
  }

  # very bad sample missingness
  samplemis = missingness(metabolites@data[, , type], by="row", exclude_features = derived_feature_excl)
  if (!all(is.na(samplemis$missingness_w_exclusions))) {
    sample_ids <- samplemis[missingness_w_exclusions >= 0.8, sample_id]
  } else {
    sample_ids <- samplemis[missingness >= 0.8, sample_id]
  }
  metabolites <- update_exclusions(metabolites, type=type, code = excl_codes[["extreme_sample_missingness"]], sample_ids = sample_ids)

  # very bad feature missingness
  featuremis  <- missingness(data = get_data(metabolites, type), by="column")
  feature_ids <- featuremis[missingness >= 0.8, feature_id]
  metabolites <- update_exclusions(metabolites, type=type, code = excl_codes[["extreme_feature_missingness"]], feature_ids = feature_ids)

  # re-estimate sample missingness
  filtered_data <- get_data(metabolites, type, apply_exclusions=TRUE)
  samplemis <- missingness(filtered_data, by="row", exclude_features = derived_feature_excl)
  if (!all(is.na(samplemis$missingness_w_exclusions))) {
    sample_ids <- samplemis[missingness_w_exclusions >= metabolites@sample_missingness, sample_id]
  } else {
    sample_ids <- samplemis[missingness >= metabolites@sample_missingness, sample_id]
  }
  metabolites <- update_exclusions(metabolites, type=type, code = excl_codes[["user_defined_sample_missingness"]], sample_ids = sample_ids)

  # re-estimate feature missingness
  filtered_data <- get_data(metabolites, type, apply_exclusions=TRUE)
  featuremis    <- missingness(filtered_data, by="column")
  feature_ids   <- featuremis[missingness >= metabolites@feature_missingness, feature_id]
  metabolites   <- update_exclusions(metabolites, type=type, code = excl_codes[["user_defined_feature_missingness"]], feature_ids = feature_ids)

  # total peak area
  filtered_data <- get_data(metabolites, type, apply_exclusions=TRUE)
  tpa <- total_peak_area(filtered_data, features_exclude = derived_feature_excl)
  tpa[, `:=`(sdev = sd(tpa_total),
             mean = mean(tpa_total))]
  tpa[, `:=`(UL   = mean + sdev * metabolites@total_peak_area_sd,
             LL   = mean - sdev * metabolites@total_peak_area_sd)]
  sample_ids <- tpa[!data.table::between(tpa_total, LL, UL), sample_id]
  metabolites   <- update_exclusions(metabolites, type=type, code = excl_codes[["user_defined_sample_totalpeakarea"]], sample_ids = sample_ids)

  # PCA data
  filtered_data <- get_data(metabolites, type, apply_exclusions=TRUE)
  if (metabolites@outlier_treatment != "leave_be") {
    omat <- outlier_detection(data = filtered_data, nsd = metabolites@outlier_udist, meansd = FALSE, by="column")
    # indices <- which(omat == 1, arr.ind = TRUE)
    # rn <- rownames(omat)[indices[, 1]]
    # cn <- colnames(omat)[indices[, 2]]
    if(metabolites@outlier_treatment == "turn_NA") {

      filtered_data[omat == 1] <- NA_real_

    } else if (outlier_treatment == "winsorize") {

      for(i in 1:ncol(filtered_data)) {
        ## identify any outliers
        w = which(omat[,i] == 1)
        if(length(w) > 0) {
          ## estimate the 'Q' quantile value from all non-outlier samples
          quantile_value = quantile(filtered_data[-w,i], probs = c(metabolites@winsorize_quantile), na.rm = TRUE)
          ## set outliers to the quantile value
          filtered_data[w,i] = quantile_value
        }
      }
      cat(paste0("\t\t- Outliers were winsorized to the ", winsorize_quantile * 100 ," quantile of remaining (non outlying) values.\n") )
    }
  }

  # add the data


  # re-identify feature independence and PC outliers
  #featuresumstats = feature_summary(wdata = pcadata, sammis = samplemis, tree_cut_height = tree_cut_height, outlier_udist = outlier_udist, feature_names_2_exclude = derived_colnames_2_exclude)
    #     } else {
    #       featuresumstats = feature.sum.stats( wdata = pcadata, sammis = samplemis, tree_cut_height = tree_cut_height, outlier_udist = outlier_udist,  feature_names_2_exclude = NA)
    #     }




  metabolites

}


#     ###########################
#     ### 12) re-identify feature independence and PC outliers
#     ###########################
#     cat( paste0("\t\t- QCstep: re-identify independent features through correlation analysis and dendrogram clustering.\n") )
#     cat( paste0("\t\t\t- using currently QCd data.\n") )
#
#     ###########################
#     ## re-estimate independent features using the qc-data to this point
#     ###########################
#     if( !is.na(derived_colnames_2_exclude[1]) ){
#       featuresumstats = feature.sum.stats( wdata = pcadata, sammis = samplemis, tree_cut_height = tree_cut_height, outlier_udist = outlier_udist, feature_names_2_exclude = derived_colnames_2_exclude)
#     } else {
#       featuresumstats = feature.sum.stats( wdata = pcadata, sammis = samplemis, tree_cut_height = tree_cut_height, outlier_udist = outlier_udist,  feature_names_2_exclude = NA)
#     }
#
#     ###########################
#     ## extract independent feature list
#     ###########################
#     w = which(featuresumstats$table$independent_features_binary == 1)
#     ind_feature_names = rownames(featuresumstats$table)[w]
#     cat( paste0("\t\t\t* ", length(ind_feature_names), " independent features identified.\n") )
#
#     ###########################
#     ## identify PC outliers
#     ###########################
#     cat( paste0("\t\t- QCstep: Perform Principle Componenet Analysis of currently QC'd data.\n") )
#     PCs_outliers = pc.and.outliers(metabolitedata =  pcadata,
#                                    indfeature_names = ind_feature_names )
#
#     ###########################
#     ## extract PCs 1-2 or 1-number of Acceleration factor PCs
#     ###########################
#     af = as.numeric( PCs_outliers[[3]] )
#     if(af<2){
#       pcs = PCs_outliers[[1]][, 1:2]
#     } else {
#       pcs = PCs_outliers[[1]][, 1:af]
#     }
#
#     ###########################
#     ## perform exclusion on top PCs to ID outliers
#     ###########################
#     cat( paste0("\t\t- QCstep: Identify PC 1-",af," outliers >= +/-", PC_out_SD , "SD of the mean.\n") )
#     if( is.na(PC_out_SD) == FALSE){
#       outliers = outlier.matrix(pcs, nsd = PC_out_SD, meansd = TRUE)
#       outliers = apply(outliers, 1, sum)
#       w = which(outliers>0)
#
#       exclusion_data[6,1] = length(w)
#
#       if(length(w)>0){
#         cat( paste0("\t\t\t* ", length(w), " samples excluded as PC outliers.\n") )
#         wdata = wdata[-w, ]
#       } else {
#         cat( paste0("\t\t\t* 0 samples excluded as PC outliers.\n") )
#       }
#     } else {
#       cat( paste0("\t\t\tYou have chosen NOT to apply a QC-filter on individuals based on their PC eigenvectors.\n") )
#       cat( paste0("\t\t\tPC_outlier_udist in the parameter file was set to NA.\n") )
#     }
#
#     ###########################
#     ## 13) put the exclusion features back
#     ###########################
#     if( exists("exdata") ){
#       cat( paste0("\t\t- QCstep: placing the initially extracted exclusion features back into the data frame.\n") )
#       ## match sample ids
#       m = match(rownames(wdata), rownames(exdata))
#       wdata = cbind(wdata, exdata[m, ])
#     }
#
#     ##
#     return( list( wdata = wdata, featuresumstats = featuresumstats, pca = PCs_outliers, exclusion_data = exclusion_data) )
#   }


