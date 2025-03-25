#' principal component analysis
#'
#' This function performs two principal component analysis. In the first, missing data is imputed to the median. In the second a probablistic PCA is run to account for the missingness.
#' Subsequent to the derivation of the PC, the median imputed PC data is used to identify the number of informative or "significant" PC by (1) an acceleration analysis, and (2) a parrallel analysis.
#' Finally the number of sample outliers are determined at 3, 4, and 5 standard deviations from the mean on the top PCs as determined by the acceleration factor analysis.
#'
#' @param metabolites an object of class Metabolites
#' @param type character, type/source of data to use
#'
#' @importFrom stats prcomp
#' @importFrom pcaMethods ppca
#' @importFrom nFactors parallel nScree
#'
#' @return a an object of class Metabolites with slots `@acceleration_factor`, `@n_parallel`, `@var_exp`, `@pcs`, and `@prob_pcs` filled.
#' @export
#'
pc_and_outliers <- new_generic("pc_and_outliers", c("metabolites"), function(metabolites, type="raw") { S7_dispatch() })
#' @name pc_and_outliers
method(pc_and_outliers, Metabolites) <- function(metabolites, type="raw") {

  if (!(type %in% dimnames(metabolites@data)[[3]])) {
    error_msg <- paste("Error: '", type, "' is not a valid type. Valid options are: ", paste(dimnames(metabolites@data)[[3]], collapse = ", "), ".", sep = "")
    stop(error_msg)
  }

  if (!(type %in% dimnames(metabolites@feature_summary)[[3]])) {
    error_msg <- paste("Error: @feature_summary slot does not contain '", type, "'. Have you run feature_summary(obj, type='", type, "') on this data?", sep = "")
    stop(error_msg)
  }

  ## set data
  indf         <- names(which(metabolites@feature_summary["independent_features_binary", , type] == 1))
  pcadata      <- get_data(metabolites, type = type, apply_exclusions = TRUE)[, indf]
  prob_pcadata <- get_data(metabolites, type = type, apply_exclusions = TRUE)[, indf]

  ##############################
  ## impute missingness as medians
  ##############################
  pcadata = median_impute(data = pcadata)

  ##############################
  ## z-transformation
  ##############################
  pcadata = apply(pcadata, 2, function(x){
    ( x - mean(x, na.rm = TRUE) ) / sd(x, na.rm = TRUE)
  })
  ##############################
  ## perform PCA
  ##############################
  mypca = stats::prcomp(pcadata, center = FALSE, scale = FALSE)
  varexp = summary(mypca)[[6]][2, ]

  ###############################
  ## perform probablisitic PCA
  ###############################
  prob_pcadata = apply(prob_pcadata, 2, function(x){
    ( x - mean(x, na.rm = TRUE) ) / sd(x, na.rm = TRUE)
  })
  ######
  prob_mypca = pcaMethods::ppca(prob_pcadata, method="ppca", nPcs = 10, seed = 1234 , maxIterations = 1000)
  colnames(prob_mypca@scores) <- paste0("prob_pc", 1:ncol(prob_mypca@scores))

  ##############################
  ## find number of sig PCs
  ##############################
  ev <- eigen(cor(pcadata)) # get eigenvalues
  ap <- nFactors::parallel(subject=nrow(pcadata),var=ncol(pcadata),
                           rep=100,cent=.05)
  nS <- nFactors::nScree(x=ev$values, aparallel=ap$eigen$qevpea)
  accelerationfactor = as.numeric( nS[[1]][["naf"]] )
  if(accelerationfactor < 2) { accelerationfactor = 2 }
  nsig_parrallel = nS[[1]][["nparallel"]]
  ####

  ## identify outliers
  Omat3 = outlier_detection(mypca$x[, 1:accelerationfactor], nsd = 3, by = "column")
  colnames(Omat3) = paste0(colnames(Omat3), "_3_sd_outlier")
  PCout = cbind(mypca$x[,1:10], Omat3)
  ####
  Omat4 = outlier_detection(mypca$x[, 1:accelerationfactor], nsd = 4, by = "column")
  colnames(Omat4) = paste0(colnames(Omat4), "_4_sd_outlier")
  PCout = cbind(PCout, Omat4)
  ####
  Omat5 = outlier_detection(mypca$x[, 1:accelerationfactor], nsd = 5, by = "column")
  colnames(Omat5) = paste0(colnames(Omat5), "_5_sd_outlier")
  PCout = cbind(PCout, Omat5)
  colnames(PCout) <- tolower(colnames(PCout))

  ## add to object
  metabolites@acceleration_factor[[type]] <- accelerationfactor
  metabolites@n_parallel[[type]] <- nsig_parrallel

  # var exp
  metabolites@var_exp <- add_layer(metabolites@var_exp, varexp, type, force=TRUE)
  metabolites@pcs <- add_layer(metabolites@pcs, PCout, type, force=TRUE)
  metabolites@prob_pcs <- add_layer(metabolites@prob_pcs, prob_mypca@scores, type, force=TRUE)

  return(metabolites)
}


median_impute = function( data ){
  rn = rownames(data)
  cn = colnames(data)
  out = sapply( 1:ncol(data), function(i){
    x = as.numeric( data[, i]  )
    m = median(x, na.rm = TRUE)
    w = which(is.na(x))
    x[w] = m
    return( t(x) )
  })
  rownames(out) = rn
  colnames(out) = cn
  return(out)
}


