#' principal component analysis
#'
#' This function performs two principal component analysis. In the first, missing data is imputed to the median. In the second a probablistic PCA is run to account for the missingness.
#' Subsequent to the derivation of the PC, the median imputed PC data is used to identify the number of informative or "significant" PC by (1) an acceleration analysis, and (2) a parrallel analysis.
#' Finally the number of sample outliers are determined at 3, 4, and 5 standard deviations from the mean on the top PCs as determined by the acceleration factor analysis.
#'
#' @param metabolitedata the metabolite data matrix. samples in row, metabolites in columns
#' @param indfeature_names a vector of independent feature names | column names.
#' @param outliers  defaulted to TRUE, a TRUE|FALSE binary flagging if you would like outliers identified.
#'
#' @keywords PCA probalistic PCA
#'
#' @importFrom stats prcomp
#' @importFrom pcaMethods ppca
#' @importFrom nFactors parallel nScree
#'
#' @return a list object of length five, with (1) a data frame of PC loadings, (2) a vector of variance explained estimates for each PC, (3) an estimate of the number of informative or top PCs determined by the acceleration factor analysis, (4) an estimate of the number of informative or top PCs determined by parrallel analysis, (5) a data frame of the probablisitic PC loadings
#'
#' @export
#'
pc_and_outliers <- new_generic("pc_and_outliers", c("metabolites"), function(metabolites) { S7_dispatch() })
#' @name pc_and_outliers
method(pc_and_outliers, Metabolites) <- function(metabolites) {

  # testing
  outliers = TRUE
  metabolites@features[, independent_features_binary := sample(0:1, .N, replace = T)]


  stopifnot("`independent_features_binary` column not found in @features, please run feature_summary() on your metabolites object first" = "independent_features_binary" %in% names(metabolites@features))

  ## set data
  indf         <- metabolites@features[independent_features_binary==1, feature_names]
  pcadata      <- metabolites@data[, indf, "raw"]
  prob_pcadata <- metabolites@data[, indf, "raw"]

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

  ##############################
  ## find number of sig PCs
  ##############################
  ev <- eigen(cor(pcadata)) # get eigenvalues
  ap <- nFactors::parallel(subject=nrow(pcadata),var=ncol(pcadata),
                           rep=100,cent=.05)
  nS <- nFactors::nScree(x=ev$values, aparallel=ap$eigen$qevpea)
  accelerationfactor = as.numeric( nS[[1]][2] )
  if(accelerationfactor < 2) { accelerationfactor = 2 }
  nsig_parrallel = nS[[1]][3]
  ####

  ## identify outliers
  Omat3 = outlier_detection(mypca$x[, 1:accelerationfactor], nsd = 3, by = "column")
  colnames(Omat3) = paste0(colnames(Omat3), "_3_SD_outlier")
  PCout = cbind(mypca$x[,1:10], Omat3)
  ####
  Omat4 = outlier_detection(mypca$x[, 1:accelerationfactor], nsd = 4, by = "column")
  colnames(Omat4) = paste0(colnames(Omat4), "_4_SD_outlier")
  PCout = cbind(PCout, Omat4)
  ####
  Omat5 = outlier_detection(mypca$x[, 1:accelerationfactor], nsd = 5, by = "column")
  colnames(Omat5) = paste0(colnames(Omat5), "_5_SD_outlier")
  PCout = cbind(PCout, Omat5)


  ## build list to return
  dataout = list(pcs = PCout, varexp = varexp, accelerationfactor = accelerationfactor, nsig_parrallel = nsig_parrallel, prob_pca = prob_mypca@scores  )

  return(dataout)
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


