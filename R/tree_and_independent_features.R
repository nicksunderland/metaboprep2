#' @title Identify Independent Features in a Numeric Matrix
#' @description
#' This function identifies independent features using Spearman's rho correlation distances, and a dendrogram
#' tree cut step.
#' @param data matrix, the metabolite data matrix. samples in row, metabolites in columns
#' @param minimum_samplesize numeric,
#' @param tree_cut_height the tree cut height. A value of 0.2 (1-Spearman's rho) is equivalent to saying that features with a rho >= 0.8 are NOT independent.
#' @param exclude_features character, a vector of feature|metabolite names to exclude from this analysis. This might be features heavily present|absent like Xenobiotics or variables derived from two or more variable already in the dataset.
#'
#' @keywords independent features
#'
#' @return a list object of (1) an hclust object, (2) independent features, (3) a data frame of feature ids, k-cluster identifiers, and a binary identifier of independent features
#'
#' @export
#'
#' @examples
#' ## define a covariance matrix
#' cmat = matrix(1, 4, 4 )
#' cmat[1,] = c(1, 0.7, 0.4, 0.2)
#' cmat[2,] = c(0.7, 1, 0.2, 0.05)
#' cmat[3,] = c(0.4, 0.2, 1, 0.375)
#' cmat[4,] = c(0.2, 0.05, 0.375,1)
#'
#' ## simulate the data (multivariable random normal)
#' set.seed(1110)
#' ex_data = MASS::mvrnorm(n = 500, mu = c(5, 45, 25, 15), Sigma = cmat )
#' rownames(ex_data) = paste0("ind", 1:nrow(ex_data))
#' colnames(ex_data) = paste0("var", 1:ncol(ex_data))
#'
#' ## run function to identify independent variables at a tree cut height
#' ## of 0.5 which is equivalent to clustering variables with a Spearman's
#' ## rho > 0.5 or (1 - tree_cut_height)
#' ind = tree_and_independent_features(ex_data, tree_cut_height = 0.5)
#'
tree_and_independent_features = function(data, minimum_samplesize = 50, tree_cut_height = 0.5, exclude_features = NA){

  #################
  ## feature ids
  #################
  if( !all(is.na(exclude_features)) || length(exclude_features)>0 ){
    ## all feature IDs
    all_fids = colnames(data)
    ## feature ids with feature exclusions
    fids = all_fids[!all_fids %in% exclude_features]
    ##
    data = data[, fids]
  } else {
    fids = colnames(data)
  }


  ## identify features with no variance
  rowvar0 <- which( apply(data, 2, function(x) var(x,na.rm=T)==0) )
  if(length(rowvar0) > 0){
    cat(paste0("\t\t\t- ", length(rowvar0) ," features excluded from analysis for 0 variance.\n"))
    data = data[, -rowvar0]
  }

  #######################
  ## missingness filter
  ## based on N, sample size
  #######################
  N = apply(data, 2, function(x){ sum(!is.na(x)) })

  f_remove = which( N < minimum_samplesize )

  if(length(f_remove) > 0){
    cat(paste0("\t\t\t- ", length(f_remove) ," features excluded from analysis for (>20% missingness) n smaller than ", minimum_samplesize ,".\n"))
    data = data[, -f_remove]
  }

  ##############################################
  ## Feature filtering
  ## Pairwise sample size (n) evaluation.
  ## Each pair of comparisons must have a complete observation 'n' >= minimum_samplesize
  ## Greedy sampling
  ##############################################

  # f_remove = greedy.pairwise.n.filter(wdata = wdata, minN = minimum_samplesize)
  # if(length(f_remove) > 0){
  #   w = which(colnames(wdata) %in% f_remove)
  #   cat(paste0("\t\t\t- ", length(f_remove) ," features excluded from analysis for n smaller than ", minimum_samplesize ," in pairwise evaluations.\n"))
  #   wdata = wdata[, -w]
  # }


  #######################
  # make tree
  #######################
  cor_matrix <- stats::cor(data, method="spearman", use = "pairwise.complete.obs")
  dist_matrix <- stats::as.dist(1-abs(cor_matrix))
  stree <- stats::hclust(dist_matrix, method = "complete")

  # restrict based on cut off
  k = stats::cutree(stree, h = tree_cut_height)

  # restrict so as to keep the feature with the least missingness
  k_group = table(k)

  # strictly independent features
  w     <- which(k_group == 1)
  ind_k <- names( k_group[w] )
  w     <- which(k %in% ind_k)
  ind   <- names(k)[w]

  ## indentify feature with least missingness in each group
  w = which(k_group > 1); k_group =  names( k_group[w] )
  ind2 = sapply(k_group, function(x){
    w = which(k %in% x); n = names( k[w] )
    o = sort(N[n], decreasing = FALSE)
    out = names(o)[1]
    return(out)
  })
  independent_features = paste( c(ind, ind2) )


  ########################################
  ## set up returning vectors (out data)
  ########################################
  ## if features where excluded at start
  ## redefine fids to bring them back in
  if(exists("all_fids")){
    fids = all_fids
  }
  ## k cluster ids
  m = match(fids, names(k))
  k = k[m]
  names(k) = fids

  ## list of independent features (0|1, 1 = yes)
  w = which(fids %in% independent_features)
  independent_features_binary = rep(0, length(fids))
  independent_features_binary[w] = 1

  out = data.frame( feature_ID = fids, k = k,
                    independent_features_binary = independent_features_binary)

  dataout = list(speartree = stree,
                 independent_features = independent_features,
                 dataframe = out)

  return(dataout)
}


