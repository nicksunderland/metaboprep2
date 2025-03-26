check_features_table <- function(tbl) {

  required_columns <- c("feature_id", "platform", "pathway", "derived_feature") #, "feature_names", "comp_id", "kegg", "hmdb")

  missing_columns <- setdiff(required_columns, names(tbl))

  check <- if (length(missing_columns) > 0) FALSE else TRUE

  return(check)
}

check_samples_table <- function(tbl) {

  required_columns <- c("sample_id")

  missing_columns <- setdiff(required_columns, names(tbl))

  check <- if (length(missing_columns) > 0) FALSE else TRUE

  return(check)

}
