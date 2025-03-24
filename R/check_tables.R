check_features_table <- function(tbl) {

  required_columns <- c("feature_names", "derived_feature", "metabolite_id", "comp_id", "platform", "pathway", "kegg", "hmdb")

  missing_columns <- setdiff(required_columns, names(tbl))

  check <- if (length(missing_columns) > 0) FALSE else TRUE

  return(check)
}

# check_features_table <- function(tbl) {
#
#   required_columns <- c("feature_names", "derived_feature", "metabolite_id", "comp_id", "platform", "pathway", "kegg", "hmdb")
#
#   missing_columns <- setdiff(required_columns, names(tbl))
#
#   if (length(missing_columns) > 0) {
#     stop(paste("Missing features table columns:", paste(missing_columns, collapse = ", ")), call. = FALSE)
#   } else {
#     return(TRUE)
#   }
#
# }
