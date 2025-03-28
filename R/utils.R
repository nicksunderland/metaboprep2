#' @title List Available Data Formats
#' @description Scans the package source files for functions starting with "read_" to determine supported data formats.
#' @return A named character vector of available data formats.
#' @export
available_data_formats <- function() {
  r_dir <- system.file(package = "metaboprep2")
  if (grepl("inst$", r_dir)) {
    r_dir <- sub("/inst","",r_dir) # local package development
  }
  namespace_file <- list.files(r_dir, pattern = "NAMESPACE", full.names = TRUE)
  lines <- readLines(namespace_file, warn = FALSE)
  read_functions <- unique(grep("^export\\(read_[a-zA-Z0-9_]+\\)", lines, value = TRUE))
  formats <- sub("^export\\(read_([a-zA-Z0-9_]+)\\)", "\\1", read_functions)
  return(formats)
}



#' @title List Available Report Templates
#' @description Scans the package source files for available report templates to write to.
#' @return A character vector of available report templates
#' @export
available_report_templates <- function() {
  r_dir <- system.file(package = "metaboprep2")
  templates <- list.dirs(file.path(r_dir, "rmarkdown", "templates"), recursive = FALSE, full.names = FALSE)
  return(templates)
}


clean_names <- function(names) {
  n <- gsub(" ", "_", names)
  n <- gsub("[-\\.]", "", n)
  tolower(n)
}


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




