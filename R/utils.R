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


clean_names <- function(names) {
  n <- gsub(" ", "_", names)
  n <- gsub("[-\\.]", "", n)
  tolower(n)
}

