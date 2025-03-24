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

add_layer <- function(current, layer, layer_name) {

  # Convert vector to 3D array (length, 1, 1)
  if (is.vector(layer)) {
    layer <- array(layer,
                   dim = c(length(layer), 1, 1),
                   dimnames = list(names(layer), "value", layer_name))
  } else if (is.matrix(layer)) {
    # Convert matrix to 3D array (rows, cols, 1)
    layer <- array(layer,
                   dim = c(nrow(layer), ncol(layer), 1),
                   dimnames = list(rownames(layer), colnames(layer), layer_name))
  }

  # If current is empty, initialize it directly with the first layer (no NAs)
  if (length(current) == 0) {
    return(layer)
  }

  # Convert a matrix to a 3D array if needed
  if (is.matrix(current)) {
    current <- array(current,
                     dim = c(nrow(current), ncol(current), 1),
                     dimnames = list(rownames(current), colnames(current), "initial"))
  }

  # Ensure row and column names match before appending
  if (!identical(rownames(current), rownames(layer)) ||
      !identical(colnames(current), colnames(layer))) {
    stop("Error: Row names and column names must match before joining layers.")
  }

  # Check if the layer_name already exists in the current array depth
  existing_layer_index <- which(dimnames(current)[[3]] == layer_name)

  if (length(existing_layer_index) > 0) {
    # If the layer exists, overwrite it
    current[, , existing_layer_index] <- layer
  } else {
    # If the layer does not exist, append the new layer
    current <- array(c(current, layer),
                     dim = c(dim(layer)[1], dim(layer)[2], dim(current)[3] + 1),
                     dimnames = list(rownames(layer),
                                     colnames(layer),
                                     c(dimnames(current)[[3]], layer_name)))
  }

  return(current)
}

