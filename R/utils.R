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


add_layer <- function(current, layer, layer_name, force=FALSE) {

  if (is.vector(layer)) {
    layer <- array(layer,
                   dim = c(length(layer), 1, 1),
                   dimnames = list(names(layer), "value", layer_name))
  } else if (is.matrix(layer)) {
    layer <- array(layer,
                   dim = c(nrow(layer), ncol(layer), 1),
                   dimnames = list(rownames(layer), colnames(layer), layer_name))
  } else if (inherits(layer, "SparseMatrix")) {
    layer <- SparseArray::SparseArray(array(layer,
                                            dim = c(nrow(layer), ncol(layer), 1),
                                            dimnames = list(rownames(layer), colnames(layer), layer_name)))
  }

  # If current is empty, initialize it directly with the first layer (no NAs)
  if (length(current) == 0) {
    return(layer)
  }

  # Ensure row and column names match before appending
  if (!identical(rownames(current), rownames(layer)) ||
      !identical(colnames(current), colnames(layer))) {

    # force together anyway (ok for things like PCs or varexp)
    if (force) {
      all_rows <- union(rownames(current), rownames(layer))
      all_cols <- union(colnames(current), colnames(layer))
      current_expanded <- array(NA_character_,
                                dim = c(length(all_rows), length(all_cols), dim(current)[3]),
                                dimnames = list(all_rows, all_cols, dimnames(current)[[3]])
      )
      current_expanded[rownames(current), colnames(current), ] <- current
      layer_expanded <- array(NA_character_,
                              dim = c(length(all_rows), length(all_cols), 1),
                              dimnames = list(all_rows, all_cols, layer_name)
      )
      layer_expanded[rownames(layer), colnames(layer), ] <- layer
      current <- current_expanded
      layer <- layer_expanded

    } else {
      stop("Error: Row names and column names must match before joining layers.")
    }
  }

  # Check if the layer_name already exists in the current array depth
  existing_layer_index <- which(dimnames(current)[[3]] == layer_name)

  # If the layer exists, overwrite it
  if (length(existing_layer_index) > 0) {

    current[, , existing_layer_index] <- layer

  # If it doesnt append layer
  } else {

    if (is.array(current)) {

      current <- array(c(current, layer),
                       dim = c(dim(layer)[1], dim(layer)[2], dim(current)[3] + 1),
                       dimnames = list(rownames(layer),
                                       colnames(layer),
                                       c(dimnames(current)[[3]], layer_name)))

    } else if (inherits(current, "SparseArray")) {

      current <- SparseArray::abind(current, layer)

    } else {

      stop("why is current not an array or sparse array")

    }

  }



  return(current)
}

