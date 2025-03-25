# Silence R CMD check
globalVariables(c("comp_id", "feature_names", "derived_feature"), package = "metaboprep2")

#' @title Read Metabolon Data (format 1)
#' @param filepath character, commercial Metabolon excel sheet with extension .xls or .xlsx
#' @returns list,  list(data = 3D matrix, samples = samples data.table, features = features data.table)
#'
#' @examples
#' filepath <- system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
#' read_metabolon_v1(filepath)
#'
#'
#' @importFrom readxl excel_sheets read_xlsx
#' @importFrom data.table setnames as.data.table
#' @export
read_metabolon_v1 <- function(filepath) {

  if(!grepl("(?i)\\.(xls|xlsx)$", filepath)){
    stop(paste0("Expected a commercial Metabolon excel sheet with extension .xls or .xlsx\n"), call.=FALSE)
  }

  sheets <- readxl::excel_sheets(filepath)
  sheets <- sheets[sheets %in% c("OrigScale", "ScaledImp")]
  std_names <- list("OrigScale" = "raw",
                    "ScaledImp" = "scaled")
  names(sheets) <- std_names[sheets]

  data_list <- list()
  features <- NULL
  batch <- NULL
  comp_ids <- NULL

  for (i in seq_along(sheets)) {

    raw <- suppressMessages(
      readxl::read_xlsx(filepath, sheet=sheets[i], col_names=FALSE, na=c("","NA","NDEF","TAG")) |> data.table::as.data.table()
    )

    data_header_row  <- which(raw[, 1] == "metabolite_id")
    batch_header_col <- which(raw[1, ] == "SAMPLE NAME")

    if (is.null(comp_ids)) {
      comp_id_col <- which(raw[data_header_row, ] == "COMP_ID")
      comp_ids    <- raw[(data_header_row+1L):nrow(raw), .SD, .SDcols=c(1,comp_id_col)]
      cnames      <- clean_names(unlist(raw[data_header_row, .SD, .SDcols=c(1,comp_id_col)]))
      data.table::setnames(comp_ids, cnames)
    }

    if (is.null(features)) {
      features <- raw[(data_header_row+1L):nrow(raw), 1L:batch_header_col]
      cnames   <- unlist(raw[data_header_row, 1L:batch_header_col])
      cnames[grep("Group", cnames)] <- "hmdb"
      data.table::setnames(features, clean_names(cnames))
      features[, feature_names := paste0("compid_", comp_id)]
      features[, derived_feature := FALSE]
      data.table::setcolorder(features, "feature_names")
    }

    if (is.null(batch)) {
      samples <- t(raw[1L:(data_header_row-1L), (batch_header_col+1L):ncol(raw)])
      samples <- data.table::as.data.table(samples)
      cnames  <- raw[1L:(data_header_row-1L), ][[batch_header_col]]
      data.table::setnames(samples, clean_names(cnames))
    }

    data           <- t(raw[(data_header_row+1L):nrow(raw), (batch_header_col+1L):ncol(raw)][, lapply(.SD, function(x) as.numeric(gsub(",","",x)))])
    metab_id_cname <- clean_names(unname(unlist(raw[(data_header_row), 1])))
    metab_ids      <- data.table::setnames(raw[(data_header_row+1L):nrow(raw), 1], metab_id_cname)
    metab_ids      <- metab_ids[comp_ids, on=metab_id_cname]
    colnames(data) <- paste0("compid_", metab_ids[, comp_id])
    rownames(data) <- raw[data_header_row, (batch_header_col+1L):ncol(raw)]


    data_list[[names(sheets)[i]]] <- data
  }

  ref_mat <- data_list[[1]]
  for (i in seq_along(data_list)) {
    if (!identical(rownames(data_list[[1]]), rownames(data_list[[i]]))) {
      stop(paste("Row names of", sheets[i], "matrix do not match across sheets!", call.=FALSE))
    }
    if (!identical(colnames(data_list[[1]]), colnames(data_list[[i]]))) {
      stop(paste("Column names of", sheets[i], "matrix do not match across sheets!", call.=FALSE))
      }
  }
  data <- array(unlist(data_list, use.names = FALSE),
                dim = c(nrow(data_list[[1]]), ncol(data_list[[1]]), length(data_list)),
                dimnames = list(rownames(data_list[[1]]), colnames(data_list[[1]]), names(sheets)))

  return(list(data       = data,
              samples    = samples,
              features   = features))
}
