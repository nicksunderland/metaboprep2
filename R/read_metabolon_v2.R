# Silence R CMD check
globalVariables(c("derived_feature", "feature_id", "pathway", "sample_id"), package = "metaboprep2")

#' @title Read Metabolon Data (format 2)
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
read_metabolon_v2 <- function(filepath) {

  # testing
  if (FALSE) {
    filepath <- system.file("extdata", "metabolon_v2_example.xlsx", package = "metaboprep2")
  }

  if(!grepl("(?i)\\.(xls|xlsx)$", filepath)){
    stop(paste0("Expected a commercial Metabolon excel sheet with extension .xls or .xlsx\n"), call.=FALSE)
  }

  sheets <- readxl::excel_sheets(filepath)
  sheets <- sheets[sheets %in% c("OrigScale", "ScaledImpData")]
  std_names <- list("OrigScale" = "raw",
                    "ScaledImpData" = "scaled")
  names(sheets) <- std_names[sheets]

  data_list <- list()
  features <- NULL
  samples <- NULL
  comp_ids <- NULL

  for (i in seq_along(sheets)) {

    raw <- suppressMessages(
      readxl::read_xlsx(filepath, sheet=sheets[i], col_names=FALSE, na=c("","NA","NDEF","TAG")) |> data.table::as.data.table()
    )

    data_header_row  <- which(!is.na(raw[, 1]))[1]
    batch_header_col <- which(!is.na(raw[1, ]))[1]

    if (is.null(features)) {
      features <- raw[(data_header_row+1L):nrow(raw), 1L:batch_header_col]
      cnames   <- unlist(raw[data_header_row, 1L:batch_header_col])
      cnames[grep("(?i)hmdb", cnames)] <- "hmdb"
      data.table::setnames(features, clean_names(cnames))

      # must have columns
      comp_id_col  <- grep("(?i)COMP.*?ID", names(features), value = TRUE)[1]
      pathway_col  <- grep("(?i)super.*?pathway", names(features), value = TRUE)[1]
      platform_col <- grep("(?i)platform", names(features), value = TRUE)[1]
      data.table::setnames(features, c(comp_id_col, pathway_col, platform_col), c("comp_id", "pathway", "platform"))

      # annotate
      features <- annotate_features(features,
                                    fixed_match_cols = list(hmdb_id = "hmdb", kegg_id = "kegg", name = "biochemical", pubchem_id = "pubchem",  chem_id = "chemical_id", comp_id = "comp_id"),
                                    fuzzy_match_cols = list(name = NULL)) # dont fuzzy match

      # must have columns
      features[, `:=`(feature_id = paste0("comp_id_", comp_id),
                      pathway  = pathway,
                      platform = platform,
                      derived_feature = anno_derived_feature)]
    }

    if (is.null(samples)) {
      samples <- t(raw[1L:(data_header_row-1L), (batch_header_col+1L):ncol(raw)])
      samples <- data.table::as.data.table(samples)
      cnames  <- raw[1L:(data_header_row-1L), ][[batch_header_col]]
      data.table::setnames(samples, clean_names(cnames))
      samples[, sample_id := get(names(samples)[grepl("(?i)sample.*?name", names(samples))][1])]
      data.table::setcolorder(samples, c("sample_id"))
    }

    data           <- t(raw[(data_header_row+1L):nrow(raw), (batch_header_col+1L):ncol(raw)][, lapply(.SD, function(x) as.numeric(gsub(",","",x)))])
    colnames(data) <- features$feature_id
    rownames(data) <- samples$sample_id

    data_list[[names(sheets)[i]]] <- data
  }

  data <- array(unlist(data_list, use.names = FALSE),
                dim = c(nrow(data_list[[1]]), ncol(data_list[[1]]), length(data_list)),
                dimnames = list(rownames(data_list[[1]]), colnames(data_list[[1]]), names(sheets)))

  return(list(data       = data,
              samples    = samples,
              features   = features))
}
