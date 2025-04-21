# Silence R CMD check
globalVariables(c("feature_id", "pathway", "sample_id", "derived_feature", "Excel column name",
                  "Biomarker name", "Unit", "Group", "Subgroup", "anno_derived_feature", "anno_main_class"), package = "metaboprep2")

#' @title Read Nightingale Data (format 2)
#' @param filepath character, commercial Nightingale excel sheet with extension .xls or .xlsx
#' @returns list,  list(data = 3D matrix, samples = samples data.table, features = features data.table)
#'
#' @examples
#' filepath <- system.file("extdata", "nightingale_v2_example.xlsx", package = "metaboprep2")
#' read_nightingale_v2(filepath)
#'
#' @importFrom readxl excel_sheets read_xlsx
#' @importFrom data.table setnames as.data.table
#'
#' @export
#'
read_nightingale_v2 <- function(filepath) {

  # testing
  if (FALSE) {
    filepath <- system.file("extdata", "nightingale_v2_example.xlsx", package = "metaboprep2")
  }

  # check excel file
  if(!grepl("(?i)\\.(xls|xlsx)$", filepath)){
    stop(paste0("Expected a commercial Nightingale excel sheet with extension .xls or .xlsx\n"), call.=FALSE)
  }

  # check sheet names
  sheets <- readxl::excel_sheets(filepath)
  expected_sheets <- c("Worksheet")
  if (!all(expected_sheets %in% sheets) ) {
    warning(
      paste0("anticipated excel sheet not found:\n",
             " - expected: ", paste("`", expected_sheets, "`", sep="", collapse=", "), "\n",
             " - observed: ", paste("`", sheets, "`", sep="", collapse=", "), "\n",
             "...proceeding with first sheet", sheets[1]
      )
    )
  }
  sheet <- sheets[1]

  # get raw data
  raw <- suppressMessages(
    readxl::read_xlsx(filepath, sheet=sheet, col_names=FALSE, na=c("","NA","NDEF","TAG")) |> data.table::as.data.table()
  )

  # get data positioning
  top_corner <- as.matrix(raw[1:ifelse(nrow(raw)<20,nrow(raw),20), 1:ifelse(nrow(raw)<10,nrow(raw),10)])
  head_inds  <- which(top_corner == "sampleid", arr.ind = TRUE)
  data_inds  <- which(top_corner == "success %", arr.ind = TRUE) + c(1,1)

  # get the samples
  raw_sample_ids <- unname(unlist(raw[data_inds[1L,"row"]:nrow(raw), head_inds[1L,"col"]:head_inds[1L,"col"]]))
  samples <- raw[data_inds[1L,"row"]:nrow(raw), 1:(head_inds[1L,"col"]-1)][, lapply(.SD, function(x) as.integer(!is.na(x)))]
  names(samples) <- unlist(raw[(head_inds[1L,"row"]-1):(head_inds[1L,"row"]-1), 1:(head_inds[1L,"col"]-1)])
  samples[, sample_id := raw_sample_ids]
  data.table::setcolorder(samples, "sample_id")

  # get the features
  raw_feature_ids <- unname(unlist(raw[head_inds[1L,"row"]:head_inds[1L,"row"], (head_inds[1L,"col"]+1):ncol(raw)]))
  features <- data.table::data.table(
    feature_id = raw_feature_ids
  )

  features <- annotate_features(features,
                                fixed_match_cols = list(name = "feature_id"),
                                fuzzy_match_cols = list(name = NULL)) # dont fuzzy match

  # must have columns
  features[, `:=`(pathway  = anno_main_class,
                  platform = NA_character_,
                  derived_feature = anno_derived_feature)]


    # get the data
  data <- as.matrix(raw[data_inds[1L,"row"]:nrow(raw), data_inds[1L,"col"]:ncol(raw)][, lapply(.SD, function(x) as.numeric(gsub(",","",x)))])
  data <- array(data,
                dim = c(nrow(data), ncol(data), 1),
                dimnames = list(raw_sample_ids, features$feature_id, "raw"))

  # return
  return(list(data       = data,
              samples    = samples,
              features   = features))
}
