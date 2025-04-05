# Silence R CMD check
globalVariables(c("feature_id", "pathway", "sample_id", "derived_feature", "Excel column name",
                  "Biomarker name", "Unit", "Group", "Subgroup"), package = "metaboprep2")

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

  # # get the feature annotations data
  # features <- suppressMessages(
  #   readxl::read_xlsx(filepath, sheet=sheets[["feature_annotations"]], na=c("","NA","NDEF","TAG")) |> data.table::as.data.table()
  # )[, list(feature_id      = clean_names(`Excel column name`),
  #          feature_name    = `Biomarker name`,
  #          feature_unit    = `Unit`,
  #          platform        = NA_character_,
  #          pathway         = `Group`,
  #          sub_pathway     = `Subgroup`,
  #          derived_feature = grepl("(?i)ratio|%", `Unit`))]
  #
  # # get the sample annotations data
  # samp_annot <- suppressMessages(
  #   readxl::read_xlsx(filepath, sheet=sheets[["sample_tags"]], na=c("","NA","NDEF","TAG")) |> data.table::as.data.table()
  # )
  #
  # # get sample data positioning
  # top_corner <- as.matrix(samp_annot[1:ifelse(nrow(samp_annot)<20,nrow(samp_annot),20), 1:ifelse(nrow(samp_annot)<10,nrow(samp_annot),10)])
  # head_inds  <- which(top_corner == "Sample id", arr.ind = TRUE) + c(-1,1)
  # data_inds  <- which(top_corner == "Sample id", arr.ind = TRUE) + c(3,0)
  #
  # # extract sample data
  # samples <- samp_annot[data_inds[1L,"row"]:nrow(samp_annot), data_inds[1L,"col"]:ncol(samp_annot)] |> data.table::as.data.table()
  # names(samples) <- c("sample_id", unlist(samp_annot[head_inds[1L,"row"]:head_inds[1L,"row"], head_inds[1L,"col"]:ncol(samp_annot)]))

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
    ng_id = raw_feature_ids
  )
  features <- annotate_features(features, id_col="ng_id")

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
