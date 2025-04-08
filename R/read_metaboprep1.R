# Silence R CMD check
globalVariables(c("CHEMICAL_NAME", "PC_outlier_SD", "PLATFORM", "SUB_PATHWAY",
                  "SUPER_PATHWAY", "feature_missingness", "feature_names",
                  "project", "raw_data", "sample_missingness", "total_peak_area_SD",
                  "tree_cut_height"), package = "metaboprep2")

#' @title Read and Assemble Metaboprep 1 Metabolite Dataset
#'
#' @description
#' Reads metabolomics and clinical data from specified file paths, merges relevant metadata,
#' and returns a list containing a post-QC metabolomics data array, sample metadata, and feature metadata.
#'
#' @param rdata_path Character. Path to the RData file output from the Metaboprep 1 pipeline.
#' @param log_path Character. Path to the Metaboprep 1 pipeline log file output.
#' @param clinical_path Character. Path to the clinical metadata file (e.g. age, sex, weight).
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`data`}{A 3D array of metabolite abundances with dimensions `[samples x features x 1]`, labeled `"post_qc"`.}
#'   \item{`samples`}{A `data.table` of sample metadata joined with clinical data.}
#'   \item{`features`}{A `data.table` of metabolite feature metadata.}
#' }
#'
#' @examples
#' \dontrun{
#' read_ieu_bbs(
#'   rdata_path = "path/to/metabolites.txt",
#'   log_path = "path/to/logfile.txt",
#'   clinical_path = "path/to/clinical.csv"
#' )
#' }
#'
#' @importFrom data.table fread
#' @export
#'
read_metaboprep1 <- function(rdata_path, log_path, clinical_path) {

  # testing
  if (FALSE) {
    rdata_path     <- file.path(Sys.getenv("BBS_METABOLITE_DIR"), "ReportData.Rdata")
    log_path   <- file.path(Sys.getenv("BBS_METABOLITE_DIR"), "bbs_only_2024_06_25_logfile.txt")
    clinical_path <- file.path(Sys.getenv("BBS_CLINICAL_DIR"),   "03_sample_clinical_data_all.csv")
  }

  # read the clinical data
  clinical_cols <- c(
    study_id         = "studyId",
    timepoint        = "timepoint",
    client_sample_id = "Client.Sample.ID.",
    sex              = "sex",
    age              = "age",
    ethnicity        = "ethnicity",
    height           = "hgt_m",
    weight           = "wgt",
    weight_op        = "wgt_op_day",
    bmi              = "bmikgm2",
    op_date          = "opdate",
    baseline_to_op   = "baseline_to_op",
    op_to_end        = "op_to_end",
    baseline_to_end  = "baseline_to_end"
  )
  clinical <- data.table::fread(clinical_path, select = unname(clinical_cols), col.names = names(clinical_cols))

  # read th emetaboprep1 pipeline data
  pipeline <- load(rdata_path)

  # read the log file
  log              <- readLines(log_path)
  iqr_outlier_line <- log[grepl("interquartile range unit distance", log)]
  treatment_line   <- log[grepl("PCA & PCA only, should be:", log)]

  # get the parameter file details
  config_props <- list(
    project_name              = project,
    created                   = as.Date(gsub("_","-",sub(".+?([0-9]+_[0-9]+_[0-9]+)_logfile.txt", "\\1", basename(log_path)))),
    feature_missingness       = feature_missingness,
    sample_missingness        = sample_missingness,
    total_peak_area_sd        = total_peak_area_SD,
    outlier_udist             = as.numeric(sub(".+? sample an outlier to be: (.+)", "\\1", iqr_outlier_line)),
    outlier_treatment         = sub(".+? PCA & PCA only, should be: (.+)", "\\1", treatment_line),
    winsorize_quantile        = 1, # fixed in current pipeline
    tree_cut_height           = tree_cut_height,
    pc_outlier_sd             = PC_outlier_SD,
    derived_var_exclusion     = TRUE,
    xenobiotics_var_exclusion = TRUE,
    batch_normalise           = FALSE # metaboprep1 pipeline normalises data prior to saving the raw data object, so don't do again
  )

  # get the features data
  features <- raw_data$feature_data |> data.table::as.data.table()
  features <- features[, list(feature_id      = feature_names,
                              pathway         = SUPER_PATHWAY,
                              sub_pathway     = SUB_PATHWAY,
                              platform        = clean_names(PLATFORM),
                              chemical_name   = CHEMICAL_NAME,
                              derived_feature = FALSE)] # will need to adjust if when reading Nightingale

  # get the samples data (including the platform columns)
  samples <- raw_data$sample_data |> data.table::as.data.table()
  data.table::setnames(samples, names(samples), clean_names(names(samples)))
  sample_cols <- c(sample_id        = "parent_sample_name",
                   client_sample_id = "client_sample_id",
                   study            = "source_study",
                   stats::setNames(unique(features$platform), unique(features$platform)))
  samples <- samples[, .SD, .SDcols = sample_cols]
  data.table::setnames(samples, names(samples), names(sample_cols))

  # join the clinical to samples
  samples <- samples[clinical, on = "client_sample_id", nomatch = NULL]

  # get the raw data
  data <- raw_data$metabolite_data[samples$sample_id, features$feature_id]
  data <- array(data,
                dim = c(nrow(data), ncol(data), 1),
                dimnames = list(rownames(data), colnames(data), "raw"))

  # return structure
  out <- list(
    data            = data,
    features        = features,
    samples         = samples,
    feature_summary = numeric(),
    sample_summary  = numeric(),
    pcs             = numeric(),
    prob_pcs        = numeric(),
    var_exp         = numeric(),
    feature_tree    = list(),
    config          = config_props
  )

  return(out)
}
