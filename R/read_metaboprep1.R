# Silence R CMD check
globalVariables(c("CHEMICAL_NAME", "PC_outlier_SD", "PLATFORM", "SUB_PATHWAY", "CHEM_ID", "PUBCHEM",
                  "KEGG", "HMDB", "INCHIKEY", "SMILES", "SUPER_PATHWAY",
                  "feature_missingness", "feature_names", "project", "raw_data",
                  "sample_missingness", "total_peak_area_SD",
                  "tree_cut_height", "barcode", "combined_id", "derived_features",
                  "i.barcode2", "label.no.units", "studyId", "study_id", "subclass",
                  "timepoint", "timepoint2"), package = "metaboprep2")

#' @title Read and Assemble Metaboprep 1 Metabolite Dataset
#'
#' @description
#' Reads metabolomics and clinical data from specified file paths, merges relevant metadata,
#' and returns a list containing a post-QC metabolomics data array, sample metadata, and feature metadata.
#'
#' @param rdata_path Character. Path to the RData file output from the Metaboprep 1 pipeline.
#' @param log_path Character. Path to the Metaboprep 1 pipeline log file output.
#' @param clinical_path Character. Path to the clinical metadata file (e.g. age, sex, weight).
#' @param manifest Character. Path to the NMR manifest prerelease description
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
#' @note
#' Will probably need to rename this to be a specific BBS study import function due to the clinical annotating
#'
read_metaboprep1 <- function(rdata_path, log_path, clinical_path, manifest=NULL) {

  # testing
  if (FALSE) {
    # rdata_path    <- file.path(Sys.getenv("BBS_NIGHTINGALE_DIR"),  "ReportData.Rdata")
    # log_path      <- file.path(Sys.getenv("BBS_NIGHTINGALE_DIR"),  "bbs_primary_nmr_2023_11_20_logfile.txt")
    # clinical_path <- file.path(Sys.getenv("BBS_CLINICAL_DIR"),   "03_sample_clinical_data_all.csv")
    # manifest      <- file.path(Sys.getenv("BBS_MANIFEST_DIR"),   "2023Q2_nmr_annotated_manifest_2023-08-14.csv")
    # rdata_path    <- file.path(Sys.getenv("BBS_METABOLON_DIR"), "ReportData.Rdata")
    # log_path      <- file.path(Sys.getenv("BBS_METABOLON_DIR"), "bbs_only_2024_06_25_logfile.txt")



  }

  # read the clinical data
  clinical_cols <- c(
    study_id         = "studyId",
    timepoint        = "timepoint",
    client_sample_id = "Client.Sample.ID.",
    study            = "source_study",
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
    baseline_to_end  = "baseline_to_end",
    site             = "site",
    time_in_freezer  = "timeinfreezer"
  )
  clinical <- data.table::fread(clinical_path, select = unname(clinical_cols), col.names = names(clinical_cols))

  # read the metaboprep1 pipeline data
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
  features_raw <- raw_data$feature_data |> data.table::as.data.table()
  if (platform=="Nightingale") {
    features     <- features_raw[, list(feature_id      = tolower(gsub("_","", feature_names)),
                                        pathway         = as.character(class),
                                        sub_pathway     = as.character(subclass),
                                        platform        = NA_character_,
                                        chemical_name   = as.character(feature_names),
                                        pubchem_id      = NA_character_,
                                        chem_id         = NA_character_,
                                        kegg            = NA_character_,
                                        hmdb            = NA_character_,
                                        inchikey        = NA_character_,
                                        smiles          = NA_character_,
                                        derived_feature = derived_features=="yes",
                                        comp_id         = NA_character_)]
  } else if (platform=="Metabolon") {
    features     <- features_raw[, list(feature_id      = tolower(gsub("_","", feature_names)),
                                        pathway         = as.character(SUPER_PATHWAY),
                                        sub_pathway     = as.character(SUB_PATHWAY),
                                        platform        = clean_names(PLATFORM),
                                        chemical_name   = as.character(CHEMICAL_NAME),
                                        pubchem_id      = as.character(PUBCHEM),
                                        chem_id         = as.character(CHEM_ID),
                                        kegg            = as.character(KEGG),
                                        hmdb            = as.character(HMDB),
                                        inchikey        = as.character(INCHIKEY),
                                        smiles          = as.character(SMILES),
                                        derived_feature = FALSE,
                                        comp_id         = sub(".*?([0-9]+).*", "\\1", feature_names))]

  } else {
    stop("Only Nightingale and Metabolon inputs supported from Metaboprep1 currently.")
  }

  # annotate from MetaboAnalystR and internal custom annotations
  features <- annotate_features(features,
                                fixed_match_cols = list(hmdb_id = "hmdb", kegg_id = "kegg", name = "chemical_name", name = "feature_id", smiles = "smiles", inchi_key = "inchikey", pubchem_id = "pubchem_id", chem_id = "chem_id", comp_id = "comp_id",  chebi_id = NULL, metlin_id = NULL),
                                fuzzy_match_cols = list(name = NULL)) # dont fuzzy match

  # get the samples data (including the platform columns)
  samples_raw <- raw_data$sample_data |> data.table::as.data.table()
  data.table::setnames(samples_raw, names(samples_raw), clean_names(names(samples_raw)))
  if (platform=="Nightingale") {

    sample_cols <- c(sample_id = "sample_id")
    samples <- samples_raw[, .SD, .SDcols = sample_cols]

    # get the NMR manifest
    stopifnot("The NMR manifest file must be provided if reading Nightingale NMR prerelease" = !is.null(manifest))
    manif <- data.table::fread(manifest)
    manif[, `:=`(barcode2   = sub("\\{\\}", "_", barcode),
                 timepoint2 = data.table::fcase(grepl("36 months", timepoint), "end",
                                                grepl("oneyear", timepoint),   "end",
                                                grepl("(?i)baseline", timepoint),  "Baseline"))]
    manif[, combined_id := paste0(studyId, "_", timepoint2)]
    data.table::setcolorder(manif, c("barcode2", "timepoint2","combined_id"))

    clinical[, combined_id := paste0(study_id, "_", timepoint)]
    clinical[manif, barcode := i.barcode2, on="combined_id"]
    data.table::setcolorder(clinical, c("barcode"))

    # join the clinical to samples
    samples_clin <- samples[clinical, on = c("sample_id"="barcode"), nomatch = NULL]

    # missing sample clinical data
    miss <- setdiff(samples$sample_id, samples_clin$sample_id)
    if (length(miss) > 0) {
      msg <- paste0(length(miss), " missing clinical data for samples: ", paste0(miss, collapse=", "))
      warning(msg)
    }

    samples <- samples_clin

  } else if (platform=="Metabolon") {

    sample_cols <- c(sample_id        = "parent_sample_name",
                     client_sample_id = "client_sample_id",
                     stats::setNames(unique(features$platform), unique(features$platform)))
    samples <- samples_raw[, .SD, .SDcols = sample_cols]
    data.table::setnames(samples, names(samples), names(sample_cols))

    # join the clinical to samples
    samples <- samples[clinical, on = "client_sample_id", nomatch = NULL]

  } else {
    stop("Only Nightingale and Metabolon inputs supported from Metaboprep1 currently.")
  }


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
