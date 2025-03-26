#
# add_import_config <- function(
#     config_name,
#     file_type,
#     sample_id,
#     feature_id,
#     top_left_cell              = list(),
#     orientation                = list(),
#     sheets                     = list(),
#     sample_parent_id           = NULL,
#     sample_client_id           = NULL,
#     sample_box_id              = NULL,
#     sample_pair                = NULL,
#     sample_lot                 = NULL,
#     sample_volume_extracted    = NULL,
#     feature_comp_id            = NULL,
#     feature_chemical_id        = NULL,
#     feature_ri                 = NULL,
#     feature_cas                = NULL,
#     feature_pubchem            = NULL,
#     feature_platform           = NULL,
#     feature_pathway            = NULL,
#     feature_sub_pathway        = NULL,
#     feature_kegg               = NULL,
#     feature_hmdb               = NULL,
#     feature_mass               = NULL,
#     feature_pathway_sort_order = NULL,
#     feature_biochemical        = NULL,
#     overwrite                  = FALSE) {
#
#   # get these before anything else created
#   call_args    <- as.list(environment())
#
#   # existing config
#   config <- yaml::read_yaml("/Users/xx20081/git/metaboprep2/data-raw/config.yml")
#   if (is.null(config)) config <- list(import_configs = list())
#
#   # check args
#   file_type <- match.arg(file_type, choices = c("excel", "flat"))
#   stopifnot("`config_name` is already in the config.yml but `overwrite=FALSE`" = !(config_name %in% names(config$import_configs) && !overwrite))
#   stopifnot("`sheets` cannot be length zero if file_type='excel'" = file_type!="excel" | file_type=="excel" && length(sheets)>0)
#   stopifnot("`top_left_cell` cannot be length zero if file_type='excel'" = file_type!="excel" | file_type=="excel" && length(top_left_cell)>0)
#   stopifnot("`orientation` cannot be length zero if file_type='excel'" = file_type!="excel" | file_type=="excel" && length(orientation)>0)
#   stopifnot("`sheets`, `top_left_cell`, `orientation` must be equal length and length > 0 if file_type='excel'" = file_type!="excel" | (file_type=="excel" && length(unique(sapply(list(sheets, top_left_cell, orientation), length)))==1 && unique(sapply(list(sheets, top_left_cell, orientation), length)) > 0))
#   stopifnot("`orientation` must be a list of lists with names 'samples' and 'features' if file_type='excel'" = file_type!="excel" | file_type=="excel" && all(sapply(orientation, function(x) length(setdiff(c("samples", "features"), names(x)))==0)))
#   stopifnot("`top_left_cell` must be a list of lists with names 'samples' and 'features' if file_type='excel'" = file_type!="excel" | file_type=="excel" && all(sapply(top_left_cell, function(x) length(setdiff(c("samples", "features"), names(x)))==0)))
#
#
#   # get the input arguments for the sample data and feature data column names
#   sample_cols  <- call_args[grepl("^sample", names(call_args))]
#   feature_cols <- call_args[grepl("^feature", names(call_args))]
#   other_args   <- call_args[!grepl("^(config_name|feature|sample|overwrite)", names(call_args))]
#
#   # list to write out
#   config$import_configs[[config_name]] <- other_args
#   config$import_configs[[config_name]][["sample_columns"]] <- sample_cols
#   config$import_configs[[config_name]][["feature_columns"]] <- feature_cols
#
#   yaml::write_yaml(config, "/Users/xx20081/git/metaboprep2/data-raw/config.yml")
#
# }
#
# add_import_config(config_name   = "metabolon_v2",
#                   file_type     = "excel",
#                   sheets        = list("OrigScale","ScaledImp"),
#                   top_left_cell = list(list(samples="A2", features="A3"), list(samples="A3", features="A4")),
#                   orientation   = list(list(samples="v", features="h"), list(samples="v", features="h")),
#                   sample_id   = "sid",
#                   feature_id  = "fid")
#
