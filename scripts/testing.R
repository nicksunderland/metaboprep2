library(devtools)
#library(metaboprep2)
load_all()

metabolites <- Metabolites(project_name = "MYPROJECT",
                           format = "metabolon_v1",
                           filepath = "/Users/xx20081/git/metaboprep/example_data/excel/metabolon_v1_example.xlsx")
metabolites <- import_data(metabolites)
metabolites


d <- metabolite_qc(metabolites)
d


d <- pc_and_outliers(metabolites)
d

d <- feature_summary(metabolites)
d

d <- sample_summary(metabolites)
d

d <- batch_normalisation(metabolites)
d


d <- get_data(metabolites, "OrigScale")
d

write_data(metabolites, "/Users/xx20081/git/metaboprep/example_data/excel")


f <- get_data(metabolites, "OrigScale", metabolite_ids = c("compid_203"), sample_ids = c("ind100"), as_df=T)
f



foo <- metaboprep::read.in.metabolon("metabolon_v1_example.xlsx",
                                     "/Users/xx20081/git/metaboprep/example_data/excel",
                                     "test_project")




Empty <- new_class("Empty",
                   properties = list(
                     x = new_property(class_numeric, default = 0),
                     y = new_property(class_character, default = ""),
                     z = new_property(class_logical, default = NA)
                   ),
                   constructor = function(x, ...) {
                     object <- new_object(S7::S7_object(),
                                          x = x)

                     return(object)
                   }
)
Empty(x=1)
