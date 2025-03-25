library(devtools)
#library(metaboprep2)
load_all()

metabolites <- Metabolites(project_name = "MYPROJECT",
                           format = "metabolon_v1",
                           total_peak_area_sd = Inf,
                           filepath = "/Users/xx20081/git/metaboprep/example_data/excel/metabolon_v1_example.xlsx")
metabolites <- import_data(metabolites)
metabolites

metabolites <- feature_summary(metabolites, type="raw")
metabolites
plot(metabolites@feature_tree$raw)

metabolites <- metabolite_qc(metabolites, source="raw", destination="post_qc")
metabolites
plot(metabolites@feature_tree$post_qc)


