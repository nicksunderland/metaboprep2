library(devtools)
#library(metaboprep2)
load_all()

metabolites <- Metabolites(project_name = "MYPROJECT",
                           format = "metabolon_v2",
                           filepath = '/Users/xx20081/Desktop/metabolomics/UNBR-0201-17ML+ Client Data Tables(220203).xlsx') # "/Users/xx20081/git/metaboprep2/inst/extdata/metabolon_v2_example.xlsx")
# system.file("extdata", "metabolon_v1_example.xlsx", package = "metaboprep2")
metabolites <- import_data(metabolites)
metabolites

metabolites <- metabolite_qc(metabolites, source="raw", destination="post_qc")
metabolites


generate_report(metabolites,
                output_dir="/Users/xx20081/git/metaboprep2/inst/rmarkdown/templates/qc_report/skeleton",
                output_filename=NULL,
                format="pdf",
                template="qc_report")





# metabolites <- feature_summary(metabolites, type="raw")
# metabolites
# plot(metabolites@feature_tree$raw)




foo<- metabolites@exclusions[,,"post_qc"]

m <- as.matrix(foo)
highlight_matrix <- apply(as.matrix(foo), c(1, 2), function(x) grepl("6", x))

library(ggplot2)
library(reshape2)

# Convert logical matrix to a data frame for ggplot2
df <- reshape2::melt(highlight_matrix)
colnames(df) <- c("Row", "Column", "Contains_2")

# Plot using ggplot2 (color based on Contains_2)
ggplot(df, aes(x = Column, y = Row, fill = Contains_2)) +
  geom_tile() +
  scale_fill_manual(values = c("white", "blue")) + # white for FALSE, blue for TRUE (highlighting "2")
  theme_minimal() +
  labs(title = "Matrix Highlighting '2'", x = "Column", y = "Row")
