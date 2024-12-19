setwd(".")
source("R/utility.R")
source("R/labelling.R")

library(data.table)


# Laod the DepMap CRIPR Cas9 effect score data
df_cripr <- read.csv('https://figshare.com/ndownloader/files/43346616')
#df_cripr <- read.csv('data/CRISPRGeneEffect.csv')
rownames(df_cripr) <- df_cripr$X
df_cripr$X <- NULL
df_cripr <- as.data.frame(t(df_cripr))  # transpose... gen by ACH

# Laod the DepMap ModelID
df_map <- read.csv('https://ndownloader.figshare.com/files/43746708')

# select the tissue/disease
line_group <- "OncotreeLineage"
tissue_list <- list("Kidney", "Lung", "CNS/Brain")
df1 <- filter_crispr_by_model(df_cripr, df_map, minlines=10, line_group=line_group)
df_nonan <- delrows_with_nan_percentage(df1, perc=95.0)
# 2-class labelling by making mode over all cell lines in different contexts
cell_lines <- select_cell_lines(df_nonan, df_map, tissue_list, nested=FALSE)
labeled_df <- labelling(df_nonan, columns = cell_lines, labelnames = list("2" = "NE", "1" = "E"),
                        n_classes = 2, algorithm = "otsu", mode = "flat-multi", verbose = TRUE)
# 3-class labelling by making mode over all cell lines in different contexts
labeled_df3 <- labelling(df_nonan, columns = cell_lines, labelnames = list("3" = "sNE", "2" = "aE", "1" = "E"), n_classes = 2, algorithm = "otsu", mode = "two-by-two", verbose = TRUE)
# 2-class labelling by making the mode over the modes on multiple contexts
cell_lines <- select_cell_lines(df_nonan, df_map, tissue_list, nested=TRUE)
labeled_df_nested <- labelling(df_nonan, columns = cell_lines, labelnames = list("2" = "NE", "1" = "E"), n_classes = 2, algorithm = "otsu", mode = "flat-multi", verbose = TRUE)
