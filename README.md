# HELP
The Help Framework for Essetial Gene detection and Prediction

![Alt text of the image](https://github.com/giordamaug/HELP/blob/main/docs/images/GraphicAbstract_IG.jpg)

# Documentation
The HELP software documentation is available at: https://giordamaug.github.io/myprojects/HELP/

# Installation
``` 
#install.packages("devtools")
devtools::install_github("giordamaug/HELPr")
```
# Usage

### Load package
```
library(HELPr)

# Laod the DepMap CRISPR Cas9 effect scores
df_cripr <- read.csv('data/CRISPRGeneEffect.csv')
rownames(df_cripr) <- df_cripr$X
df_cripr$X <- NULL
df_cripr <- as.data.frame(t(df_cripr))  # transpose... gen by ACH

# Laod the DepMap ModelID
df_map <- read.csv('data/Model.csv')

# select the tissue/disease
line_group <- "OncotreeLineage"
tissue_list <- list("Kidney", "Lung", "CNS/Brain")
df1 <- filter_crispr_by_model(df_cripr, df_map, minlines=10, line_group=line_group)
df_nonan <- delrows_with_nan_percentage(df1, perc=95.0)
# 2-class labelling by making mode over all cell lines in different contexts
cell_lines <- select_cell_lines(df_nonan, df_map, tissue_list, nested=FALSE)
labeled_df <- labelling(df_nonan, columns = cell_lines, labelnames = list("2" = "NE", "1" = "E"),
                        n_classes = 2, algorithm = "otsu", mode = "flat-multi", verbose = TRUE)
```

# Credits
The HELP Framework was developed by the Computational Data Science group of High Performance Computing and Networking Institute of National Research Council of Italy (ICAR-CNR).

# Cite
If you use want to reference this software, please use the DOI: doi/10.5281/zenodo.10964743 

[![DOI](https://zenodo.org/badge/753478555.svg)](https://zenodo.org/doi/10.5281/zenodo.10964743)

If you want to cite the work in which this software is first used and described, 
please cite the following article:

```
@article {Granata2024.04.16.589691,
	author = {Ilaria Granata and Lucia Maddalena and Mario Manzo and Mario  Rosario Guarracino and Maurizio Giordano},
	title = {HELP: A computational framework for labelling and predicting human context-specific essential genes},
	elocation-id = {2024.04.16.589691},
	year = {2024},
	doi = {10.1101/2024.04.16.589691},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/04/20/2024.04.16.589691},
	eprint = {https://www.biorxiv.org/content/early/2024/04/20/2024.04.16.589691.full.pdf},
	journal = {bioRxiv}
}
```

