# Unconscious WM meta-analysis

This is the repository for the paper "Unconscious Visual Working Memory. A critical review and
Bayesian meta-analysis". Details about the project can be found on OSF https://osf.io/4bfnm/.

# Analysis Folder

This folder contains the full `.Rproj` in order to reproduce all the analysis, figures, tables and supplementary materials. The folder structure is the following:

- `data/` contains the raw and cleaned data
- `figures/` contains figures generated for the paper
- `functions/` custom functions used in the analysis
- `mod/` and `mod/sensitivity_analysis` are the folders for fitted models and sensitivity analysis. These folders are empty because of the size. The content of this folder can be reproduced using the `1_models.R` and `2_sensitivity_analysis.R` script
- `objects/` contains the results of post-processing of the `mod` folder such as models tables and other relevant metrics
- `renv/` this folder contains a reproducible R environment in order to have the appropriate version for each relevant package
- `supplementary_materials/` contains the `.Rmd` file to reproduce the supplementary materials `.pdf`
- `tables/` contains all tables in `.pdf` and `.docx` format

The scripts are numerated and allow to reproduce the entire analysis:

- `0_preprocessing.R`: takes the raw data and create the cleaned version for the analysis
- `1_models.R`: for fitting all the relevant models
- `2_sensitivity_analysis.R`: for computing the prior and LOO sensitivity analysis
- `3_post_processing_models.R`: for creating models tables and summary statistics
- `4_post_processing_sensitivity.R`: for creating sensitivity analysis tables and summary statistics
- `5_creating_figures.R`: to reproduce figures of the paper
- `6_creating_tables.R`: to reproduce tables of the paper
