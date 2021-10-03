## -----------------------------------------------------------------------------
## Script: Post-processing Sensitivity
##
## Author: Filippo Gambarota
##
## -----------------------------------------------------------------------------

# This script compute all relevant measures and diagnostic from sensitivity analysis

# Packages ----------------------------------------------------------------

library(tidyverse)
library(tidybayes)
library(brms)
library(broomExtra)

# Data ----------------------------------------------------------------

meta_clean <- read_rds("data/cleaned/meta_analysis_data.rds")

mod_list <- load_models()

# Functions ----------------------------------------------------------------

source("functions/post_processing_functions.R")

clean_sens <- function(data){
  data %>% 
    map(., get_params_posterior) %>% 
    map(., compute_posterior_summary) %>% 
    bind_rows(., .id = "mod") %>% 
    clean_param_names(., col_param = .param)
}

# Sensitivity Prior -------------------------------------------------------

# Average

sens_prior_avg <- read_rds("mod/sensitivity_analysis/sens_prior_avg.rds")
sens_prior_avg <- add_to_list(sens_prior_avg, mod_list$fit0, name = "prior_2")

get_params_posterior(sens_prior_avg$prior_0.1)

sens_prior_avg_clean <- sens_prior_avg %>% 
  clean_sens() %>%
  mutate(prior = parse_number(mod)) %>% 
  filter(.param == "Average")

saveRDS(sens_prior_avg_clean, "objects/sens_prior_avg_clean.rds")

# Tau

sens_prior_tau <- read_rds("mod/sensitivity_analysis/sens_prior_tau.rds")

sens_prior_tau_clean <- sens_prior_tau %>% 
  clean_sens() %>%
  mutate(prior = parse_number(mod)) %>% 
  filter(startsWith(.param, "Tau"))

saveRDS(sens_prior_tau_clean, "objects/sens_prior_tau_clean.rds")

# Loo Sensitivity ---------------------------------------------------------

# Paper

sens_fit_paper <- read_rds("mod/sensitivity_analysis/sens_fit_paper.rds")

sens_fit_paper_clean <- sens_fit_paper %>% 
  clean_sens()

saveRDS(sens_fit_paper_clean, "objects/sens_fit_paper_clean.rds")

# Study

sens_fit_study <- read_rds("mod/sensitivity_analysis/sens_fit_study.rds")

sens_fit_study_clean <- sens_fit_study %>% 
  clean_sens()

saveRDS(sens_fit_study_clean, "objects/sens_fit_study_clean.rds")