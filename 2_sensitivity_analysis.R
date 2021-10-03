## -----------------------------------------------------------------------------
## Script: Sensitivity Analysis Models
##
## Author: Filippo Gambarota
##
## -----------------------------------------------------------------------------

# This script compute the sensitivity analysis for priors and the leave-one-out
# analysis.

seed <- 2021

set.seed(seed)

# Packages ----------------------------------------------------------------

library(brms)
library(tidyverse)

# Loading Data ------------------------------------------------------------

meta_clean <- read_rds("data/cleaned/meta_analysis_data.rds")

# Loading Models ----------------------------------------------------------

fit0 <- read_rds("mod/fit0.rds")

############ LOO Sensitivity Analysis

# LOO Sensitivity Analysis - Paper Level - Overall Fit --------------------

sens_data_paper <- get_data_without_id(meta_clean, grouping_factor = "study_meta")

sens_fit_paper <- map(sens_data_paper, function(dat) refit_model(fit0, newdata = dat, seed = seed))

# LOO Sensitivity Analysis - Study Level - Overall Fit --------------------

sens_data_study <- get_data_without_id(meta_clean, grouping_factor = "study_id")

sens_fit_study <- map(sens_data_study, function(dat) refit_model(fit0, newdata = dat, seed = seed))

# Saving ------------------------------------------------------------------

saveRDS(sens_fit_paper, file = "mod/sensitivity_analysis/sens_fit_paper.rds")
saveRDS(sens_fit_study, file = "mod/sensitivity_analysis/sens_fit_study.rds")

############ Prior Sensitivity Analysis

# Average Effect Sensitivity ----------------------------------------------

avg_prior_sd <- c(0.1, 0.5, 1, 5, 10, 30)

avg_prior_list <- map(avg_prior_sd, function(x) create_prior(type = "b", x))
names(avg_prior_list) <- paste0("prior_", avg_prior_sd)
sens_prior_avg <- map(avg_prior_list, function(prior) refit_model(fit0, prior = prior, seed = seed))

# Tau Effect Sensitivity --------------------------------------------------

tau_prior_sd <- c(0.1, 0.5, 1, 5, 10, 30)

tau_prior_list <- map(tau_prior_sd, function(x) create_prior(type = "tau", x))
names(tau_prior_list) <- paste0("prior_", tau_prior_sd)

sens_prior_tau <- map(tau_prior_list, function(prior) refit_model(fit0, prior = prior, seed = seed))

# Saving ------------------------------------------------------------------

saveRDS(sens_prior_tau, file = "mod/sensitivity_analysis/sens_prior_tau.rds")
saveRDS(sens_prior_avg, file = "mod/sensitivity_analysis/sens_prior_avg.rds")