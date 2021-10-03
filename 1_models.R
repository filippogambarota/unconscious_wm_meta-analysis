## -----------------------------------------------------------------------------
## Script: Fitting Models and Information Criteria
##
## Author: Filippo Gambarota
##
## -----------------------------------------------------------------------------

# This script fit all relevant models using the brms package

# Environment -------------------------------------------------------------

seed <- 2021
iter <- 15000
chains <- 6

set.seed(seed)

# Packages ---------------------------------------------------------------

library(brms)
library(metafor)
library(tidyverse)
library(tidybayes)

# Loading Data ------------------------------------------------------------

#source("meta_analysis_preprocessing.R") # preprocessing data

meta_clean <- read_rds("data/cleaned/meta_analysis_data.rds")

# General fitting function for multiple models to be less error prones

brm_fun <- function(eff_formula, data, priors, filename, path = "mod",
                    niter = iter, nchains = chains, 
                    nseed = seed){
  
  meta_formula <- "eff_size_g|se(eff_size_se_g) ~"
  ranef_formula <- "(1|study_meta/substudy)"
  
  filepath <- file.path(path, paste0(filename, ".rds"))
  
  fit_formula <- formula(paste(meta_formula, eff_formula, "+", ranef_formula))
  
  fit <- brm(fit_formula,
      data = data,
      prior = priors,
      iter = niter,
      cores = parallel::detectCores(),
      chains = nchains,
      sample_prior = T,
      save_pars = save_pars(all = TRUE), 
      control = list(adapt_delta = 0.99),
      file = filepath,
      seed = nseed)
  
  return(fit)
  
}

#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#
#------------------------- META-ANALYSIS ----------------------------------#

# Prior List --------------------------------------------------------------

prior_list <- list(
  mod0_prior = c(
    prior(normal(0, 2), class = "b", coef = "Intercept"),
    prior(cauchy(0, 1), class = "sd")
  ),
  trialsample_prior = c(
    prior(normal(0, 2), class = "b", coef = "Intercept"),
    prior(normal(0, 0.01), class = "b", coef = "ntrials0"),
    prior(normal(0, 0.01), class = "b", coef = "sample_size0"),
    prior(cauchy(0, 1), class = "sd")
  ),
  blind_prior = c(
    prior(normal(0, 2), class = "b"),
    prior(cauchy(0, 1), class = "sd")
  ),
  paradigm_prior = c(
    prior(normal(0, 2), class = "b"),
    prior(cauchy(0, 1), class = "sd")
  ),
  author_prior = c(
    prior(normal(0, 2), class = "b",),
    prior(cauchy(0, 1), class = "sd")
  ),
  targetdur_bin_prior = c(
    prior(normal(0, 2), class = "b"),
    prior(cauchy(0, 1), class = "sd")
  )
)

# Overall Model -----------------------------------------------------------

# Overall Effect

fit0 <- brm_fun("0 + Intercept", 
                data = meta_clean, 
                filename = "fit0",
                priors = prior_list$mod0_prior)

# No grey literature

fit0_nogray <- brm_fun("0 + Intercept", 
                       data = meta_clean %>% filter(published == "Published"),
                       filename = "fit0_nogray",
                       priors = prior_list$mod0_prior)

#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#
#------------------------ METAREGRESSION ----------------------------------#

# Trials and Sample Size Model --------------------------------------------

fit_trial_sample <- brm_fun("0 + Intercept + ntrials0 + sample_size0",
                            priors = prior_list$trialsample_prior,
                            filename = "fit_trial_sample",
                            data = meta_clean)

# Blinding Paradigm Model -------------------------------------------------

# collapse Nakano sandwitch masking

meta_clean$blinding_paradigm_m <- ifelse(meta_clean$blinding_paradigm == "Sandwitch Masking",
                                         "Backward Masking", meta_clean$blinding_paradigm)

fit_blind <- brm_fun("0 + blinding_paradigm_m",
                     priors = prior_list$blind_prior,
                     filename = "fit_blind",
                     data = meta_clean)

# WM Paradigm Model -------------------------------------------------------

fit_paradigm <- brm_fun("0 + paradigm",
                     priors = prior_list$blind_prior,
                     filename = "fit_paradigm",
                     data = meta_clean)

# Author Model ------------------------------------------------------------

fit_author <- brm_fun("0 + first_author_mod",
                      priors = prior_list$author_prior,
                      filename = "fit_author", 
                      data = meta_clean)

# Target Duration Bin -----------------------------------------------------

meta_clean %>% 
  mutate(target_duration_bin = case_when(target_duration <= 50 ~ "16-50ms",
                                         target_duration > 50 & target_duration <= 500 ~ "150-500ms",
                                         target_duration > 500 ~ "3000ms")) -> meta_clean

fit_targetdur_bin <- brm_fun("0 + target_duration_bin",
                      priors = prior_list$author_prior,
                      filename = "fit_targetdur_bin", 
                      data = meta_clean)

# Publication Bias --------------------------------------------------------

# Egger Test

# For the egger test we use metafor::rma() for fitting the model and then
# the regtest function for the funnel plot asymmetry

fit_rma <- rma(eff_size_g, eff_size_variance_g, data = meta_clean)

fit_egger <- regtest(fit_rma, model = "rma", predictor = "sei")

saveRDS(fit_egger, file = "objects/fit0_egger.rds")
saveRDS(fit_rma, file = "objects/fit0_rma.rds")

# Add Information criteria and re-save models -----------------------------

# Cite this https://link.springer.com/article/10.1007/s11222-016-9696-4
# Here for terminology https://avehtari.github.io/modelselection/CV-FAQ.html

mod_list <- list.files("mod", pattern = ".rds", full.names = T)

map(mod_list, function(i) add_ic_and_save(i, criteria = c("waic", "loo")))

# Exploratory Analysis ----------------------------------------------------

meta_clean$implicit <- factor(meta_clean$implicit)

prior_impl <- c(
  prior(normal(0, 2), coef = "implicityes"),
  prior(normal(0, 2), coef = "Intercept"),
  prior(cauchy(0, 1), class = "sd")
)

fit_impl <- brm(eff_size_g|se(eff_size_se_g) ~ 0 + Intercept + implicit + (1|study_meta/substudy),
                data = meta_clean,
                iter = iter,
                cores = parallel::detectCores(),
                chains = chains,
                sample_prior = T,
                prior = prior_impl,
                save_pars = save_pars(all = TRUE), 
                control = list(adapt_delta = 0.99),
                file = "mod/fit_impl.rds",
                seed = seed)
