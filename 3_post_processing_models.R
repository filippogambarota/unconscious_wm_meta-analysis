## -----------------------------------------------------------------------------
## Script: Post-processing Models
##
## Author: Filippo Gambarota
##
## -----------------------------------------------------------------------------

# This script compute all relevant measures and diagnostic from fitted models.

# Packages ----------------------------------------------------------------

library(tidyverse)
library(tidybayes)
library(brms)
library(bayestestR)
library(broomExtra)
library(pwr)

# Loading Data ------------------------------------------------------------

meta_clean <- read_rds("data/cleaned/meta_analysis_data.rds")

mod_list <- load_models()

# General Settings --------------------------------------------------------

ROPE <- c(-Inf, 0.1)
set.seed(2021)

# ROPE Metaregression -----------------------------------------------------

# For the rope approach see:
#  - https://cran.r-project.org/web/packages/bayestestR/vignettes/region_of_practical_equivalence.html
#  - Kruschke, J. K., & Liddell, T. M. (2018). The bayesian new statistics: Hypothesis testing, estimation, 
#    meta-analysis, and power analysis from a bayesian perspective. Psychonomic Bulletin & Review, 25(1), 
#    178–206. (this for the -1, 1 range) 

# For the ROPE approach applied to meta-analysis https://journals.sagepub.com/doi/pdf/10.1177/2515245918771304
# in particular see the supplementary materials https://osf.io/fchdr/

# Check this for a meta-analysis on binary data 
# http://doingbayesiandataanalysis.blogspot.com/2011/08/extrasensory-perception-esp-bayesian.html#comment-form

rope_table <- mod_list %>% 
  map(., get_params_posterior) %>% 
  map(., get_rope_table) %>% 
  bind_rows(., .id = "mod")

saveRDS(rope_table, "objects/rope_table_new.rds")

# Computing Model Weights -------------------------------------------------

mod_to_keep <- c("fit_author", "fit_blind", "fit_paradigm",
                 "fit_targetdur_bin", "fit_trial_sample", "fit0")

mod_list_sel <- mod_list[names(mod_list) %in% mod_to_keep]

loo_weights <- model_weights(
  mod_list_sel$fit_author,
  mod_list_sel$fit_blind,
  mod_list_sel$fit_paradigm,
  mod_list_sel$fit_targetdur_bin,
  mod_list_sel$fit_trial_sample, 
  mod_list_sel$fit0,
  weights = "loo")

loo_tab <- map(mod_list_sel, function(x) {
  temp <- data.frame(loo(x)$estimates)
  temp$index <- rownames(temp)
  rownames(temp) <- NULL
  return(temp)
})

loo_tab <- bind_rows(loo_tab, .id = "mod")

loo_tab$weights <- rep(loo_weights, each = 3)

loo_tab <- loo_tab %>% 
  filter(index == "elpd_loo") %>% 
  rename("Model" = mod,
         "Weight" = weights) %>% 
  clean_mod_names(., col_mod = Model)

saveRDS(loo_tab, "objects/loo_tab.rds")

# R2 ----------------------------------------------------------------------

r2_tab <- tibble(
  mod = names(mod_list),
  r2_2 = map_dbl(mod_list, function(x) meta_r2(mod_list$fit0, x)$rsquared_2),
  r2_3 = map_dbl(mod_list, function(x) meta_r2(mod_list$fit0, x)$rsquared_3)
)

saveRDS(r2_tab, "objects/r2_tab.rds")

# Mod Tab -----------------------------------------------------------------

mod_tab <- mod_list %>% 
  map(., get_model_table) %>% 
  bind_rows(., .id = "mod") %>% 
  mutate(term = ifelse(mod == "fit_impl" & term == "b_(Intercept)",
                       "Explicit", term)) %>% 
  clean_param_names(., term) %>% 
  clean_mod_names(., mod) %>% 
  select(-component)

saveRDS(mod_tab, "objects/mod_tab.rds")

# ICC ---------------------------------------------------------------------

mod_list %>% 
  map(., meta_icc) %>% 
  bind_rows(., .id = "mod") %>%
  saveRDS(., "objects/icc_tab.rds")

# Trials Table ------------------------------------------------------------

summ_trials <- meta_clean %>% 
  drop_na(ntrials) %>% # removing Pan (2014)
  summarise(
    .mean = mean(ntrials),
    .sd = sd(ntrials),
    .min = min(ntrials),
    .max = max(ntrials))
  

blind_trials <- meta_clean %>% 
  drop_na(ntrials) %>% # removing Pan (2014)
  mutate(new_blind = ifelse(blinding_paradigm == "Sandwitch Masking",
                            "Backward Masking", blinding_paradigm),
         drop_percentage = (n_trials_valid - ntrials) / n_trials_valid * 100) %>% 
  group_by(new_blind) %>% 
  summarise(
    .mean = mean(ntrials),
    .mean_drop = mean(drop_percentage))

trials_tables <- list(
  summ_trials = summ_trials,
  blind_trials = blind_trials
)

saveRDS(trials_tables, "objects/trials_tables.rds")


# Prediction Interval -----------------------------------------------------

pred_interval <- bind_rows(
  get_pred_interval(mod_list$fit0),
  get_pred_interval(mod_list$fit0_nogray),
) %>% 
  mutate(mod = c("Overall model", "Overall model (only published literature)"))
  
saveRDS(pred_interval, "objects/pred_interval.rds")


# Power Contour -----------------------------------------------------------

p <- seq(0.5, 1, 0.1)
ntrials <- seq(10, 200, 1)
chance <- 0.5

sim_trials <- expand_grid(
  p,
  ntrials
)

sim_trials$h <- pwr::ES.h(sim_trials$p, chance)

get_power <- function(h, n){
  pwr::pwr.p.test(h, n)$power
}

sim_trials <- mutate(sim_trials, power = map2_dbl(h, ntrials, get_power))

saveRDS(sim_trials, "objects/sim_trials.rds")

# Sensitivity Analysis ----------------------------------------------------

clean_sens <- function(data){
  data %>% 
    map(., get_params_posterior) %>% 
    map(., compute_posterior_summary) %>% 
    bind_rows(., .id = "mod") %>% 
    clean_param_names(., col_param = .param)
}