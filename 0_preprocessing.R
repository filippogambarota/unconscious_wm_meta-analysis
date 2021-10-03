## -----------------------------------------------------------------------------
## Script: Pre-processing
##
## Author: Filippo Gambarota
##
## -----------------------------------------------------------------------------

# This script create the meta-analysis table for computing all the analyses

# Packages ----------------------------------------------------------------

library(tidyverse)
library(stringr)

# Importing Data ----------------------------------------------------------

meta <- read_csv("data/raw/meta_table_raw.csv", na = "NA")

# Dataset Cleaning --------------------------------------------------------

meta_clean <- meta %>% 
  mutate(subresult = paste0(substudy, "_", subresult), # recode subresult variable
         retention_interval = unlist(retention_interval),
         year = unlist(year),
         study_id = paste(first_author, substudy), # study index for plots
         study_meta = paste(first_author, year),
         st_error = ifelse(is.na(st_error), # All papers
                           yes = (mean_acc - chance_level)/t_value,
                           no = st_error),
         st_error = ifelse(is.na(st_error), # Nakano Paper
                           yes = st_deviation/sqrt(sample_size),
                           no = st_error),
         st_deviation = st_error * sqrt(sample_size),
         
         # Compute Effect Size 
         
         eff_size = ifelse(is.na(eff_size),
                           yes = (mean_acc - chance_level)/st_deviation,
                           no = eff_size)) %>% 
  group_by(substudy) %>% 
  mutate(ntrials = sum(n_trials_1), # calculate mean trials for each substudy
         nresult = n(), # number of results for each study
         mean_acc = mean(mean_acc, na.rm = T), # calculate mean accuracy for each substudy
         eff_size_variance = (1/sample_size)+eff_size^2/(2*sample_size)) %>% # Calculate Effect Size Variance according to Wilson
  ungroup()

# Compute Hedges g --------------------------------------------------------

# Hedges g is the corrected version of the cohen d for small sample size. See Borenstein (2009, pp. 27)

meta_clean %>% 
  mutate(hedges_corr =  1-(3/(4*(sample_size - 1) - 1)), # correction factor
         eff_size_g = eff_size * hedges_corr, # adjusting effect size
         eff_size_variance_g = hedges_corr^2 * eff_size_variance) -> meta_clean # adjusting variance

# Compute the effect size standard error for the bayesian model -----------

meta_clean %>% 
  mutate(eff_size_se_g = sqrt(eff_size_variance_g),
         eff_size_se = sqrt(eff_size_variance),
         weight_se = 1/eff_size_se,
         weight_var = 1/eff_size_variance_g) -> meta_clean

# Computing Pooled Variance -----------------------------------------------

# These functions compute the formula for pooling different subresults within the same study
# See Borenstein (2009, pp 228)

boren_method12 <- function(x, r = 0.5){
  (1/length(x))^2 * (sum(x) + (2*r*prod(sqrt(x))))
}

meta_clean %>% 
  filter(nresult < 3) %>% 
  group_by(substudy) %>% 
  mutate(eff_size_g = mean(eff_size_g)) %>% # combine effect size
  group_modify(~ mutate(.,
                        eff_size_variance_g = if_else(.$collapse == "YES",
                                                         true = boren_method12(eff_size_variance_g), # combine effect size variance
                                                         false = eff_size_variance_g)), .keep = F) %>% 
  ungroup() -> meta_12res

meta_clean <- meta_12res

# Combine Final Dataset ---------------------------------------------------

meta_clean %>% 
  select(study, substudy, first_author, first_author_mod, study_id, study_meta, published, year,
         paradigm, blinding_paradigm, target_duration, target_type, task_type,
         treshold_stimulus, aware_question, retention_interval, implicit,
         sample_size, ntrials, n_trials_valid, mean_acc, chance_level, acc_measure,
         eff_size_g, eff_size_variance_g) %>% 
  distinct() %>% 
  mutate(tot_trials = sample_size * ntrials, # estimated number of total trials
         ntrials0 = ntrials - mean(ntrials, na.rm = T), # centering trials
         sample_size0 = sample_size - mean(sample_size), # centering sample size
         eff_size_se_g = sqrt(eff_size_variance_g)) %>% # standard error for brms
  arrange(study) -> meta_clean

# Saving Data -------------------------------------------------------------

saveRDS(meta_clean, file = "data/cleaned/meta_analysis_data.rds")
write_csv(meta_clean, file = "data/cleaned/meta_analysis_data.csv")