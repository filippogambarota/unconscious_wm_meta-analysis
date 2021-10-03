###### Modelling and post-processing functions

# meta_r2 -----------------------------------------------------------------

# Compute a R2 value at level-two and level three (using the formula 
# provided by Cheung (2014) using the tau estimation from a null model 
# and the model with a predictor. As default the function use the posterior
# median of taus distribution.

meta_r2 <- function(mod_null, mod){
  taus_mod_null <- get_tau(mod_null)
  taus_mod <- get_tau(mod)
  
  # Computing Medians
  taus_mod_null <- apply(taus_mod_null, 2, median, simplify = F)
  taus_mod <- apply(taus_mod, 2, median, simplify = F)
  rsquared_2 <- 1 - taus_mod$tau2/taus_mod_null$tau2
  rsquared_3 <- 1 - taus_mod$tau3/taus_mod_null$tau3
  
  # Returning a list
  out <- list(rsquared_2 = rsquared_2, rsquared_3 = rsquared_3)
  return(out)
}

# meta_icc ----------------------------------------------------------------

# Calculating ICC at level-two and level-three using the formula provided 
# by Cheung 2014

meta_icc <- function(mod, summarise = FALSE){
  out <- get_tau(mod) # extract posterior distribution taus
  out$icc2 <- with(out, tau2 / (tau2 + tau3))
  out$icc3 <- with(out, tau3 / (tau2 + tau3))
  
  return(out)
}

# get_tau -----------------------------------------------------------------

# Utility function to extract posterior distributions of taus from a fitted
# model

get_tau <- function(mod){
  tau3 <- posterior_samples(mod, "sd_study_meta__Intercept")[[1]]^2
  tau2 <- posterior_samples(mod, "sd_study_meta:substudy__Intercept")[[1]]^2
  out <- tibble::tibble(tau2, tau3)
  return(out)
}

get_tau_tidy <- function(mod){
  spread_draws(mod, `sd_.*`, regex = TRUE) %>% # this return only taus
    pivot_longer(starts_with("sd"), names_to = ".param", values_to = ".value")
  
}

# get_mu_se ---------------------------------------------------------------

# Return the average effect standard error from a given model

get_mu_se <- function(mod){
  mu_se <- fixef(mod)[2]
  return(mu_se)
}

get_mu <- function(mod){
  mu <- fixef(mod)[1]
  return(mu)
}

# get_pred_interval -------------------------------------------------------

# return the prediction interval given a model and the width using the metafor
# function check https://www.metafor-project.org/doku.php/faq#for_random-effects_models_fitt

get_pred_interval <- function(mod, width = 0.89){
  taus <- get_tau(mod)
  taus <- apply(taus, 2, median)
  mu_se <- get_mu_se(mod)
  mu <- get_mu(mod)
  sqrt_var <- sqrt(mu_se^2 + sum(taus^2))
  z_val <- abs(qnorm((1-width)/2))
  out <- tibble(
    mu = mu, 
    mu_se = mu_se, 
    lower = mu - z_val*sqrt_var,
    upper = mu + z_val*sqrt_var,
    width = paste0(width*100, "%")
  )
  return(out)
}

# add_ic_and_save ---------------------------------------------------------

# Load a model, add information criteria and save again the model

add_ic_and_save <- function(file, criteria = c("loo", "waic")){
  mod <- readr::read_rds(file)
  mod <- add_criterion(mod, criteria, moment_match = T)
  cli::cli_alert_success(paste("Criteria added to model", basename(file), "and saved!"))
}

#### Utils for cleaning names

# clean_mod_names ---------------------------------------------------------

# these functions prepares the tables names for a nice output. The syntax is written
# using tidy-evaluation

clean_mod_names <- function(tab, col_mod){
  tab %>% 
    mutate(!!enquo(col_mod) := case_when(!!enquo(col_mod) == "fit_author" ~ "~ Author",
                                     !!enquo(col_mod) == "fit_blind" ~ " ~ Blinding Paradigm",
                                     !!enquo(col_mod) == "fit_paradigm" ~ " ~ WM Paradigm",
                                     !!enquo(col_mod) == "fit_targetdur_bin" ~ " ~ Target Duration",
                                     !!enquo(col_mod) == "fit_trial_sample" ~ " ~ Trials + Sample Size",
                                     !!enquo(col_mod) == "fit0" ~ "Overall Model",
                                     !!enquo(col_mod) == "fit0_nogray" ~ "Overall Model (only published)",
                                     !!enquo(col_mod) == "fit_impl" ~ "~ Instructions",
                                      TRUE ~ !!enquo(col_mod)))
}

clean_param_names <- function(tab, col_param, col_mod = NULL){
  tab %>% 
    mutate(
      !!enquo(col_param) := stringr::str_remove(!!enquo(col_param), "b_first_author_mod"),
      !!enquo(col_param) := stringr::str_remove(!!enquo(col_param), "b_blinding_paradigm_m"),
      !!enquo(col_param) := stringr::str_remove(!!enquo(col_param), "b_paradigm"),
      !!enquo(col_param) := stringr::str_remove(!!enquo(col_param), "b_target_duration_bin"),
      !!enquo(col_param) := case_when(!!enquo(col_param) == "AttentionalBlink" ~ "AB",
                                      !!enquo(col_param) == "b_Intercept" ~ "Average",
                                      !!enquo(col_param) == "BackwardMasking" ~ "BM",
                                      !!enquo(col_param) == "ContinuousFlashSuppression" ~"CFS",
                                      !!enquo(col_param) == "MetacontrastMasking" ~ "MM",
                                      !!enquo(col_param) == "150M500ms" ~ "150-500ms",
                                      !!enquo(col_param) == "16M50ms" ~ "16-50ms",
                                      !!enquo(col_param) == "sd_study_meta__(Intercept)" ~ "Tau Paper",
                                      !!enquo(col_param) == "sd_study_meta:substudy__(Intercept)" ~ "Tau Study",
                                      !!enquo(col_param) == "sd_study_meta__Intercept" ~ "Tau Paper",
                                      !!enquo(col_param) == "sd_study_meta:substudy__Intercept" ~ "Tau Study",
                                      !!enquo(col_param) == "b_Intercept" ~ "Average",
                                      !!enquo(col_param) == "b_(Intercept)" ~ "Average",
                                      !!enquo(col_param) == "b_ntrials0" ~ "Trials",
                                      !!enquo(col_param) == "b_sample_size0" ~ "Sample Size",
                                      !!enquo(col_param) == "b_eff_size_se_g" ~ "Standard Error",
                                      !!enquo(col_param) == "b_implicityes" ~ "Implicit - Explicit",
                                      TRUE ~ !!enquo(col_param))
    )
}

# Loo Cleaning ------------------------------------------------------------

# these functions prepares the tables for the loo-sensitivity analysis

loo_clean <- function(mod_list, conf.method = "HPDinterval", conf.level = 0.89) {
  
  temp_list <- lapply(mod_list, function(x) {
    broomExtra::tidy(x, parameters = c("^sd_", "^b_"),
                     robust = T, rhat = T, ess = T, conf.method = conf.method, conf.level = conf.level)
  })
  
  mod_tab <- bind_rows(temp_list, .id = "mod") %>% 
    mutate(term = ifelse(startsWith(term, "b_"), "b",
                         ifelse(startsWith(term, "sd_study_meta__(Intercept)"),
                                "tau_study", "tau_substudy")))
  
  return(mod_tab)
  
}

loo_post_clean <- function(mod_list) {
  
  temp_list <- lapply(mod_list, function(x) {
    spread_draws(x, b_Intercept, sd_study_meta__Intercept)
  })
  
  mod_tab <- bind_rows(temp_list, .id = "mod") %>% 
    rename_at(vars(starts_with("sd_")), ~ "tau") %>% 
    rename_at(vars(starts_with("b_")), ~ "b")
  
  return(mod_tab)
  
}

# Model Summary -----------------------------------------------------------

# create a model summary with all relevant information

get_model_table <- function(mod, conf_method = "HPDinterval", conf_level = 0.89){
  
  temp <- broom.mixed::tidy(mod, parameters = c("^sd_", "^b_"), 
                            conf.level = conf_level,
                            conf.method = conf_method)
  temp$conf_level <- conf_level
  
  # Check if the model is a two or a three-level model
  
  if(length(summary(mod)$random) > 1){
    tau1 <- summary(mod)$random[[1]][, c("Rhat", "Bulk_ESS", "Tail_ESS")]
    tau2 <- summary(mod)$random[[2]][, c("Rhat", "Bulk_ESS", "Tail_ESS")]
    tau = rbind(tau1, tau2)
    
  }else{
    tau <- summary(mod)$random[[1]][, c("Rhat", "Bulk_ESS", "Tail_ESS")]
  }
  
  # Get diagnostic

  fix <- summary(mod)$fixed[, c("Rhat", "Bulk_ESS", "Tail_ESS")]
  
  diagn <- rbind(tau, fix)
  
  out <- tibble(cbind(temp, diagn))
  
  return(out)
}

# Posterior Analysis ------------------------------------------------------

get_params_posterior <- function(mod, only_betas = FALSE){
  if(only_betas){
    post <- spread_draws(mod, `b_.*`, regex = TRUE)
  }else{
    post <- spread_draws(mod, `(b_.*)|(sd_.*)`, regex = TRUE)
  }
  post %>% 
    ungroup() %>% 
    # this create a long format of all betas
    pivot_longer(4:ncol(.), names_to = ".param", values_to = ".value")
}

get_rope_table <- function(params_posterior, ROPE, ci = 0.95, ci_method = "HDI"){
  
  params_posterior %>% 
    group_by(.param) %>%
    summarise(median_hdi(.value, .width = 0.89),
              .se = sd(.value),
              within_rope = rope(.value, range = c(-Inf, 0.1))$ROPE_Percentage,
              outside_rope = 1-within_rope) %>% 
    rename(".value" = y,
           ".lower" = ymin,
           ".upper" = ymax)
}

compute_posterior_summary <- function(posterior, ci = 0.89){
  posterior %>% 
    group_by(.param) %>%
    summarise(median_hdi(.value, .width = 0.89),
              .se = sd(.value),
              median_hdi(.value, .width = ci)) %>% 
    rename(".value" = y,
           ".lower" = ymin,
           ".upper" = ymax) %>% 
    select(.param, .value, everything())
}

# Sensitivity Analysis ----------------------------------------------------

# get_data_without_id -----------------------------------------------------

# Given a dataset and a grouping factor, the function return a list of datasets
# with one level of grouping factor removed. In addition append the full dataset
# for comparison

get_data_without_id <- function(data, grouping_factor){
  id_to_remove <- unique(data[[grouping_factor]]) # get id
  
  # remove all ids once
  out <- map(id_to_remove, function(id) {
    data %>% filter(!!ensym(grouping_factor) != id)
  })
  
  # append the full dataset
  out[[length(out) + 1]] <- data
  
  names(out) <- c(id_to_remove, "Overall")
  
  return(out)
}

prior_est_change <- function(mod, mod_list){
  
  ref_mod <- fixef(fit0, robust = T)[, "Estimate"]
  
  est = unlist(map(bf_fit0, function(x) fixef(x, robust = T)[, "Estimate"]))
  prior = names(est)
  
  sens_mod <- tibble(est, prior)
  
  sens_mod$change <- round(((est - ref_mod)/ref_mod) * 100, 2)
  
  return(sens_mod)
}


# refit_model -------------------------------------------------------------

# This function is a wrapper of brsm::update() for setting useful default
# options

refit_model <- function(mod, newdata = NULL, prior = NULL, seed){
  update(mod, newdata = newdata, cores = 4, prior = prior, seed = seed)
}

# get_prior_list ----------------------------------------------------------

# This funcion return a list of prior objects for tau and the avg effect
# from a list of sds

create_prior <- function(type, sd){
  
  if(type == "tau"){
    dist = "cauchy"
    class = "sd"
    prior_chr <- sprintf("%s(0, %s)", dist, sd) 
    out <- prior_string(prior_chr, class = class)
  }else{
    dist = "normal"
    class = "b"
    coef = "Intercept"
    prior_chr <- sprintf("%s(0, %s)", dist, sd) 
    out <- prior_string(prior_chr, class = class, coef = coef)
  }

  return(out)
}