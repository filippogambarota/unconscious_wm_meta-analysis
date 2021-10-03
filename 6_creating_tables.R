## -----------------------------------------------------------------------------
## Script: Creating Tables
##
## Author: Filippo Gambarota
##
## -----------------------------------------------------------------------------

# This script create tables included in the paper. The idea is to use the
# flextable package and create word documents with each table. The script also
# create a pdf version of all tables.

# Packages ----------------------------------------------------------------

library(flextable)
library(tidyverse)
library(officer)
library(broomExtra)
library(tidybayes)

# Loading Data ------------------------------------------------------------

meta_clean <- read_rds("data/cleaned/meta_analysis_data.rds")

# Loading Models ----------------------------------------------------------

mod_list <- load_models()

# Meta-analysis table -----------------------------------------------------

# Setup table

sect_properties <- prop_section(
  page_size = page_size(orient = "landscape"),
  type = "continuous",
  page_margins = page_mar()
)

# Create

meta_clean %>% 
  select(first_author, year, published,
         paradigm, blinding_paradigm, target_duration, target_type, retention_interval,
         sample_size, ntrials, aware_question,
         mean_acc, eff_size_g) -> tab

col_names <- c("Authors", "Year", "Publication Status", "WM Paradigm",
               "Blinding Paradigm", "Target Duration (ms)", "Target", "WM Delay (ms)",
               "Sample Size", "Trials", "Awareness Question", "Accuracy", "Hedges's g")

colnames(tab) <- col_names

meta_clean %>%
  select(first_author, year, published, sample_size, ntrials,
         paradigm, blinding_paradigm, acc_measure,
         target_type, task_type, eff_size_g, eff_size_se_g) %>% 
  mutate(author_year = paste0(first_author, " (", year, ")")) %>% 
  select(-first_author, -year) %>% 
  select(author_year, everything()) %>% 
  mutate(blinding_paradigm = case_when(blinding_paradigm == "Backward Masking" ~ "BM",
                                       blinding_paradigm == "Metacontrast Masking" ~ "MM",
                                       blinding_paradigm == "Continuous Flash Suppression" ~ "CFS",
                                       blinding_paradigm == "Sandwitch Masking" ~ "SM",
                                       blinding_paradigm == "Attentional Blink" ~ "AB")) %>% 
  rename("Study" = author_year,
         "Published" = published,
         "Sample Size" = sample_size,
         "Mean Trials" = ntrials,
         "WM Paradigm" = paradigm,
         "Blinding" = blinding_paradigm,
         "Acc Measure" = acc_measure,
         "Target" = target_type,
         "WM Task" = task_type,
         "Hedges's g" = eff_size_g,
         "Hedges's g SE" = eff_size_se_g) %>% 
  arrange(Study) %>% 
  flextable() %>% 
  colformat_double(digits = 2) %>% 
  theme_vanilla() %>% 
  autofit() %>% 
  line_spacing(space = 2) %>% 
  align(align = "center", part = "all") %>% 
  fontsize(size = 8, part = "body") %>% 
  fontsize(size = 9, part = "header") %>% 
  merge_v("Study") %>% 
  footnote(i = 1, j = c("WM Paradigm", "Blinding"),
           value = as_paragraph(
             c("CDT = Change Detection Task, DD = Delayed Detection, DMS = Delayed Match-to-sample Task, DET = Delayed Estimation Task",
               "CFS = Continuous Flash Suppression, AB = Attentional Blink, BM = Backward Masking, SM = Sandwich Masking, MM = Metacontrast Masking")
           ),
           ref_symbols = c("a", "b"),
           part = "header") -> studies_table_paper
  
save_as_docx(studies_table_paper, path = "tables/studies_table_paper.docx", pr_section = sect_properties)

# Mod Table ---------------------------------------------------------

# Definining table properties

sect_properties <- prop_section(
  page_size = page_size(orient = "landscape"),
  type = "continuous",
  page_margins = page_mar()
)

# Loading Table

rope_tab <- read_rds("objects/rope_table.rds")

mod_to_keep <- c("fit_author",
                 "fit_blind",
                 "fit_paradigm",
                 "fit_targetdur_bin",
                 "fit_trial_sample",
                 "fit0",
                 "fit0_nogray",
                 "fit_impl")

rope_tab <- rope_tab %>% 
  filter(mod %in% mod_to_keep)

rope_tab %>% 
  clean_param_names(.param) %>% 
  clean_mod_names(mod) %>% 
  mutate(within_rope = ifelse(startsWith(.param, "Tau"), "--", round(within_rope, 2)),
         outside_rope = ifelse(startsWith(.param, "Tau"), "--", round(outside_rope, 2)),
         within_rope = ifelse(mod == "~ Instructions", "--", within_rope),
         outside_rope = ifelse(mod == "~ Instructions", "--", outside_rope)) %>% 
  select(-.width, -.point, -.interval) %>% 
  rename("Model" = mod,
         "Parameter" = .param,
         "Effect" = .value,
         "SE" = .se,
         "Low" = .lower,
         "High" = .upper, 
         "% Within" = within_rope,
         "% Outside"= outside_rope)  %>% 
  arrange(Model) %>%
  select(Model, Parameter, Effect, SE, everything()) %>% 
  flextable() %>% 
  add_header(Model = "", 
             Parameter = "",
             Effect = "",
             SE = "",
             Low = "89% HPDI",
             High = "89% HPDI",
             `% Within` = "ROPE",
             `% Outside` = "ROPE") %>% 
  merge_h(part = "header") %>% 
  merge_v(1) %>% 
  theme_vanilla() %>% 
  colformat_double(digits = 2) %>% 
  align(align = "center", part = "header") %>% 
  align(align = "center", part = "body") %>% 
  autofit() -> rope_table_paper

save_as_docx(rope_table_paper, path = "tables/rope_table_paper.docx")

# LOO Table ---------------------------------------------------------------

loo_tab <- read_rds("objects/loo_tab.rds")
r2_tab <- read_rds("objects/r2_tab.rds")

r2_tab %>% 
  clean_mod_names(mod) %>% 
  rename("Model" = mod) %>% 
  right_join(., loo_tab, by = "Model") %>% 
  select(Model, Estimate, SE, Weight, r2_2, r2_3) %>% 
  mutate(r2_2 = ifelse(r2_2 < 0, 0, round(r2_2, 2)),
         r2_3 = ifelse(r2_3 < 0, 0, round(r2_3, 2))) %>% 
  flextable() %>%
  flextable::compose(part = "header", j = "r2_2", 
          value = as_paragraph("R",  as_sup("2"), as_sub("level-2"))) %>% 
  flextable::compose(part = "header", j = "r2_3", 
                     value = as_paragraph("R", as_sup("2"), as_sub("level-3"))) %>% 
  theme_vanilla() %>% 
  colformat_double(digits = 2) %>% 
  align(align = "center", part = "header",) %>% 
  align(align = "center", part = "body", j = c(2,3,4,5,6)) %>% 
  autofit() -> loo_tab_paper

save_as_docx(loo_tab_paper, path = "tables/loo_tab.docx")

# Saving Tables -----------------------------------------------------------

tables <- list(
  studies_table_paper = studies_table_paper,
  rope_table_paper = rope_table_paper,
  loo_tab_paper = loo_tab_paper
)

saveRDS(tables, "objects/tables_report.rds")

# Rendering Table ---------------------------------------------------------

# This create a PDF version of table included in the papers

pagedown::chrome_print("tables/tables.Rmd")