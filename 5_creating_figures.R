## -----------------------------------------------------------------------------
## Script: Creating Figures
##
## Author: Filippo Gambarota
##
## -----------------------------------------------------------------------------

# This script create all the relevant figures for the papers and save as ".png"
# eps or tiff or ".rds" files to use in the supplementary materials.

# Packages ----------------------------------------------------------------

library(tidyverse)
library(gridExtra)
library(latex2exp)
library(tidybayes)
library(metafor)
library(PRISMAstatement)
library(cowplot)
library(ggthemes)
library(metaviz)
library(pdftools)
library(brms)

# Loading Data ------------------------------------------------------------

ROPE <- c(-Inf, 0.1)

meta_clean <- read_rds("data/cleaned/meta_analysis_data.rds")

# Loading Models

mod_list <- load_models()

# Functions ---------------------------------------------------------------

source("functions/plotting_functions.R")

save_plot <- function(plot, name, height, width){
  
  # Word file
  ggsave(plot = plot, 
         filename = paste0("figures/", name, ".png"), 
         dpi = "retina", 
         device = "png", 
         height = height, 
         width = width)
  # Upload
  ggsave(plot = plot, 
         filename = paste0("figures/", name, ".tiff"), 
         dpi = "retina", 
         device = "tiff", 
         height = height, 
         width = width)
}

# Priors Plot -------------------------------------------------------------

ggplot(NULL, aes(c(0,1))) + 
  geom_area(stat = "function", fun = extraDistr::dhcauchy, args = list(0.5),
            fill = "lightblue", alpha = 0.7, color = "black", size = 0.5) +
  xlab("Value") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 15)) +
  ggtitle(TeX("$\\tau$ Prior Distribution")) -> tau_prior

ggplot(NULL, aes(c(-7,7))) + 
  geom_area(stat = "function", fun = dnorm, args = list(0, 2),
            fill = "lightblue", alpha = 0.7, color = "black", size = 0.5) +
  xlab("Value") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 15)) +
  ggtitle(TeX("$\\mu$ Prior Distribution")) -> mu_prior

prior_plot <- cowplot::plot_grid(mu_prior, tau_prior, align = "v")

saveRDS(prior_plot, file = "objects/prior_plot.rds")

# Fit0 forest plot --------------------------------------------------------

fit0_list <- prep_forest(mod_list$fit0, meta_clean)

fit0_forest <- plot_forest(fit0_list, size_fac = 1.3)

saveRDS(fit0_forest, file = "objects/fit0_forest.rds")

save_plot(fit0_forest, "fit0_forest", height = 16, width = 13)

# Fit0 Publication Bias --------------------------------------------------------

meta_clean_metafor <- meta_clean %>% 
  select(study_meta, substudy, eff_size_g, eff_size_variance_g) %>% 
  rename("yi" = eff_size_g,
         "vi" = eff_size_variance_g)

# Fitting the model with metafor for using the funnel function

fit0_metafor <- rma.mv(yi = yi, vi, random = ~1|study_meta/substudy,
                       data = meta_clean_metafor, method = "REML")

fit0_funnel <- funnel(fit0_metafor)

meta_clean %>% 
  select(study_meta, substudy, eff_size_g, eff_size_se_g, published) %>% 
  mutate(published = ifelse(published == "Published", "P", "U")) %>% 
  as.data.frame() -> temp_funnel

fit0_funnel <- metaviz::viz_funnel(temp_funnel[, 3:4],
                                   sig_contours = F,
                                   detail_level = 10,
                                   point_size = 3,
                                   group_legend = T,
                                   xlab = latex2exp::TeX("$\\Hedges's \\; g"),
                                   ylab = latex2exp::TeX("$\\Standard \\; Error"))

fit0_funnel +
  geom_point(data = temp_funnel, aes(x = eff_size_g, y = eff_size_se_g, color = published), size = 4) +
  cowplot::theme_minimal_hgrid() +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none") +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 0.7,
            linetype = "dashed",
            alpha = 0.7) -> fit0_funnel_egger

save_plot(fit0_funnel_egger, "fit0_pub_bias", height = 8, width = 8)
saveRDS(fit0_funnel_egger, file = "objects/fit0_funnel_egger.rds")

# PRISMA ------------------------------------------------------------------

prisma <- PRISMAstatement::prisma(
  found = 1038,
  found_other = 2,
  no_dupes = 757,
  screened = 757,
  screen_exclusions = 703,
  full_text = 54, 
  full_text_exclusions = 41,
  quantitative = 13,
  qualitative = 13,
  font_size = 13,
  width = 200,
  height = 200
)

PRISMAstatement:::prisma_pdf(prisma, filename = "figures/prisma_flowchart.pdf")
pdftools::pdf_convert("figures/prisma_flowchart.pdf", 
                      format = "tiff", 
                      dpi = 300,
                      filenames = "figures/prisma_flowchart.tiff")

# ROPE Plot ---------------------------------------------------------------

mod_sel <- c("fit_author", "fit_blind", "fit_paradigm", "fit_targetdur_bin", "fit_impl")

mod_list_sel <- mod_list[names(mod_list) %in% mod_sel]

post_list <- mod_list_sel %>% 
  map(., get_params_posterior, only_betas = TRUE) %>% 
  map(., clean_param_names, .param)

# Tweaks specific stuff

post_list$fit_targetdur_bin <- post_list$fit_targetdur_bin %>% 
  mutate(.param = factor(.param, levels = c("16-50ms", "150-500ms", "3000ms")))

post_list$fit_impl <- post_list$fit_impl %>% 
  mutate(.param = ifelse(.param == "Implicit - Explicit", "Impl - Expl", .param))

# Clean the fit impl

post_list$fit_impl <- post_list$fit_impl %>%
  mutate(.param = ifelse(.param == "Average", "Explicit", .param)) %>% 
  pivot_wider(names_from = ".param", values_from = .value) %>%
  mutate("Implicit" = Explicit + `Impl - Expl`) %>%
  pivot_longer(4:ncol(.), names_to = ".param", values_to = ".value") %>%
  mutate(.fct_order = ifelse(.param == "Implicit", 1,
                             ifelse(.param == "Implicit", 2, 3)),
         .param = fct_reorder(.param, .fct_order))

post_plot <- map(post_list, meta_reg_plot_h, ROPE, 1.1)

rope_plot <- plot_grid(
  post_plot$fit_blind,
  post_plot$fit_impl + remove_axes(axes = "y"),
  post_plot$fit_paradigm,
  post_plot$fit_targetdur_bin + remove_axes(axes = "y"),
  labels = "AUTO",
  label_size = 20
)

author_plot <- meta_reg_plot_h(post_list$fit_author, ROPE, fac_size = 1.1) +
  scale_x_discrete(guide = guide_axis(angle = 45))

saveRDS(rope_plot, file = "objects/rope_plot.rds")
saveRDS(author_plot, file = "objects/author_plot.rds")

save_plot(rope_plot, "rope_plot", height = 14, width = 14)
save_plot(author_plot, "author_plot", height = 10, width = 12)

# trial_figure_example -------------------------------------------------------

# This require inkscape installed and activated as command line tool

system("inkscape --export-type='eps' figures/example_trial.svg")