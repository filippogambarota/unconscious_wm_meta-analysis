# prep_forest -------------------------------------------------------------

# This function takes the model and the dataset and returns a list of objects
# for plotting the overall forest plot.

prep_forest <- function(mod, data){
  
  meta_weight <- data %>% 
    arrange(first_author) %>% 
    mutate(paper = paste0(first_author, " (", year, ")")) %>% 
    select(paper, eff_size_g, eff_size_variance_g) %>% 
    mutate(weight = 1/eff_size_variance_g)
  
  post <- add_fitted_draws(data %>% select(first_author, year, study_meta, substudy, eff_size_se_g), mod)
  
  pred_interval <- get_pred_interval(mod)
  
  avg_summ <- spread_draws(mod, b_Intercept) %>% 
    summarise(median_hdci(b_Intercept, .width = 0.89)) %>% 
    mutate(paper = "Average") %>% 
    select(paper, everything())
  
  avg <- spread_draws(mod, b_Intercept) %>% 
    mutate(paper = "Average") %>% 
    select(paper, everything())
  
  post %>% 
    ungroup() %>% 
    mutate(paper = paste0(first_author, " (", year, ")")) %>% 
    group_by(paper) %>% 
    summarise(median_hdci(.value, .width = 0.89)) %>% 
    bind_rows(avg_summ, .) %>% 
    mutate_if(.predicate = is.numeric, round, 2) %>% 
    mutate(published = if_else(startsWith(paper, "Barton") | startsWith(paper, "Taglia"),
                               "U", "P"),
           col = if_else(paper == "Average", "col1", "col2"),
           paper = as_factor(paper)) -> post
  
  out <- list("avg_summ" = avg_summ, "avg" = avg, "post" = post, "weight" = meta_weight, pred_interval = pred_interval)
  return(out)
}

# plot_forest -------------------------------------------------------------

# This is the actual plotting function that takes the list of the prep_forest()
# function and return the ggplot object

plot_forest <- function(meta_list, size_fac = 1){
  
  ggplot() +
    geom_vline(xintercept = meta_list$avg_summ$y, # mean
               col = "darkgrey", size = 0.9, linetype = "dashed") +
    geom_segment(aes(x = meta_list$pred_interval$lower, 
                     xend = meta_list$pred_interval$upper,
                     y = 0.2, 
                     yend = 0.2),
                 size = 3,
                 color = "darkgray") +
    geom_pointrange(data = meta_list$post %>% filter(!paper == "Average"),
                    aes(x = y, y = paper, xmin = ymin, xmax = ymax,
                        col = published, shape = published), size = 1.2) +
    stat_halfeye(data = meta_list$avg, aes(x = b_Intercept, y = paper), fill = "lightblue", size = 4) +
    coord_cartesian(xlim = c(-0.5, 2.3)) +
    geom_label(data = meta_list$post %>% filter(!paper == "Average"),
               aes(y = paper, label = glue::glue("{y} [{ymin}, {ymax}]"), x = Inf),
               hjust = "inward", size = 5, fill = "white", color = "black", label.size = NA,
               family = "TT Times New Roman") + 
    geom_label(data = meta_list$post %>% filter(paper == "Average"),
               aes(y = paper, label = glue::glue("{y} [{ymin}, {ymax}]"), x = Inf),
               hjust = "inward", size = 5, fill = "white", color = "black", label.size = NA, 
               family = "TT Times New Roman") +
    geom_point(data = meta_list$weight, 
               aes(x = eff_size_g, y = paper, size = weight),
               show.legend=FALSE,
               position = position_nudge(y = -.35),
               alpha = 0.7, color = "darkgray") +
    scale_y_discrete(breaks = levels(meta_list$post$paper),
                     limits = c(levels(meta_list$post$paper)[1],
                                "skip",
                                levels(meta_list$post$paper)[-1])) +
    xlim(c(-0.5, 4)) +
    scale_color_manual(values = c("#e41a1c", "#377eb8")) +
    scale_shape_manual(values = c(15,16)) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 20*size_fac),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 15*size_fac, color = "black"),
          axis.text.x = element_text(size = 15*size_fac),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=10*size_fac),
          legend.position = "bottom",
          legend.direction = "horizontal",
          text = element_text(family = "TT Times New Roman")) +
    xlab(latex2exp::TeX("Hedges' $g$"))
}


# meta_reg_plot_h -----------------------------------------------------------

# this is the plotting function for metaregression and ROPE analysis. The
# function takes the posterior distrubutions as raw and summary and the desired
# ROPE range and return the ggplot object

meta_reg_plot_h <- function(post, ROPE, fac_size = 1){
  post %>% 
    ggplot(aes(x = .param, y = .value, fill = stat(y > ROPE[2]))) +
    geom_hline(yintercept = ROPE[2], col = "red", size = 1) +
    stat_halfeye(.width = c(0.89), size = 5) +
    cowplot::theme_minimal_hgrid() +
    coord_cartesian(ylim = c(-2,2)) +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 25*fac_size),
          axis.text.x = element_text(size = 20*fac_size),
          axis.text.y = element_text(size = 20*fac_size),
          axis.title.x = element_blank(),
          plot.margin=unit(c(1,1,1.5,1.2),"cm"),
          text = element_text(family = "TT Times New Roman")) +
    ylab("Hedges's g") +
    scale_fill_manual(values = c("lightgrey", "lightblue"))
}


# remove_axes -------------------------------------------------------------

# wrapper for removing axes from composite plots

remove_axes <- function(axes = "x"){
  if(axes == "x"){
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  }else{
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank())
  }
}

# remove_title -------------------------------------------------------------

# wrapper for removing title from composite plots

remove_title <- function(axes = "x"){
  if(axes == "x"){
    theme(axis.title.x = element_blank())
  }else{
    theme(axis.title.y = element_blank())
  }
}