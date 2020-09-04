library(tidyverse)
library(cowplot)
library(metR)
library(reshape2)
library(scales)
library(viridis)

#==============================================================================
# Read in data
#==============================================================================

# Contour plots
contour_path <- "figures/data/fig2_contour_1k_scenario.tsv.gz"
contour_data <- suppressMessages(read_tsv(contour_path)) %>% 
  filter(generation_alpha == 0.064)

#==============================================================================
# Specify plotting theme information
#==============================================================================

# Fonts
font <- "sans" # Main figure font
titlefont <- font # Font for axis titles etc. (if different)
fontsize_base <- 12 # Basic figure font size
fontscale_title <- 1.5 # Default axis-title scale relative to regular font
fontscale_main <- 1.5 # Default plot-title scale
fontscale_label <- 2 # Default subfigure-label scale (A, B, etc)
fontscale_legend <- 1

# Theme
theme_base <-   theme_bw() + theme(
  legend.position = "bottom",
  axis.text = element_text(size = fontsize_base, family = font, 
                           colour = "black"),
  axis.title.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(r=3,unit="mm")),
  axis.title.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t=5,unit="mm")),
  legend.text = element_text(size = fontsize_base * fontscale_legend, 
                             family = font, colour = "black"),
  legend.title = element_text(size = fontsize_base * fontscale_legend,
                              family = font, colour = "black",
                              face = "bold", vjust=0.5, 
                              margin=margin(r=3, unit="mm")),
  plot.title = element_text(size = fontsize_base * fontscale_main,
                            family = titlefont, colour="black",
                            hjust = 0.5, vjust = 0.5,
                            margin = margin(b=4, unit="mm"),
                            face = "bold"),
  plot.margin = margin(t=1.5, l=0.5, r=0.5, b = 0.5, unit="cm"),
  strip.background = element_blank(),
  strip.text.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t = 3, b = 3, unit="mm")),
  strip.text.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(l=3, r=3, unit="mm")),
  legend.justification = "center",
  panel.border = element_blank(),
  axis.line = element_line(size=0.5, colour="black")
)

#==============================================================================
# Labelling functions
#==============================================================================

label_backtrace <- function(x){
  ifelse(x==0, "Forward tracing only",
         "Bidirectional tracing")
  # "Forward and reverse\ntracing")
}

label_limits <- function(x){
  ifelse(as.numeric(x) == Inf, "No trace limit", paste0(x, "-day trace limit"))
}

label_limits_manual <- function(x){
  ifelse(as.numeric(x) == Inf, "No manual limit", paste0(x, "-day manual limit"))
}


label_vars <- function(x){
  ifelse(x=="p_data_sharing_auto", "% of cases sharing data",
         "% of cases with\ntrace-enabled smartphones")
}

label_pc <- function(x) round(as.numeric(x)*100)

#==============================================================================
# Define plotting functions
#==============================================================================

#------------------------------------------------------------------------------
# Contour plot palettes
#------------------------------------------------------------------------------

col_bidir <- "#d95f02"
col_fwd   <- "#1b9e77"
step_dense <- 0.05
step_sparse <- 0.1

make_palette <- function(colour_high, step_size, data = NULL,
                         contour_var = NULL, bt_dist = NULL, 
                         palette_max = NULL, palette_min = NULL,
                         reverse_scale = FALSE){
  # Calibrate range of palette
  if (is.null(palette_max) | is.null(palette_min)){
    data_range <- data %>% filter(backtrace_distance == bt_dist) %>%
      pull(contour_var) %>% unique
  }
  if (is.null(palette_max)){
    palette_max <- data_range %>% max %>%
      (function(x) ceiling(x/step_size)*step_size)
  }
  if (is.null(palette_min)){
    palette_min <- data_range %>% min %>%
      (function(x) floor(x/step_size)*step_size)
  }
  palette_val <- seq(palette_min, palette_max, step_size)
  # Calculate hue values
  palette_fn  <- scales::seq_gradient_pal(low="white", high=colour_high, "Lab")
  palette_hue <- palette_fn((seq(0,1,length.out=length(palette_val))))
  if (reverse_scale) palette_hue <- rev(palette_hue)
  return(tibble(level = palette_val, colour = palette_hue))
}

# Control rate
contour_palette_ctrl_bidir <- make_palette(col_bidir, step_sparse,
                                           palette_min = 0, palette_max = 1)
contour_palette_ctrl_fwd <- make_palette(col_fwd, step_sparse,
                                         palette_min = 0, palette_max = 1)

# R_eff
contour_palette_reff_bidir <- make_palette(col_bidir, step_dense,
                                           data = contour_data, bt_dist = Inf,
                                           reverse_scale = TRUE,
                                           contour_var = "effective_r0_mean")
contour_palette_reff_fwd <- make_palette(col_fwd, step_dense,
                                         data = contour_data, bt_dist = 0,
                                         reverse_scale = TRUE,
                                         contour_var = "effective_r0_mean")

#------------------------------------------------------------------------------
# Line plots
#------------------------------------------------------------------------------

ctrl_plot <- function(data, err_width = 0.05, pt_size = 2,
                      coord_ratio = 1){
  ggplot(data, aes(x=p_traced_auto, y=p_controlled,
                   ymin=p_controlled_lower, ymax=p_controlled_upper,
                   colour=factor(backtrace_distance),
                   linetype = factor(contact_limit_auto, levels=c(2,Inf)),
                   shape = factor(contact_limit_auto, levels=c(2,Inf)))) +
    geom_errorbar(width=err_width) + geom_line() + geom_point(size=pt_size) +
    coord_fixed(ratio = coord_ratio) + 
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_x_continuous(name = "Probability of trace\nsuccess (%)",
                       breaks = seq(0,1,0.2), labels = label_pc) +
    scale_colour_brewer(palette = "Dark2", labels=label_backtrace,
                        name="Trace type") +
    scale_linetype_discrete(name = NULL, labels = label_limits) +
    scale_shape_discrete(name = NULL, labels = label_limits) +
    guides(colour=guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base + theme(
      # legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      # legend.spacing.y = unit(legend_vspace, "mm"),
    )
}

reff_plot <- function(data, pt_size = 2, coord_ratio = 1){
  ggplot(data, aes(x=p_traced_auto, y=effective_r0_mean,
                   colour=factor(backtrace_distance),
                   linetype = factor(contact_limit_auto, levels=c(2,Inf)),
                   shape = factor(contact_limit_auto, levels=c(2,Inf)))) +
    geom_line() + geom_point(size=pt_size) +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    geom_hline(yintercept = 2.5, linetype = "dotted", size=1, colour = "red") +
    annotate("text", x=1, y=2.4, hjust=1, vjust = 1,
             colour = "red", size=4.6, label=expression(italic(R)[0])) +
    coord_fixed(ratio = coord_ratio/2.5) + 
    scale_y_continuous(name = expression(paste("Mean ",italic(R)[eff])), limits = c(0,NA),
                       breaks = seq(0,100,0.5)) +
    scale_x_continuous(name = "Probability of trace\nsuccess (%)",
                       breaks = seq(0,1,0.2), labels = label_pc) +
    scale_colour_brewer(palette = "Dark2", labels=label_backtrace,
                        name="Trace type") +
    scale_linetype_discrete(name = NULL, labels = label_limits) +
    scale_shape_discrete(name = NULL, labels = label_limits) +
    guides(colour=guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base + theme(
      # legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      # legend.spacing.y = unit(legend_vspace, "mm"),
    )
}

#------------------------------------------------------------------------------
# Contour plots
#------------------------------------------------------------------------------

contour_plot_reff <- function(data, inc_manual, bt_dist,
                              palette_tib = contour_palette_reff_bidir,
                              contour_skip = 1, step_size = 0.1){
  # Prepare data
  data_filtered <- data %>% filter((p_traced_manual != 0) == inc_manual,
                                   backtrace_distance == bt_dist)
  data_indexed <- data_filtered %>%
    group_by(p_smartphone_overall) %>% arrange(p_smartphone_overall) %>% 
    mutate(index_smartphone = group_indices()) %>%
    group_by(p_data_sharing_auto) %>% arrange(p_data_sharing_auto) %>% 
    mutate(index_sharing = group_indices())
  # Smooth data
  smoothed_reff <- sapply(1:max(data_indexed$index_smartphone), function(p)
    sapply(1:max(data_indexed$index_sharing), function(s)
      data_indexed %>% filter(abs(index_smartphone-p) <= 1,
                              abs(index_sharing-s) <= 1) %>% 
        pull(effective_r0_mean) %>% mean)) %>%
    melt() %>% as_tibble %>%
    rename(index_smartphone = Var2, index_sharing = Var1, r_eff_smoothed = value)
  data_smoothed <- inner_join(data_indexed, smoothed_reff, 
                              by=c("index_sharing", "index_smartphone"))
  # Prepare palette
  min_level <- data_smoothed %>% pull(r_eff_smoothed) %>% 
    (function(x) floor(x/step_size)*step_size) %>% min
  palette <- palette_tib %>% filter(level >= min_level) %>% pull(colour)
  # Make plot
  g <- ggplot(data_smoothed, aes(x=p_smartphone_overall, y=p_data_sharing_auto,
                                 z = r_eff_smoothed)) +
    geom_contour_filled(colour="black", breaks=seq(0,10,step_size)) +
    geom_text_contour(colour="black", stroke=0.2, breaks=seq(0,10,step_size),
                      skip=contour_skip) +
    scale_x_continuous(name = "% of cases with\nchirping smartphones",
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_y_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                       labels = function(x) round(x*100)) +
    scale_fill_manual(values = palette) +
    coord_fixed() +
    theme_base + theme(legend.position = "none")
  return(g)
}

contour_plot_ctrl <- function(data, inc_manual, bt_dist,
                              palette_tib = contour_palette_ctrl_bidir,
                              contour_skip = 0, step_size = 0.1){
  # Prepare data
  data_filtered <- data %>% filter((p_traced_manual != 0) == inc_manual,
                                   backtrace_distance == bt_dist)
  data_indexed <- data_filtered %>%
    group_by(p_smartphone_overall) %>% arrange(p_smartphone_overall) %>% 
    mutate(index_smartphone = group_indices()) %>%
    group_by(p_data_sharing_auto) %>% arrange(p_data_sharing_auto) %>% 
    mutate(index_sharing = group_indices())
  # Smooth data
  smoothed_ctrl <- sapply(1:max(data_indexed$index_smartphone), function(p)
    sapply(1:max(data_indexed$index_sharing), function(s)
      data_indexed %>% filter(abs(index_smartphone-p) <= 1,
                              abs(index_sharing-s) <= 1) %>% 
        pull(p_controlled) %>% mean)) %>%
    melt() %>% as_tibble %>%
    rename(index_smartphone = Var2, index_sharing = Var1, p_ctrl_smoothed = value)
  data_smoothed <- inner_join(data_indexed, smoothed_ctrl, 
                              by=c("index_sharing", "index_smartphone")) %>%
    mutate(pc_ctrl_smoothed = p_ctrl_smoothed*100)
  # Prepare palette
  min_level <- data_smoothed %>% pull(p_ctrl_smoothed) %>% 
    (function(x) floor(x/step_size)*step_size) %>% min
  palette <- palette_tib %>% filter(level >= min_level) %>% pull(colour)
  # Make plot
  step_size <- step_size * 100
  g <- ggplot(data_smoothed, aes(x=p_smartphone_overall, y=p_data_sharing_auto,
                                 z = pc_ctrl_smoothed)) +
    geom_contour_filled(colour="black", breaks=seq(0,100,step_size)) +
    geom_text_contour(colour="black", stroke=0.2, breaks=seq(0,100,step_size),
                      skip=contour_skip) +
    scale_x_continuous(name = "% of cases with\nchirping smartphones",
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_y_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                       labels = function(x) round(x*100)) +
    scale_fill_manual(values = palette) +
    coord_fixed() +
    theme_base + theme(legend.position = "none")
  return(g)
}

#==============================================================================
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Contour plots
#------------------------------------------------------------------------------

contour_bidir_reff_2day <- contour_data %>% filter(contact_limit_manual == 2) %>% 
  contour_plot_reff(TRUE, Inf, contour_skip = 1, step_size = 0.05)
contour_bidir_ctrl_2day <- contour_data %>% filter(contact_limit_manual == 2) %>% 
  contour_plot_ctrl(TRUE, Inf, contour_skip = 0)
contour_bidir_grid_2day <- plot_grid(
  contour_bidir_ctrl_2day + ggtitle("% of outbreaks\ncontrolled"),
  contour_bidir_reff_2day + ggtitle("Mean effective\nreproduction number"),
  nrow = 1, align = "hv", axis = "l")

contour_bidir_reff_7day <- contour_data %>% filter(contact_limit_manual == 7) %>% 
  contour_plot_reff(TRUE, Inf, contour_skip = 0, step_size = 0.05)
contour_bidir_ctrl_7day <- contour_data %>% filter(contact_limit_manual == 7) %>% 
  contour_plot_ctrl(TRUE, Inf, contour_skip = 0)
contour_bidir_grid_7day <- plot_grid(
  contour_bidir_ctrl_7day + ggtitle("% of outbreaks\ncontrolled"),
  contour_bidir_reff_7day + ggtitle("Mean effective\nreproduction number"),
                  nrow = 1, align = "hv", axis = "l")

contour_bidir_grid_7day_annot <- plot_grid(
  contour_bidir_ctrl_7day + ggtitle("% of outbreaks\ncontrolled") +
    geom_vline(xintercept=0.8, colour="red", linetype = "dashed", size=0.8) +
    geom_hline(yintercept=0.9, colour="red", linetype = "dashed", size = 0.8),
  contour_bidir_reff_7day + ggtitle("Mean effective\nreproduction number") +
    geom_vline(xintercept=0.8, colour="red", linetype = "dashed", size = 0.8) +
    geom_hline(yintercept=0.9, colour="red", linetype = "dashed", size = 0.8),
  nrow = 1, align = "hv", axis = "l")

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio, device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/presentation_hybrid"

# 48% pre-symptomatic
plot_scale <- 13
plot_ratio <- 21/plot_scale 
save_fig(path_prefix, "_contour_2day.png", contour_bidir_grid_2day, plot_scale,
         plot_ratio, device="png")
save_fig(path_prefix, "_contour_7day.png", contour_bidir_grid_7day, plot_scale,
         plot_ratio, device="png")
save_fig(path_prefix, "_contour_7day_annot.png", contour_bidir_grid_7day_annot,
         plot_scale, plot_ratio, device="png")