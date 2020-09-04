library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Universal coverage
limits_path <- "figure_data/si/si_limits_baseline_1k_scenario.tsv.gz"
limits_data <- suppressMessages(read_tsv(limits_path))


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
         "Forward and reverse\ntracing")
}

label_ttype <- function(x){
  ifelse(x==0, "Manual tracing only",
         "Manual+automated\ntracing")
}

label_pc <- function(x) round(as.numeric(x)*100)

label_trace <- function(x){
  paste0(round(as.numeric(x)*100), "% contacts traced")
}

#==============================================================================
# Make plots
#==============================================================================

control_plot <- function(data, x_var_name, x_axis_name, x_breaks, x_lab,
                         lt_var_name = NULL, 
                         lt_var_levels = NULL, lt_labeller = NULL,
                         err_width = 0.05, pt_size = 2, coord_ratio = 1){
  # Prepare data
  data[["x"]] <- data[[x_var_name]]
  if (!is.null(lt_var_name)) data[["lt"]] <- data[[lt_var_name]]
  # Make base plot
  if (is.null(lt_var_name)){
    g <- ggplot(data, aes(x=x, y=p_controlled, 
                          ymin=p_controlled_lower, ymax=p_controlled_upper,
                          colour = factor(backtrace_distance)))
  } else {
    g <- ggplot(data, aes(x=x, y=p_controlled, 
                          ymin=p_controlled_lower, ymax=p_controlled_upper,
                          colour = factor(backtrace_distance),
                          linetype = factor(lt, levels=lt_var_levels),
                          shape = factor(lt, levels=lt_var_levels)))
  }
  # Add geoms
  g <- g + geom_errorbar(width=err_width) + geom_line() + 
    geom_point(size=pt_size) +
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_colour_brewer(palette = "Dark2", labels=label_backtrace,
                        name=NULL) +
    scale_x_continuous(name = x_axis_name, breaks = x_breaks, labels = x_lab
    ) +
    coord_fixed(ratio = coord_ratio) + 
    theme_base + theme(
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.spacing.y = unit(1, "mm"),
    )
  if (!is.null(lt_var_name)){
    g <- g + scale_linetype_discrete(name = NULL, labels = lt_labeller) +
      scale_shape_discrete(name = NULL, labels = lt_labeller)
  }
  return(g)
}

r_eff_plot <- function(data, x_var_name, x_axis_name, x_breaks, x_lab,
                       col_var_name, col_var_levels, col_labeller,
                       lt_var_name = NULL, lt_var_levels = NULL, lt_labeller = NULL,
                       err_width = 0.05, pt_size = 2, coord_ratio = 0.4,
                       legend_pos = c(0.01, 0.21), legend_vspace = -1,
                       y_title_short = TRUE, y_break_max = 10,
                       y_limits = c(0,NA), col_palette = "Dark2",
                       r0_base = 2.5, r0_lab_x = 0.98, r0_lab_y = 2.37,
                       r0_lab_size = 4.6, r0_col = "red"){
  # Prepare data
  data[["x"]] <- data[[x_var_name]]
  data[["col"]] <- data[[col_var_name]]
  if (!is.null(lt_var_name)) data[["lt"]] <- data[[lt_var_name]]
  y_title <- "Mean effective reproduction number"
  # Make base plot
  if (is.null(lt_var_name)){
    g <- ggplot(data, aes(x=x, y=effective_r0_mean, 
                          colour = factor(col, levels = col_var_levels)))
  } else {
    g <- ggplot(data, aes(x=x, y=effective_r0_mean, 
                          colour = factor(col, levels = col_var_levels),
                          linetype = factor(lt, levels=lt_var_levels),
                          shape = factor(lt, levels=lt_var_levels)))
  }
  # Add geoms etc
  g <- g + geom_line() + geom_point(size=pt_size) +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    scale_y_continuous(name = y_title, breaks = seq(0,y_break_max,0.5),
                       limits = y_limits) +
    
    scale_colour_brewer(palette = col_palette, labels=col_labeller,
                        name=NULL) +
    scale_x_continuous(name = x_axis_name, breaks = x_breaks, labels = x_lab) +
    coord_fixed(ratio = coord_ratio) + 
    theme_base + theme(
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.spacing.y = unit(1, "mm"),
    )
  if (!is.null(lt_var_name)){
    g <- g + scale_linetype_discrete(name = NULL, labels = lt_labeller) +
      scale_shape_discrete(name = NULL, labels = lt_labeller)
  }
  if (!is.null(r0_base)){
    g <- g + geom_hline(yintercept = r0_base, linetype = "dotted", size=1,
                        colour = "red") +
      annotate("text", x=r0_lab_x, y=r0_lab_y, hjust=0.5, vjust = 0.5,
               colour = r0_col, size=r0_lab_size, label=expression(R[0]))
  }
  return(g)
}

#------------------------------------------------------------------------------
# Limit plots
#------------------------------------------------------------------------------

# Manual only
limit_plot_ctrl_manual <- limits_data %>% filter(p_traced_manual >= 0.5,
                                             p_traced_auto == 0) %>%
  control_plot("contact_limit_manual",
                           "Trace limit (days)", seq(0,10,1), function(x) x, 
                           coord_ratio = 4) + 
  facet_wrap(~p_traced_manual, labeller = labeller(p_traced_manual = label_trace)) +
  theme(legend.position = "bottom")

limit_plot_reff_manual <- limits_data %>% filter(p_traced_manual >= 0.5,
                                           p_traced_auto == 0) %>%
  r_eff_plot("contact_limit_manual", "Trace limit (days)", 
             seq(0,10,1), function(x) x, "backtrace_distance", c(0, Inf), 
             label_backtrace,
             coord_ratio = 3.2*4/7, r0_lab_x = 6.8, r0_lab_y = 2.3) + 
  facet_wrap(~p_traced_manual, labeller = labeller(p_traced_manual = label_trace)) +
  theme(legend.position = "bottom")

limit_plot_manual <- plot_grid(limit_plot_ctrl_manual, limit_plot_reff_manual,
                               labels = "auto", nrow = 2, ncol = 1, align = "hv",
                               axis = "l", label_size = fontsize_base * fontscale_label,
                               label_fontfamily = titlefont, label_colour = "black")

# Combined
limit_plot_ctrl_combo <- limits_data %>% filter(p_traced_manual >= 0.5,
                                             p_traced_auto < 0) %>%
  control_plot("contact_limit_manual",
               "Trace limit (days)", seq(0,10,1), function(x) x, 
               coord_ratio = 4) + 
  facet_wrap(~p_traced_manual, labeller = labeller(p_traced_manual = label_trace)) +
  theme(legend.position = "bottom")

limit_plot_reff_combo <- limits_data %>% filter(p_traced_manual >= 0.5,
                                           p_traced_auto < 0) %>%
  r_eff_plot("contact_limit_manual", "Trace limit (days)", 
             seq(0,10,1), function(x) x, "backtrace_distance", c(0, Inf), 
             label_backtrace,
             coord_ratio = 3.2*4/7, r0_lab_x = 6.8, r0_lab_y = 2.3) + 
  facet_wrap(~p_traced_manual, labeller = labeller(p_traced_manual = label_trace)) +
  theme(legend.position = "bottom")

limit_plot_combo <- plot_grid(limit_plot_ctrl_combo, limit_plot_reff_combo,
                               labels = "auto", nrow = 2, ncol = 1, align = "hv",
                               axis = "l", label_size = fontsize_base * fontscale_label,
                               label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = "png", width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures/si/limits_"

# 48% pre-symptomatic
row_height <- 8
save_fig(path_prefix, "manual.png", limit_plot_manual, row_height*5, 1.4/2)
save_fig(path_prefix, "combo.png", limit_plot_combo, row_height*5, 1.4/2)


