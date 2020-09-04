library(tidyverse)
library(cowplot)
library(metR)
library(reshape2)
library(scales)
library(viridis)

#==============================================================================
# Read in data
#==============================================================================

# Manual tracing (with and without automated)
data_path <- "figures/data/fig2_univ_env_1k_scenario.tsv.gz"
data <- suppressMessages(read_tsv(data_path)) %>%
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
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Base functions
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

combo_plot <- function(data, err_width = 0.05, pt_size = 2, coord_ratio = 1){
  #' Make control and R_eff plots for presentation
  ctrl <- ctrl_plot(data = data, err_width = err_width,
                    pt_size = pt_size, coord_ratio = coord_ratio)
  reff <- reff_plot(data = data, pt_size = pt_size,
                    coord_ratio = coord_ratio)
  grid <- plot_grid(ctrl + theme(legend.position = "none"),
                    reff + theme(legend.position = "none"), 
                    nrow = 1, align = "hv", axis = "l")
  legend_a <- get_legend(ctrl + theme(legend.justification = "center"))
  out <- plot_grid(grid, legend_a, ncol = 1, rel_heights = c(1,0.1))
  return(out)
}

#------------------------------------------------------------------------------
# Subfigure plots
#------------------------------------------------------------------------------

# 7-day bidirectional
digital_7day_bidir <- data %>% combo_plot

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio, device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/presentation_digital"

# 48% pre-symptomatic
row_height <- 14
save_fig(path_prefix, "_7day_bidir.png", digital_7day_bidir,
         row_height, 1.5, device="png")