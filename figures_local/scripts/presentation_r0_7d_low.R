library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

# Median scenario
median_path <- "figures_local/data/fig3_baseline_explore_1k_scenario.tsv.gz"
median_data <- suppressMessages(read_tsv(median_path)) %>%
  filter(p_traced_auto == 0)

median_digital_path <- "figures_local/data/fig4_baseline_morer0_1k_scenario.tsv.gz"
median_data_digital <- suppressMessages(read_tsv(median_digital_path)) %>%
  filter(p_traced_manual == 0 | p_traced_manual == p_traced_auto)

median_data <- bind_rows(median_data, median_data_digital) %>%
mutate(trace_type = ifelse(p_traced_auto == 0,
    ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
    ifelse(p_traced_manual == 0, "Digital only", "Manual + digital"))) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels=c(Inf,0)))

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
  ifelse(x==0, "Forward only",
         "Bidirectional")
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
  ggplot(data, aes(x=r0_base, y=p_controlled,
                   ymin=p_controlled_lower, ymax=p_controlled_upper,
                   colour=trace_type,
                   linetype = backtrace_distance,
                   shape = backtrace_distance)) +
    geom_errorbar(width=err_width) + geom_line() + geom_point(size=pt_size) +
    coord_fixed(ratio = coord_ratio * (max(data$r0_base)-min(data$r0_base))) + 
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_x_continuous(name = expression(R[0]), breaks = seq(0,100,0.5)) +
    scale_colour_brewer(palette = "Set1", name="Trace type") +
    scale_linetype_discrete(name = NULL, labels = label_backtrace) +
    scale_shape_discrete(name = NULL, labels = label_backtrace) +
    guides(colour=guide_legend(ncol = 2, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base + theme(
      # legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      # legend.spacing.y = unit(legend_vspace, "mm"),
    )
}

reff_plot <- function(data, pt_size = 2, coord_ratio = 1){
  ggplot(data, aes(x=r0_base, y=effective_r0_mean,
                   colour=trace_type,
                   linetype = backtrace_distance,
                   shape = backtrace_distance)) +
    geom_line() + geom_point(size=pt_size) +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    coord_fixed(ratio = coord_ratio * (max(data$r0_base)-min(data$r0_base))/max(data$r0_base)) + 
    scale_y_continuous(name = expression(paste("Mean ",italic(R)[eff])),
                       breaks = seq(0,100,0.5), limits = c(0, max(data$r0_base))) +
    scale_x_continuous(name = expression(R[0]), breaks = seq(0,100,0.5)) +
    scale_colour_brewer(palette = "Set1", name="Trace type") +
    scale_linetype_discrete(name = NULL, labels = label_backtrace) +
    scale_shape_discrete(name = NULL, labels = label_backtrace) +
    guides(colour=guide_legend(ncol = 2, order = 1),
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
# Make different figure versions
#------------------------------------------------------------------------------

# 7-day trace limit, 90% tracing, high uptake
data_7d_090_low <- median_data %>%
  filter(contact_limit_manual == 7, p_traced_manual != 0.8,
         p_traced_auto != 0.8, trace_neg_symptomatic == TRUE,
         (p_traced_auto == 0 | p_smartphone_overall == 0.53))
r0_7d_090_low <- data_7d_090_low %>% combo_plot()


#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio, device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/presentation_r0"

# 48% pre-symptomatic
plot_scale <- 13
plot_ratio <- 21/plot_scale 
save_fig(path_prefix, "_7day_090_low.png", r0_7d_090_low,
         plot_scale, plot_ratio, device="png")
