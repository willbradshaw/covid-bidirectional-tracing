library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Universal coverage
data_path <- "figure_data/si/data_window_pessimistic_500_scenario.tsv.gz"
data <- suppressMessages(read_tsv(data_path)) %>%
  filter(data_limit_auto < Inf) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
    ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
    ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")))


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

label_limits <- function(x){
  ifelse(as.numeric(x) == Inf, "No data limit", paste0(x, "-day data limit"))
}

label_vars <- function(x){
  ifelse(x=="p_data_sharing_auto", "% of cases sharing data",
         "% of cases with\ntrace-enabled smartphones")
}

label_pc <- function(x) round(as.numeric(x)*100)

#==============================================================================
# Make plots
#==============================================================================
ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

control_plot <- data %>%
  ggplot(aes(x=r0_base, y=p_controlled,
             ymin=p_controlled_lower, ymax=p_controlled_upper,
             colour=factor(trace_type,
                           levels = ttype_levels),
             linetype=factor(data_limit_auto, levels = c(14,28)),
             shape=factor(data_limit_auto, levels = c(14,28)))) +
  geom_errorbar(width=0.05) + geom_line() + geom_point(size=2) +
  facet_grid(backtrace_distance~., labeller=labeller(backtrace_distance=label_backtrace))+
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_colour_brewer(palette = "Set1", name="Trace type") +
  scale_x_continuous(name = expression(R[0]~"(pessimistic scenario)")) +
  scale_linetype_discrete(name = "Data-retention limit (days)") +
  scale_shape_discrete(name = "Data-retention limit (days)") +
  coord_fixed(ratio = 3) + 
  guides(colour=guide_legend(ncol = 2, order = 1),
         linetype = guide_legend(ncol = 1, order = 2),
         shape = guide_legend(ncol = 1, order = 2)) +
  theme_base + theme(
    legend.position = "bottom",
  )

reff_plot <- data %>%
  ggplot(aes(x=r0_base, y=effective_r0_mean,
             colour=factor(trace_type,
                           levels = ttype_levels),
             linetype=factor(data_limit_auto, levels = c(14,28)),
             shape=factor(data_limit_auto, levels = c(14,28)))) +
  geom_line() + geom_point(size=2) +
  facet_grid(backtrace_distance~., labeller=labeller(backtrace_distance=label_backtrace))+
  scale_y_continuous(name = "Mean effective\nreproduction number", limits = c(0,4),
                     breaks = seq(0,10,0.5)) +
  scale_colour_brewer(palette = "Set1", name="Trace type") +
  scale_x_continuous(name = expression(R[0]~"(pessimistic scenario)")) +
  scale_linetype_discrete(name = "Data-retention limit (days)") +
  scale_shape_discrete(name = "Data-retention limit (days)") +
  coord_fixed(ratio = 1) + 
  guides(colour=guide_legend(ncol = 2, order = 1),
         linetype = guide_legend(ncol = 1, order = 2),
         shape = guide_legend(ncol = 1, order = 2)) +
  theme_base + theme(
    legend.position = "bottom",
  )

# Combine together
pessim_plot <- plot_grid(control_plot, reff_plot, 
                         labels = "auto", nrow = 1, ncol = 2, align = "hv",
                         axis = "l", label_size = fontsize_base * fontscale_label,
                         label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = "svg", width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures/si/data_limits_pessim"
row_height <- 12
save_fig(path_prefix, ".svg", pessim_plot, row_height*2, 1.5)
