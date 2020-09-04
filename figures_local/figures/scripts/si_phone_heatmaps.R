library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

data_path <- "figure_data/si/si_phone_heatmap_500_scenario.tsv.gz"
data <- suppressMessages(read_tsv(data_path))


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
# Make plots
#==============================================================================

gradient_low <- "white"

# Forward-only tracing
heatmap_fwd <- data %>% filter(backtrace_distance == 0) %>%
  ggplot(aes(y=p_data_sharing_auto, x=p_smartphone_overall, fill=p_controlled)) +
  geom_tile() +
  geom_text(aes(y=p_data_sharing_auto, x=p_smartphone_overall, label=round(p_controlled*100)),
            colour="black", size=3)+
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_y_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_fill_gradient(low=gradient_low, high="#7570b3", name = "% of outbreaks\ncontrolled",
                      labels = function(x) round(x*100), limits=c(0,1)) +
  coord_fixed(ratio=1) +
  theme_base + theme(
    legend.position = "bottom",
  ) +
  ggtitle("Forward tracing only")

# Bidirectional tracing
heatmap_bi <- data %>% filter(backtrace_distance == Inf) %>%
  ggplot(aes(y=p_data_sharing_auto, x=p_smartphone_overall, fill=p_controlled)) +
  geom_tile() +
  geom_text(aes(y=p_data_sharing_auto, x=p_smartphone_overall, label=round(p_controlled*100)),
            colour="black", size=3)+
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_y_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_fill_gradient(low=gradient_low, high="#7570b3", name = "% of outbreaks\ncontrolled",
                      labels = function(x) round(x*100), limits=c(0,1)) +
  coord_fixed(ratio=1) +
  theme_base + theme(
    legend.position = "bottom",
  ) +
  ggtitle("Bidirectional tracing")

# Combine
heatmaps <- plot_grid(heatmap_fwd, heatmap_bi,
                  labels = "auto", nrow=1, ncol=2,
                  align="hv", axis="l",
                  label_size = fontsize_base * fontscale_label,
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

path_prefix <- "figures/si/phone_heatmaps"
row_height <- 14

save_fig(path_prefix, ".png", heatmaps, row_height*1.1, 2 * 0.8)



