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
  filter(generation_alpha == 0.064, contact_limit_manual == 7)

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

contour_palette_viridis <- tibble(level = seq(0,1,0.1), colour = viridis(11))
contour_palette_bidir   <- tibble(level = seq(0,1,0.1), 
                                  colour = scales::seq_gradient_pal(low="white", high="#d95f02", "Lab")(seq(0,1,length.out=11)))
contour_palette_fwd   <- tibble(level = seq(0,1,0.1), 
                                colour = scales::seq_gradient_pal(low="white", high="#1b9e77", "Lab")(seq(0,1,length.out=11)))

contour_plot <- function(data, inc_manual, bt_dist,
                         palette_tib = contour_palette_bidir){
  # Prepare data
  data_filtered <- data %>% filter((p_traced_manual != 0) == inc_manual,
                                   backtrace_distance == bt_dist)
  data_indexed <- data_filtered %>%
    group_by(p_smartphone_overall) %>% arrange(p_smartphone_overall) %>% 
    mutate(index_smartphone = group_indices()) %>%
    group_by(p_data_sharing_auto) %>% arrange(p_data_sharing_auto) %>% 
    mutate(index_sharing = group_indices())
  # Smooth data
  smoothed_control <- sapply(1:max(data_indexed$index_smartphone), function(p)
    sapply(1:max(data_indexed$index_sharing), function(s)
      data_indexed %>% filter(abs(index_smartphone-p) <= 1,
                              abs(index_sharing-s) <= 1) %>% 
        pull(p_controlled) %>% mean)) %>%
    melt() %>% as_tibble %>%
    rename(index_smartphone = Var2, index_sharing = Var1, p_controlled_smoothed = value)
  data_smoothed <- inner_join(data_indexed, smoothed_control, 
                              by=c("index_sharing", "index_smartphone"))
  # Prepare palette
  min_level <- data_smoothed %>% pull(p_controlled_smoothed) %>% 
    (function(x) floor(x*10)/10) %>% min
  palette <- palette_tib %>% filter(level >= min_level) %>% pull(colour)
  # Make plot
  g <- ggplot(data_smoothed, aes(x=p_smartphone_overall, y=p_data_sharing_auto,
                                 z = round(p_controlled_smoothed*100))) +
    geom_contour_filled(colour="black", breaks=seq(0,100,10)) +
    geom_text_contour(colour="black", stroke=0.2, breaks=seq(0,100,10),
                      check_overlap = TRUE) +
    scale_x_continuous(name = "% of cases with\nchirping smartphones",
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_y_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                       labels = function(x) round(x*100)) +
    scale_fill_manual(values = palette) +
    coord_fixed() +
    theme_base + theme(legend.position = "none",
                       # panel.background = element_rect(fill = NA),
                       # panel.ontop = TRUE
    )
  return(g)
}

#------------------------------------------------------------------------------
# Make contour plots
#------------------------------------------------------------------------------

contour_plot_digital_fwd <- contour_plot(contour_data, FALSE, 0,
                                         contour_palette_fwd) +
  ggtitle("Forward-only digital\ntracing")
contour_plot_combo_fwd     <- contour_plot(contour_data, TRUE, 0,
                                           contour_palette_fwd) +
  ggtitle("Forward-only hybrid\ntracing (7-day manual limit)")

contour_plot_grid <- plot_grid(contour_plot_digital_fwd,
                               contour_plot_combo_fwd,
                               labels = "auto", nrow = 1, ncol = 2, align = "hv",
                               axis = "l", label_size = fontsize_base * fontscale_label,
                               label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio, device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/si_contour_fwd"

row_height <- 13
save_fig(path_prefix, ".png", contour_plot_grid, row_height, 2 * 0.8, device="png")