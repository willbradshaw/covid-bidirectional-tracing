library(tidyverse)
library(cowplot)
library(metR)
library(reshape2)
library(scales)
library(viridis)

#==============================================================================
# Read in data
#==============================================================================

asym_path <- "figures_local/data/si_asym_contour_500_scenario.tsv.gz"
asym_data <- suppressMessages(read_tsv(asym_path)) %>%
  filter(backtrace_distance == 0 | contact_limit_manual == 7)

#==============================================================================
# Process data
#==============================================================================

#------------------------------------------------------------------------------
# Processing functions
#------------------------------------------------------------------------------

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

label_data <- function(data, p_traced = 0.9){
  data %>% filter(p_traced_auto %in% c(-1, 0, p_traced),
                  p_traced_manual %in% c(-1, 0, p_traced)) %>%
    mutate(bt_type = ifelse(backtrace_distance == 0, "Forward-only",
                            "Bidirectional"),
           trace_window = paste0(contact_limit_manual, "-day manual window"),
           row_label = paste(bt_type, trace_window, sep=",\n"),
           trace_type = ifelse(p_traced_auto == 0,
             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital"))) %>%
    mutate(trace_type = factor(trace_type, levels = ttype_levels))
}

index_data <- function(data, key_x = "p_asymptomatic", key_y = "rel_r0_asymptomatic"){
  # Compute index levels
  levels_x <- sort(unique(data[[key_x]]))
  levels_y <- sort(unique(data[[key_y]]))
  # Assign indices
  data_indexed <- data %>%
    mutate(index_x = match(!!sym(key_x), levels_x),
           index_y = match(!!sym(key_y), levels_y))
  return(data_indexed)
}

smooth_data_single <- function(data){
  indices_x <- sort(unique(data$index_x))
  indices_y <- sort(unique(data$index_y))
  index_tab <- expand_grid(index_x = indices_x, index_y = indices_y)
  smoothed_control <- sapply(1:nrow(index_tab), function(n)
    data %>% filter(abs(index_x - index_tab$index_x[n]) <= 1,
                    abs(index_y - index_tab$index_y[n]) <= 1) %>%
      pull(p_controlled) %>% mean)
  smoothed_r_eff <- sapply(1:nrow(index_tab), function(n)
    data %>% filter(abs(index_x - index_tab$index_x[n]) <= 1,
                    abs(index_y - index_tab$index_y[n]) <= 1) %>%
      pull(effective_r0_mean) %>% mean)
  index_tab_smoothed <- mutate(index_tab,
                               p_controlled_smoothed = smoothed_control,
                               r_eff_smoothed = smoothed_r_eff)
  data_smoothed <- inner_join(data, index_tab_smoothed,
                              by=c("index_x", "index_y"))
  return(data_smoothed)
}

smooth_data <- function(data, group_vars = c("trace_type", "row_label")){
  data_split <- data %>% group_by_at(group_vars) %>% group_split()
  data_smoothed_split <- lapply(data_split, smooth_data_single)
  return(bind_rows(data_smoothed_split) %>% ungroup %>%
           mutate(pc_controlled_smoothed = p_controlled_smoothed*100))
}

#------------------------------------------------------------------------------
# Apply to data
#------------------------------------------------------------------------------

# Label data
asym_data_labelled <- label_data(asym_data)

# Index data
asym_data_indexed <- index_data(asym_data_labelled)

# Smooth data
asym_data_smoothed <- smooth_data(asym_data_indexed)
  

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
# Configure palettes
#==============================================================================

col_bidir <- "#d95f02"
col_fwd   <- "#1b9e77"
step_dense <- 0.05
step_sparse <- 0.1

make_palette_ctrl <- function(col_high, step_size){
  # Calibrate range of palette
  palette_max <- 1
  palette_min <- 0
  palette_val <- seq(palette_min, palette_max, step_size)
  # Calculate hue values
  palette_fn  <- scales::seq_gradient_pal(low="white", high=col_high, "Lab")
  palette_hue <- palette_fn((seq(0,1,length.out=length(palette_val)))) %>% rev
  return(tibble(level = palette_val, colour = palette_hue))
}

bidir_palette_dense <- make_palette_ctrl(col_bidir, step_dense)
bidir_palette_sparse <- make_palette_ctrl(col_bidir, step_sparse)
fwd_palette_dense <- make_palette_ctrl(col_fwd, step_dense)
fwd_palette_sparse <- make_palette_ctrl(col_fwd, step_sparse)


#==============================================================================
# Make plots 
#==============================================================================

contour_plot_base <- function(data, z_key, contour_breaks,
                              label_map = aes(label = ..level..),
                              coord_ratio = 1, lwd = 0.8){
  data %>%
    ggplot(aes_string(x="p_asymptomatic", y="rel_r0_asymptomatic",
                      z=z_key)) +
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    # geom_hline(yintercept=c(0.53, 0.8), linetype = "dashed", colour="red", size=lwd) +
    # geom_vline(xintercept = 2.5, linetype = "dashed", colour = "red", size = lwd) +
    geom_text_contour(colour="black", stroke = 0.2, breaks = contour_breaks,
                      check_overlap = TRUE, mapping = label_map) +
    facet_grid(row_label ~ trace_type) +
    scale_x_continuous(name = "% of asymptomatic carriers", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_y_continuous(name = expression(paste("Relative ",italic(R)[0]," of asymptomatic carriers")),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    coord_fixed(ratio = coord_ratio) +
    theme_base + theme(
      legend.position = "none",
    )
}

#------------------------------------------------------------------------------
# Absolute values
#------------------------------------------------------------------------------

contour_plot_ctrl <- function(data){
  contour_plot_base(asym_data_smoothed, "pc_controlled_smoothed", seq(-10, 110, 10))
}

contour_plot_reff <- function(data){
  contour_plot_base(asym_data_smoothed, "r_eff_smoothed", seq(0,10,0.2))
  
}

asym_contour_ctrl <- contour_plot_ctrl(asym_data_smoothed)
asym_contour_reff <- contour_plot_reff(asym_data_smoothed)


#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, 
                     plot_scale,
                     plot_ratio, 
                     device = "png"){
  ggsave(filename=paste0(path_prefix, path_suffix, ".", device), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

row_height <- 8.5
plot_scale <- row_height * 3
plot_ratio <- 4/3
plot_prefix <- "figures_local/img/presentation_asym_contour_090_7d_"

save_fig(plot_prefix, "ctrl", asym_contour_ctrl, plot_scale, plot_ratio)
save_fig(plot_prefix, "reff", asym_contour_reff, plot_scale, plot_ratio)

