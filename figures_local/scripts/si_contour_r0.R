library(tidyverse)
library(cowplot)
library(metR)
library(reshape2)
library(scales)
library(viridis)

#==============================================================================
# Read in data
#==============================================================================

#------------------------------------------------------------------------------
# Hybrid tracing
#------------------------------------------------------------------------------

# Median scenario
median_path_hybrid <- "figures_local/data/si_contour_baseline_500_scenario.tsv.gz"
median_data_hybrid <- suppressMessages(read_tsv(median_path_hybrid)) %>%
  mutate(scenario_type = "Median scenario")

# Optimistic scenario
optim_path_hybrid <- "figures_local/data/si_contour_optimistic_500_scenario.tsv.gz"
optim_data_hybrid <- suppressMessages(read_tsv(optim_path_hybrid)) %>%
  mutate(scenario_type = "Optimistic scenario")

# Pessimistic scenario
pessim_path_hybrid <- "figures_local/data/si_contour_pessimistic_500_scenario.tsv.gz"
pessim_data_hybrid <- suppressMessages(read_tsv(pessim_path_hybrid)) %>%
  mutate(scenario_type = "Pessimistic scenario")

# Combine
data_hybrid <- bind_rows(median_data_hybrid, optim_data_hybrid,
                         pessim_data_hybrid) %>% mutate(trace_type = "hybrid")

#------------------------------------------------------------------------------
# Manual tracing
#------------------------------------------------------------------------------

# Median scenario
median_path_manual <- "figures_local/data/si_contour_baseline_manual_500_scenario.tsv.gz"
median_data_manual <- suppressMessages(read_tsv(median_path_manual)) %>%
  mutate(scenario_type = "Median scenario")

# Optimistic scenario
optim_path_manual <- "figures_local/data/si_contour_optimistic_manual_500_scenario.tsv.gz"
optim_data_manual <- suppressMessages(read_tsv(optim_path_manual)) %>%
  mutate(scenario_type = "Optimistic scenario")

# Pessimistic scenario
pessim_path_manual <- "figures_local/data/si_contour_pessimistic_manual_500_scenario.tsv.gz"
pessim_data_manual <- suppressMessages(read_tsv(pessim_path_manual)) %>%
  mutate(scenario_type = "Pessimistic scenario")

# Combine
data_manual <- bind_rows(median_data_manual, optim_data_manual,
                         pessim_data_manual) %>% mutate(trace_type = "manual")
data <- bind_rows(data_hybrid, data_manual)

#==============================================================================
# Process data
#==============================================================================

#------------------------------------------------------------------------------
# Processing functions
#------------------------------------------------------------------------------

label_data <- function(data, p_traced){
  data %>% filter(p_traced_manual == p_traced) %>%
    mutate(bt_type = ifelse(backtrace_distance == 0, "Forward-only",
                            "Bidirectional"),
           trace_window = paste0(contact_limit_manual, "-day manual window"),
           row_label = paste(bt_type, trace_window, sep=",\n"))
}

index_data <- function(data, key_x = "r0_base", key_y = "p_smartphone_overall"){
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

smooth_data <- function(data, group_vars = c("scenario_type", "trace_type",
                                             "row_label")){
  data_split <- data %>% group_by_at(group_vars) %>% group_split()
  data_smoothed_split <- lapply(data_split, smooth_data_single)
  return(bind_rows(data_smoothed_split) %>% ungroup)
}

diff_data_backtrace <- function(data, group_vars = c("scenario_type", "contact_limit_manual",
                                           "r0_base", "p_smartphone_overall",
                                           "trace_type")){
  spread_control <- data %>% group_by_at(c(group_vars, "index_x", "index_y")) %>%
    select(pc_controlled_smoothed, backtrace_distance) %>%
    spread(backtrace_distance, pc_controlled_smoothed, sep="_", drop=TRUE) %>%
    mutate(pc_controlled_smoothed_diff = backtrace_distance_Inf - backtrace_distance_0) %>%
    select(-matches("backtrace_distance"))
  spread_r_eff <- data %>% group_by_at(c(group_vars, "index_x", "index_y")) %>%
    select(r_eff_smoothed, backtrace_distance) %>%
    spread(backtrace_distance, r_eff_smoothed, sep="_", drop=TRUE) %>%
    mutate(r_eff_smoothed_diff = backtrace_distance_Inf - backtrace_distance_0) %>%
    select(-matches("backtrace_distance"))
  data_diff <- inner_join(spread_control, spread_r_eff,
                          by=c(group_vars, "index_x", "index_y")) %>%
    mutate(row_label = paste0(contact_limit_manual, "-day manual window"))
  return(data_diff)
}

diff_data_hybrid <- function(data, group_vars = c("scenario_type", "contact_limit_manual",
                                                  "r0_base", "p_smartphone_overall",
                                                  "row_label")){
  spread_control <- data %>% group_by_at(c(group_vars, "index_x", "index_y")) %>%
    select(pc_controlled_smoothed, trace_type) %>%
    spread(trace_type, pc_controlled_smoothed, sep="_", drop=TRUE) %>%
    mutate(pc_controlled_smoothed_diff = trace_type_hybrid - trace_type_manual) %>%
    select(-matches("trace_type"))
  spread_r_eff <- data %>% group_by_at(c(group_vars, "index_x", "index_y")) %>%
    select(r_eff_smoothed, trace_type) %>%
    spread(trace_type, r_eff_smoothed, sep="_", drop=TRUE) %>%
    mutate(r_eff_smoothed_diff = trace_type_hybrid - trace_type_manual) %>%
    select(-matches("trace_type"))
  data_diff <- inner_join(spread_control, spread_r_eff,
                          by=c(group_vars, "index_x", "index_y"))
  return(data_diff)
}


#------------------------------------------------------------------------------
# Hybrid tracing
#------------------------------------------------------------------------------

# Label data
data_090_labelled <- label_data(data, 0.9)
data_080_labelled <- label_data(data, 0.8)

# Index data
data_090_indexed <- index_data(data_090_labelled)
data_080_indexed <- index_data(data_080_labelled)

# Smooth data
data_090_smoothed <- smooth_data(data_090_indexed) %>%
  mutate(pc_controlled_smoothed = p_controlled_smoothed * 100)
data_080_smoothed <- smooth_data(data_080_indexed) %>%
  mutate(pc_controlled_smoothed = p_controlled_smoothed * 100)

# Diff data (bidirectional vs forward-only)
data_090_diff <- diff_data_backtrace(data_090_smoothed)
data_080_diff <- diff_data_backtrace(data_080_smoothed)

# Diff data (manual vs hybrid)
data_090_diff_hybrid <- diff_data_hybrid(data_090_smoothed)
data_080_diff_hybrid <- diff_data_hybrid(data_080_smoothed)


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
                              coord_ratio = 3, lwd = 0.8){
  data %>%
    ggplot(aes_string(x="r0_base", y="p_smartphone_overall",
                      z=z_key)) +
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    geom_hline(yintercept=c(0.53, 0.8), linetype = "dashed", colour="red", size=lwd) +
    geom_vline(xintercept = 2.5, linetype = "dashed", colour = "red", size = lwd) +
    geom_text_contour(colour="black", stroke = 0.2, breaks = contour_breaks,
                      check_overlap = TRUE, mapping = label_map) +
    facet_grid(row_label ~ scenario_type) +
    scale_x_continuous(name = expression(italic(R)[0]),
                       breaks = seq(0,10,0.5)) +
    scale_y_continuous(name = "% of cases with chirping smartphones",
                       breaks = seq(0,1,0.2),
                       labels = function(x) round(x*100)) +
    coord_fixed(ratio = coord_ratio) +
    theme_base + theme(
      legend.position = "none",
    )
}

#------------------------------------------------------------------------------
# Absolute values
#------------------------------------------------------------------------------

contour_plot_facet_viridis <- function(data, ttype = "hybrid",
                                    z_key = "pc_controlled_smoothed",
                                    contour_breaks = seq(0, 110, 10)){
  contour_plot_base(data = filter(data, trace_type == ttype),
                    z_key = z_key, contour_breaks = contour_breaks)
}

contour_plot_ctrl <- function(data) contour_plot_facet_viridis(data)
contour_plot_reff <- function(data) contour_plot_facet_viridis(data, z_key = "r_eff_smoothed",
                                                               contour_breaks = seq(0,10,0.2))

# 90% contacts traced
contour_r0_090_ctrl <- contour_plot_ctrl(data_090_smoothed)
contour_r0_090_reff <- contour_plot_reff(data_090_smoothed)

# 80% contacts traced
contour_r0_080_ctrl <- contour_plot_ctrl(data_080_smoothed)
contour_r0_080_reff <- contour_plot_reff(data_080_smoothed)

#------------------------------------------------------------------------------
# Fwd vs bidirectional
#------------------------------------------------------------------------------

contour_plot_facet_diff <- function(data, ttype = "hybrid",
                                    z_key = "pc_controlled_smoothed_diff",
                                    contour_breaks = seq(-10, 100, 10)){
  contour_plot_base(data = filter(data, trace_type == ttype),
                    z_key = z_key, contour_breaks = contour_breaks,
                    label_map = aes(label=paste0(ifelse(..level.. >= 0, "+", ""), ..level..)),
                    coord_ratio = 3)
}

diff_plot_ctrl <- function(data) contour_plot_facet_diff(data)
diff_plot_reff <- function(data) contour_plot_facet_diff(data, z_key = "r_eff_smoothed_diff",
                                                               contour_breaks = seq(0.1,-10,-0.1))

# 90% contacts traced
diff_r0_090_ctrl <- diff_plot_ctrl(data_090_diff)
diff_r0_090_reff <- diff_plot_reff(data_090_diff)

# 80% contacts traced
diff_r0_080_ctrl <- diff_plot_ctrl(data_080_diff)
diff_r0_080_reff <- diff_plot_reff(data_080_diff)

#------------------------------------------------------------------------------
# Manual vs hybrid
#------------------------------------------------------------------------------

contour_plot_hybrid_diff <- function(data, z_key = "pc_controlled_smoothed_diff",
                                     contour_breaks = seq(-10, 100, 10)){
  contour_plot_base(data = data, z_key = z_key, contour_breaks = contour_breaks,
                    label_map = aes(label=paste0(ifelse(..level.. >= 0, "+", ""), ..level..)),
                    coord_ratio = 3)
}

diff_hybrid_plot_ctrl <- function(data) contour_plot_hybrid_diff(data)
diff_hybrid_plot_reff <- function(data) contour_plot_hybrid_diff(data, z_key = "r_eff_smoothed_diff",
                                                         contour_breaks = seq(0.1,-10,-0.1))

# 90% contacts traced
diff_r0_090_hybrid_ctrl <- diff_hybrid_plot_ctrl(data_090_diff_hybrid)
diff_r0_090_hybrid_reff <- diff_hybrid_plot_reff(data_090_diff_hybrid)

# 80% contacts traced
diff_r0_080_hybrid_ctrl <- diff_hybrid_plot_ctrl(data_080_diff_hybrid)
diff_r0_080_hybrid_reff <- diff_hybrid_plot_reff(data_080_diff_hybrid)

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio, device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/si_contour_r0_"

row_height <- 10

save_fig(path_prefix, "090_ctrl.png", contour_r0_090_ctrl,
         row_height * 4, 0.8, device="png")
save_fig(path_prefix, "090_reff.png", contour_r0_090_reff,
         row_height * 4, 0.8, device="png")
save_fig(path_prefix, "080_ctrl.png", contour_r0_080_ctrl,
         row_height * 4, 0.8, device="png")
save_fig(path_prefix, "080_reff.png", contour_r0_080_reff,
         row_height * 4, 0.8, device="png")

path_prefix <- "figures_local/img/si_contour_r0_diff_"

save_fig(path_prefix, "090_ctrl.png", diff_r0_090_ctrl,
         row_height * 2, 1.4, device="png")
save_fig(path_prefix, "090_reff.png", diff_r0_090_reff,
         row_height * 2, 1.4, device="png")
save_fig(path_prefix, "080_ctrl.png", diff_r0_080_ctrl,
         row_height * 2, 1.4, device="png")
save_fig(path_prefix, "080_reff.png", diff_r0_080_reff,
         row_height * 2, 1.4, device="png")

path_prefix <- "figures_local/img/si_contour_r0_diff_hybrid_"

save_fig(path_prefix, "090_ctrl.png", diff_r0_090_hybrid_ctrl,
         row_height * 4, 0.8, device="png")
save_fig(path_prefix, "090_reff.png", diff_r0_090_hybrid_reff,
         row_height * 4, 0.8, device="png")
save_fig(path_prefix, "080_ctrl.png", diff_r0_080_hybrid_ctrl,
         row_height * 4, 0.8, device="png")
save_fig(path_prefix, "080_reff.png", diff_r0_080_hybrid_reff,
         row_height * 4, 0.8, device="png")
