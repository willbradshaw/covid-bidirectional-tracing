#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

cat("done.\n")

#==============================================================================
# Read in data
#==============================================================================

cat("\nReading in data...")

median_path <- "figures_local/data/si_window_median_1k_scenario.tsv.gz"
median_data <- suppressMessages(read_tsv(median_path)) %>%
  mutate(scenario = "Median scenario")

optim_path <- "figures_local/data/si_window_optimistic_1k_scenario.tsv.gz"
optim_data <- suppressMessages(read_tsv(optim_path)) %>%
  mutate(scenario = "Optimistic scenario")

pessim_path <- "figures_local/data/si_window_pessimistic_1k_scenario.tsv.gz"
pessim_data <- suppressMessages(read_tsv(pessim_path)) %>%
  mutate(scenario = "Pessimistic scenario")

cat("done.\n")

#==============================================================================
# Process data
#==============================================================================

cat("\nProcessing data...")

# Combine scenarios
data <- bind_rows(median_data, optim_data, pessim_data)

# Label variables
# Contour variables: p_traced_manual, contact_limit_manual
# Other variables: p_smartphone_overall, p_traced_auto, backtrace_distance, scenario
data_labelled <- data %>%
  mutate(trace_type = ifelse(backtrace_distance == 0, "Forward-only tracing",
                             "Bidirectional tracing"),
         uptake = paste(ifelse(p_smartphone_overall == 0.8, "high", "low"),
                        "digital uptake"),
         label = ifelse(p_traced_auto == 0, trace_type,
                        paste0(trace_type, ",\n", uptake)))

# Prepare contour data
data_indexed <- index_data(data_labelled, key_x = "contact_limit_manual", 
                           key_y = "p_traced_manual")
data_smoothed <- smooth_data(data_indexed, 
                             group_vars = c("scenario", "trace_type",
                                            "uptake", "label")) %>%
  mutate(pc_controlled_smoothed = p_controlled_smoothed*100)

cat("done.\n")

#==============================================================================
# Prepare palettes
#==============================================================================

cat("\nPreparing palettes...")

# R_eff, bidirectional
bidir_palette_reff_dense <- data_smoothed %>%
  filter(backtrace_distance > 0) %>%
  make_contour_palette_reff(col_bidir, step_dense, ., Inf)
bidir_palette_reff_sparse <- data_smoothed %>%
  filter(backtrace_distance > 0) %>%
  make_contour_palette_reff(col_bidir, step_sparse, ., Inf)

# R_eff, forward-only
fwd_palette_reff_dense <- data_smoothed %>%
  filter(backtrace_distance == 0) %>%
  make_contour_palette_reff(col_fwd, step_dense, ., 0)
fwd_palette_reff_sparse <- data_smoothed %>%
  filter(backtrace_distance == 0) %>%
  make_contour_palette_reff(col_fwd, step_sparse, ., 0)

# Control
bidir_palette_ctrl_dense <- make_contour_palette_ctrl(col_bidir, step_dense)
bidir_palette_ctrl_sparse <- make_contour_palette_ctrl(col_bidir, step_sparse)
fwd_palette_ctrl_dense <- make_contour_palette_ctrl(col_fwd, step_dense)
fwd_palette_ctrl_sparse <- make_contour_palette_ctrl(col_fwd, step_sparse)

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

cat("\nGenerating plots...")

#------------------------------------------------------------------------------
# Plotting function
#------------------------------------------------------------------------------

contour_window <- function(data, contour_fn = contour_plot_reff_abs,
                           palette = bidir_palette_reff_sparse,
                           contour_breaks = seq(0,10,step_sparse),
                           upper = FALSE, lower = FALSE,
                           join_scale = 20){
  g <- data %>%
    contour_fn("contact_limit_manual", "p_traced_manual", palette = palette,
               breaks = contour_breaks, coord_ratio = 10) +
    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
    geom_vline(xintercept = c(2,7), linetype = "dashed", colour = "red") +
    scale_x_continuous(name = "Trace limit (days)", limits = c(0,10),
                       breaks = seq(0,10,2)) +
    scale_y_continuous(name = "Probability of trace\nsuccess (%)",
                       limits = c(0,1), breaks = seq(0,1,0.2),
                       labels = label_pc) +
    facet_grid(label~scenario)
  if(upper) {
    g <- g + theme(axis.text.x = element_blank(),
                           axis.title.x = element_blank(),
                           plot.margin = margin(b = -join_scale))
  }
  if (lower){
    g <- g + theme(strip.text.x = element_blank(),
                   plot.margin = margin(t = -join_scale))
  }
  return(g)
}

#------------------------------------------------------------------------------
# Manual tracing
#------------------------------------------------------------------------------

# Data
manual_data <- data_smoothed %>% filter(p_traced_auto == 0,
                                        p_smartphone_overall == 0.8)

# R_eff
reff_contour_manual_fwd <- manual_data %>% filter(backtrace_distance == 0) %>%
  contour_window(palette = fwd_palette_reff_sparse, 
                 contour_breaks = seq(0,10,step_sparse), upper = TRUE)
reff_contour_manual_bidir <- manual_data %>% filter(backtrace_distance == Inf) %>%
  contour_window(palette = bidir_palette_reff_sparse, 
                 contour_breaks = seq(0,10,step_sparse), lower = TRUE)
reff_contour_manual <- plot_grid(reff_contour_manual_fwd, 
                                 reff_contour_manual_bidir, nrow = 2)

# Control
ctrl_contour_manual_fwd <-  manual_data %>% filter(backtrace_distance == 0) %>%
  contour_window(contour_plot_ctrl_abs, palette = fwd_palette_ctrl_dense,
                 contour_breaks = seq(0,100,step_dense*100), upper = TRUE)
ctrl_contour_manual_bidir <-  manual_data %>% filter(backtrace_distance == Inf) %>%
  contour_window(contour_plot_ctrl_abs, palette = bidir_palette_ctrl_dense,
                 contour_breaks = seq(0,100,step_dense*100), lower = TRUE)
ctrl_contour_manual <- plot_grid(ctrl_contour_manual_fwd, 
                                 ctrl_contour_manual_bidir, nrow = 2)

# Overall grid
contour_manual <- plot_grid(reff_contour_manual, ctrl_contour_manual, nrow = 2,
                            labels = "auto", align = "hv",
                            axis = "l", label_size = fontsize_base * fontscale_label,
                            label_fontfamily = titlefont, label_colour = "black")


#------------------------------------------------------------------------------
# Hybrid tracing
#------------------------------------------------------------------------------

# Data
hybrid_data <- data_smoothed %>% filter(p_traced_auto != 0)

# R_eff
reff_contour_hybrid_fwd <- hybrid_data %>% filter(backtrace_distance == 0) %>%
  contour_window(palette = fwd_palette_reff_sparse, 
                 contour_breaks = seq(0,10,step_sparse), upper = TRUE)
reff_contour_hybrid_bidir <- hybrid_data %>% filter(backtrace_distance == Inf) %>%
  contour_window(palette = bidir_palette_reff_sparse, 
                 contour_breaks = seq(0,10,step_sparse), lower = TRUE)
reff_contour_hybrid <- plot_grid(reff_contour_hybrid_fwd, 
                                 reff_contour_hybrid_bidir, nrow = 2)

# Control
ctrl_contour_hybrid_fwd <-  hybrid_data %>% filter(backtrace_distance == 0) %>%
  contour_window(contour_plot_ctrl_abs, palette = fwd_palette_ctrl_dense,
                 contour_breaks = seq(0,100,step_dense*100), upper = TRUE)
ctrl_contour_hybrid_bidir <-  hybrid_data %>% filter(backtrace_distance == Inf) %>%
  contour_window(contour_plot_ctrl_abs, palette = bidir_palette_ctrl_dense,
                 contour_breaks = seq(0,100,step_dense*100), lower = TRUE)
ctrl_contour_hybrid <- plot_grid(ctrl_contour_hybrid_fwd, 
                                 ctrl_contour_hybrid_bidir, nrow = 2)

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/si_window"
plot_scale_cm <- 11
plot_scale_in <- plot_scale_cm/2.54

# Manual
cowplot::save_plot(filename=paste0(path_prefix, "_manual", ".png"),
       plot=contour_manual, ncol=2.5, nrow=3.75, base_height=plot_scale_in,
       base_asp = 1)

# Hybrid R_eff
cowplot::save_plot(filename=paste0(path_prefix, "_hybrid_reff", ".png"),
                   plot=reff_contour_hybrid, ncol=2.5, nrow=3.25,
                   base_height=plot_scale_in, base_asp = 1)

# Hybrid control
cowplot::save_plot(filename=paste0(path_prefix, "_hybrid_ctrl", ".png"),
                   plot=ctrl_contour_hybrid, ncol=2.5, nrow=3.25,
                   base_height=plot_scale_in, base_asp = 1)

