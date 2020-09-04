source("figures_local/scripts/format_plots.R")

#==============================================================================
# Read in data
#==============================================================================

data_path <- "figures_local/data/si_contour_uptake_lowasc_1k.tsv.gz"
data <- suppressMessages(read_tsv(data_path))

#==============================================================================
# Process data
#==============================================================================

# Label variables
# Contour variables: p_smartphone_overall, p_data_sharing_auto
# Other variables: contact_limit, backtrace_distance, scenario
data_labelled <- data %>%
  mutate(trace_type = ifelse(backtrace_distance == 0, "Forward-only tracing",
                             "Bidirectional tracing"),
         limit = paste0(contact_limit_manual, "-day manual limit"))

# Prepare contour data
data_indexed <- index_data(data_labelled, key_x = "p_smartphone_overall", 
                           key_y = "p_data_sharing_auto")
data_smoothed <- smooth_data(data_indexed, 
                             group_vars = c("trace_type", "limit")) %>%
  mutate(pc_controlled_smoothed = p_controlled_smoothed*100)

#==============================================================================
# Prepare palettes
#==============================================================================

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

#==============================================================================
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Plotting function
#------------------------------------------------------------------------------

contour_uptake <- function(data, contour_fn = contour_plot_reff_abs,
                           palette = bidir_palette_reff_sparse,
                           contour_breaks = seq(0,10,step_sparse),
                           upper = FALSE, lower = FALSE,
                           join_scale = 20){
  g <- data %>%
    contour_fn("p_smartphone_overall", "p_data_sharing_auto", palette = palette,
               breaks = contour_breaks, coord_ratio = 1) +
    # geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
    # geom_vline(xintercept = c(0.53, 0.8), linetype = "dashed", colour = "red") +
    scale_x_continuous(name = "% of cases with\nchirping smartphones",
                       limits = c(0,1), breaks = seq(0,1,0.2),
                       labels = label_pc) +
    scale_y_continuous(name = "% of cases sharing data",
                       limits = c(0,1), breaks = seq(0,1,0.2),
                       labels = label_pc) +
    facet_grid(trace_type~limit)
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

# R_eff
reff_contour_fwd <- data_smoothed %>% filter(backtrace_distance == 0) %>%
  contour_uptake(palette = fwd_palette_reff_sparse, 
                 contour_breaks = seq(0,10,step_sparse), upper = TRUE)
reff_contour_bidir <- data_smoothed %>% filter(backtrace_distance == Inf) %>%
  contour_uptake(palette = bidir_palette_reff_sparse, 
                 contour_breaks = seq(0,10,step_sparse), lower = TRUE)
reff_contour <- plot_grid(reff_contour_fwd,
                          reff_contour_bidir, nrow = 2)

# Control
ctrl_contour_fwd <-  data_smoothed %>% filter(backtrace_distance == 0) %>%
  contour_uptake(contour_plot_ctrl_abs, palette = fwd_palette_ctrl_dense,
                 contour_breaks = seq(0,100,step_dense*100), upper = TRUE)
ctrl_contour_bidir <-  data_smoothed %>% filter(backtrace_distance == Inf) %>%
  contour_uptake(contour_plot_ctrl_abs, palette = bidir_palette_ctrl_dense,
                 contour_breaks = seq(0,100,step_dense*100), lower = TRUE)
ctrl_contour <- plot_grid(ctrl_contour_fwd,
                          ctrl_contour_bidir, nrow = 2)

# Specific subfigure: bidirectional, R_eff, 2-day
reff_bidir_2day_raw <- data_smoothed %>%
  filter(backtrace_distance == Inf, contact_limit_manual == 2) %>%
  contour_uptake(palette = bidir_palette_reff_sparse, 
                 contour_breaks = seq(0,10,step_sparse))
reff_bidir_2day <- reff_bidir_2day_raw +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank())


#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/si_uptake_lowasc"
plot_scale_cm <- 11
plot_scale_in <- plot_scale_cm/2.54

# R_eff
cowplot::save_plot(filename=paste0(path_prefix, "_reff", ".png"),
                   plot=reff_contour, ncol=2.2, nrow=2.5,
                   base_height=plot_scale_in, base_asp = 1)

# Control
cowplot::save_plot(filename=paste0(path_prefix, "_ctrl", ".png"),
                   plot=ctrl_contour, ncol=2, nrow=2,
                   base_height=plot_scale_in, base_asp = 1)

# R_eff single
cowplot::save_plot(filename=paste0(path_prefix, "_reff_bidir_2day", ".png"),
                   plot=reff_bidir_2day, ncol=1.3, nrow=1.3,
                   base_height=plot_scale_in, base_asp = 1)
