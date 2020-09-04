source("figures_local/scripts/format_plots.R")

#==============================================================================
# Read in data
#==============================================================================

# Contour plots
contour_path <- "figures/data/fig2_contour_1k_scenario.tsv.gz"
contour_data <- suppressMessages(read_tsv(contour_path)) %>% 
  filter(generation_alpha == 0.064) %>% 
  prep_contour_data(0.9, "p_smartphone_overall", "p_data_sharing_auto",
                    c("p_traced_manual", "backtrace_distance", "contact_limit_manual"))


#==============================================================================
# Make plots
#==============================================================================

# Palettes
bidir_palette_sparse <- make_contour_palette_ctrl(col_bidir, step_sparse)
fwd_palette_sparse <- make_contour_palette_ctrl(col_fwd, step_sparse)

# Labelling
contour_label <- "% outbreaks controlled, given\n90% prob. of trace success"
label_contour <- function(contour_plot, size = 4, label=contour_label){
  contour_plot + annotate("label", x=0.02, y=0.02, label=label,
                          hjust = 0, vjust = 0, size = size,
                          label.size=NA, fill=alpha("white", 0.8))
}

#------------------------------------------------------------------------------
# Forward-only
#------------------------------------------------------------------------------

ctrl_contour_digital_fwd <- (contour_data %>%
  filter(p_traced_manual == 0, backtrace_distance == 0,
         contact_limit_manual == 7) %>%
  contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_sparse) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  ggtitle("Forward-only digital\ntracing")) %>%
  label_contour

ctrl_contour_hybrid_2day_fwd <- (contour_data %>%
  filter(p_traced_manual != 0, backtrace_distance == 0,
         contact_limit_manual == 2) %>%
  contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_sparse) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  ggtitle("Forward-only hybrid\ntracing (2-day manual limit)")) %>%
  label_contour

ctrl_contour_hybrid_7day_fwd <- (contour_data %>%
  filter(p_traced_manual != 0, backtrace_distance == 0,
         contact_limit_manual == 7) %>%
  contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_sparse) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
    ggtitle("Forward-only hybrid\ntracing (7-day manual limit)")) %>%
  label_contour

#------------------------------------------------------------------------------
# Bidirectional
#------------------------------------------------------------------------------

ctrl_contour_digital_bidir <- (contour_data %>%
                               filter(p_traced_manual == 0, backtrace_distance == Inf,
                                      contact_limit_manual == 7) %>%
                               contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                                                     palette = bidir_palette_sparse) +
                               scale_x_continuous(name = "% of cases with\nchirping smartphones",
                                                  labels = label_pc, breaks = seq(0,1,0.2)) +
                               scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                                                  breaks = seq(0,1,0.2)) +
                                 ggtitle("Bidirectional digital\ntracing")) %>%
  label_contour

ctrl_contour_hybrid_2day_bidir <- (contour_data %>%
                                   filter(p_traced_manual != 0, backtrace_distance == Inf,
                                          contact_limit_manual == 2) %>%
                                   contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                                                         palette = bidir_palette_sparse) +
                                   scale_x_continuous(name = "% of cases with\nchirping smartphones",
                                                      labels = label_pc, breaks = seq(0,1,0.2)) +
                                   scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                                                      breaks = seq(0,1,0.2)) +
                                     ggtitle("Bidirectional hybrid\ntracing (2-day manual limit)")) %>%
  label_contour

ctrl_contour_hybrid_7day_bidir <- (contour_data %>%
                                   filter(p_traced_manual != 0, backtrace_distance == Inf,
                                          contact_limit_manual == 7) %>%
                                   contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                                                         palette = bidir_palette_sparse) +
                                   scale_x_continuous(name = "% of cases with\nchirping smartphones",
                                                      labels = label_pc, breaks = seq(0,1,0.2)) +
                                   scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                                                      breaks = seq(0,1,0.2)) +
                                     ggtitle("Bidirectional hybrid\ntracing (7-day manual limit)")) %>%
  label_contour

#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

grid_fig2 <- plot_grid(ctrl_contour_digital_fwd, ctrl_contour_hybrid_2day_fwd,
                       ctrl_contour_hybrid_7day_fwd, ctrl_contour_digital_bidir,
                       ctrl_contour_hybrid_2day_bidir, ctrl_contour_hybrid_7day_bidir,
                       labels = "auto", nrow = 2, ncol = 3, align = "hv",
                       axis = "l", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/si_contour_ctrl"

save_plot(path_prefix, "", grid_fig2, 26, 1.2)
