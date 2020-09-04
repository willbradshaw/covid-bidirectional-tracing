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

#------------------------------------------------------------------------------
# R_eff contour plots
#------------------------------------------------------------------------------

# Palettes
fwd_palette_dense <- contour_data %>%
  filter(backtrace_distance == 0, contact_limit_manual == 7) %>%
  make_contour_palette_reff(col_fwd, step_dense, ., 0)
fwd_palette_sparse <- contour_data %>%
  filter(backtrace_distance == 0, contact_limit_manual == 7) %>%
  make_contour_palette_reff(col_fwd, step_sparse, ., 0)

# Digital fwdectional
reff_contour_digital <- contour_data %>%
  filter(p_traced_manual == 0, backtrace_distance == 0,
         contact_limit_manual == 7) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_sparse, breaks = seq(0,10,step_sparse)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))
reff_contour_hybrid <- contour_data %>%
  filter(p_traced_manual != 0, backtrace_distance == 0,
         contact_limit_manual == 7) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_dense, breaks = seq(0,10,step_dense)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))

# Add labels
contour_label <- "Effective reprod. number, given\n90% prob. of trace success"
label_contour <- function(contour_plot, size = 3.8, label=contour_label){
  contour_plot + annotate("label", x=0.02, y=0.02, label=label,
                          hjust = 0, vjust = 0, size = size,
                          label.size=NA, fill=alpha("white", 0.8))
}
reff_contour_digital_labelled <- label_contour(reff_contour_digital)
reff_contour_hybrid_labelled <- label_contour(reff_contour_hybrid)


#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

grid_fig2 <- plot_grid(reff_contour_digital_labelled +
                         ggtitle("Forward-only digital\ntracing"),
                       reff_contour_hybrid_labelled + 
                         ggtitle("Forward-only hybrid\ntracing (7-day manual window)"),
                       labels = "auto", nrow = 1, ncol = 2, align = "hv",
                       axis = "l", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/si_contour_fwd_reff"

save_plot(path_prefix, "", grid_fig2, 14, 1.7)
