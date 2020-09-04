source("figures_local/scripts/format_plots.R")

#==============================================================================
# Read in data
#==============================================================================

# Digital tracing, universal coverage
univ_path <- "figures/data/fig2_univ_env_1k_scenario.tsv.gz"
univ_data <- suppressMessages(read_tsv(univ_path)) %>%
  filter(generation_alpha == 0.064) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(0,Inf)),
         contact_limit_auto = factor(contact_limit_auto, levels=c(Inf,2)))

# Digital tracing, high coverage
digital_path <- "figures_local/data/contour_line_1000_scenario.tsv.gz"
digital_data <- suppressMessages(read_tsv(digital_path)) %>%
  filter(generation_alpha == 0.064, p_smartphone_overall == 0.8,
         contact_limit_auto == Inf) %>%
  mutate(trace_type = "Digital tracing only", p_traced = p_traced_auto,
         backtrace_distance = factor(backtrace_distance, levels=c(0,Inf)))

# Manual tracing (with and without automated)
manual_path <- "figures/data/fig2_manual_equiv_1000_scenario.tsv.gz"
manual_data <- suppressMessages(read_tsv(manual_path)) %>%
  filter(generation_alpha == 0.064) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(0,Inf)),
         contact_limit_manual = factor(contact_limit_manual, levels=c(7,2)),
         p_traced = p_traced_manual,
         trace_type = ifelse(p_traced_auto == 0, "Manual tracing only", 
                           "Manual + digital tracing"),
       p_traced = p_traced_manual)

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
# Single-strategy R_eff plots
#------------------------------------------------------------------------------

vspace = -1.5

reff_manual <- ggplot(manual_data %>% filter(p_traced_auto == 0), 
                      aes(x=p_traced_manual, y=effective_r0_mean,
                               colour=backtrace_distance,
                               linetype=contact_limit_manual,
                               shape=contact_limit_manual)) %>%
  format_reff(r0=2.5) %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

reff_hybrid <- ggplot(manual_data %>% filter(p_traced_auto != 0), 
                      aes(x=p_traced_manual, y=effective_r0_mean,
                          colour=backtrace_distance,
                          linetype=contact_limit_manual,
                          shape=contact_limit_manual)) %>%
  format_reff(r0=2.5) %>% colour_backtrace %>% 
  linetype_window(label = label_limits_manual) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

reff_univ <- ggplot(univ_data, 
                      aes(x=p_traced_auto, y=effective_r0_mean,
                          colour=backtrace_distance,
                          linetype=contact_limit_auto,
                          shape=contact_limit_auto)) %>%
  format_reff(r0=2.5) %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

#------------------------------------------------------------------------------
# Comparative control plot
#------------------------------------------------------------------------------

comp_data <- bind_rows(manual_data %>% filter(contact_limit_manual == 7) %>%
                         select(-contact_limit_manual),
                       digital_data) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(Inf, 0)),
         trace_type = factor(trace_type, 
                             levels=c("Manual tracing only", "Manual + digital tracing",
                                      "Digital tracing only")))

ctrl_comp <- ggplot(comp_data,
                    aes(x=p_traced, y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                            colour=trace_type,
                            linetype=backtrace_distance,
                            shape=backtrace_distance)) %>%
  format_ctrl %>% colour_ttype(ncol=1, name=NULL) %>% 
  linetype_backtrace(label=label_backtrace) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 1), yjust = 1, xjust = 0, vspace=-1)

#------------------------------------------------------------------------------
# R_eff contour plots
#------------------------------------------------------------------------------

# Palettes
bidir_palette_dense <- contour_data %>%
  filter(backtrace_distance == Inf, contact_limit_manual == 7) %>%
  make_contour_palette_reff(col_bidir, step_dense, ., Inf)
bidir_palette_sparse <- contour_data %>%
  filter(backtrace_distance == Inf, contact_limit_manual == 7) %>%
  make_contour_palette_reff(col_bidir, step_sparse, ., Inf)

# Digital bidirectional
reff_contour_digital <- contour_data %>%
  filter(p_traced_manual == 0, backtrace_distance == Inf,
         contact_limit_manual == 7) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_sparse, breaks = seq(0,10,step_sparse)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))

reff_contour_hybrid <- contour_data %>%
  filter(p_traced_manual != 0, backtrace_distance == Inf,
         contact_limit_manual == 7) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_dense, breaks = seq(0,10,step_dense)) +
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

grid_fig2 <- plot_grid(reff_manual + ggtitle("Manual tracing only"),
                       reff_univ + ggtitle("Digital tracing only\n(universal coverage)"),
                       reff_contour_digital_labelled +
                         ggtitle("Bidirectional digital\ntracing (partial coverage)"),
                       reff_hybrid + ggtitle("Manual + digital\n(hybrid) tracing"),
                       ctrl_comp + ggtitle("% outbreaks controlled"),
                       reff_contour_hybrid_labelled + 
                         ggtitle("Bidirectional hybrid\ntracing (partial coverage)"),
                       labels = "auto", nrow = 2, ncol = 3, align = "hv",
                       axis = "l", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/fig2_reff"

save_plot(path_prefix, "", grid_fig2, 26, 1.2)