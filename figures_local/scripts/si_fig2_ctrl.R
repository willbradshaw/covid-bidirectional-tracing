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

#==============================================================================
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Single-strategy R_eff plots
#------------------------------------------------------------------------------

vspace = -0.5

ctrl_manual <- ggplot(manual_data %>% filter(p_traced_auto == 0), 
                      aes(x=p_traced_manual, y=p_controlled,
                          ymin=p_controlled_lower, ymax=p_controlled_upper,
                               colour=backtrace_distance,
                               linetype=contact_limit_manual,
                               shape=contact_limit_manual)) %>%
  format_ctrl() %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 1), yjust = 1, xjust = 0, vspace=vspace)

ctrl_hybrid <- ggplot(manual_data %>% filter(p_traced_auto != 0), 
                      aes(x=p_traced_manual, y=p_controlled,
                          ymin=p_controlled_lower, ymax=p_controlled_upper,
                          colour=backtrace_distance,
                          linetype=contact_limit_manual,
                          shape=contact_limit_manual)) %>%
  format_ctrl %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 1), yjust = 1, xjust = 0, vspace=vspace)

ctrl_univ <- ggplot(univ_data, 
                      aes(x=p_traced_auto, y=p_controlled,
                          ymin=p_controlled_lower, ymax=p_controlled_upper,
                          colour=backtrace_distance,
                          linetype=contact_limit_auto,
                          shape=contact_limit_auto)) %>%
  format_ctrl %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 1), yjust = 1, xjust = 0, vspace=vspace)

#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

grid_out <- plot_grid(ctrl_manual + ggtitle("Manual tracing only"),
                       ctrl_univ + ggtitle("Digital tracing only\n(universal coverage)"),
                       ctrl_hybrid + ggtitle("Manual + digital\n(hybrid) tracing"),
                       labels = "auto", nrow = 1, ncol = 3, align = "hv",
                       axis = "l", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/si_fig2_ctrl"

save_plot(path_prefix, "", grid_out, 14, 2.3)