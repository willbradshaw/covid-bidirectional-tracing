source("figures_local/scripts/format_plots.R")

#==============================================================================
# Read in data
#==============================================================================

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

# Manual tracing vs ascertainment
asc_path <- "figures_local/data/si_ascertainment_baseline_1000_scenario.tsv.gz"
asc_data <- suppressMessages(read_tsv(asc_path)) %>%
  filter(p_smartphone_overall == 0.8) %>%
  filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital"))) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels=c(Inf,0)))

# Manual tracing vs sensitivity
sens_path <- "figures_local/data/si_sensitivity_1000_scenario.tsv.gz"
sens_data <- suppressMessages(read_tsv(sens_path)) %>%
  filter(p_smartphone_overall == 0.8) %>%
  filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital"))) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels=c(Inf,0)))

#==============================================================================
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Trace probability
#------------------------------------------------------------------------------

ctrl_plot_asc <- function(data){
  ggplot(data, aes(x=p_ident_sym, y=p_controlled,
                   ymin=p_controlled_lower, ymax=p_controlled_upper,
                   colour=trace_type,
                   linetype=backtrace_distance,
                   shape=backtrace_distance)) %>%
    format_ctrl %>% colour_ttype %>% linetype_backtrace %>% x_ascertainment %>%
    theme_external
}

reff_plot_asc <- function(data){
  ggplot(data, aes(x=p_ident_sym, y=effective_r0_mean,
                   colour=trace_type,
                   linetype=backtrace_distance,
                   shape=backtrace_distance)) %>%
    format_reff(r0=2.5) %>% colour_ttype %>% linetype_backtrace %>% x_ascertainment %>%
    theme_external
}


ctrl_asc <- ctrl_plot_asc(asc_data)
reff_asc <- reff_plot_asc(asc_data)

#------------------------------------------------------------------------------
# Test sensitivity
#------------------------------------------------------------------------------

ctrl_plot_sens <- function(data, tr_neg_sym = TRUE){
  data %>% filter(trace_neg_symptomatic == tr_neg_sym) %>%
  ggplot(aes(x=test_sensitivity, y=p_controlled,
                   ymin=p_controlled_lower, ymax=p_controlled_upper,
             colour=trace_type,
             linetype=backtrace_distance,
             shape=backtrace_distance)) %>%
    format_ctrl %>% colour_ttype %>% linetype_backtrace %>% 
    x_sensitivity %>% theme_external
}

reff_plot_sens <- function(data, tr_neg_sym = TRUE){
  data %>% filter(trace_neg_symptomatic == tr_neg_sym) %>%
    ggplot(aes(x=test_sensitivity, y=effective_r0_mean,
               colour=trace_type,
               linetype=backtrace_distance,
               shape=backtrace_distance)) %>%
    format_reff(r0=2.5) %>% colour_ttype %>% linetype_backtrace %>% 
    x_sensitivity %>% theme_external
}

ctrl_sens_notest <- ctrl_plot_sens(sens_data)
reff_sens_notest <- reff_plot_sens(sens_data)

ctrl_sens_test <- ctrl_plot_sens(sens_data, FALSE)
reff_sens_test <- reff_plot_sens(sens_data, FALSE)


#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

grid_reff <- plot_grid(reff_asc + theme(legend.position = "none") + ggtitle("Ascertainment"),
                       reff_sens_notest + theme(legend.position = "none") + ggtitle("Sensitivity\n(Test not required)"),
                       reff_sens_test + theme(legend.position = "none") + ggtitle("Sensitivity\n(Test required)"),
                      labels = "auto", nrow = 1, ncol = 3, align = "hv",
                      axis = "l", label_size = fontsize_base * fontscale_label,
                      label_fontfamily = titlefont, label_colour = "black")

legend_a <- get_legend(reff_asc + theme(legend.justification = "center"))
grid_out <- plot_grid(grid_reff, legend_a, ncol = 1, rel_heights = c(1, 0.1))

#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/fig2.5"

save_plot(path_prefix, "", grid_out, 13, 2)