source("figures_local/scripts/format_plots.R")

#==============================================================================
# Read in data
#==============================================================================

# Optimistic scenario
optim_path <- "figures_local/data/r0_optimistic_lowasc_1k_scenario.tsv.gz"
optim_data <- suppressMessages(read_tsv(optim_path)) %>%
  filter(generation_alpha == 0.397) %>%
  mutate(scenario_type = "Optimistic scenario")

# Pessimistic scenario
pessim_path <- "figures_local/data/r0_pessimistic_lowasc_1k_scenario.tsv.gz"
pessim_data <- suppressMessages(read_tsv(pessim_path)) %>% 
  filter(generation_alpha == -0.095) %>%
  mutate(scenario_type = "Pessimistic scenario")

# Median scenario
median_path <- "figures_local/data/r0_median_lowasc_1k_scenario.tsv.gz"
median_data <- suppressMessages(read_tsv(median_path)) %>%
  filter(generation_alpha == 0.064) %>%
  mutate(scenario_type = "Median scenario")

#==============================================================================
# Process data
#==============================================================================

# Prepare combination
comb_data <- bind_rows(optim_data, pessim_data, median_data) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")),
  ) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels = c(Inf, 1, 0)))

# Filter to preferred values
p_traced <- 0.9
trace_neg_sym <- FALSE
trace_limit <- 7
uptake <- 0.53

comb_data_filtered <- comb_data %>%
  filter(p_traced_manual %in% c(0, -1, p_traced),
         p_traced_auto %in% c(0, -1, p_traced),
         trace_neg_symptomatic == trace_neg_sym,
         contact_limit_manual == trace_limit,
         p_smartphone_overall == uptake)
  

#==============================================================================
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Base plots
#------------------------------------------------------------------------------

# Control rate
ctrl_plot <- ggplot(comb_data_filtered, aes(x=r0_base, y=p_controlled,
                                            ymin=p_controlled_lower,
                                            ymax=p_controlled_upper,
                                            colour=trace_type,
                                            linetype=backtrace_distance,
                                            shape=backtrace_distance)) %>%
  format_ctrl(coord_ratio = 2, ylab = "% of outbreaks\ncontrolled") %>% x_r0 %>% 
  colour_ttype %>% linetype_backtrace(label = label_backtrace_gen) %>%
  theme_external
ctrl_plot <- ctrl_plot + facet_grid(.~scenario_type) +
  theme(strip.text.x = element_text(face = "bold"))

# R_eff
reff_plot <- ggplot(comb_data_filtered, aes(x=r0_base, y=effective_r0_mean,
                                            colour=trace_type,
                                            linetype=backtrace_distance,
                                            shape=backtrace_distance)) %>%
  format_reff(coord_ratio = (3/4)/1.5) %>% x_r0 %>% colour_ttype %>% 
  linetype_backtrace(label = label_backtrace_gen) %>%
  theme_external
reff_plot <- reff_plot + facet_grid(.~scenario_type) +
  theme(strip.text.x = element_text(face = "bold"))

#------------------------------------------------------------------------------
# Make grid
#------------------------------------------------------------------------------

ctrl_plot_top <- ctrl_plot + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  legend.position = "none",
  plot.margin = unit(c(-3.5, 1, -4.5, 0), "cm"),
)

reff_plot_btm <- reff_plot + theme(
  strip.text.x = element_blank(),
  strip.background.x = element_blank(),
  strip.background = element_blank(),
  plot.margin = unit(c(0, 1, 0, 0), "cm"),
)

# legend_out <- get_legend(reff_plot + 
#                            theme(legend.justification = "center",
#                                  legend.margin = margin(b=0.5,t=-0.5)))
  
out_plot <- plot_grid(ctrl_plot_top, reff_plot_btm,
                      nrow = 2, rel_heights = c(1,1.1), align = "v",
                      axis = "l")

#==============================================================================
# Save output
#==============================================================================

path_prefix <- "figures_local/img/r0_horiz_lowasc_lowuptake"
save_plot(path_prefix, "", out_plot, plot_scale = 15*1.2, plot_ratio = 1.5)