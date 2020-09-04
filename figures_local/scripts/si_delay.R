source("figures_local/scripts/format_plots.R")

#==============================================================================
# Read in data
#==============================================================================

delay_path <- "figures_local/data/si_delay_1k_scenario.tsv.gz"
delay_data <- suppressMessages(read_tsv(delay_path))

#==============================================================================
# Process data
#==============================================================================

# Trace types
delay_data_processed <- delay_data %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", 
                                    "Manual + digital")),
         backtrace_distance = factor(backtrace_distance, levels = c(Inf, 0)))

# Delay times
delay_times <- delay_data_processed %>% group_by(delay_time) %>% summarise %>%
  mutate(delay_time_mean = sapply(delay_time, function(d) mean(eval(parse(text=d))(1e6)))) %>%
  mutate(delay_time_mean = round(delay_time_mean, 2))

delay_data_means <- inner_join(delay_data_processed, delay_times, by="delay_time")

#==============================================================================
# Make plots
#==============================================================================

# Control plots
ctrl_delay <- ggplot(delay_data_means, 
                    aes(x=delay_time_mean, y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                        colour=trace_type,
                        linetype=backtrace_distance,
                        shape=backtrace_distance)) %>%
  format_ctrl(coord_ratio = 5.1) %>% colour_ttype(ncol=2, name="Trace type") %>% 
  linetype_backtrace(label=label_backtrace) %>%
  theme_external()
ctrl_delay <- ctrl_delay +
  geom_vline(xintercept = 3.83, colour = "#888888", linetype = "dashed",
             size = 0.7) +
  scale_x_continuous(name="Mean isolation delay (days)", limits=c(-0.05, 5.05),
                     breaks = seq(0,5,1))

# R_eff
reff_delay <- ggplot(delay_data_means, 
                     aes(x=delay_time_mean, y=effective_r0_mean,
                         colour=trace_type,
                         linetype=backtrace_distance,
                         shape=backtrace_distance)) %>%
  format_reff(r0 = 2.5, coord_ratio = 5.1, r0_lab_x = 0.5) %>% 
  colour_ttype(ncol=2, name="Trace type") %>% 
  linetype_backtrace(label=label_backtrace) %>%
  theme_external()
reff_delay <- reff_delay +
  geom_vline(xintercept = 3.83, colour = "#888888", linetype = "dashed",
             size = 0.7) +
  scale_x_continuous(name="Mean isolation delay (days)", limits=c(-0.05, 5.05),
                     breaks = seq(0,5,1))


legend <- get_legend(reff_delay)

delay_plot <- plot_grid(ctrl_delay + theme(legend.position = "none"),
                        reff_delay + theme(legend.position = "none"),
                        labels = "auto", nrow = 1, ncol = 2, align = "hv",
                        axis = "l", label_size = fontsize_base * fontscale_label,
                        label_fontfamily = titlefont, label_colour = "black")
delay_plot_legend <- plot_grid(delay_plot, legend, nrow = 2,
                               rel_heights = c(1, 0.15))


#==============================================================================
# Save output
#==============================================================================

plot_scale <- 14
plot_prefix <- "figures_local/img/si_delay"

save_plot(plot_prefix, "", delay_plot_legend, plot_scale = plot_scale,
          plot_ratio = 1.5)