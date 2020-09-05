#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Specify parameters
ascertainment_main <- snakemake@params[["ascertainment_main"]]
uptake_main <- snakemake@params[["uptake_main"]]
window_main <- snakemake@params[["window_main"]]
trace_neg_main <- snakemake@params[["trace_neg_main"]]
legend_fill <- snakemake@params[["legend_fill"]]
plot_scale_cm <- snakemake@params[["panel_scale"]]
# ascertainment_main <- 0.5
# uptake_main <- 0.8
# window_main <- 6
# trace_neg_main <- FALSE
# ctrl_upper <- TRUE
# legend_fill <- "#F6F6F6"
# plot_scale_cm <- 9

plot_scale_in <- plot_scale_cm/2.54

# Specify input paths
optim_path <- snakemake@input[["optimistic"]]
pessim_path <- snakemake@input[["pessimistic"]]
median_path <- snakemake@input[["median"]]
# optim_path <- "data/main_r0_optimistic_1k_scenario.tsv.gz"
# pessim_path <- "data/main_r0_pessimistic_1k_scenario.tsv.gz"
# median_path <- "data/main_r0_median_1k_scenario.tsv.gz"

# Specify output paths
output_path <- snakemake@output[[1]]
# output_path <- "output_files/dev_si_rel_reff.png"

# Modify theme
theme_base <- theme_base + theme(
  plot.margin = margin(l = 1, r = 1, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm")),
  legend.box.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.key = element_rect(colour = NA, fill = NA),
  strip.text.x = element_text(margin=margin(b=0.4, unit="cm")),
  legend.margin = margin(r = 0.5, t = 0.2, b = 0.2, l = 0.2, unit = "cm"),
  aspect.ratio = 1
)

cat("done.\n")

#==============================================================================
# Read in data
#==============================================================================

cat("\nReading in data...")

# Optimistic scenario
optim_data <- suppressMessages(read_tsv(optim_path)) %>%
  mutate(scenario_type = "Optimistic scenario")

# Pessimistic scenario
pessim_data <- suppressMessages(read_tsv(pessim_path)) %>% 
  mutate(scenario_type = "Pessimistic scenario")

# Median scenario
median_data <- suppressMessages(read_tsv(median_path)) %>%
  mutate(scenario_type = "Median scenario")

cat("done.\n")

#==============================================================================
# Process data
#==============================================================================

cat("\nProcessing data...")

# Combine and filter datasets
data_comb <- bind_rows(optim_data, pessim_data, median_data) %>%
  filter(p_ident_sym == ascertainment_main,
         p_smartphone_overall == uptake_main,
         #contact_limit_manual == window_main,
         trace_neg_symptomatic == trace_neg_main) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "none", "manual"),
                             ifelse(p_traced_manual == 0, "digital", "hybrid")),
  ) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels = c(0, Inf)),
         contact_limit_manual = factor(contact_limit_manual, levels = c(6,2)))

# Spread by trace type and compute relative R_eff
data_spread <- data_comb %>% spread(trace_type, effective_r0_mean, sep = "_") %>%
  group_by(scenario_type, r0_base, backtrace_distance, contact_limit_manual) %>%
  summarise(trace_type_none = sum(trace_type_none, na.rm = TRUE),
            trace_type_manual = sum(trace_type_manual, na.rm = TRUE),
            trace_type_digital = sum(trace_type_digital, na.rm = TRUE),
            trace_type_hybrid = sum(trace_type_hybrid, na.rm = TRUE))

# Compute relative R_eff
data_rel <- data_spread %>%
  mutate(`Hybrid vs manual` = trace_type_hybrid/trace_type_manual,
         `Hybrid vs no tracing` = trace_type_hybrid/trace_type_none,
         `Manual vs no tracing` = trace_type_manual/trace_type_none) %>%
  select(-starts_with("trace_type")) %>%
  gather(comparison, rel_reff, -scenario_type, -r0_base, -backtrace_distance,
         -contact_limit_manual)

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

cat("\nGenerating plots...")

g_rel_reff_raw <- ggplot(data_rel, aes(x=r0_base, y=rel_reff,
                                       colour = backtrace_distance,
                                       linetype = contact_limit_manual, 
                                       shape = contact_limit_manual)) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = expression(Relative~italic(R)[eff])) +
  facet_grid(comparison~scenario_type)

g_rel_reff <- g_rel_reff_raw %>%
  x_r0 %>% colour_backtrace %>% linetype_window(label = label_limits_manual) %>%
  theme_external

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

cowplot::save_plot(filename=output_path, plot=g_rel_reff,
                   ncol = 3.1, nrow = 3.2, base_height = plot_scale_in,
                   base_asp = 1)
cat("done.\n")
