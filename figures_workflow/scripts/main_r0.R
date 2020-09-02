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
ctrl_upper <- snakemake@params[["ctrl_upper"]]
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
main_path <- snakemake@output[["main"]]
ascertainment_path <- snakemake@output[["si_ascertainment"]]
uptake_path <- snakemake@output[["si_uptake"]]
window_path <- snakemake@output[["si_window"]]
test_path <- snakemake@output[["si_test"]]
# main_path <- "output_files/dev_main_r0.png"
# ascertainment_path <- "output_files/dev_si_r0_ascertainment.png"
# uptake_path <- "output_files/dev_si_r0_uptake.png"
# window_path <- "output_files/dev_si_r0_window.png"
# test_path <- "output_files/dev_si_r0_test.png"

# Modify theme
theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.2, r = 0.2, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm")),
  legend.box.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.key = element_rect(colour = NA, fill = NA),
  strip.text.x = element_text(face = "bold", margin=margin(b=0.4, unit="cm")),
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

# Prepare combination
comb_data <- bind_rows(optim_data, pessim_data, median_data) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")),
  ) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels = c(Inf, 0)))

# Specify filtering function

data_main <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         p_smartphone_overall == uptake_main,
         contact_limit_manual == window_main,
         trace_neg_symptomatic == trace_neg_main)

data_ascertainment <- comb_data %>% 
  filter(p_ident_sym != ascertainment_main,
         p_smartphone_overall == uptake_main,
         contact_limit_manual == window_main,
         trace_neg_symptomatic == trace_neg_main)

data_uptake <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         p_smartphone_overall != uptake_main,
         contact_limit_manual == window_main,
         trace_neg_symptomatic == trace_neg_main)

data_window <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         p_smartphone_overall == uptake_main,
         contact_limit_manual != window_main,
         trace_neg_symptomatic == trace_neg_main)

data_test <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         p_smartphone_overall == uptake_main,
         contact_limit_manual == window_main,
         trace_neg_symptomatic != trace_neg_main)

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

cat("\nGenerating plots...")

#------------------------------------------------------------------------------
# Plotting functions
#------------------------------------------------------------------------------

ctrl_plot <- function(data, coord_ratio_base = 1.5){
  r0_min <- min(data$r0_base)
  r0_max <- max(data$r0_base)
  g <- ggplot(data, aes(x=r0_base, y=p_controlled,
                        ymin=p_controlled_lower,
                        ymax=p_controlled_upper,
                        colour=trace_type,
                        linetype=backtrace_distance,
                        shape=backtrace_distance)) %>%
    format_ctrl(coord_ratio = (r0_max-r0_min)/coord_ratio_base, 
                ylab = "% outbreaks controlled") %>% 
    x_r0 %>% 
    colour_ttype %>% linetype_backtrace(label = label_backtrace_brief) %>%
    theme_external
  g <- g + facet_grid(.~scenario_type)
  return(g)
}

reff_plot <- function(data, coord_ratio_base = 1.5){
  r0_min <- min(data$r0_base)
  r0_max <- max(data$r0_base)
  g <- ggplot(data, aes(x=r0_base, y=effective_r0_mean,
                        colour=trace_type,
                        linetype=backtrace_distance,
                        shape=backtrace_distance)) %>%
    format_reff(coord_ratio = (r0_max-r0_min)/r0_max/coord_ratio_base) %>% 
    x_r0 %>% colour_ttype %>% 
    linetype_backtrace(label = label_backtrace_brief) %>%
    theme_external
  g <- g + facet_grid(.~scenario_type)
  return(g)
}

grid_plot <- function(data, coord_ratio_base = 1.5,
                      ctrl_upper = TRUE, join_scale = 20){
  ctrl <- ctrl_plot(data, coord_ratio_base)
  reff <- reff_plot(data, coord_ratio_base)
  if (ctrl_upper){
    ctrl <- ctrl %>% strip_upper(join_scale)
    reff <- reff %>% strip_lower(join_scale)
  } else {
    ctrl <- ctrl %>% strip_lower(join_scale)
    reff <- reff %>% strip_upper(join_scale)
  }
  grid_out <- plot_grid(ctrl, reff, nrow = 2, align = "vh", axis = "l")
  return(grid_out)
}

#------------------------------------------------------------------------------
# Make plots
#------------------------------------------------------------------------------

grid_main <- grid_plot(data_main)
grid_ascertainment <- grid_plot(data_ascertainment)
grid_uptake <- grid_plot(data_uptake)
grid_window <- grid_plot(data_window)
grid_test <- grid_plot(data_test)

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

save_fig <- function(path, graph){
  cowplot::save_plot(filename=path, plot=graph,
                     ncol = 3, nrow = 2.25, base_height = plot_scale_in,
                     base_asp = 1.2)
}

save_fig(main_path, grid_main)
save_fig(ascertainment_path, grid_ascertainment)
save_fig(uptake_path, grid_uptake)
save_fig(window_path, grid_window)
save_fig(test_path, grid_test)

cat("done.\n")