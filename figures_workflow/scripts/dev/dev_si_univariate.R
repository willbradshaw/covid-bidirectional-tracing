#==============================================================================
# Preamble
#==============================================================================

# logfile <- file(snakemake@log[[1]], open = "wt")
# sink(logfile ,type = "output")
# sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Specify parameters
#legend_fill <- snakemake@params[["legend_fill"]]
#plot_scale_cm <- snakemake@params[["panel_scale"]]
legend_fill <- "#F6F6F6"
plot_scale_cm <- 10

plot_scale_in <- plot_scale_cm/2.54

# Specify input paths
# input_paths <- list(
#   environmental = snakemake@input[["environmental"]],
#   overdispersion = snakemake@input[["overdispersion"]],
#   presymptomatic = snakemake@input[["presymptomatic"]],
#   sensitivity = snakemake@input[["sensitivity"]],
#   ascertainment = snakemake@input[["ascertainment"]],
#   isolation = snakemake@input[["isolation"]]
# )
input_paths <- list(
  environmental = "data/si_univariate_environmental_500_scenario.tsv.gz",
  overdispersion = "data/si_univariate_overdispersion_500_scenario.tsv.gz",
  presymptomatic = "data/si_univariate_presymptomatic_500_scenario.tsv.gz",
  sensitivity = "data/si_univariate_sensitivity_500_scenario.tsv.gz",
  ascertainment = "data/si_univariate_ascertainment_500_scenario.tsv.gz",
  isolation = "data/si_univariate_isolation_500_scenario.tsv.gz"
)

input_paths

# Specify output paths
# output_paths <- list(
#   environmental = snakemake@output[["environmental"]],
#   overdispersion = snakemake@output[["overdispersion"]],
#   presymptomatic = snakemake@output[["presymptomatic"]],
#   sensitivity = snakemake@output[["sensitivity"]],
#   ascertainment = snakemake@output[["ascertainment"]],
#   isolation = snakemake@output[["isolation"]]
# )
output_paths <- list(
  environmental = "output_files/dev_si_univariate_environmental.png",
  overdispersion = "output_files/dev_si_univariate_overdispersion.png",
  presymptomatic = "output_files/dev_si_univariate_presymptomatic.png",
  sensitivity = "output_files/dev_si_univariate_sensitivity.png",
  ascertainment = "output_files/dev_si_univariate_ascertainment.png",
  isolation = "output_files/dev_si_univariate_isolation.png"
)

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
  aspect.ratio = 3/4,
  panel.spacing.y = unit(0.4, "cm"),
)

cat("done.\n")

#==============================================================================
# Read in data
#==============================================================================

cat("\nReading in data...")

data <- lapply(input_paths, function(p) suppressMessages(read_tsv(p)))

cat("done.\n")

#==============================================================================
# Process data
#==============================================================================

cat("\nProcessing data...")

# Add trace types
data_processed <- lapply(data, function(d)
  mutate(d, trace_type = ifelse(p_traced_auto == 0,
                                ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                                ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")),
         trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels = c(Inf, 0))))
       
# Infer presymptomatic transmission
alphas <- c(-1.38, -0.73, -0.325, 0, 0.325, 0.73, 1.38, 0.064, 0.397, -0.095)
p_presym <- c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.48, 0.38, 0.53)
a_presym <- tibble(generation_alpha = alphas, p_presymptomatic = p_presym)
data_processed$presymptomatic <- inner_join(data_processed$presymptomatic,
                                            a_presym, by = "generation_alpha")

# Infer average isolation delays
delay_times <- data_processed$isolation %>% group_by(delay_time) %>% summarise %>%
  mutate(delay_time_mean = sapply(delay_time, function(d) mean(eval(parse(text=d))(1e6)))) %>%
  mutate(delay_time_mean = round(delay_time_mean, 2))
data_processed$isolation <- inner_join(data_processed$isolation, delay_times,
                                       by = "delay_time")

# TODO: Any variable-specific operations go here

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

cat("\nGenerating plots...")

#------------------------------------------------------------------------------
# Plotting functions
#------------------------------------------------------------------------------

ctrl_plot <- function(data, xvar, xscale, ylim = 1,
                      ybreaks = seq(0,1,0.2)){
  g <- ggplot(data, aes_string(x=xvar, y="p_controlled",
                               ymin="p_controlled_lower",
                               ymax="p_controlled_upper",
                               colour="trace_type",
                               linetype="backtrace_distance",
                               shape="backtrace_distance")) + xscale
  h <- g %>% format_ctrl(ylab = "% outbreaks controlled", 
                         ylimits = c(0, ylim), ybreaks = ybreaks,
                         coord_ratio = NULL) %>% 
    colour_ttype %>% linetype_backtrace(label = label_backtrace_brief) %>%
    theme_external
  return(h)
}

reff_plot <- function(data, xvar, xscale, r0_lab_x){
  r0_min <- min(data$r0_base)
  r0_max <- max(data$r0_base)
  g <- ggplot(data, aes_string(x=xvar, y="effective_r0_mean",
                               colour="trace_type",
                               linetype="backtrace_distance",
                               shape="backtrace_distance")) + xscale
  h <- g %>% format_reff(r0 = r0_max, r0_lab_size = fontsize_base * 5/14,
                         r0_lab_x = r0_lab_x, coord_ratio = NULL) %>%
    colour_ttype %>% 
    linetype_backtrace(label = label_backtrace_brief) %>%
    theme_external
  return(h)
}

grid_plot_base <- function(reff, ctrl, legend_height = 0.2){
  grid_base <- plot_grid(reff + theme(legend.position = "none"),
                         ctrl + theme(legend.position = "none"),
                         ncol = 2, align = "vh", axis = "l",
                         label_size = fontsize_base * fontscale_label,
                         labels = "auto",
                         label_fontfamily = titlefont, label_fontface = "bold")
  grid_out <- plot_grid(grid_base, get_legend(ctrl), nrow = 2,
                        rel_heights = c(1, legend_height))
  return(grid_out)
}

grid_plot <- function(data, xvar, xscale, ylim_ctrl, ybreaks_ctrl,
                      r0_lab_x_reff, legend_height = 0.2){
  ctrl <- ctrl_plot(data, xvar, xscale, ylim_ctrl, ybreaks_ctrl)
  reff <- reff_plot(data, xvar, xscale, r0_lab_x_reff)
  return(grid_plot_base(reff, ctrl, legend_height))
}

label_policy <- function(x){
  ifelse(x, "Test not required", "Test required")
}

#------------------------------------------------------------------------------
# Make plots
#------------------------------------------------------------------------------

# x-scales
xscale_environmental  <- scale_x_continuous(name = "% environmental transmission", 
                                            breaks = seq(0,1,0.1), labels = label_pc)
xscale_ascertainment  <- scale_x_continuous(name = "% symptomatic cases identified", 
                                            breaks = seq(0,1,0.2), labels = label_pc)
xscale_presymptomatic <- scale_x_continuous(name = "% presymptomatic transmission", 
                                            breaks = seq(0,1,0.2), labels = label_pc)
xscale_isolation      <- scale_x_continuous(name="Mean isolation delay (days)",
                                            breaks = seq(0,5,1))
xscale_sensitivity    <- scale_x_continuous(name = "Test sensitivity (%)", 
                                            breaks = seq(0,1,0.2), labels = label_pc)
xscale_overdispersion <- scale_x_log10(name = "Overdispersion (k)")

# Basic plots
plot_environmental  <- grid_plot(data_processed$environmental, "p_environmental",
                                 xscale_environmental, 0.5, seq(0,1,0.1), 0.5)
plot_ascertainment  <- grid_plot(data_processed$ascertainment, "p_ident_sym",
                                 xscale_ascertainment, 0.5, seq(0,1,0.1), 1)
plot_presymptomatic <- grid_plot(data_processed$presymptomatic, "p_presymptomatic",
                                 xscale_presymptomatic, 0.5, seq(0,1,0.1), 0.25)
plot_isolation      <- grid_plot(data_processed$isolation, "delay_time_mean",
                                 xscale_isolation, 0.65, seq(0,1,0.2), 0.2)
plot_overdispersion <- grid_plot(data_processed$overdispersion, "dispersion",
                                 xscale_overdispersion, 1, seq(0,1,0.2), 0.013)

# Sensitivity (faceted by testing strategy)
reff_sensitivity <- reff_plot(data_processed$sensitivity, "test_sensitivity", 
                              xscale_sensitivity, 1) +
  facet_grid(trace_neg_symptomatic~., labeller=labeller(trace_neg_symptomatic = label_policy))
ctrl_sensitivity <- ctrl_plot(data_processed$sensitivity, "test_sensitivity",
                              xscale_sensitivity, 0.5, seq(0,1,0.1)) +
  facet_grid(trace_neg_symptomatic~., labeller=labeller(trace_neg_symptomatic = label_policy))
plot_sensitivity <- grid_plot_base(reff_sensitivity, ctrl_sensitivity, 0.1)

plots <- list(
  environmental = plot_environmental,
  overdispersion = plot_overdispersion,
  presymptomatic = plot_presymptomatic,
  sensitivity = plot_sensitivity,
  ascertainment = plot_ascertainment,
  isolation = plot_isolation
)

cat("done.\n")

# TODO: Fix labels & horizontal panel spacing

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

save_fig <- function(path, graph, label){
  nrow <- ifelse(label == "sensitivity", 1.7, 1.3)
  cowplot::save_plot(filename=path, plot=graph,
                     ncol = 2.4, nrow = nrow, base_height = plot_scale_in,
                     base_asp = 1)
}

lapply(names(data_processed), function(n)
  save_fig(output_paths[[n]], plots[[n]], n))

cat("done.\n")