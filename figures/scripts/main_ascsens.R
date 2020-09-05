#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Specify input paths
asc_path <- snakemake@input[["ascertainment"]]
sens_path <- snakemake@input[["sensitivity"]]
# asc_path <- "data/main_ascsens_ascertainment_1k_scenario.tsv.gz"
# sens_path <- "data/main_ascsens_sensitivity_1k_scenario.tsv.gz"

# Specify parameters
uptake_main <- snakemake@params[["uptake_main"]]
window_main <- snakemake@params[["window_main"]]
legend_fill <- snakemake@params[["legend_fill"]]
plot_scale_cm <- snakemake@params[["panel_scale"]]
# uptake_main <- 0.8 # TODO: Switch to lower uptake
# window_main <- 6
# plot_scale_cm <- 13
# legend_fill <- "#F6F6F6"
plot_scale_in <- plot_scale_cm/2.54

# Specify output path(s)
main_path <- snakemake@output[[1]]
# main_path <- "output_files/dev_main_ascsens.png"

theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.2, r = 0.2, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm")),
  legend.box.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.key = element_rect(colour = NA, fill = NA),
)

cat("done.\n")


#==============================================================================
# Read in data
#==============================================================================

cat("\nReading in data...")

# Manual tracing vs ascertainment
asc_data <- suppressMessages(read_tsv(asc_path)) %>%
  filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital"))) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels=c(Inf,0)))

# Manual tracing vs sensitivity
sens_data <- suppressMessages(read_tsv(sens_path)) %>%
  filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital"))) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels=c(Inf,0)))

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Ascertainment
#------------------------------------------------------------------------------

cat("\nGenerating ascertainment plots...")

reff_plot_asc <- function(data){
  ggplot(data, aes(x=p_ident_sym, y=effective_r0_mean,
                   colour=trace_type,
                   linetype=backtrace_distance,
                   shape=backtrace_distance)) %>%
    format_reff(r0=2.5) %>% colour_ttype %>% linetype_backtrace %>% x_ascertainment %>%
    theme_external
}

asc_data_main <- asc_data %>% 
  filter(p_smartphone_overall == uptake_main,
         contact_limit_manual == window_main)

reff_asc_main <- reff_plot_asc(asc_data_main)

cat("done.\n")

#------------------------------------------------------------------------------
# Test sensitivity
#------------------------------------------------------------------------------
cat("\nGenerating sensitivity plots...")

reff_plot_sens <- function(data, tr_neg_sym = TRUE){
  data %>% filter(trace_neg_symptomatic == tr_neg_sym) %>%
    ggplot(aes(x=test_sensitivity, y=effective_r0_mean,
               colour=trace_type,
               linetype=backtrace_distance,
               shape=backtrace_distance)) %>%
    format_reff(r0=2.5) %>% colour_ttype %>% linetype_backtrace %>% 
    x_sensitivity %>% theme_external
}

sens_data_main <- sens_data %>% 
  filter(p_smartphone_overall == uptake_main,
         contact_limit_manual == window_main)

reff_sens_notest_main <- reff_plot_sens(sens_data_main)
reff_sens_test_main <- reff_plot_sens(sens_data_main, FALSE)

cat("done.\n")

#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

cat("\nAssembling final plot...")

grid_reff <- plot_grid(reff_asc_main + theme(legend.position = "none") + ggtitle("Ascertainment"),
                       reff_sens_test_main + theme(legend.position = "none") + ggtitle("Sensitivity\n(Test required)"),
                       reff_sens_notest_main + theme(legend.position = "none") + ggtitle("Sensitivity\n(Test not required)"),
                      labels = "auto", nrow = 1, ncol = 3, align = "hv",
                      axis = "l", label_size = fontsize_base * fontscale_label,
                      label_fontfamily = titlefont, label_colour = "black")

legend_a <- get_legend(reff_asc_main + theme(legend.justification = "center"))
grid_out <- plot_grid(grid_reff, legend_a, ncol = 1, rel_heights = c(1, 0.1))

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

cowplot::save_plot(filename=main_path, plot=grid_out,
                   ncol = 3, nrow = 1.1, base_height = plot_scale_in,
                   base_asp = 0.8)
cat("done.\n")
