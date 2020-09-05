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
trace_neg_main <- snakemake@params[["trace_neg_main"]]
legend_fill <- snakemake@params[["legend_fill"]]
plot_scale_cm <- snakemake@params[["panel_scale"]]
# ascertainment_main <- 0.5
# trace_neg_main <- FALSE
# legend_fill <- "#F6F6F6"
# plot_scale_cm <- 7

plot_scale_in <- plot_scale_cm/2.54
uptake_high <- 0.8
uptake_low <- 0.53
blab_buffer_default = 0.20
blab_colour_default = "#377eb8"
bar_width_default <- 0.6
aspect_default <- 1/1.5
blab_label_default <- "‡"

# Specify input paths
optim_path <- snakemake@input[["optimistic"]]
pessim_path <- snakemake@input[["pessimistic"]]
median_path <- snakemake@input[["median"]]
# optim_path <- "data/main_summary_optimistic_1k_scenario.tsv.gz"
# pessim_path <- "data/main_summary_pessimistic_1k_scenario.tsv.gz"
# median_path <- "data/main_summary_median_1k_scenario.tsv.gz"

# Specify output paths
main_path <- snakemake@output[["main"]]
ascertainment_path <- snakemake@output[["si_ascertainment"]]
test_path <- snakemake@output[["si_test"]]
# main_path <- "output_files/dev_main_summary.png"
# ascertainment_path <- "output_files/dev_si_summary_ascertainment.png"
# test_path <- "output_files/dev_si_summary_test.png"

# Modify theme
theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.4, r = 0.4, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.3, unit = "cm")),
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle=45,hjust=1,vjust=1,
                             size = fontsize_base * fontscale_mid,
                             colour = "black"),
  legend.box.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.background = element_rect(fill = legend_fill, linetype = "blank"),
  legend.key = element_rect(colour = NA, fill = NA),
  legend.title = element_blank(),
  legend.text = element_text(margin = margin(l=0.3, r = 0.4, unit = "cm")),
  legend.spacing.x = unit(0, "cm"),
  legend.margin = margin(l=0, r=0, t=0, b=0, "cm"),
  legend.box.margin = margin(t=0.4, l=0.4, r = -.02, b = 0.4, unit="cm"),
  legend.box.spacing = unit(0.6, "cm"),
  strip.text.x = element_text(face = "bold", margin=margin(b=0.4, unit="cm")),
  aspect.ratio = aspect_default,
  panel.spacing.y = unit(0.65, "cm")
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

data <- bind_rows(optim_data, pessim_data, median_data)

cat("done.\n")

#==============================================================================
# Process data
#==============================================================================

cat("\nProcessing data...")

# Label datasets
data_processed <- data %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual"),
                             ifelse(p_traced_manual == 0, "Digital", "Hybrid")),
         uptake = ifelse(p_traced_auto == 0, "",
                         ifelse(p_smartphone_overall == uptake_high, ", high-uptake", ", low-uptake")),
         trace_window = ifelse(p_traced_manual == 0, "", paste0(", ", contact_limit_manual, "-day")),
         label = paste0(trace_type, uptake, trace_window),
         backtrace_distance = factor(backtrace_distance, levels = c(0, Inf)),
         r0_base_lab = format(r0_base, nsmall = 1)
         )

# Filter combinations
data_main <- data_processed %>% 
  filter(p_ident_sym == ascertainment_main,
         trace_neg_symptomatic == trace_neg_main)
data_ascertainment <- data_processed %>% 
  filter(p_ident_sym != ascertainment_main,
         trace_neg_symptomatic == trace_neg_main)
data_test <- data_processed %>% 
  filter(p_ident_sym == ascertainment_main,
         trace_neg_symptomatic != trace_neg_main)

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

cat("\nGenerating plots...")

#------------------------------------------------------------------------------
# Plotting functions
#------------------------------------------------------------------------------

label_levels <- c("No tracing", "Manual, 2-day",
                  "Manual, 6-day",
                  "Digital, low-uptake",
                  "Digital, high-uptake",
                  "Hybrid, low-uptake, 2-day",
                  "Hybrid, low-uptake, 6-day",
                  "Hybrid, high-uptake, 2-day",
                  "Hybrid, high-uptake, 6-day")

barplot_reff <- function(data, bar_width = bar_width_default){
  g <- ggplot(data, aes(x=factor(label, levels=label_levels),
                        fill=backtrace_distance)) +
    geom_col(aes(y=effective_r0_mean), 
             position = "dodge", width=bar_width) +
    geom_hline(yintercept = 1, linetype = "dotted", size = 1, colour = "black") +
    geom_hline(aes(yintercept = r0_base), linetype = "dotted", size = 1, colour = "red") +
    scale_y_continuous(name = "Mean effective reproduction number",
                       limits = c(0,NA),
                       breaks = seq(0,10,0.5),
                       expand = c(0,0)) +
    facet_grid(r0_base_lab~scenario_type, scales = "fixed",
               labeller = label_bquote(cols = .(scenario_type),
                                       rows = italic(R)[0]~"="~.(r0_base_lab))) +
    theme_base
  return(g %>% fill_backtrace(ncol = 2))
}

label_current_practice <- function(g, data, yvar, blab_buffer, blab_colour,
                                   blab_label = "‡"){
  data[["y_raw"]] <- data[[yvar]]
  blab_tab <- data %>% filter(label %in% c("No tracing", "Manual, 2-day")) %>%
    group_by(label, scenario_type, r0_base_lab) %>%
    filter(y_raw == max(y_raw)) %>%
    filter(row_number() == 1) %>%
    mutate(y=y_raw+blab_buffer, colour = blab_colour)
  h <- g + geom_text(aes(x=label, y=y, colour=colour),
                     label = blab_label, data = blab_tab,
                     size = fontsize_base * 5/14) +
    scale_colour_identity(guide = "legend", name = NULL,
                          labels = "Current practice   ")
  return(h)
}

out_plot <- function(data, bar_width = bar_width_default,
                     blab_buffer = blab_buffer_default,
                     blab_colour = blab_colour_default,
                     blab_label = blab_label_default){
  barplot_reff(data, bar_width) %>%
    label_current_practice(data, "effective_r0_mean", blab_buffer,
                           blab_colour, blab_label)
}

#------------------------------------------------------------------------------
# Make plots
#------------------------------------------------------------------------------

grid_main <- out_plot(data_main)
grid_ascertainment <- out_plot(data_ascertainment)
grid_test <- out_plot(data_test)

cat("done.\n")

#==============================================================================
# Save outputs
#==============================================================================

cat("\nSaving output...")

save_fig <- function(path, graph){
  cowplot::save_plot(filename=path, plot=graph,
                     ncol = 3.1, nrow = 4, base_height = plot_scale_in,
                     base_asp = 1/aspect_default)
}

save_fig(main_path, grid_main)
save_fig(ascertainment_path, grid_ascertainment)
save_fig(test_path, grid_test)

cat("done.\n")

