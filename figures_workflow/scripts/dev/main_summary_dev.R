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
#ascertainment_main <- snakemake@params[["ascertainment_main"]]
#trace_neg_main <- snakemake@params[["trace_neg_main"]]
#r0_main <- snakemake@params[["r0_main"]]
#ctrl_upper <- snakemake@params[["ctrl_upper"]]
#legend_fill <- snakemake@params[["legend_fill"]]
#plot_scale_cm <- snakemake@params[["panel_scale"]]

ascertainment_main <- 0.5
trace_neg_main <- FALSE
ctrl_upper <- TRUE
legend_fill <- "#F6F6F6"
plot_scale_cm <- 9

ctrl_upper_default <- ctrl_upper
plot_scale_in <- plot_scale_cm/2.54
uptake_high <- 0.8
uptake_low <- 0.53
blab_buffer_ctrl = 0.08
blab_buffer_reff = 0.20
blab_colour_default = "#377eb8"
err_width_default <- 0.2
bar_width_default <- 0.6
aspect_default <- 1/1.5
r0_x_default = 9
r0_ydiff_default <- 0.1
ydiff_default <- 0.1
blab_label_default <- "‡"
grid_heights <- c(1,1.7)

# Specify input paths
# optim_path <- snakemake@input[["optimistic"]]
# pessim_path <- snakemake@input[["pessimistic"]]
# median_path <- snakemake@input[["median"]]
optim_path <- "data/main_summary_optimistic_1k_scenario.tsv.gz"
pessim_path <- "data/main_summary_pessimistic_1k_scenario.tsv.gz"
median_path <- "data/main_summary_median_1k_scenario.tsv.gz"

# Specify output paths
#main_path <- snakemake@output[["main"]]
# r0_20_path <- snakemake@output[["r0_20_path"]]
# r0_30_path <- snakemake@output[["r0_30_path"]]
#ascertainment_path <- snakemake@output[["si_ascertainment"]]
#test_path <- snakemake@output[["si_test"]]
main_path <- "output_files/dev_main_summary.png"
r0_20_path <- "output_files/dev_si_summary_r0_20.png"
r0_30_path <- "output_files/dev_si_summary_r0_30.png"
ascertainment_path <- "output_files/dev_si_summary_ascertainment.png"
test_path <- "output_files/dev_si_summary_test.png"

# Modify theme
theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.4, r = 0.4, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
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

# Combine datasets
comb_data <- bind_rows(optim_data, pessim_data, median_data) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual"),
                             ifelse(p_traced_manual == 0, "Digital", "Hybrid")),
         uptake = ifelse(p_traced_auto == 0, "",
                         ifelse(p_smartphone_overall == uptake_high, ", high-uptake", ", low-uptake")),
         trace_window = ifelse(p_traced_manual == 0, "", paste0(", ", contact_limit_manual, "-day")),
         label = paste0(trace_type, uptake, trace_window),
         backtrace_distance = factor(backtrace_distance, levels = c(0, Inf))
         )

# Filter combination
data_main <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         r0_base == 2.5,
         trace_neg_symptomatic == trace_neg_main)
data_r0_20 <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         r0_base == 2.0,
         trace_neg_symptomatic == trace_neg_main)
data_r0_30 <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         r0_base == 3.0,
         trace_neg_symptomatic == trace_neg_main)
data_ascertainment <- comb_data %>% 
  filter(p_ident_sym != ascertainment_main,
         r0_base == 2.5,
         trace_neg_symptomatic == trace_neg_main)
data_test <- comb_data %>% 
  filter(p_ident_sym == ascertainment_main,
         r0_base == 2.5,
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

barplot_ctrl <- function(data, err_width = err_width_default,
                         bar_width = bar_width_default,
                         aspect_ratio = aspect_default){
  g <- ggplot(data, aes(x=factor(label, levels=label_levels),
                        fill=backtrace_distance)) +
    geom_col(aes(y=p_controlled), position="dodge", width=bar_width) +
    facet_grid(.~scenario_type) +
    geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                  width=err_width, position = position_dodge(width=bar_width)) + 
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    theme_base + theme(aspect.ratio = aspect_ratio)
  return(g %>% fill_backtrace(ncol = 2))
}

barplot_reff <- function(data, bar_width = bar_width_default,
                         aspect_ratio = aspect_default,
                         r0_lab_x = r0_x_default,
                         r0_lab_ydiff = r0_ydiff_default,
                         ylim_diff = ydiff_default){
  r0 <- data %>% group_by(r0_base) %>% summarise %>% pull(r0_base)
  g <- ggplot(data, aes(x=factor(label, levels=label_levels),
                        fill=backtrace_distance)) +
    geom_col(aes(y=effective_r0_mean), 
             position = "dodge", width=bar_width) +
    geom_hline(yintercept = 1, linetype = "dotted", size = 1, colour = "black") +
    geom_hline(yintercept = r0, linetype = "dotted", size = 1, colour = "red") +
    annotate("text", x = r0_lab_x, y = r0 - r0_lab_ydiff, hjust = 1, vjust = 1,
             colour = "red", size = fontsize_base * 5/14,
             label = expression(R[0])) +
    scale_y_continuous(name = "Mean effective\nreproduction number",
                       limits = c(0,r0 + ylim_diff),
                       breaks = seq(0,10,0.5)) +
    facet_grid(.~scenario_type) +
    theme_base + theme(aspect.ratio = aspect_ratio)
  return(g %>% fill_backtrace(ncol = 2))
}

label_current_practice <- function(g, data, yvar, blab_buffer, blab_colour,
                                   blab_label = "‡"){
  data[["y_raw"]] <- data[[yvar]]
  blab_tab <- data %>% filter(label %in% c("No tracing", "Manual, 2-day")) %>%
    group_by(label, scenario_type) %>%
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

grid_plot <- function(data, err_width = err_width_default,
                      bar_width = bar_width_default,
                      aspect_ratio = aspect_default,
                      r0_lab_x = r0_x_default,
                      r0_lab_ydiff = r0_ydiff_default,
                      ylim_diff = ydiff_default,
                      buffer_ctrl = blab_buffer_ctrl,
                      buffer_reff = blab_buffer_reff,
                      blab_colour = blab_colour_default,
                      blab_label = blab_label_default,
                      ctrl_upper = ctrl_upper_default,
                      join_scale = 20,
                      rel_heights = grid_heights){
  ctrl <- barplot_ctrl(data, err_width, bar_width, aspect_ratio) %>%
    label_current_practice(data, "p_controlled_upper", buffer_ctrl,
                           blab_colour, blab_label)
  reff <- barplot_reff(data, bar_width, aspect_ratio,
                       r0_lab_x, r0_lab_ydiff, ylim_diff) %>%
    label_current_practice(data, "effective_r0_mean", buffer_reff,
                           blab_colour, blab_label)
  if (ctrl_upper){
    ctrl <- ctrl %>% strip_upper(join_scale)
    reff <- reff %>% strip_lower(join_scale)
    grid_out <- plot_grid(ctrl, reff, ncol = 1, nrow = 2, 
                          align = "vh", axis = "l",
                          rel_heights = grid_heights)
  } else {
    ctrl <- ctrl %>% strip_lower(join_scale)
    reff <- reff %>% strip_upper(join_scale)
    grid_out <- plot_grid(reff, ctrl, ncol = 1, nrow = 2, 
                          align = "vh", axis = "l",
                          rel_heights = grid_heights)
  }
  return(grid_out)
}

#------------------------------------------------------------------------------
# Make plots
#------------------------------------------------------------------------------

grid_main <- grid_plot(data_main)
grid_r0_20 <- grid_plot(data_r0_20)
grid_r0_30 <- grid_plot(data_r0_30)
grid_ascertainment <- grid_plot(data_ascertainment)
grid_test <- grid_plot(data_test)

cat("done.\n")

#==============================================================================
# Save outputs
#==============================================================================

cat("\nSaving output...")

save_fig <- function(path, graph){
  cowplot::save_plot(filename=path, plot=graph,
                     ncol = 3, nrow = 2.65, base_height = plot_scale_in,
                     base_asp = 1.2)
}

save_fig(main_path, grid_main)
save_fig(r0_20_path, grid_r0_20)
save_fig(r0_30_path, grid_r0_30)
save_fig(ascertainment_path, grid_ascertainment)
save_fig(test_path, grid_test)

cat("done.\n")

