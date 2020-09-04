logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Input paths
median_path <- snakemake@input[["median"]]
optim_path <- snakemake@input[["optimistic"]]
pessim_path <- snakemake@input[["pessimistic"]]
# median_path <- "data/si_r0_contour_median_1k_scenario.tsv.gz"
# optim_path <- "data/si_r0_contour_optimistic_1k_scenario.tsv.gz"
# pessim_path <- "data/si_r0_contour_pessimistic_1k_scenario.tsv.gz"

# Run parameters
ascertainment <- snakemake@params[["ascertainment"]]
# ascertainment <- 0.5

# Output paths & parameters
plot_scale_cm <- snakemake@params[["panel_scale"]]
abs_reff_path <- snakemake@output[["abs_reff"]]
abs_ctrl_path <- snakemake@output[["abs_ctrl"]]
bidir_reff_path <- snakemake@output[["bidir_reff"]]
bidir_ctrl_path <- snakemake@output[["bidir_ctrl"]]
hybrid_reff_path <- snakemake@output[["hybrid_reff"]]
hybrid_ctrl_path <- snakemake@output[["hybrid_ctrl"]]
hybrid_path <- snakemake@output[["hybrid"]]
# plot_scale_cm <- 9
# abs_reff_path <- "output_files/dev_si_r0-contour-abs-reff.png"
# abs_ctrl_path <- "output_files/dev_si_r0-contour-abs-ctrl.png"
# bidir_reff_path <- "output_files/dev_si_r0-contour-bidir-reff.png"
# bidir_ctrl_path <- "output_files/dev_si_r0-contour-bidir-ctrl.png"
# hybrid_reff_path <- "output_files/dev_si_r0-contour-hybrid-reff.png"
# hybrid_ctrl_path <- "output_files/dev_si_r0-contour-hybrid-ctrl.png"

plot_scale_in <- plot_scale_cm/2.54

# Modify theme
theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.2, r = 0.2, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm")),
  aspect.ratio = 1,
  legend.position = "none",
)

cat("done.\n")


#==============================================================================
# Read in data
#==============================================================================

cat("\nPreparing data...")

median_data <- suppressMessages(read_tsv(median_path)) %>%
  mutate(scenario_type = "Median scenario")

optim_data <- suppressMessages(read_tsv(optim_path)) %>%
  mutate(scenario_type = "Optimistic scenario")

pessim_data <- suppressMessages(read_tsv(pessim_path)) %>%
  mutate(scenario_type = "Pessimistic scenario")

data <- bind_rows(median_data, optim_data, pessim_data)

#==============================================================================
# Process data
#==============================================================================

# Label variables
# Contour variables: p_smartphone_overall, r0_base
# Other variables: scenario_type, backtrace_distance, contact_limit_manual, p_traced_auto

levels_col <- c("Median scenario",
                "Optimistic scenario",
                "Pessimistic scenario")

data_labelled <- data %>%
  filter(p_ident_sym == ascertainment) %>%
  mutate(bidir_type = ifelse(backtrace_distance == 0, "Forward-only tracing",
                             "Bidirectional tracing"),
         window = paste0(contact_limit_manual, "-day manual window"),
         label_col = scenario_type,
         label_row = paste0(bidir_type, ",\n", window)) %>%
  mutate(label_col = factor(label_col, levels = levels_col))

# Prepare contour data
data_indexed <- index_data(data_labelled, key_x = "r0_base", 
                           key_y = "p_smartphone_overall")
data_smoothed <- smooth_data(data_indexed, 
                             group_vars = c("label_col", "label_row",
                                            "p_traced_auto", "bidir_type",
                                            "window", "scenario_type")) %>%
  mutate(pc_controlled_smoothed = p_controlled_smoothed*100)

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

cat("\nPreparing plots...")

#------------------------------------------------------------------------------
# Plotting functions
#------------------------------------------------------------------------------

# Auxiliary functions

contour_r0_base <- function(g){
  g + geom_vline(xintercept = 2.5, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = c(0.53, 0.8), linetype = "dashed", colour = "red") +
    scale_y_continuous(name = "% of cases with chirping smartphones",
                       limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = label_pc) +
    scale_x_continuous(name = expression(italic(R)[0]),
                       limits = c(1,4), breaks = seq(0,10,0.5)) +
    facet_grid(label_row~label_col) + theme_base
}

contour_r0_reff_base <- function(data){
  g <- data %>% ggplot(aes(x=r0_base, y=p_smartphone_overall, z=r_eff_smoothed))
  return(g)
}

contour_r0_ctrl_base <- function(data){
  g <- data %>% ggplot(aes(x=r0_base, y=p_smartphone_overall, z=pc_controlled_smoothed))
  return(g)
}

contour_r0_rel_base <- function(g, contour_breaks){
  g <- g +   
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    geom_text_contour(colour="black", stroke = 0.1, breaks = contour_breaks,
                      check_overlap = TRUE, 
                      mapping = aes(label=paste0(ifelse(..level.. >= 0, "+", ""), ..level..)),
                      size = fontsize_base * 5/14)
  return(contour_r0_base(g))
}

contour_r0_abs_base <- function(g, contour_breaks){
  g <- g +   
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    geom_text_contour(colour="black", stroke = 0.1, breaks = contour_breaks,
                      check_overlap = TRUE, mapping = aes(label = ..level..),
                      size = fontsize_base * 5/14)
  return(contour_r0_base(g))
}

# Absolute values

contour_r0_reff_abs <- function(data, step = 0.2){
  g <- contour_r0_reff_base(data)
  h <- contour_r0_abs_base(g, seq(0,10,step))
  return(h)
}

contour_r0_ctrl_abs <- function(data, step = 10){
  g <- contour_r0_ctrl_base(data)
  h <- contour_r0_abs_base(g, seq(0,110,step))
  return(h)
}

# Relative values

contour_r0_reff_rel <- function(data, step = 0.1){
  g <- contour_r0_reff_base(data)
  h <- contour_r0_rel_base(g, round(seq(-10,10,step), 1))
  return(h)
}

contour_r0_ctrl_rel <- function(data, step = 10){
  g <- contour_r0_ctrl_base(data)
  h <- contour_r0_rel_base(g, round(seq(-110,110,step)))
  return(h)
}


#------------------------------------------------------------------------------
# Absolute values
#------------------------------------------------------------------------------

data_abs <- data_smoothed %>% filter(p_traced_auto != 0)

reff_contour_abs <- contour_r0_reff_abs(data_abs)
ctrl_contour_abs <- contour_r0_ctrl_abs(data_abs)

#------------------------------------------------------------------------------
# Forward vs bidirectional
#------------------------------------------------------------------------------

data_comp_bidir <- data_smoothed %>% 
  filter(p_traced_auto != 0) %>%
  group_by(r0_base, p_smartphone_overall, p_traced_auto, window,
           scenario_type, label_col) %>%
  summarise(r_eff_smoothed = r_eff_smoothed[backtrace_distance == Inf] -
              r_eff_smoothed[backtrace_distance != Inf],
            pc_controlled_smoothed = pc_controlled_smoothed[backtrace_distance == Inf] -
              pc_controlled_smoothed[backtrace_distance != Inf]) %>%
  mutate(label_row = window)

reff_contour_rel_bidir <- contour_r0_reff_rel(data_comp_bidir)
ctrl_contour_rel_bidir <- contour_r0_ctrl_rel(data_comp_bidir)

#------------------------------------------------------------------------------
# Hybrid vs manual
#------------------------------------------------------------------------------

data_comp_hybrid <- data_smoothed %>% 
  group_by(r0_base, p_smartphone_overall, backtrace_distance, window,
           scenario_type, label_row, label_col) %>%
  summarise(r_eff_smoothed = r_eff_smoothed[p_traced_auto != 0] -
              r_eff_smoothed[p_traced_auto == 0],
            pc_controlled_smoothed = pc_controlled_smoothed[p_traced_auto != 0] -
              pc_controlled_smoothed[p_traced_auto == 0])

reff_contour_rel_hybrid <- contour_r0_reff_rel(data_comp_hybrid, 0.1)
ctrl_contour_rel_hybrid <- contour_r0_ctrl_rel(data_comp_hybrid)

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

# Absolute values
cat("\n\tAbsolute values...")
cowplot::save_plot(filename=abs_reff_path,
                   plot=reff_contour_abs, nrow=4.5, ncol=3.5, 
                   base_height=plot_scale_in, base_asp = 1)
cowplot::save_plot(filename=abs_ctrl_path,
                   plot=ctrl_contour_abs, nrow=4.5, ncol=3.5, 
                   base_height=plot_scale_in, base_asp = 1)

cat("\n\tForward-only vs bidirectional...")
cowplot::save_plot(filename=bidir_reff_path,
                   plot=reff_contour_rel_bidir, nrow=2.5, ncol=3.5, 
                   base_height=plot_scale_in, base_asp = 1)
cowplot::save_plot(filename=bidir_ctrl_path,
                   plot=ctrl_contour_rel_bidir, nrow=2.5, ncol=3.5, 
                   base_height=plot_scale_in, base_asp = 1)

cat("\n\tHybrid vs manual...")
cowplot::save_plot(filename=hybrid_reff_path,
                   plot=reff_contour_rel_hybrid, nrow=4.5, ncol=3.5, 
                   base_height=plot_scale_in, base_asp = 1)
cowplot::save_plot(filename=hybrid_ctrl_path,
                   plot=ctrl_contour_rel_hybrid, nrow=4.5, ncol=3.5, 
                   base_height=plot_scale_in, base_asp = 1)

cat("\n...done.\n")
