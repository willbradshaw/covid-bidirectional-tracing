# logfile <- file(snakemake@log[[1]], open = "wt")
# sink(logfile ,type = "output")
# sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Input paths
# input_path <- snakemake@input[[1]]
input_path <- "data/si_asymptomatic_500_scenario.tsv.gz"

# Output paths & parameters
#plot_scale_cm <- snakemake@params[["panel_scale"]]
# out_path_reff_high <- snakemake@output[["reff_high"]]
# out_path_ctrl_high <- snakemake@output[["ctrl_high"]]
# out_path_reff_low <- snakemake@output[["reff_low"]]
# out_path_ctrl_low <- snakemake@output[["ctrl_low"]]
plot_scale_cm <- 7
out_path_reff_high <- "output_files/dev_si_asymptomatic_reff_highuptake.png"
out_path_ctrl_high <- "output_files/dev_si_asymptomatic_ctrl_highuptake.png"
out_path_reff_low <- "output_files/dev_si_asymptomatic_reff_lowuptake.png"
out_path_ctrl_low <- "output_files/dev_si_asymptomatic_ctrl_lowuptake.png"

plot_scale_in <- plot_scale_cm/2.54
uptake_high <- 0.8
uptake_low  <- 0.53

# Modify theme
theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.2, r = 0.2, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm")),
  aspect.ratio = 1,
  legend.position = "none",
  panel.spacing = unit(0.5, "cm"),
)

cat("done.\n")


#==============================================================================
# Read in data
#==============================================================================

cat("\nPreparing data...")

data <- suppressMessages(read_tsv(input_path))

#==============================================================================
# Process data
#==============================================================================

# Label variables
# Contour variables: p_asymptomatic, rel_r0_asymptomatic
# Other variables: scenario_type, backtrace_distance, contact_limit_manual, p_traced_auto

data_labelled <- data %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")),
         bidir_type = ifelse(backtrace_distance == 0, "Forward-only tracing",
                             "Bidirectional tracing"),
         window = paste0(contact_limit_manual, "-day manual window"),
         label_col = factor(trace_type, levels = ttype_levels),
         label_row = paste0(bidir_type, ",\n", window))

# Prepare contour data
data_indexed <- index_data(data_labelled, key_x = "p_asymptomatic", 
                           key_y = "rel_r0_asymptomatic")
data_smoothed <- smooth_data(data_indexed, 
                             group_vars = c("label_col", "label_row",
                                            "trace_type", "bidir_type", "window",
                                            "p_smartphone_overall")) %>%
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

contour_asym_base <- function(g){
  g + scale_x_continuous(name = "% of asymptomatic carriers", limits = c(0,1),
                         breaks = seq(0,1,0.2), labels = label_pc) +
    scale_y_continuous(name = expression(paste("Relative ",italic(R)[0]," of asymptomatic carriers (%)")),
                       limits = c(0,1), breaks = seq(0,1,0.2), 
                       labels = label_pc) +
    facet_grid(label_row~label_col) + theme_base
}

contour_asym_reff_base <- function(data){
  g <- data %>% ggplot(aes(x=p_asymptomatic, y=rel_r0_asymptomatic, z=r_eff_smoothed))
  return(g)
}

contour_asym_ctrl_base <- function(data){
  g <- data %>% ggplot(aes(x=p_asymptomatic, y=rel_r0_asymptomatic, z=pc_controlled_smoothed))
  return(g)
}

contour_asym_abs_base <- function(g, contour_breaks){
  g <- g +   
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    geom_text_contour(colour="black", stroke = 0.1, breaks = contour_breaks,
                      check_overlap = TRUE, mapping = aes(label = ..level..),
                      size = fontsize_base * 5/14)
  return(contour_asym_base(g))
}

# Absolute values

contour_asym_reff <- function(data, step = 0.2){
  g <- contour_asym_reff_base(data)
  h <- contour_asym_abs_base(g, seq(0,10,step))
  return(h)
}

contour_asym_ctrl <- function(data, step = 10){
  g <- contour_asym_ctrl_base(data)
  h <- contour_asym_abs_base(g, seq(-10,110,step))
  return(h)
}

#------------------------------------------------------------------------------
# Plots
#------------------------------------------------------------------------------

# Filter by uptake
data_high <- data_smoothed %>% filter(p_smartphone_overall == uptake_high)
data_low <- data_smoothed %>% filter(p_smartphone_overall == uptake_low)

# R_eff plots
reff_contour_high <- contour_asym_reff(data_high)
reff_contour_low  <- contour_asym_reff(data_low)

# Ctrl plots
ctrl_contour_high <- contour_asym_ctrl(data_high)
ctrl_contour_low  <- contour_asym_ctrl(data_low)

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

save_fig <- function(path, graph){
  cowplot::save_plot(filename=path, plot=graph,
                     ncol = 4.5, nrow = 4.5, base_height = plot_scale_in,
                     base_asp = 1)
}

save_fig(out_path_reff_high, reff_contour_high)
save_fig(out_path_reff_low, reff_contour_low)
save_fig(out_path_ctrl_high, ctrl_contour_high)
save_fig(out_path_ctrl_low, ctrl_contour_low)

cat("done.\n")
