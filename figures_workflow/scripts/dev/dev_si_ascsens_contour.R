# logfile <- file(snakemake@log[[1]], open = "wt")
# sink(logfile ,type = "output")
# sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.2, r = 0.2, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm")),
  aspect.ratio = 1,
  legend.position = "none",
)

# Input paths
#input_path <- snakemake@input[[1]]
input_path <- "data/si_ascsens_contour_1k_scenario.tsv.gz"

# Output paths & parameters
#plot_scale_cm <- snakemake@params[["panel_scale"]]
#manual_path <- snakemake@output[["manual"]]
#hybrid_path <- snakemake@output[["hybrid"]]
plot_scale_cm <- 9
plot_scale_in <- plot_scale_cm/2.54
abs_path <- "output_files/dev_si_ascsens_contour_abs.png"
rel_path <- "output_files/dev_si_ascsens_contour_rel.png"

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
# Contour variables: p_ident_sym, test_sensitivity
# Other variables: p_smartphone_overall, p_traced_auto, backtrace_distance, trace_neg_symptomatic

levels_row <- c("Manual tracing",
                "Hybrid tracing,\nlow digital uptake",
                "Hybrid tracing,\nhigh digital uptake" )
levels_col <- c("Bidirectional tracing,\ntest required",
                "Bidirectional tracing,\ntest not required",
                "Forward-only tracing,\ntest required",
                "Forward-only tracing,\ntest not required")

data_labelled <- data %>%
  filter(p_traced_auto != 0 | p_smartphone_overall == 0.8) %>%
  mutate(trace_type = ifelse(p_traced_auto == 0, "Manual tracing",
                             "Hybrid tracing"),
         bidir_type = ifelse(backtrace_distance == 0, "Forward-only tracing",
                             "Bidirectional tracing"),
         uptake = paste(ifelse(p_smartphone_overall == 0.8, "high", "low"),
                        "digital uptake"),
         test_req = ifelse(trace_neg_symptomatic, "test not required",
                           "test required"),
         label_col = paste0(bidir_type, ",\n", test_req),
         label_row = ifelse(p_traced_auto == 0, trace_type,
                        paste0(trace_type, ",\n", uptake))) %>%
  mutate(label_col = factor(label_col, levels = levels_col),
         label_row = factor(label_row, levels = levels_row))

# Prepare contour data
data_indexed <- index_data(data_labelled, key_x = "p_ident_sym", 
                           key_y = "test_sensitivity")
data_smoothed <- smooth_data(data_indexed, 
                             group_vars = c("label_col", "label_row",
                                            "trace_type", "bidir_type",
                                            "uptake", "test_req")) %>%
  mutate(pc_controlled_smoothed = p_controlled_smoothed*100)

cat("done.\n")

#==============================================================================
# Make plots
#==============================================================================

cat("\nPreparing plots...")

#------------------------------------------------------------------------------
# Plotting function
#------------------------------------------------------------------------------

contour_ascsens_abs <- function(data, palette = bidir_palette_reff_sparse,
                           contour_breaks = seq(0,10,step_sparse)){
  g <- data %>%
    ggplot(aes(x=test_sensitivity, y=p_ident_sym, z=r_eff_smoothed)) +
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    geom_text_contour(colour="black", stroke = 0.2, breaks = contour_breaks,
                      check_overlap = TRUE, mapping = aes(label = ..level..),
                      size = fontsize_base * 5/14) +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
    scale_y_continuous(name = "% symptomatic cases identified", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = label_pc) +
    scale_x_continuous(name = "Test sensitivity (%)",
                       limits = c(0,1), breaks = seq(0,1,0.2),
                       labels = label_pc) +
    facet_grid(label_row~label_col) + theme_base
  return(g)
}

contour_ascsens_rel <- function(data, palette = bidir_palette_reff_sparse,
                                contour_breaks = round(seq(-10,10,step_sparse),1)){
  g <- data %>%
    ggplot(aes(x=test_sensitivity, y=p_ident_sym, z=r_eff_smoothed)) +
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    geom_text_contour(colour="black", stroke = 0.2, breaks = contour_breaks,
                      check_overlap = TRUE, size = fontsize_base * 5/14,
                      mapping = aes(label=paste0(ifelse(..level.. >= 0, "+", ""), ..level..))) +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
    scale_y_continuous(name = "% symptomatic cases identified", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = label_pc) +
    scale_x_continuous(name = "Test sensitivity (%)",
                       limits = c(0,1), breaks = seq(0,1,0.2),
                       labels = label_pc) +
    facet_grid(label_row~label_col) + theme_base
  return(g)
}



#------------------------------------------------------------------------------
# Absolute values
#------------------------------------------------------------------------------

reff_contour_abs <- data_smoothed %>% contour_ascsens_abs()

#------------------------------------------------------------------------------
# Forward vs bidirectional
#------------------------------------------------------------------------------

data_comp_bidir <- data_smoothed %>% group_by(trace_type, uptake, test_req,
                                              test_sensitivity, p_ident_sym,
                                              label_row) %>%
  summarise(r_eff_smoothed = r_eff_smoothed[bidir_type == "Bidirectional tracing"] - 
              r_eff_smoothed[bidir_type != "Bidirectional tracing"]) %>%
  mutate(label_col = str_to_sentence(test_req),
         label_col = factor(label_col, 
                            levels = c("Test required","Test not required")))

reff_contour_rel <- contour_ascsens_rel(data_comp_bidir)


cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

# Manual
cat("\n\tAbsolute values...")
cowplot::save_plot(filename=abs_path,
                   plot=reff_contour_abs, ncol=4.5, nrow=3.5, 
                   base_height=plot_scale_in, base_asp = 1)

cat("\n\tRelative values...")
cowplot::save_plot(filename=rel_path,
                   plot=reff_contour_rel, ncol=2.5, nrow=3.5, 
                   base_height=plot_scale_in, base_asp = 1)

cat("\n...done.\n")
