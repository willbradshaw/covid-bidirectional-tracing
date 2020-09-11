#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

theme_base <- theme_base + theme(
  plot.margin = margin(l = 0.2, r = 0.2, b = -0.3,
                       t = -0.2, unit = "cm"),
  axis.title.y = element_text(margin = margin(r = 0.15, unit = "cm")),
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm"))
)

# Input paths
manual_path <- snakemake@input[["manual"]]
digital_path <- snakemake@input[["digital"]]
uptake_path <- snakemake@input[["uptake"]]
window_path <- snakemake@input[["window"]]
# manual_path <- "data/main_strategies_manual_1k_scenario.tsv.gz"
# digital_path <- "data/main_strategies_digital_1k_scenario.tsv.gz"
# uptake_path <- "data/main_strategies_uptake_1k_scenario.tsv.gz"
# window_path <- "data/si_window_median_1k_scenario.tsv.gz"

# Output paths & parameters
plot_scale_main_cm <- snakemake@params[["panel_scale_main"]]
plot_scale_si_cm <- snakemake@params[["panel_scale_si"]]
main_path <- snakemake@output[["main"]]
si_uptake_reff_path <- snakemake@output[["si_uptake_reff"]]
si_digital_univ_path <- snakemake@output[["si_digital_univ"]]
si_hybrid_delta_path <- snakemake@output[["si_hybrid_delta"]]
# plot_scale_main_cm <- 13
# plot_scale_si_cm <- 11
# main_path <- "output_files/dev_main_strategies.png"
# si_uptake_reff_path <- "output_files/dev_si_uptake_reff.png"
# si_digital_univ_path <- "output_files/dev_si_digital_univ.png"
# si_hybrid_delta_path <- "output_files/dev_si_hybrid_delta.png"
plot_scale_main_in <- plot_scale_main_cm/2.54
plot_scale_si_in <- plot_scale_si_cm/2.54

cat("done.\n")

#==============================================================================
# Read in data
#==============================================================================

cat("\nReading in data...")

# Manual tracing (with and without automated)
manual_data <- suppressMessages(read_tsv(manual_path)) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(0,Inf)),
         contact_limit_manual = factor(contact_limit_manual, levels=c(6,2)),
         trace_type = ifelse(p_traced_auto == 0, "Manual tracing only", 
                             "Manual + digital tracing"),
         p_traced = p_traced_manual)

# Digital tracing alone (universal or high coverage)
digital_data <- suppressMessages(read_tsv(digital_path)) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(0,Inf)),
         contact_limit_auto = factor(contact_limit_auto, levels=c(Inf,2)),
         trace_type = "Digital tracing only",
         p_traced = p_traced_auto)

# Uptake contour plots
uptake_data <- suppressMessages(read_tsv(uptake_path)) %>% 
  filter(p_traced_manual == 0.9 | contact_limit_manual == 6) %>%
  prep_contour_data(0.9, "p_smartphone_overall", "p_data_sharing_auto",
                    c("p_traced_manual", "backtrace_distance", "contact_limit_manual"))

# Window contour plots

window_data <- suppressMessages(read_tsv(window_path)) %>% 
  filter(p_traced_auto == 0, backtrace_distance == Inf, p_smartphone_overall == 0.8) %>%
  index_data(key_x = "contact_limit_manual", key_y = "p_traced_manual") %>%
  smooth_data(group_vars = c())

# Relative R_eff data (manual vs hybrid)
manual_only_data <- manual_data %>%
  filter(p_traced_auto == 0, p_smartphone_overall == 0.8) %>%
  select(-p_smartphone_overall, -p_traced_auto) %>%
  rename(effective_r0_manual = effective_r0_mean)

rel_data_hybrid <- manual_data %>%
  filter(p_traced_auto != 0) %>%
  inner_join(manual_only_data, 
             by = c("contact_limit_manual", "p_traced_manual", "backtrace_distance")) %>%
  mutate(effective_r0_rel = effective_r0_mean - effective_r0_manual)

# Relative R_eff data (bidir vs current practice)
current_practice <- manual_data %>%
  filter(p_traced_auto == 0, p_smartphone_overall == 0.8, contact_limit_manual == 2) %>%
  select(-p_traced_auto, -p_smartphone_overall, -contact_limit_manual) %>%
  rename(effective_r0_current = effective_r0_mean) %>%
  mutate(effective_r0_current_rel = effective_r0_current - effective_r0_current[p_traced_manual == 0])
current_practice_fwd <- current_practice %>%
  filter(backtrace_distance == 0) %>% select(-backtrace_distance)
current_practice_bidir <- current_practice %>%
  filter(backtrace_distance == Inf) %>% select(-backtrace_distance)

rel_data_bidir_fwd <- manual_data %>%
  filter(backtrace_distance == Inf) %>%
  inner_join(current_practice_fwd, by="p_traced_manual") %>%
  mutate(effective_r0_rel_sub = effective_r0_mean - effective_r0_current,
         effective_r0_rel_div = effective_r0_rel_sub / effective_r0_current_rel)
rel_data_bidir_bidir <- manual_data %>%
  filter(backtrace_distance == Inf) %>%
  inner_join(current_practice_bidir, by="p_traced_manual") %>%
  mutate(effective_r0_rel = effective_r0_mean - effective_r0_current)
  
cat("done.\n")

#==============================================================================
# Make main-text figure
#==============================================================================

#------------------------------------------------------------------------------
# Single-strategy R_eff plots
#------------------------------------------------------------------------------

cat("\nGenerating single-strategy R_eff plots...")

vspace <- -1.5

cat("\n\t1")

reff_manual <- manual_data %>%
  filter(p_traced_auto == 0, p_smartphone_overall == 0.8) %>%
  ggplot(aes(x=p_traced_manual, y=effective_r0_mean,
             colour=backtrace_distance, 
             linetype=contact_limit_manual,
             shape=contact_limit_manual)) %>%
  format_reff(r0=2.5) %>% colour_backtrace %>% 
  linetype_window(label = label_windows) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

cat("\n\t2")

hybrid_data <- manual_data %>% filter(p_traced_auto != 0) %>%
  mutate(bidir_type = ifelse(backtrace_distance == 0, "Forward tracing only",
                             "Bidirectional tracing"),
         uptake_pc = paste0(round(p_smartphone_overall * 100), "%"),
         uptake_type = ifelse(backtrace_distance == 0 & p_smartphone_overall == 0.53,
                              paste0(" (", uptake_pc, "\nsmartphone coverage)"),
                              paste0(" (", uptake_pc, ")")),
         label = paste0(bidir_type, uptake_type)) %>%
  mutate(label = factor(label, levels = c("Forward tracing only (53%\nsmartphone coverage)",
                                          "Forward tracing only (80%)",
                                          "Bidirectional tracing (53%)",
                                          "Bidirectional tracing (80%)")))
reff_hybrid_raw <- hybrid_data %>%
  ggplot(aes(x=p_traced_manual, y=effective_r0_mean,
             colour=label, 
             linetype=contact_limit_manual,
             shape=contact_limit_manual)) %>%
  format_reff(r0=2.5) %>% #colour_backtrace %>% 
  linetype_window(label = label_windows_manual) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)
reff_hybrid <- reff_hybrid_raw + 
  scale_colour_manual(name = NULL, values = c("#1B9E77", "#7570B3", "#D95F02", "#E7298A")) +
  guides(colour=guide_legend(ncol = 1, order = 1)) + theme(
    legend.background = element_rect(fill = alpha("white", 0.8))
  )

cat("\n...done.\n")

#------------------------------------------------------------------------------
# Relative R_eff line plots
#------------------------------------------------------------------------------

# Hybrid - manual

rel_data_hybrid_labelled <- rel_data_hybrid %>%
  mutate(bidir_type = ifelse(backtrace_distance == 0, "Forward tracing only",
                             "Bidirectional tracing"),
         uptake_pc = paste0(round(p_smartphone_overall * 100), "%"),
         uptake_type = ifelse(backtrace_distance == 0 & p_smartphone_overall == 0.53,
                              paste0(" (", uptake_pc, "\nsmartphone coverage)"),
                              paste0(" (", uptake_pc, ")")),
         label = paste0(bidir_type, uptake_type)) %>%
  mutate(label = factor(label, levels = c("Forward tracing only (53%\nsmartphone coverage)",
                                          "Forward tracing only (80%)",
                                          "Bidirectional tracing (53%)",
                                          "Bidirectional tracing (80%)")))

reff_rel_hybrid_raw <- rel_data_hybrid_labelled %>%
  ggplot(aes(x=p_traced_manual, y=effective_r0_rel,
             colour = label,
             linetype=contact_limit_manual,
             shape=contact_limit_manual)) %>%
  linetype_window(label = label_windows_manual) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)
reff_rel_hybrid <- reff_rel_hybrid_raw +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "black", size = 1) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = expression(paste("Δ",italic(R)[eff], " (hybrid – manual)")),
                     breaks = seq(-10,10,0.2), limits = c(-0.7, 0.05)) +
  scale_colour_manual(name = NULL, values = c("#1B9E77", "#7570B3", "#D95F02", "#E7298A")) +
  guides(colour=guide_legend(ncol = 1, order = 1)) + theme(
    legend.background = element_rect(fill = alpha("white", 0.8))
  ) +
  theme(aspect.ratio = 1)

#------------------------------------------------------------------------------
# R_eff contour plots
#------------------------------------------------------------------------------

cat("\nGenerating R_eff contour plots...")

# Palettes
bidir_palette_dense <- uptake_data %>%
  filter(backtrace_distance == Inf, contact_limit_manual == 6) %>%
  make_contour_palette_reff(col_bidir, step_dense, ., Inf)
bidir_palette_sparse <- uptake_data %>%
  filter(backtrace_distance == Inf, contact_limit_manual == 6) %>%
  make_contour_palette_reff(col_bidir, step_sparse, ., Inf)

# Digital bidirectional
reff_contour_digital <- uptake_data %>%
  filter(p_traced_manual == 0, backtrace_distance == Inf,
         contact_limit_manual == 6) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_sparse, breaks = seq(0,10,step_sparse)) +
  geom_vline(xintercept = c(0.53, 0.8), linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))

reff_contour_hybrid_6day <- uptake_data %>%
  filter(p_traced_manual != 0, backtrace_distance == Inf,
         contact_limit_manual == 6) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_sparse, breaks = seq(0,10,step_sparse)) +
  geom_vline(xintercept = c(0.53, 0.8), linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))

reff_contour_hybrid_2day <- uptake_data %>%
  filter(p_traced_manual != 0, backtrace_distance == Inf,
         contact_limit_manual == 2) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_sparse, breaks = seq(0,10,step_sparse)) +
  geom_vline(xintercept = c(0.53, 0.8), linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))

# Add labels
contour_label <- expression(paste(italic(R)[eff], "  given 90% trace success"))#"Effective reprod. number, given\n90% prob. of trace success"
label_contour <- function(contour_plot, size = fontsize_base * 5/14,
                          label=contour_label){
  contour_plot + annotate("label", x=0.02, y=0.02, label=label,
                          hjust = 0, vjust = 0, size = size,
                          label.size=NA, fill=alpha("white", 0.8))
}
reff_contour_digital_labelled <- label_contour(reff_contour_digital)
reff_contour_hybrid_6day_labelled <- label_contour(reff_contour_hybrid_6day)
reff_contour_hybrid_2day_labelled <- label_contour(reff_contour_hybrid_2day)

cat("done.\n")

#------------------------------------------------------------------------------
# Window contour plot
#------------------------------------------------------------------------------

cat("\nGenerating trace-window contour plots...")

window_contour <- window_data %>%
  contour_plot_reff_abs("p_traced_manual", "contact_limit_manual",
             palette = bidir_palette_sparse,
             breaks =  seq(0,10,step_sparse)) +
  scale_y_continuous(name = "Manual tracing window (days)", limits = c(0,10),
                     breaks = seq(0,10,2)) +
  scale_x_continuous(name = "Probability of trace\nsuccess (%)",
                     limits = c(0,1), breaks = seq(0,1,0.2),
                     labels = label_pc) +
  theme(aspect.ratio = 1)

window_contour_labelled <- window_contour +
  annotate("label", x=0.02, y=0.2, label=expression(italic(R)[eff]),
           hjust = 0, vjust = 0, size = fontsize_base * 5/14,
           label.size=NA, fill=alpha("white", 0.8))

cat("done.\n")

#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

cat("\nAssembling main figure...\n")

grid_main <- plot_grid(reff_manual + ggtitle("Manual tracing only"),
                       window_contour_labelled + 
                         ggtitle("Effect of tracing window\n (bidir. manual tracing)"),
                       reff_contour_digital_labelled +
                         ggtitle("Effect of digital uptake\n(bidir. digital tracing)"),
                       reff_hybrid + ggtitle("Manual + digital\n(hybrid) tracing"),
                       # reff_rel +
                       #   ggtitle("Relative benefit of\nhybrid tracing"),
                       reff_contour_hybrid_2day_labelled + 
                         ggtitle("Effect of digital uptake\n(bidir. hybrid tracing,\n2-day manual window)"),
                       reff_contour_hybrid_6day_labelled +
                         ggtitle("Effect of digital uptake\n(bidir. hybrid tracing,\n6-day manual window)"),
                       labels = "auto", nrow = 2, ncol = 3, align = "hv",
                       axis = "l", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

cat("\n...done.\n")

#==============================================================================
# Supplementary versions
#==============================================================================

#------------------------------------------------------------------------------
# Single-strategy R_eff plot for universal digital tracing
#------------------------------------------------------------------------------

cat("\nGenerating universal digital plot...")

si_reff_univ <- digital_data %>%
  filter(p_smartphone_overall == 1.0, p_data_sharing_auto == 1.0) %>%
  ggplot(aes(x=p_traced_auto, y=effective_r0_mean,
             colour=backtrace_distance,
             linetype=contact_limit_auto,
             shape=contact_limit_auto)) %>%
  format_reff(r0=2.5) %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

cat("done.\n")

#------------------------------------------------------------------------------
# Generalised contour plots
#------------------------------------------------------------------------------

cat("\nGenerating supplementary contour plots...")

# Format data
cat("\n\tPreparing data...")
si_uptake_data <- uptake_data %>%
  mutate(row_lab = paste(bt_type, "tracing"),
         col_lab = ifelse(trace_type == "Digital only", "Digital tracing",
                          paste0("Hybrid tracing\n(", trace_window_manual, ")")))

# Define missing palettes
cat("\n\tPreparing palettes...")
bidir_palette_dense_ctrl <- make_contour_palette_ctrl(col_bidir, step_dense)
bidir_palette_sparse_ctrl <- make_contour_palette_ctrl(col_bidir, step_sparse)
fwd_palette_dense_ctrl <- make_contour_palette_ctrl(col_fwd, step_dense)
fwd_palette_sparse_ctrl <- make_contour_palette_ctrl(col_fwd, step_sparse)
fwd_palette_dense <- uptake_data %>%
  filter(backtrace_distance == 0, contact_limit_manual == 6) %>%
  make_contour_palette_reff(col_fwd, step_dense, ., 0)
fwd_palette_sparse <- uptake_data %>%
  filter(backtrace_distance == 0, contact_limit_manual == 6) %>%
  make_contour_palette_reff(col_fwd, step_sparse, ., 0)


# R_eff SI plot
cat("\n\tMaking R_eff contour plots...")
si_uptake_reff_bidir <- si_uptake_data %>% filter(bt_type == "Bidirectional") %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_sparse, breaks = seq(0,10,step_sparse)) +
  geom_vline(xintercept = c(0.53, 0.8), linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  scale_x_continuous(name = "% of cases with chirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  facet_grid(row_lab~col_lab) +
  theme_base + theme(legend.position="none")
si_uptake_reff_fwd <- si_uptake_data %>% filter(bt_type == "Forward-only") %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_sparse, breaks = seq(0,10,step_sparse)) +
  geom_vline(xintercept = c(0.53, 0.8), linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  scale_x_continuous(name = "% of cases with chirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  facet_grid(row_lab~col_lab) +
  theme_base + theme(legend.position="none")
si_uptake_reff <- plot_grid(si_uptake_reff_fwd %>% strip_upper(),
                            si_uptake_reff_bidir %>% strip_lower(),
                            nrow = 2, align = "v", axis= "l")

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

cat("\n\tMain figure...")
cowplot::save_plot(filename=main_path, plot=grid_main,
                   ncol = 3, nrow = 2.1, base_height = plot_scale_main_in,
                   base_asp = 0.8)

cat("\n\tSI R_eff contour plots...")
cowplot::save_plot(filename=si_uptake_reff_path, plot=si_uptake_reff,
                   ncol = 3, nrow = 2, base_height = plot_scale_si_in,
                   base_asp = 0.83)

cat("\n\tDigital R_eff plots...")
cowplot::save_plot(filename=si_digital_univ_path, plot=si_reff_univ,
                   ncol = 1.2, nrow = 1.2, base_height = plot_scale_si_in,
                   base_asp = 0.83)

cat("\n\tHybrid/manual relative R_eff plots...")
cowplot::save_plot(filename=si_hybrid_delta_path, plot=reff_rel_hybrid,
                   ncol = 1.2, nrow = 1.2, base_height = plot_scale_si_in,
                   base_asp = 0.83)


cat("done.\n")