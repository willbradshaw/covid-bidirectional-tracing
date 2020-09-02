#==============================================================================
# Preamble
#==============================================================================

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
  axis.title.x = element_text(margin = margin(t = 0.2, b = 0, unit = "cm"))
)

# Input paths
#manual_path <- snakemake@input[["manual"]]
#digital_path <- snakemake@input[["digital"]]
#uptake_path <- snakemake@input[["uptake"]]
#window_path <- snakemake@input[["window"]]
manual_path <- "data/main_strategies_manual_1k_scenario.tsv.gz"
digital_path <- "data/main_strategies_digital_1k_scenario.tsv.gz"
uptake_path <- "data/main_strategies_uptake_1k_scenario.tsv.gz"
window_path <- "data/si_window_median_1k_scenario.tsv.gz"


# Output paths & parameters
#plot_scale_cm <- snakemake@params[["panel_scale"]]
#main_path <- snakemake@output[["main"]]
#si_uptake_reff_path <- snakemake@output[["si_uptake_reff"]]
plot_scale_cm <- 13
plot_scale_in <- plot_scale_cm/2.54
main_path <- "output_files/dev_main_strategies.png"
si_uptake_reff_path <- "output_files/dev_si_uptake_reff.png"

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
  filter(p_traced_auto != 0, p_smartphone_overall == 0.53,
         backtrace_distance == Inf) %>%
  index_data(key_x = "contact_limit_manual", key_y = "p_traced_manual") %>%
  smooth_data(group_vars = c())

cat("done.\n")

#==============================================================================
# Make main-text figure
#==============================================================================
;
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
  format_reff(r0=2.5) %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
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
  linetype_window(label = label_limits_manual) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)
reff_hybrid <- reff_hybrid_raw + 
  scale_colour_manual(name = NULL, values = c("#1B9E77", "#7570B3", "#D95F02", "#E7298A")) +
  guides(colour=guide_legend(ncol = 1, order = 1)) + theme(
    legend.background = element_rect(fill = alpha("white", 0.8))
  )

cat("\n\t3")

reff_univ <- digital_data %>%
  filter(p_smartphone_overall == 1.0, p_data_sharing_auto == 1.0) %>%
  ggplot(aes(x=p_traced_auto, y=effective_r0_mean,
             colour=backtrace_distance,
             linetype=contact_limit_auto,
             shape=contact_limit_auto)) %>%
  format_reff(r0=2.5) %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

cat("\n...done.\n")

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
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))

reff_contour_hybrid <- uptake_data %>%
  filter(p_traced_manual != 0, backtrace_distance == Inf,
         contact_limit_manual == 6) %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_dense, breaks = seq(0,10,step_dense)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2))

# Add labels
contour_label <- "Effective reprod. number, given\n90% prob. of trace success"
label_contour <- function(contour_plot, size = fontsize_base * 5/14,
                          label=contour_label){
  contour_plot + annotate("label", x=0.02, y=0.98, label=label,
                          hjust = 0, vjust = 1, size = size,
                          label.size=NA, fill=alpha("white", 0.8))
}
reff_contour_digital_labelled <- label_contour(reff_contour_digital)
reff_contour_hybrid_labelled <- label_contour(reff_contour_hybrid)

cat("done.\n")

#------------------------------------------------------------------------------
# Window contour plot
#------------------------------------------------------------------------------

cat("\nGenerating trace-window contour plots...")

window_contour <- window_data %>%
  contour_plot_reff_abs("p_traced_manual", "contact_limit_manual",
             palette = bidir_palette_sparse,
             breaks =  seq(0,10,step_sparse)) +
  scale_y_continuous(name = "Manual trace limit (days)", limits = c(0,10),
                     breaks = seq(0,10,2)) +
  scale_x_continuous(name = "Probability of trace\nsuccess (%)",
                     limits = c(0,1), breaks = seq(0,1,0.2),
                     labels = label_pc) +
  theme(aspect.ratio = 1)

label_window <- "Effective reprod. number, given\n53% smartphone coverage"
window_contour_labelled <- window_contour +
  annotate("label", x=0.02, y=9.8, label=label_window,
           hjust = 0, vjust = 1, size = fontsize_base * 5/14,
           label.size=NA, fill=alpha("white", 0.8))

cat("done.\n")

#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

cat("\nAssembling main figure...\n")

grid_main <- plot_grid(reff_manual + ggtitle("Manual tracing only"),
                       reff_univ + ggtitle("Digital tracing only\n(universal uptake)"),
                       reff_contour_digital_labelled +
                         ggtitle("Effect of digital uptake\n(bidir. digital tracing)"),
                       reff_hybrid + ggtitle("Manual + digital\n(hybrid) tracing"),
                       # ctrl_comp + ggtitle("% outbreaks controlled"),
                       window_contour_labelled + 
                         ggtitle("Effect of tracing window\n (bidir. hybrid tracing)"),
                       reff_contour_hybrid_labelled + 
                         ggtitle("Effect of digital uptake\n(bidir. hybrid tracing)"),
                       labels = "auto", nrow = 2, ncol = 3, align = "hv",
                       axis = "l", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

cat("\n...done.\n")

#==============================================================================
# Supplementary versions
#==============================================================================

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
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  facet_grid(row_lab~col_lab) +
  theme_base + theme(legend.position="none")
si_uptake_reff_fwd <- si_uptake_data %>% filter(bt_type == "Forward-only") %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_sparse, breaks = seq(0,10,step_sparse)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
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
                   ncol = 3, nrow = 2, base_height = plot_scale_in,
                   base_asp = 0.8)

cat("\n\tSI R_eff contour plots...")
cowplot::save_plot(filename=si_uptake_reff_path, plot=si_uptake_reff,
                   ncol = 3, nrow = 2, base_height = plot_scale_in,
                   base_asp = 0.83)

cat("done.\n")