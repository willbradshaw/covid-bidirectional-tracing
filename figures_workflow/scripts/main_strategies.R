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
cat("done.\n")

#==============================================================================
# Read in data
#==============================================================================

cat("\nReading in data...")

# Manual tracing (with and without automated)
manual_path <- snakemake@input[["manual"]]
#manual_path <- "data/main_strategies_manual_1k_scenario.tsv.gz"
manual_data <- suppressMessages(read_tsv(manual_path)) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(0,Inf)),
         contact_limit_manual = factor(contact_limit_manual, levels=c(6,2)),
         trace_type = ifelse(p_traced_auto == 0, "Manual tracing only", 
                             "Manual + digital tracing"),
         p_traced = p_traced_manual)

# Digital tracing alone (universal or high coverage)
digital_path <- snakemake@input[["digital"]]
#digital_path <- "data/main_strategies_digital_1k_scenario.tsv.gz"
digital_data <- suppressMessages(read_tsv(digital_path)) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(0,Inf)),
         contact_limit_auto = factor(contact_limit_auto, levels=c(Inf,2)),
         trace_type = "Digital tracing only",
         p_traced = p_traced_auto)

# Uptake contour plots
uptake_path <- snakemake@input[["uptake"]]
#uptake_path <- "data/main_strategies_uptake_1k_scenario.tsv.gz"
uptake_data <- suppressMessages(read_tsv(uptake_path)) %>% 
  filter(p_traced_manual == 0.9 | contact_limit_manual == 6) %>%
  prep_contour_data(0.9, "p_smartphone_overall", "p_data_sharing_auto",
                    c("p_traced_manual", "backtrace_distance", "contact_limit_manual"))

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
  format_reff(r0=2.5) %>% colour_backtrace %>% linetype_window %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

cat("\n\t2")

reff_hybrid <- manual_data %>%
  filter(p_traced_auto != 0, p_smartphone_overall == 0.8) %>%
  ggplot(aes(x=p_traced_manual, y=effective_r0_mean,
             colour=backtrace_distance, 
             linetype=contact_limit_manual,
             shape=contact_limit_manual)) %>%
  format_reff(r0=2.5) %>% colour_backtrace %>% 
  linetype_window(label = label_limits_manual) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 0.01), yjust = 0, xjust = 0, vspace=vspace)

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
# Comparative control plot
#------------------------------------------------------------------------------

cat("\nGenerating comparative control plot...")

cat("\n\tAssembling data...")
digital_data_prep <- digital_data %>% filter(p_data_sharing_auto == 0.9,
                                             contact_limit_auto == Inf) %>%
  select(-contact_limit_auto)
manual_data_prep <- manual_data %>% filter(contact_limit_manual == 6) %>%
  select(-contact_limit_manual)
comp_data <- bind_rows(manual_data_prep, digital_data_prep) %>%
  filter(p_smartphone_overall == 0.8) %>%
  mutate(backtrace_distance = factor(backtrace_distance, levels=c(Inf, 0)),
         trace_type = factor(trace_type, 
                             levels=c("Manual tracing only", "Manual + digital tracing",
                                      "Digital tracing only")))

cat("\n\tMaking plot...")
ctrl_comp <- ggplot(comp_data,
                    aes(x=p_traced, y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                        colour=trace_type,
                        linetype=backtrace_distance,
                        shape=backtrace_distance)) %>%
  format_ctrl(ylimits = c(0,0.5), coord_ratio = 2) %>% 
  colour_ttype(ncol=1, name=NULL) %>% 
  linetype_backtrace(label=label_backtrace) %>% x_ptrace %>%
  theme_internal(pos = c(0.02, 1), yjust = 1, xjust = 0, vspace=-1)

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
label_contour <- function(contour_plot, size = 3.8, label=contour_label){
  contour_plot + annotate("label", x=0.02, y=0.02, label=label,
                          hjust = 0, vjust = 0, size = size,
                          label.size=NA, fill=alpha("white", 0.8))
}
reff_contour_digital_labelled <- label_contour(reff_contour_digital)
reff_contour_hybrid_labelled <- label_contour(reff_contour_hybrid)

cat("done.\n")

#------------------------------------------------------------------------------
# Grid plot
#------------------------------------------------------------------------------

cat("\nAssembling main figure...\n")

grid_main <- plot_grid(reff_manual + ggtitle("Manual tracing only"),
                       reff_univ + ggtitle("Digital tracing only\n(universal coverage)"),
                       reff_contour_digital_labelled +
                         ggtitle("Bidirectional digital\ntracing (partial coverage)"),
                       reff_hybrid + ggtitle("Manual + digital\n(hybrid) tracing"),
                       ctrl_comp + ggtitle("% outbreaks controlled"),
                       reff_contour_hybrid_labelled + 
                         ggtitle("Bidirectional hybrid\ntracing (partial coverage)"),
                       labels = "auto", nrow = 2, ncol = 3, align = "hv",
                       axis = "l", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

cat("\n...done.\n")

#==============================================================================
# Supplementary versions
#==============================================================================

#------------------------------------------------------------------------------
# Single-strategy control-rate plots
#------------------------------------------------------------------------------

cat("\nGenerating single-strategy control-rate plots...")

cat("\n\tPreparing data...")
si_ctrl_data_digital <- digital_data %>%
  filter(p_smartphone_overall == 1.0, p_data_sharing_auto == 1.0) %>%
  mutate(trace_type = "Digital tracing only\n(universal coverage)",
         contact_limit = as.character(contact_limit_auto)) %>%
  select(-contact_limit_auto)

si_ctrl_data_manual <- manual_data %>%
  filter(p_smartphone_overall == 0.8) %>%
  mutate(trace_type = ifelse(trace_type == "Manual + digital tracing",
                             "Manual + digital\n(hybrid) tracing",
                             trace_type),
         contact_limit = as.character(contact_limit_manual)) %>%
  select(-contact_limit_manual)

si_ctrl_data <- bind_rows(si_ctrl_data_digital, si_ctrl_data_manual) %>%
  mutate(contact_limit = factor(contact_limit, levels=c("Inf", "6", "2")))

cat("\n\tMaking plot...")
si_ctrl <- ggplot(si_ctrl_data,
                    aes(x=p_traced, y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                        colour=backtrace_distance, 
                        linetype=contact_limit,
                        shape=contact_limit)) %>%
  format_ctrl(ylimits = c(0,0.5), coord_ratio = 2) %>% colour_backtrace %>% 
  linetype_window %>% x_ptrace %>%
  theme_external()

si_ctrl <- si_ctrl + facet_grid(.~trace_type)

cat("\n...done.\n")

#------------------------------------------------------------------------------
# Generalised contour plots
#------------------------------------------------------------------------------

cat("\nGenerating supplementary contour plots...")

# Format data
cat("\n\tPreparing data...")
si_uptake_data <- uptake_data %>%
  mutate(row_lab = paste(bt_type, "tracing,\n"),
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
                        palette = bidir_palette_dense, breaks = seq(0,10,step_dense)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  facet_grid(row_lab~col_lab) +
  theme_base + theme(legend.position="none")
si_uptake_reff_fwd <- si_uptake_data %>% filter(bt_type == "Forward-only") %>%
  contour_plot_reff_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_dense, breaks = seq(0,10,step_dense)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  facet_grid(row_lab~col_lab) +
  theme_base + theme(legend.position="none")
si_uptake_reff <- plot_grid(si_uptake_reff_fwd %>% strip_upper(),
                            si_uptake_reff_bidir %>% strip_lower(),
                            nrow = 2, align = "v", axis= "l")

# Ctrl SI plot
cat("\n\tMaking control-rate contour plots...")
si_uptake_ctrl_bidir <- si_uptake_data %>% filter(bt_type == "Bidirectional") %>%
  contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = bidir_palette_dense_ctrl, 
                        breaks = seq(0,100,step_dense*100)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  facet_grid(row_lab~col_lab) +
  theme_base + theme(legend.position="none")
si_uptake_ctrl_fwd <- si_uptake_data %>% filter(bt_type == "Forward-only") %>%
  contour_plot_ctrl_abs("p_smartphone_overall", "p_data_sharing_auto",
                        palette = fwd_palette_dense_ctrl, 
                        breaks = seq(0,100,step_dense*100)) +
  scale_x_continuous(name = "% of cases with\nchirping smartphones",
                     labels = label_pc, breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "% of cases sharing data", labels = label_pc,
                     breaks = seq(0,1,0.2)) +
  facet_grid(row_lab~col_lab) +
  theme_base + theme(legend.position="none")
si_uptake_ctrl <- plot_grid(si_uptake_ctrl_fwd %>% strip_upper(),
                            si_uptake_ctrl_bidir %>% strip_lower(),
                            nrow = 2, align = "v", axis = "l")

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

plot_scale_cm <- snakemake@params[["panel_scale"]]
#plot_scale_cm <- 13
plot_scale_in <- plot_scale_cm/2.54

cat("\n\tMain figure...")
main_path <- snakemake@output[["main"]]
#main_path <- "output_files/dev_main_strategies.png"
cowplot::save_plot(filename=main_path, plot=grid_main,
                   ncol = 3, nrow = 2, base_height = plot_scale_in,
                   base_asp = 0.8)

cat("\n\tSI R_eff contour plots...")
si_uptake_reff_path <- snakemake@output[["si_uptake_reff"]]
#si_uptake_reff_path <- "output_files/dev_si_uptake_reff.png"
cowplot::save_plot(filename=si_uptake_reff_path, plot=si_uptake_reff,
                   ncol = 3, nrow = 2, base_height = plot_scale_in,
                   base_asp = 0.83)

cat("\n\tSI control-rate contour plots...")
si_uptake_ctrl_path <- snakemake@output[["si_uptake_ctrl"]]
#si_uptake_ctrl_path <- "output_files/dev_si_uptake_ctrl.png"
cowplot::save_plot(filename=si_uptake_ctrl_path, plot=si_uptake_ctrl,
                   ncol = 3, nrow = 2, base_height = plot_scale_in,
                   base_asp = 0.83)

cat("\n\tSI control-rate line plots...")
si_ptrace_ctrl_path <- snakemake@output[["si_ptrace_ctrl"]]
#si_ptrace_ctrl_path <- "output_files/dev_si_ptrace_ctrl.png"
cowplot::save_plot(filename=si_ptrace_ctrl_path, plot=si_ctrl,
                   ncol = 3, nrow = 1, base_height = plot_scale_in,
                   base_asp = 0.65)

cat("done.\n")