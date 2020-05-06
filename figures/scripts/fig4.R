library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Paths
optim_path <- "figures/data/fig4_optimistic_1k_scenario.tsv.gz"
pessim_path <- "figures/data/fig4_pessimistic_1k_scenario.tsv.gz"
median_path <- "figures/data/fig4_baseline_1k_scenario.tsv.gz"

# Data
optim_data <- suppressMessages(read_tsv(optim_path))
pessim_data <- suppressMessages(read_tsv(pessim_path))
median_data <- suppressMessages(read_tsv(median_path))

#==============================================================================
# Process data
#==============================================================================

process_data_fig4 <- function(data, p_traced, contact_limit, trace_neg_sym,
                              r0 = 2.5, uptake_high = 0.8, uptake_low = 0.53){
  # Basic filtering
  data_filtered <- data %>%
    filter(p_traced_manual %in% c(0, p_traced),
           p_traced_auto %in% c(0, p_traced),
           contact_limit_manual == 7, trace_neg_symptomatic == TRUE)
  if (!is.null(r0)) data_filtered <- filter(data_filtered, r0_base == r0)
  # Phone usage rate doesn't matter if no digital tracing
  data_filtered <- filter(data_filtered, 
                          !(p_traced_auto == 0 & p_smartphone_overall == uptake_low))
  # Label scenarios
  data_labelled <- data_filtered %>% mutate(
    trace_type = ifelse(p_traced_auto == 0,
                        ifelse(p_traced_manual == 0, "No tracing", "Manual tracing"),
                        ifelse(p_traced_manual == 0, "Digital tracing", "Combined tracing")),
    uptake = ifelse(p_smartphone_overall == uptake_high, "high uptake", "low uptake"),
    label = ifelse(p_traced_auto == 0, trace_type,
                   paste(trace_type, uptake, sep=", "))

  )
  return(data_labelled)
}

median_data_processed <- process_data_fig4(median_data, 0.9, 7, TRUE)
optim_data_processed <- process_data_fig4(optim_data, 0.9, 7, TRUE)
pessim_data_processed <- process_data_fig4(pessim_data, 0.9, 7, TRUE)

#==============================================================================
# Plotting themes and functions
#==============================================================================

# Fonts
font <- "sans" # Main figure font
titlefont <- font # Font for axis titles etc. (if different)
fontsize_base <- 12 # Basic figure font size
fontscale_title <- 1.5 # Default axis-title scale relative to regular font
fontscale_mid <- 1.25
fontscale_main <- 1.5 # Default plot-title scale
fontscale_label <- 2 # Default subfigure-label scale (A, B, etc)
fontscale_legend <- 1

# Theme
theme_base <-   theme_bw() + theme(
  legend.position = "bottom",
  axis.text = element_text(size = fontsize_base, family = font, 
                           colour = "black"),
  axis.title.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(r=3,unit="mm")),
  axis.title.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t=5,unit="mm")),
  legend.text = element_text(size = fontsize_base * fontscale_legend, 
                             family = font, colour = "black"),
  legend.title = element_text(size = fontsize_base * fontscale_legend,
                              family = font, colour = "black",
                              face = "bold", vjust=0.5, 
                              margin=margin(r=3, unit="mm")),
  plot.title = element_text(size = fontsize_base * fontscale_main,
                            family = titlefont, colour="black",
                            hjust = 0.5, vjust = 0.5,
                            margin = margin(b=4, unit="mm"),
                            face = "bold"),
  plot.margin = margin(t=1.5, l=0.5, r=0.5, b = 0.5, unit="cm"),
  strip.background = element_blank(),
  strip.text.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t = 3, b = 3, unit="mm")),
  strip.text.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(l=3, r=3, unit="mm")),
  legend.justification = "center",
  panel.border = element_blank(),
  axis.line = element_line(size=0.5, colour="black")
)

label_backtrace <- function(x){
  ifelse(x==0, "Forward tracing only",
         "Forward and reverse tracing")
}

label_levels <- c("No tracing", "Manual tracing",
                  "Digital tracing, low uptake",
                  "Digital tracing, high uptake",
                  "Combined tracing, low uptake",
                  "Combined tracing, high uptake")

make_barplot_fig4_ctrl <- function(data, 
                                   err_width = 0.2, bar_width=0.6,
                                   legend_pos = c(0.05,0.95), legend_vspace=0){
  g <- ggplot(data, aes(x=factor(label, levels=label_levels), y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                        fill=factor(backtrace_distance))) +
    geom_col(position="dodge", width=bar_width) + 
    geom_errorbar(width=err_width, position = position_dodge(width=bar_width)) + 
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_fill_brewer(palette = "Dark2", name=NULL, labels=label_backtrace) +
    guides(fill=guide_legend(ncol=1,order=1)) +
    coord_fixed(ratio=5) +
    theme_base + theme(
      axis.text.x = element_text(angle=45,hjust=1,vjust=1,
                                 size = fontsize_base * fontscale_mid),
      axis.title.x = element_blank(),
      legend.justification = c("left","top"),
      legend.background = element_rect(fill=alpha("white", 0.5)),
      legend.key = element_rect(fill=alpha("white", 0)),
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.position = legend_pos,
      legend.spacing.y = unit(legend_vspace, "mm"),
    )
  return(g)
}

make_barplot_fig4_reff <- function(data, r0=2.5, bar_width=0.6,
                                   legend_pos = c(0,1), legend_vspace=0,
                                   ylim_max = 3.5,
                                   r0_lab_x = 5, r0_lab_y = 2.35,
                                   r0_lab_size = 4.6, r0_col = "red"){
  g <- ggplot(data, aes(x=factor(label, levels=label_levels), y=effective_r0_mean,
                        fill=factor(backtrace_distance))) +
    geom_col(position="dodge", width=bar_width) + 
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    geom_hline(yintercept = r0, linetype = "dotted", size=1,
               colour = "red") +
    annotate("text", x=r0_lab_x, y=r0_lab_y, hjust=0.5, vjust = 0.5,
             colour = r0_col, size=r0_lab_size, label=expression(R[0])) +
    scale_y_continuous(name = "Average effective\nreproduction number", limits = c(0,ylim_max),
                       breaks = seq(0,4,0.5)) +
    scale_fill_brewer(palette = "Dark2", name=NULL, labels=label_backtrace) +
    guides(fill=guide_legend(ncol=1,order=1)) +
    coord_fixed(ratio=5/4) +
    theme_base + theme(
      axis.text.x = element_text(angle=45,hjust=1,vjust=1,
                                 size = fontsize_base * fontscale_mid),
      axis.title.x = element_blank(),
      legend.justification = c("left","top"),
      legend.background = element_rect(fill=alpha("white", 0.5)),
      legend.key = element_rect(fill=alpha("white", 0)),
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.position = legend_pos,
      legend.spacing.y = unit(legend_vspace, "mm"),
    )
  return(g)
}

make_ctrl_grid <- function(median_data, optim_data, pessim_data,
                           p_traced, contact_limit, trace_neg_sym,
                           r0 = 2.5, uptake_high = 0.8, uptake_low = 0.53,
                           err_width = 0.2, bar_width=0.6,
                           legend_pos = c(0.05,0.95), legend_vspace=0){
  data <- list(median_data, optim_data, pessim_data)
  data_processed <- lapply(data, function(d)
    process_data_fig4(d, p_traced, contact_limit, trace_neg_sym,
                      r0, uptake_high, uptake_low))
  titles <- c("Median scenario", "Optimistic scenario", "Pessimistic scenario")
  plots <- lapply(data_processed, function(d)
    make_barplot_fig4_ctrl(d, err_width, bar_width, legend_pos, legend_vspace))
  plots <- lapply(1:length(plots), function(n) plots[[n]] + ggtitle(titles[n]))
  plot_out <- plot_grid(plots[[1]], plots[[2]], plots[[3]],
                        labels = "auto", nrow = 1, ncol = 3, align = "hv",
                        axis = "l", label_size = fontsize_base * fontscale_label,
                        label_fontfamily = titlefont, label_colour = "black")
  return(plot_out)
}

make_reff_grid <- function(median_data, optim_data, pessim_data,
                           p_traced, contact_limit, trace_neg_sym,
                           r0 = 2.5, uptake_high = 0.8, uptake_low = 0.53,
                           bar_width=0.6, ylim_max = 3.5,
                           legend_pos = c(0, 1), legend_vspace=0,
                           r0_lab_x = 6, r0_lab_y = 2.3,
                           r0_lab_size = 4.6, r0_col = "red"){
  data <- list(median_data, optim_data, pessim_data)
  data_processed <- lapply(data, function(d)
    process_data_fig4(d, p_traced, contact_limit, trace_neg_sym,
                      r0, uptake_high, uptake_low))
  titles <- c("Median scenario", "Optimistic scenario", "Pessimistic scenario")
  plots <- lapply(data_processed, function(d)
    make_barplot_fig4_reff(d, r0, bar_width, legend_pos, legend_vspace,
                           r0_lab_x = r0_lab_x, r0_lab_y = r0_lab_y,
                           r0_lab_size = r0_lab_size, r0_col = r0_col,
                           ylim_max = ylim_max))
  plots <- lapply(1:length(plots), function(n) plots[[n]] + ggtitle(titles[n]))
  plot_out <- plot_grid(plots[[1]], plots[[2]], plots[[3]],
                        labels = "auto", nrow = 1, ncol = 3, align = "hv",
                        axis = "l", label_size = fontsize_base * fontscale_label,
                        label_fontfamily = titlefont, label_colour = "black")
  return(plot_out)
}


#==============================================================================
# Make plots
#==============================================================================

# R0 2.5, 90% contacts traced, 7d manual window, trace before testing
ctrl_7d_090_pre_r025 <- make_ctrl_grid(median_data, optim_data, pessim_data,
                                       0.9, 7, TRUE, 2.5)
reff_7d_090_pre_r025 <- make_reff_grid(median_data, optim_data, pessim_data,
                                       0.9, 7, TRUE, 2.5)

# R0 2.0, 90% contacts traced, 7d manual window, trace before testing
ctrl_7d_090_pre_r020 <- make_ctrl_grid(median_data, optim_data, pessim_data,
                                       0.9, 7, TRUE, 2.0)
reff_7d_090_pre_r020 <- make_reff_grid(median_data, optim_data, pessim_data,
                                       0.9, 7, TRUE, 2.0)


# 
# # R0 = 2.5, 90% contacts traced, 7d manual window, trace before testing
# median_plot_7d_090_pre_r025 <- median_data %>%
#   process_data_fig4(0.9, 7, TRUE) %>%
#   make_barplot_fig4() + ggtitle("Median parameters")
# optim_plot_7d_090_pre_r025 <- optim_data %>%
#   process_data_fig4(0.9, 7, TRUE) %>%
#   make_barplot_fig4() + ggtitle("Optimistic parameters")
# pessim_plot_7d_090_pre_r025 <- pessim_data %>%
#   process_data_fig4(0.9, 7, TRUE) %>%
#   make_barplot_fig4() + ggtitle("Pessimistic parameters")
# plot_7d_090_pre_r025 <- plot_grid(median_plot_7d_090_pre_r025,
#                                   optim_plot_7d_090_pre_r025,
#                                   pessim_plot_7d_090_pre_r025,
# 
# 
# # R0 = 2.0, 90% contacts traced, 7d manual window, trace before testing
# median_plot_7d_090_pre_r020 <- median_data %>%
#   process_data_fig4(0.9, 7, TRUE, r0=2.0) %>%
#   make_barplot_fig4()
# optim_plot_7d_090_pre_r020 <- optim_data %>%
#   process_data_fig4(0.9, 7, TRUE, r0=2.0) %>%
#   make_barplot_fig4()
# pessim_plot_7d_090_pre_r020 <- pessim_data %>%
#   process_data_fig4(0.9, 7, TRUE, r0=2.0) %>%
#   make_barplot_fig4()
# 

#==============================================================================
# Save outputs
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio, device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

# Main figure (R0 = 2.5 only)
row_height <- 16
path_prefix <- "figures/img/fig4_7d_090_pre_r025"
save_fig(path_prefix, "_ctrl.png", ctrl_7d_090_pre_r025, row_height, 3*0.6)
save_fig(path_prefix, "_reff.png", reff_7d_090_pre_r025, row_height, 3*0.6)
save_fig(path_prefix, "_ctrl.svg", ctrl_7d_090_pre_r025, row_height, 3*0.6, device="svg")
save_fig(path_prefix, "_reff.svg", reff_7d_090_pre_r025, row_height, 3*0.6, device="svg")


# R0 = 2.0
row_height <- 16
path_prefix <- "figures/img/fig4_7d_090_pre_r020"
save_fig(path_prefix, "_ctrl.png", ctrl_7d_090_pre_r020, row_height, 3*0.6)
save_fig(path_prefix, "_reff.png", reff_7d_090_pre_r020, row_height, 3*0.6)
