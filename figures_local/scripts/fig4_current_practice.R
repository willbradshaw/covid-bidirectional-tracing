library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Paths
optim_path <- "figures_local/data/fig4_optimistic_variants_1k_scenario.tsv.gz"
pessim_path <- "figures_local/data/fig4_pessimistic_variants_1k_scenario.tsv.gz"
median_path <- "figures_local/data/fig4_baseline_variants_1k_scenario.tsv.gz"

# Data
optim_data <- suppressMessages(read_tsv(optim_path))
pessim_data <- suppressMessages(read_tsv(pessim_path))
median_data <- suppressMessages(read_tsv(median_path))

#==============================================================================
# Process data
#==============================================================================

uptake_high <- 0.80
uptake_low  <- 0.53

comb_data_all <- function(median_data, optim_data, pessim_data){
  bind_rows(optim_data %>% mutate(scenario_type = "Optimistic scenario"),
            median_data %>% mutate(scenario_type = "Median scenario"),
            pessim_data %>% mutate(scenario_type = "Pessimistic scenario")) %>%
    filter((p_traced_auto == 0 | p_traced_manual == 0 | p_traced_auto == p_traced_manual)) %>%
    mutate(scenario_type = factor(scenario_type, 
                                  levels = c("Median scenario", "Optimistic scenario",
                                             "Pessimistic scenario")),
           trace_type = ifelse(p_traced_auto == 0,
                               ifelse(p_traced_manual == 0, "No tracing", "Manual"),
                               ifelse(p_traced_manual == 0, "Digital", "Hybrid")),
           uptake = ifelse(p_traced_auto == 0, "",
                           ifelse(p_smartphone_overall == uptake_high, ", high-uptake", ", low-uptake")),
           trace_window = ifelse(p_traced_manual == 0, "", paste0(", ", contact_limit_manual, "-day")),
           label = paste0(trace_type, uptake, trace_window)) %>%
    # Filter by uptake for non-digital cases and trace limit for non-manual cases
    filter((p_traced_auto != 0 | p_smartphone_overall == uptake_low),
           (p_traced_manual != 0 | contact_limit_manual == 2))
}

filter_comb <- function(data, p_traced, r0, trace_neg_sym){
  filter(data, p_traced_manual %in% c(-1, 0, p_traced),
         p_traced_auto %in% c(-1, 0, p_traced),
         r0_base == r0, trace_neg_symptomatic == trace_neg_sym)
}

# Combine datasets
comb_data <- comb_data_all(median_data, optim_data, pessim_data)

# Combine data
comb_data_090_pre_r025 <- filter_comb(comb_data, 0.9, 2.5, TRUE)
comb_data_090_pre_r020 <- filter_comb(comb_data, 0.9, 2.0, TRUE)
comb_data_090_pre_r030 <- filter_comb(comb_data, 0.9, 3.0, TRUE)
comb_data_090_post_r025 <- filter_comb(comb_data, 0.9, 2.5, FALSE)
comb_data_080_pre_r025 <- filter_comb(comb_data, 0.8, 2.5, TRUE)

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
fontscale_legend <- 1.25

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
         "Bidirectional tracing")
}


label_levels <- c("No tracing", "Manual, 2-day",
                  "Manual, 7-day",
                  "Digital, low-uptake",
                  "Digital, high-uptake",
                  "Hybrid, low-uptake, 2-day",
                  "Hybrid, low-uptake, 7-day",
                  "Hybrid, high-uptake, 2-day",
                  "Hybrid, high-uptake, 7-day")

blab_size = 6
blab_buffer_ctrl = 0.05
blab_buffer_reff = 0.17
blab_col = "#377eb8"
blab_annot = "‡"

make_barplot_fig4_base <- function(g, data, coord_ratio, y_ref = "p_controlled_upper",
                                   baseline_lab_size = blab_size, 
                                   baseline_lab_buffer = blab_buffer_ctrl,
                                   baseline_lab_col = blab_col,
                                   baseline_lab_annot = blab_annot){
  data_annot <- data %>% filter(label %in% c("No tracing", "Manual, 2-day"), 
                                backtrace_distance == Inf)
  data_annot[["y_annot"]] <- data_annot[[y_ref]] + baseline_lab_buffer
  g <- g +
    geom_text(aes(y=y_annot, x=label), data=data_annot, 
              label = baseline_lab_annot, size = baseline_lab_size,
              colour = baseline_lab_col, vjust = 0, hjust = 0.5) +
    facet_grid(.~scenario_type) +
    scale_fill_brewer(palette = "Dark2", name=NULL, labels=label_backtrace) +
    guides(fill=guide_legend(nrow=1,order=1), alpha=FALSE) +
    coord_fixed(ratio=coord_ratio) +
    theme_base + theme(
      axis.text.x = element_text(angle=45,hjust=1,vjust=1,
                                 size = fontsize_base * fontscale_mid,
                                 colour = "black"),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.margin = margin(t=0.7, unit="cm"),
      strip.text = element_text(face = "bold", size = fontsize_base * fontscale_title),
    )
  return(g)
}

make_barplot_fig4_ctrl <- function(data, err_width = 0.2, bar_width=0.6,
                                   coord_ratio = 6, 
                                   baseline_lab_size = blab_size, 
                                   baseline_lab_buffer = blab_buffer_ctrl,
                                   baseline_lab_col = blab_col,
                                   baseline_lab_annot = blab_annot){
  g <- ggplot(data, aes(x=factor(label, levels=label_levels), y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                        fill=factor(backtrace_distance))) +
    geom_col(position="dodge", width=bar_width) + 
    geom_errorbar(width=err_width, position = position_dodge(width=bar_width)) + 
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100))
  h <- make_barplot_fig4_base(g, data, coord_ratio, "p_controlled_upper", baseline_lab_size,
                              baseline_lab_buffer, baseline_lab_col, baseline_lab_annot)
  return(h)
}

make_barplot_fig4_reff <- function(data, r0=2.5, bar_width=0.6,
                                   legend_pos = c(0,1), legend_vspace=0,
                                   coord_ratio = 6/2.5,
                                   r0_lab_x = 8.5, r0_lab_ydiff = 0.15,
                                   r0_lab_size = 4.6, ylim_diff = 0.05,
                                   baseline_lab_size = blab_size, 
                                   baseline_lab_buffer = blab_buffer_reff,
                                   baseline_lab_col = blab_col,
                                   baseline_lab_annot = blab_annot){
  g <- ggplot(data, aes(x=factor(label, levels=label_levels), y=effective_r0_mean,
                        fill=factor(backtrace_distance))) +
    geom_col(position="dodge", width=bar_width) + 
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    geom_hline(yintercept = r0, linetype = "dotted", size=1,
               colour = "red") +
    annotate("text", x=r0_lab_x, y=r0-r0_lab_ydiff, hjust=0.5, vjust = 0.5,
             colour = "red", size=r0_lab_size, label=expression(R[0])) +
    scale_y_continuous(name = "Mean effective\nreproduction number",
                       limits = c(0,r0 + ylim_diff),
                       breaks = seq(0,4,0.5))
  h <- make_barplot_fig4_base(g, data, coord_ratio, "effective_r0_mean", 
                              baseline_lab_size,
                              baseline_lab_buffer, baseline_lab_col,
                              baseline_lab_annot)
  return(h)
}

make_barplot_fig4_grid <- function(data, err_width = 0.2, bar_width=0.6,
                              coord_ratio = 6, ylim_diff = 0.05,
                              baseline_lab_size = blab_size, 
                              baseline_lab_buffer_ctrl = blab_buffer_ctrl,
                              baseline_lab_buffer_reff = blab_buffer_reff,
                              baseline_lab_col = blab_col,
                              baseline_lab_annot = blab_annot, r0 = 2.5,
                              r0_lab_x = 8.5, r0_lab_ydiff = 0.15,
                              r0_lab_size = 4.6, squash = TRUE){
  ctrl_plot <- make_barplot_fig4_ctrl(data, err_width = err_width,
                                      bar_width=bar_width, coord_ratio=coord_ratio,
                                      baseline_lab_size = baseline_lab_size,
                                      baseline_lab_buffer = baseline_lab_buffer_ctrl,
                                      baseline_lab_col = baseline_lab_col,
                                      baseline_lab_annot = baseline_lab_annot)
  reff_plot <- make_barplot_fig4_reff(data, r0 = r0, ylim_diff = ylim_diff,
                                      bar_width=bar_width, coord_ratio=coord_ratio/r0,
                                      r0_lab_x = r0_lab_x, r0_lab_ydiff = r0_lab_ydiff,
                                      r0_lab_size = r0_lab_size,
                                      baseline_lab_size = baseline_lab_size,
                                      baseline_lab_buffer = baseline_lab_buffer_reff,
                                      baseline_lab_col = baseline_lab_col,
                                      baseline_lab_annot = baseline_lab_annot)
  if (squash){
    ctrl_plot <- ctrl_plot + theme(
      legend.position = "none",
      axis.text.x = element_blank()
    )
    reff_plot <- reff_plot + theme(
      strip.text.y = element_blank(),
      strip.text.x = element_blank(),
    )
    grid_plot <- plot_grid(ctrl_plot, reff_plot,
                           nrow = 2, ncol = 1, align = "hv",
                           axis = "l", rel_heights = c(1,2))
  } else {
    rel_heights <- c(1,1)
    grid_plot <- plot_grid(ctrl_plot, reff_plot,
                           labels = "auto", nrow = 2, ncol = 1, align = "hv",
                           axis = "l", label_size = fontsize_base * fontscale_label,
                           label_fontfamily = titlefont, label_colour = "black",
                           rel_heights = rel_heights
    )
  }
  return(grid_plot)
}

#==============================================================================
# Make plots
#==============================================================================

# R0 2.5, 90% contacts traced, trace before testing
grid_090_pre_r025 <- make_barplot_fig4_grid(comb_data_090_pre_r025)

# R0 3.0, 90% contacts traced, trace before testing
grid_090_pre_r030 <- make_barplot_fig4_grid(comb_data_090_pre_r030,
                                            squash = FALSE, r0 = 3.0)

# R0 2.0, 90% contacts traced, trace before testing
grid_090_pre_r020 <- make_barplot_fig4_grid(comb_data_090_pre_r020,
                                            squash = FALSE, r0 = 2.0)

# R0 2.5, 80% contacts traced, trace before testing
grid_080_pre_r025 <- make_barplot_fig4_grid(comb_data_080_pre_r025,
                                            squash = FALSE)

# R0 2.5, 90% contacts traced, trace after testing
grid_090_post_r025 <- make_barplot_fig4_grid(comb_data_090_post_r025,
                                             squash = FALSE)

#==============================================================================
# Save outputs
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, 
                     plot_scale = 34,
                     plot_ratio = 1,
                     device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix, ".", device), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

# Main figure (R0 = 2.5 only)
path_prefix <- "figures_local/img/"
path_prefix_main <- paste0(path_prefix, "fig4_090_pre_r025")
save_fig(path_prefix_main, "", grid_090_pre_r025, device = "svg")

# Supplementary variants
save_fig(path_prefix, "fig4_090_pre_030", grid_090_pre_r030)
save_fig(path_prefix, "fig4_090_pre_020", grid_090_pre_r020)
save_fig(path_prefix, "fig4_080_pre_025", grid_080_pre_r025)
save_fig(path_prefix, "fig4_090_post_025", grid_090_post_r025)

save_fig(path_prefix, "fig4_090_pre_030", grid_090_pre_r030, device="svg")
save_fig(path_prefix, "fig4_090_pre_020", grid_090_pre_r020, device="svg")




