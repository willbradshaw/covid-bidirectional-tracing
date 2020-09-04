library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Optimistic scenario
optim_path <- "figure_data/si/fig3_optimistic_explore_1k_scenario.tsv.gz"
optim_data <- suppressMessages(read_tsv(optim_path))

optim_digital_path <- "figure_data/si/fig3_baseline_explore_digital_1k_scenario.tsv.gz"
optim_data_digital <- suppressMessages(read_tsv(optim_digital_path))

optim_data <- optim_data %>% filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  bind_rows(optim_data_digital)

# pessimistic scenario
pessim_path <- "figure_data/si/fig3_pessimistic_explore_1k_scenario.tsv.gz"
pessim_data <- suppressMessages(read_tsv(pessim_path))

pessim_digital_path <- "figure_data/si/fig3_baseline_explore_digital_1k_scenario.tsv.gz"
pessim_data_digital <- suppressMessages(read_tsv(pessim_digital_path))

pessim_data <- pessim_data %>% filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  bind_rows(pessim_data_digital)

#==============================================================================
# Specify plotting theme information
#==============================================================================

# Fonts
font <- "sans" # Main figure font
titlefont <- font # Font for axis titles etc. (if different)
fontsize_base <- 12 # Basic figure font size
fontscale_title <- 1.5 # Default axis-title scale relative to regular font
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

#==============================================================================
# Labelling functions
#==============================================================================

label_backtrace <- function(x){
  ifelse(x==0, "Forward tracing only",
         "Forward and reverse\ntracing")
}

label_limits <- function(x){
  ifelse(as.numeric(x) == Inf, "No trace limit", paste0(x, "-day trace limit"))
}

label_limits_manual <- function(x){
  ifelse(as.numeric(x) == Inf, "No manual limit", paste0(x, "-day manual limit"))
}


label_vars <- function(x){
  ifelse(x=="p_data_sharing_auto", "% of cases sharing data",
         "% of cases with\ntrace-enabled smartphones")
}

label_pc <- function(x) round(as.numeric(x)*100)

#==============================================================================
# Make plots
#==============================================================================

#------------------------------------------------------------------------------
# Base functions
#------------------------------------------------------------------------------

ttype_levels <- c("Manual tracing only", "Manual + digital",
                  "Digital only", "No tracing")

add_trace_type <- function(data, ttype_exclude = NULL){
  data_ttype <- mutate(data, trace_type = ifelse(p_traced_auto == 0,
    ifelse(p_traced_manual == 0, "No tracing", "Manual tracing only"),
    ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")))
  if (!is.null(ttype_exclude)){
    data_ttype <- filter(data_ttype, !(trace_type %in% ttype_exclude))
  }
  return(data_ttype)
}

control_plot_fig3 <- function(data, lt_var_name, lt_var_levels, lt_labeller,
                              err_width = 0.05, pt_size = 2, x_axis_title = expression(R[0]),
                              coord_ratio = 2.5, legend_pos = c(1,1),
                              legend_vspace = -1, ttype_exclude = NULL,
                              col_palette = "Set1"){
  data <- add_trace_type(data, ttype_exclude)
  data[["lt"]] <- data[[lt_var_name]]
  g <- ggplot(data, aes(x=r0_base, y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                        colour=factor(trace_type,
                                      levels = ttype_levels),
                        linetype=factor(lt, levels = lt_var_levels),
                        shape=factor(lt, levels = lt_var_levels))) +
    geom_errorbar(width=err_width) + geom_line() + geom_point(size=pt_size) +
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_colour_brewer(palette = col_palette, name=NULL) +
    scale_x_continuous(name = x_axis_title) +
    scale_linetype_discrete(name = NULL, labels = lt_labeller) +
    scale_shape_discrete(name = NULL, labels = lt_labeller) +
    coord_fixed(ratio = coord_ratio) + 
    guides(colour=guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base + theme(
      legend.justification = c("right","top"),
      legend.background = element_rect(fill=alpha("white", 0.5)),
      legend.key = element_rect(fill=alpha("white", 0)),
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.position = legend_pos,
      legend.spacing.y = unit(legend_vspace, "mm"),
    )
  return(g)
}

reff_plot_fig3 <- function(data, lt_var_name, lt_var_levels, lt_labeller,
                           pt_size = 2, x_axis_title = expression(R[0]),
                           coord_ratio = 1, legend_pos = c(0.01, 1), 
                           legend_vspace = -1, ttype_exclude = NULL,
                           y_title_short = TRUE, y_break_max = 10,
                           y_limits = c(0, NA), col_palette = "Set1"){
  data <- add_trace_type(data, ttype_exclude)
  data[["lt"]] <- data[[lt_var_name]]
  y_title <- ifelse(y_title_short, "Mean effective reprod. number",
                    "Mean effective\nreproduction number")
  g <- ggplot(data, aes(x=r0_base, y=effective_r0_mean,
                        colour=factor(trace_type,
                                      levels = ttype_levels),
                        linetype=factor(lt, levels = lt_var_levels),
                        shape=factor(lt, levels = lt_var_levels))) +
    geom_line() + geom_point(size=pt_size) +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    scale_y_continuous(name = y_title, limits = y_limits,
                       breaks = seq(0,y_break_max,0.5)) +
    scale_colour_brewer(palette = col_palette, name=NULL) +
    scale_x_continuous(name = x_axis_title) +
    scale_linetype_discrete(name = NULL, labels = lt_labeller) +
    scale_shape_discrete(name = NULL, labels = lt_labeller) +
    coord_fixed(ratio = coord_ratio) + 
    guides(colour=guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base + theme(
      legend.justification = c("left","top"),
      legend.background = element_rect(fill=alpha("white", 0.5)),
      legend.key = element_rect(fill=alpha("white", 0)),
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.position = legend_pos,
      legend.spacing.y = unit(legend_vspace, "mm"),
    )
  return(g)
}

prepare_fig3_op <- function(optim_data, pessim_data, 
                                  p_traced, trace_neg_sym, contact_limit = NULL, bt_distance = NULL,
                                  ttype_exclude = NULL, pt_size = 2,
                                  err_width_ctrl = 0.05,
                                  coord_ratio_ctrl = 2.5, legend_pos_ctrl = c(1,1),
                                  coord_ratio_reff = 0.625, legend_pos_reff = c(0.01, 1),
                                  y_limits_reff = c(0, 4), y_title_short_reff = TRUE,
                                  legend_vspace = -1, col_palette = "Set1"
){
  if (!xor(is.null(contact_limit),is.null(bt_distance))){
    stop("Exactly one of contact_limit and bt_distance must be non-null.")
  }
  # Filter data by p_traced and linetype criterion
  data <- list(optim_data, pessim_data)
  data_filtered <- lapply(data, function(d) 
    filter(d, p_traced_manual %in% c(0, -1, p_traced),
           p_traced_auto %in% c(0, -1, p_traced),
           trace_neg_symptomatic == trace_neg_sym))
  if (!is.null(contact_limit)) data_filtered <- lapply(data_filtered, function(d)
    filter(d, contact_limit_manual == contact_limit))
  if (!is.null(bt_distance)) data_filtered <- lapply(data_filtered, function(d)
    filter(d, backtrace_distance == bt_distance))
  # Add and filter by trace types
  data_ttype <- lapply(data_filtered, function(d) add_trace_type(d, ttype_exclude))
  # Determine linetype configuration
  lt_var_name <- ifelse(is.null(contact_limit), "contact_limit_manual",
                        "backtrace_distance")
  lt_labeller <- ifelse(is.null(contact_limit), label_limits_manual,
                        label_backtrace)
  lt_var_levels <- ifelse(rep(is.null(contact_limit), 2), c(7,2),
                          c(Inf, 0))
  # Specify x-axis titles
  x_axis_titles <- c(expression(R[0]~"(optimistic scenario)"),
                     expression(R[0]~"(pessimistic scenario)"))
  # Make control plots
  ctrl_plots <- lapply(1:2, function(n)
    control_plot_fig3(data_ttype[[n]], lt_var_name, lt_var_levels, lt_labeller,
                      err_width_ctrl, pt_size, x_axis_titles[n], coord_ratio_ctrl,
                      legend_pos_ctrl, legend_vspace, NULL, col_palette))
  # Make R_eff plots
  reff_plots <- lapply(1:2, function(n)
    reff_plot_fig3(data_ttype[[n]], lt_var_name, lt_var_levels, lt_labeller, 
                   pt_size, x_axis_titles[n],
                   coord_ratio_reff, legend_pos_reff, legend_vspace,
                   NULL, y_title_short_reff, 100, y_limits_reff,
                   col_palette))
  # Assemble figure
  fig <- plot_grid(ctrl_plots[[1]], reff_plots[[1]],
                   ctrl_plots[[2]], reff_plots[[2]],
                   labels = "auto", nrow = 2, ncol = 2, align = "hv",
                   axis = "l", label_size = fontsize_base * fontscale_label,
                   label_fontfamily = titlefont, label_colour = "black")
  return(fig)
}

#------------------------------------------------------------------------------
# Make different figure versions
#------------------------------------------------------------------------------

ttype_exclude <- NULL

# Combinations: 2/7 day trace limit, 80/90% tracing, trace_neg_symptomatic policy

fig_7d_080_pre <- prepare_fig3_op(optim_data, pessim_data, 0.8, TRUE,
                                  7, NULL, y_title_short_reff = FALSE,
                                  ttype_exclude = ttype_exclude)

fig_2d_090_pre <- prepare_fig3_op(optim_data, pessim_data, 0.9, TRUE,
                                  2, NULL, y_title_short_reff = FALSE,
                                  ttype_exclude = ttype_exclude)

fig_2d_080_pre <- prepare_fig3_op(optim_data, pessim_data, 0.8, TRUE,
                                  2, NULL, y_title_short_reff = FALSE,
                                  ttype_exclude = ttype_exclude)

fig_7d_090_post <- prepare_fig3_op(optim_data, pessim_data, 0.9, FALSE,
                                  7, NULL, y_title_short_reff = FALSE,
                                  ttype_exclude = ttype_exclude)

fig_7d_080_post <- prepare_fig3_op(optim_data, pessim_data, 0.8, FALSE,
                                  7, NULL, y_title_short_reff = FALSE,
                                  ttype_exclude = ttype_exclude)

fig_2d_090_post <- prepare_fig3_op(optim_data, pessim_data, 0.9, FALSE,
                                  2, NULL, y_title_short_reff = FALSE,
                                  ttype_exclude = ttype_exclude)

fig_2d_080_post <- prepare_fig3_op(optim_data, pessim_data, 0.8, FALSE,
                                  2, NULL, y_title_short_reff = FALSE,
                                  ttype_exclude = ttype_exclude)

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = "png", width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures/si/op_explore_"
row_height <- 12

save_fig(path_prefix, "7d_080_pre.png", fig_7d_080_pre, row_height*2, 2 * 0.7)
save_fig(path_prefix, "2d_080_pre.png", fig_2d_080_pre, row_height*2, 2 * 0.7)
save_fig(path_prefix, "2d_090_pre.png", fig_2d_090_pre, row_height*2, 2 * 0.7)
save_fig(path_prefix, "7d_080_post.png", fig_7d_080_post, row_height*2, 2 * 0.7)
save_fig(path_prefix, "7d_090_post.png", fig_7d_090_post, row_height*2, 2 * 0.7)
save_fig(path_prefix, "2d_080_post.png", fig_2d_080_post, row_height*2, 2 * 0.7)
save_fig(path_prefix, "2d_090_post.png", fig_2d_090_post, row_height*2, 2 * 0.7)