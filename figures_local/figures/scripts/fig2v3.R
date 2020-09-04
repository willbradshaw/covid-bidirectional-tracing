library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Universal coverage

univ_path <- "figure_data/fig2/fig2_univ_env_1k_scenario.tsv.gz"
univ_data <- suppressMessages(read_tsv(univ_path)) %>%
  filter(generation_alpha == 0.064)

# Data sharing
sharing_path <- "figure_data/fig2/sharing_1k_scenario.tsv.gz"
sharing_data <- suppressMessages(read_tsv(sharing_path)) %>%
  filter(generation_alpha == 0.064)

# Smartphone usage
phone_path <- "figure_data/fig2/smartphones_1k_scenario.tsv.gz"
phone_data <- suppressMessages(read_tsv(phone_path)) %>%
  filter(generation_alpha == 0.064)

# Manual tracing (with and without automated)
manual_path <- "figure_data/fig2/fig2_manual_equiv_1000_scenario.tsv.gz"
manual_data <- suppressMessages(read_tsv(manual_path)) %>%
  filter(generation_alpha == 0.064)

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

control_plot <- function(data, x_var_name, x_axis_name, x_breaks, x_lab,
                         lt_var_name = NULL, 
                         lt_var_levels = NULL, lt_labeller = NULL,
                         err_width = 0.05, pt_size = 2, coord_ratio = 1,
                         legend_pos = c(0.05, 0.75), legend_vspace = 0){
  # Prepare data
  data[["x"]] <- data[[x_var_name]]
  if (!is.null(lt_var_name)) data[["lt"]] <- data[[lt_var_name]]
  # Make base plot
  if (is.null(lt_var_name)){
    g <- ggplot(data, aes(x=x, y=p_controlled, 
                          ymin=p_controlled_lower, ymax=p_controlled_upper,
                          colour = factor(backtrace_distance)))
  } else {
    g <- ggplot(data, aes(x=x, y=p_controlled, 
                          ymin=p_controlled_lower, ymax=p_controlled_upper,
                          colour = factor(backtrace_distance),
                          linetype = factor(lt, levels=lt_var_levels),
                          shape = factor(lt, levels=lt_var_levels)))
  }
  # Add geoms
  g <- g + geom_errorbar(width=err_width) + geom_line() + 
    geom_point(size=pt_size) +
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_colour_brewer(palette = "Dark2", labels=label_backtrace,
                        name=NULL) +
    scale_x_continuous(name = x_axis_name, breaks = x_breaks, labels = x_lab
                       ) +
    coord_fixed(ratio = coord_ratio) + 
    guides(colour=guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base + theme(
      legend.justification = c("left","center"),
      legend.background = element_rect(fill=alpha("white", 0.5)),
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.position = legend_pos,
      legend.spacing.y = unit(legend_vspace, "mm"),
    )
  if (!is.null(lt_var_name)){
    g <- g + scale_linetype_discrete(name = NULL, labels = lt_labeller) +
      scale_shape_discrete(name = NULL, labels = lt_labeller)
  }
  return(g)
}

r_eff_plot <- function(data, x_var_name, x_axis_name, x_breaks, x_lab,
                       col_var_name, col_var_levels, col_labeller,
                       lt_var_name = NULL, lt_var_levels = NULL, lt_labeller = NULL,
                       err_width = 0.05, pt_size = 2, coord_ratio = 0.4,
                       legend_pos = c(0.01, 0.21), legend_vspace = -1,
                       y_title_short = TRUE, y_break_max = 10,
                       y_limits = c(0,NA), col_palette = "Set1",
                       r0_base = 2.5, r0_lab_x = 0.98, r0_lab_y = 2.37,
                       r0_lab_size = 4.6, r0_col = "red"){
  # Prepare data
  data[["x"]] <- data[[x_var_name]]
  data[["col"]] <- data[[col_var_name]]
  if (!is.null(lt_var_name)) data[["lt"]] <- data[[lt_var_name]]
  y_title <- ifelse(y_title_short, "Mean effective reprod. number",
                    "Mean effective\nreproduction number")
  # Make base plot
  if (is.null(lt_var_name)){
    g <- ggplot(data, aes(x=x, y=effective_r0_mean, 
                          colour = factor(col, levels = col_var_levels)))
  } else {
    g <- ggplot(data, aes(x=x, y=effective_r0_mean, 
                          colour = factor(col, levels = col_var_levels),
                          linetype = factor(lt, levels=lt_var_levels),
                          shape = factor(lt, levels=lt_var_levels)))
  }
  # Add geoms etc
  g <- g + geom_line() + geom_point(size=pt_size) +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    scale_y_continuous(name = y_title, breaks = seq(0,y_break_max,0.5),
                       limits = y_limits) +
    
    scale_colour_brewer(palette = col_palette, labels=col_labeller,
                        name=NULL) +
    scale_x_continuous(name = x_axis_name, breaks = x_breaks, labels = x_lab) +
    coord_fixed(ratio = coord_ratio) + 
    guides(colour=guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base + theme(
      legend.justification = c("left","center"),
      legend.background = element_rect(fill=alpha("white", 0.5)),
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.position = legend_pos,
      legend.spacing.y = unit(legend_vspace, "mm"),
    )
  if (!is.null(lt_var_name)){
    g <- g + scale_linetype_discrete(name = NULL, labels = lt_labeller) +
      scale_shape_discrete(name = NULL, labels = lt_labeller)
  }
  if (!is.null(r0_base)){
    g <- g + geom_hline(yintercept = r0_base, linetype = "dotted", size=1,
               colour = "red") +
      annotate("text", x=r0_lab_x, y=r0_lab_y, hjust=0.5, vjust = 0.5,
               colour = r0_col, size=r0_lab_size, label=expression(R[0]))
  }
  return(g)
}

show_assumption <- function(plot, line_x, lab_x, lab_y,
                            colour = "#7570b3",
                            lab_txt = "assumption going forward",
                            line_type = "dashed", lab_size = 4.5){
  plot_out <- plot + geom_vline(xintercept = line_x, linetype = line_type,
                                colour = colour, size = 0.5) +
    annotate("text", hjust = 1, vjust = 0.5, x = lab_x, y = lab_y,
             colour = colour, label = lab_txt, size = lab_size)
  return(plot_out)
}

#------------------------------------------------------------------------------
# Subfigure plots
#------------------------------------------------------------------------------

# Universal coverage
univ_plot <- (univ_data %>%
                control_plot("p_traced_auto", "% of non-environmental\ncontacts traced",
                             seq(0,1,0.2), label_pc, "contact_limit_auto",
                             c(Inf, 2), label_limits) +
                ggtitle("Digital tracing only\n(universal coverage)")) %>%
  show_assumption(0.9, 0.88, 0.01)

# Data sharing
sharing_plot <- (sharing_data %>%
  control_plot("p_data_sharing_auto",
               "% of cases sharing data", seq(0,1,0.2), label_pc,
               legend_pos = c(0.05, 0.85)) +
  ggtitle("Digital tracing only\n(partial data sharing)")) %>%
  show_assumption(0.9, 0.88, 0.01)

# Smartphone coverage
phone_plot <- (phone_data %>%
  control_plot("p_smartphone_overall",
               "% of cases with chirping\nsmartphones",
               seq(0,1,0.2), label_pc, legend_pos = c(0.05, 0.85)) +
  ggtitle("Digital tracing only\n(partial smartphone usage)")) %>%
  show_assumption(0.8, 0.78, 0.01)

# Manual tracing only
manual_plot <- manual_data %>% filter(p_traced_auto == 0) %>%
  control_plot("p_traced_manual", "% of non-environmental\ncontacts traced",
               seq(0,1,0.2), label_pc, "contact_limit_manual",
               c(2, 7), label_limits_manual) +
  ggtitle("Manual tracing\nonly")

# Combined tracing
combo_plot <- manual_data %>% filter(p_traced_auto < 0) %>%
  control_plot("p_traced_manual", "% of non-environmental\ncontacts traced",
               seq(0,1,0.2), label_pc, "contact_limit_manual",
               c(2, 7), label_limits_manual) +
  ggtitle("Manual + digital\ntracing")
  

# R_eff plot
reff_plot <- manual_data %>% 
  mutate(trace_type = ifelse(p_traced_auto == 0, "Manual tracing only",
                             "Manual + digital")) %>%
  filter(backtrace_distance == Inf) %>%
  r_eff_plot("p_traced_manual", "% of non-environmental\ncontacts traced", seq(0,1,0.2),
             label_pc, "trace_type", c("Manual tracing only", "Manual + digital"),
             function(x) x, "contact_limit_manual", c(2,7), label_limits_manual) +
  ggtitle(expression(bold(atop("Effect on"~R[eff],"(forward + reverse tracing)"))))

# Combine together
fig2_v1 <- plot_grid(univ_plot, sharing_plot, phone_plot,
                         manual_plot, combo_plot, reff_plot,
                         labels = "auto", nrow = 2, ncol = 3, align = "hv",
                         axis = "l", label_size = fontsize_base * fontscale_label,
                         label_fontfamily = titlefont, label_colour = "black")

fig2_v2 <- plot_grid(manual_plot, univ_plot, sharing_plot, 
                         phone_plot, combo_plot, reff_plot,
                         labels = "auto", nrow = 3, ncol = 2, align = "hv",
                         axis = "l", label_size = fontsize_base * fontscale_label,
                         label_fontfamily = titlefont, label_colour = "black")

# R_eff plot for forward-only tracing (SI)
reff_plot_fwd <- manual_data %>% 
  mutate(trace_type = ifelse(p_traced_auto == 0, "Manual tracing only",
                             "Manual + digital")) %>%
  filter(backtrace_distance == 0) %>%
  r_eff_plot("p_traced_manual", "% of non-environmental\ncontacts traced", seq(0,1,0.2),
             label_pc, "trace_type", c("Manual tracing only", "Manual + digital"),
             function(x) x, "contact_limit_manual", c(2,7), label_limits_manual) +
  ggtitle(expression(bold(atop("Effect on"~R[eff],"(forward tracing only)"))))


#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = "png", width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures/fig2/fig2"

# 48% pre-symptomatic
row_height <- 14
save_fig(path_prefix, "v1.png", fig2_v1, row_height*2, 3/2 * 0.8)
save_fig(path_prefix, "v2.png", fig2_v2, row_height*3, 2/3 * 0.8)
save_fig("figures/si/reff_fwd", ".png", reff_plot_fwd, row_height, 0.8)