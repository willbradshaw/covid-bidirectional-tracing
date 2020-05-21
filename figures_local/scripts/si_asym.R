library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

asym_path <- "figures_local/data/si_asym_contour_500_scenario.tsv.gz"
asym_data <- suppressMessages(read_tsv(asym_path))

#==============================================================================
# Process data
#==============================================================================

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")
rel_r0_keep <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)

add_trace_type <- function(data, ttype_exclude = NULL){
  data_ttype <- mutate(data, trace_type = ifelse(p_traced_auto == 0,
    ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
    ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")))
  if (!is.null(ttype_exclude)){
    data_ttype <- filter(data_ttype, !(trace_type %in% ttype_exclude))
  }
  return(data_ttype)
}

asym_data_processed <- add_trace_type(asym_data) %>%
  mutate(trace_type = factor(trace_type, levels = ttype_levels),
         backtrace_distance = factor(backtrace_distance, levels = c(Inf, 0))) %>%
  filter(rel_r0_asymptomatic %in% rel_r0_keep)

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
  strip.text.x = element_text(size = fontsize_base,# * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t = 3, b = 3, unit="mm")),
  strip.text.y = element_text(size = fontsize_base,# * fontscale_title,
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

label_rel_r0 <- function(x){
  paste("Rel. R0:", x)
}

label_p_asym <- function(x){
  paste0(round(as.numeric(x)*100), "% asym. carriers")
}

#==============================================================================
# Make p_asymptomatic plots
#==============================================================================

#------------------------------------------------------------------------------
# Auxiliary functions
#------------------------------------------------------------------------------

make_control_plot <- function(data, p_traced, trace_limit, coord_ratio = 0.4){
  data %>% filter(p_traced_auto %in% c(-1, 0, p_traced),
                  p_traced_manual %in% c(-1, 0, p_traced),
                  contact_limit_manual == trace_limit) %>%
    ggplot(aes(x=p_asymptomatic, y=p_controlled,
               ymin=p_controlled_lower, ymax=p_controlled_upper,
               colour = trace_type,
               linetype = backtrace_distance, shape = backtrace_distance)) +
    geom_errorbar(width = 0.03) + geom_line() + geom_point(size = 2) +
    facet_grid(rel_r0_asymptomatic~.,
               labeller = labeller(rel_r0_asymptomatic = label_rel_r0)) +
    scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_x_continuous(name = "% of asymptomatic carriers", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
    scale_linetype_discrete(name=NULL, labels=label_backtrace) +
    scale_shape_discrete(name=NULL, labels=label_backtrace) +
    guides(colour=guide_legend(nrow=2, order=1),
           linetype=guide_legend(nrow=2, order=2),
           shape=guide_legend(nrow=2, order=2)) +
    coord_fixed(ratio=coord_ratio) +
    theme_base
}

make_reff_plot <- function(data, p_traced, trace_limit, coord_ratio = 0.4/2.5){
  data %>% filter(p_traced_auto %in% c(-1, 0, p_traced),
                  p_traced_manual %in% c(-1, 0, p_traced),
                  contact_limit_manual == trace_limit) %>%
    ggplot(aes(x=p_asymptomatic, y=effective_r0_mean,
               colour = trace_type,
               linetype = backtrace_distance, shape = backtrace_distance)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    geom_hline(yintercept = 2.5, linetype = "dotted", size=1,
               colour = "red") +
    annotate("text", x=0.98, y=2.53, hjust=0.5, vjust = 0,
             colour = "red", size=4.5, label=expression(italic(R)[0])) +
    facet_grid(rel_r0_asymptomatic~.,
               labeller = labeller(rel_r0_asymptomatic = label_rel_r0)) +
    scale_y_continuous(name = "Mean effective\nreproduction number",
                       limits=c(0,2.8), breaks = seq(0,10,0.5)) +
    scale_x_continuous(name = "% of asymptomatic carriers", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
    scale_linetype_discrete(name=NULL, labels=label_backtrace) +
    scale_shape_discrete(name=NULL, labels=label_backtrace) +
    guides(colour=guide_legend(nrow=2, order=1),
           linetype=guide_legend(nrow=2, order=2),
           shape=guide_legend(nrow=2, order=2)) +
    coord_fixed(ratio=coord_ratio) +
    theme_base
}


#------------------------------------------------------------------------------
# Make p_asymptomatic plots
#------------------------------------------------------------------------------

# Control plots
ctrl_asym_090_7d <- make_control_plot(asym_data_processed, 0.9, 7)
ctrl_asym_090_2d <- make_control_plot(asym_data_processed, 0.9, 2)

# R_eff plots
reff_asym_090_7d <- make_reff_plot(asym_data_processed, 0.9, 7)
reff_asym_090_2d <- make_reff_plot(asym_data_processed, 0.9, 2)


# Combine
asym_090_7d <- plot_grid(ctrl_asym_090_7d, reff_asym_090_7d,
                         labels = "auto", nrow = 1, ncol = 2, align = "hv",
                         axis = "l", label_size = fontsize_base * fontscale_label,
                         label_fontfamily = titlefont, label_colour = "black")

asym_090_2d <- plot_grid(ctrl_asym_090_2d, reff_asym_090_2d,
                         labels = "auto", nrow = 1, ncol = 2, align = "hv",
                         axis = "l", label_size = fontsize_base * fontscale_label,
                         label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, 
                     plot_scale,
                     plot_ratio, 
                     device = "svg"){
  ggsave(filename=paste0(path_prefix, path_suffix, ".", device), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

row_height <- 8
plot_scale <- row_height * 6
plot_ratio <- 0.8
plot_prefix <- "figures_local/img/si_asym_line_090_"

save_fig(plot_prefix, "7d", asym_090_7d, row_height * 5, plot_ratio)
save_fig(plot_prefix, "2d", asym_090_2d, row_height * 5, plot_ratio)
