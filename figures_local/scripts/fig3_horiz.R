library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Optimistic scenario
optim_path <- "figures/data/fig3_optimistic_fine_1k_scenario.tsv.gz"
optim_data <- suppressMessages(read_tsv(optim_path)) %>%
  filter(generation_alpha == 0.397) %>%
  mutate(scenario_type = "Optimistic scenario")

# Pessimistic scenario
pessim_path <- "figures/data/fig3_pessimistic_fine_1k_scenario.tsv.gz"
pessim_data <- suppressMessages(read_tsv(pessim_path)) %>% 
  filter(generation_alpha == -0.095) %>%
  mutate(scenario_type = "Pessimistic scenario")

# Median scenario
median_path <- "figures/data/fig3_baseline_fine_1k_scenario.tsv.gz"
median_data <- suppressMessages(read_tsv(median_path)) %>%
  filter(generation_alpha == 0.064) %>%
  mutate(scenario_type = "Median scenario")


#==============================================================================
# Process data
#==============================================================================

uptake_high <- 0.80
uptake_low  <- 0.53

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

prepare_data_fig3 <- function(data, trace_limit = 7, p_traced = 0.9,
                         trace_neg_sym = TRUE){
  data %>% mutate(trace_type = ifelse(p_traced_auto == 0,
      ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
      ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")),
    ) %>%
    mutate(trace_type = factor(trace_type, levels = ttype_levels),
           backtrace_distance = factor(backtrace_distance, levels = c(Inf, 0))) %>%
    filter(p_traced_manual %in% c(0, -1, p_traced),
           p_traced_auto %in% c(0, -1, p_traced),
           trace_neg_symptomatic == trace_neg_sym,
           contact_limit_manual == trace_limit)
}

comb_data_fig3 <- function(median_data, optim_data, pessim_data, 
                           trace_limit = 7, p_traced = 0.9,
                           trace_neg_sym = TRUE){
  bind_rows(optim_data %>% mutate(scenario_type = "Optimistic scenario"),
            median_data %>% mutate(scenario_type = "Median scenario"),
            pessim_data %>% mutate(scenario_type = "Pessimistic scenario")) %>%
    prepare_data_fig3(trace_limit, p_traced, trace_neg_sym)
}


# Combine datasets
comb_data <- comb_data_fig3(median_data, optim_data, pessim_data)

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
                              margin = margin(t = 3, b = 3, unit="mm"),
                              face="bold"),
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
  ifelse(x==0, "Forward only",
         "Bidirectional")
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

r0_plot_base <- function(g, pt_size, coord_ratio){
  g + geom_line() + geom_point(size=pt_size) +
    facet_grid(.~scenario_type) +
    scale_colour_brewer(palette = "Set1", name="Trace type") +
    scale_x_continuous(name = expression(R[0])) +
    scale_linetype_discrete(name = NULL, labels = label_backtrace) +
    scale_shape_discrete(name = NULL, labels = label_backtrace) +
    coord_fixed(ratio = coord_ratio) + 
    guides(colour=guide_legend(ncol = 2, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base
}

r0_plot_ctrl <- function(data, err_width, pt_size, coord_ratio){
  g <- ggplot(data, aes(x=r0_base, y=p_controlled,
                        ymin=p_controlled_lower, ymax=p_controlled_upper,
                        colour=trace_type,
                        linetype=backtrace_distance,
                        shape=backtrace_distance)) +
    geom_errorbar(width=err_width) + 
    scale_y_continuous(name = "% outbreaks controlled", limits = c(0,1),
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100))
  h <- r0_plot_base(g, pt_size, coord_ratio)
  return(h)
}

r0_plot_reff <- function(data, pt_size, coord_ratio){
  g <- ggplot(data, aes(x=r0_base, y=effective_r0_mean,
                        colour=trace_type,
                        linetype=backtrace_distance,
                        shape=backtrace_distance)) +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    scale_y_continuous(name = "Mean effective\nreproduction number",
                       limits = c(0, data %>% pull(r0_base) %>% max),
                       breaks = seq(0,100,0.5))
  h <- r0_plot_base(g, pt_size, coord_ratio)
  return(h)
}

r0_plot_grid <- function(data, err_width = 0.05, pt_size = 2,
                         coord_ratio=2, squash = TRUE){
  ctrl_plot <- r0_plot_ctrl(data, err_width, pt_size, coord_ratio)
  reff_plot <- r0_plot_reff(data, pt_size, coord_ratio/(data %>% pull(r0_base) %>% max))
  if (squash){
    ctrl_plot <- ctrl_plot + theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )
    reff_plot <- reff_plot + theme(
      strip.text.y = element_blank(),
      strip.text.x = element_blank(),
    )
    grid_plot <- plot_grid(ctrl_plot, reff_plot,
                           nrow = 2, ncol = 1, align = "hv",
                           axis = "l", rel_heights = c(1,1.3))
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

#------------------------------------------------------------------------------
# Make figure
#------------------------------------------------------------------------------

fig3_all <- r0_plot_grid(comb_data, coord_ratio = 2)

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, 
                     plot_scale = 24,
                     plot_ratio = 1.3,
                     device="png"){
  ggsave(filename=paste0(path_prefix, path_suffix, ".", device), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/fig3_horiz"
save_fig(path_prefix, "", fig3_all, device = "svg")
save_fig(path_prefix, "", fig3_all, device = "png")
