library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Optimistic scenario
median_path <- "figures_local/data/fig3_baseline_explore_1k_scenario.tsv.gz"
median_data <- suppressMessages(read_tsv(median_path))

median_digital_path <- "figures_local/data/fig3_baseline_explore_digital_1k_scenario.tsv.gz"
median_data_digital <- suppressMessages(read_tsv(median_digital_path))

median_data <- median_data %>% filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  bind_rows(median_data_digital) %>% mutate(scenario_type = "Median scenario")

# Optimistic scenario
optim_path <- "figures_local/data/fig3_optimistic_explore_1k_scenario.tsv.gz"
optim_data <- suppressMessages(read_tsv(optim_path))

optim_digital_path <- "figures_local/data/fig3_optimistic_explore_digital_1k_scenario.tsv.gz"
optim_data_digital <- suppressMessages(read_tsv(optim_digital_path))

optim_data <- optim_data %>% filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  bind_rows(optim_data_digital) %>% mutate(scenario_type = "Optimistic scenario")

# pessimistic scenario
pessim_path <- "figures_local/data/fig3_pessimistic_explore_1k_scenario.tsv.gz"
pessim_data <- suppressMessages(read_tsv(pessim_path))

pessim_digital_path <- "figures_local/data/fig3_pessimistic_explore_digital_1k_scenario.tsv.gz"
pessim_data_digital <- suppressMessages(read_tsv(pessim_digital_path))

pessim_data <- pessim_data %>% filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  bind_rows(pessim_data_digital) %>% mutate(scenario_type = "Pessimistic scenario")

# Combine data
data <- bind_rows(median_data, optim_data, pessim_data)

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
         "Bidirectional tracing")
}

label_limits <- function(x){
  ifelse(as.numeric(x) == Inf, "No trace limit", paste0(x, "-day trace limit"))
}

label_limits_manual <- function(x){
  ifelse(as.numeric(x) == Inf, "No manual limit", paste0(x, "-day manual limit"))
}

label_pc <- function(x) round(as.numeric(x)*100)

label_policy <- function(x){
  ifelse(x, "Test not required", "Test required")
}

label_trace <- function(x){
  paste0(round(as.numeric(x)*100), "% contacts traced")
}

#==============================================================================
# Process data (combined vs manual tracing only)
#==============================================================================

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

prepare_data <- function(data, p_traced, trace_neg_sym){
  data %>% mutate(trace_type = ifelse(p_traced_auto == 0,
      ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
      ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")),
    ) %>%
    mutate(trace_type = factor(trace_type, levels = ttype_levels)) %>%
    filter(p_traced_manual %in% c(0, -1, p_traced),
           p_traced_auto %in% c(0, -1, p_traced),
           trace_neg_symptomatic == trace_neg_sym)
}

process_data_vmanual <- function(data, p_traced = 0.9, trace_neg_sym = TRUE){
  # Prepare data for hybrid vs manual comparison
  data %>% prepare_data(p_traced, trace_neg_sym) %>%
    filter(trace_type %in% c("Manual only", "Manual + digital")) %>%
    group_by(r0_base, backtrace_distance, contact_limit_manual,
             trace_neg_symptomatic, scenario_type) %>%
    spread(trace_type, effective_r0_mean) %>%
    summarise(eff_r0_manual = sum(`Manual only`, na.rm = TRUE),
              eff_r0_combo  = sum(`Manual + digital`, na.rm = TRUE)) %>%
    mutate(rel_eff_r0 = eff_r0_combo/eff_r0_manual) %>% ungroup %>%
    mutate(comparison = "Hybrid vs manual")
}

process_data_vno <- function(data, p_traced = 0.9, trace_neg_sym = TRUE){
  # Prepare data for hybrid vs no tracing
  data %>% prepare_data(p_traced, trace_neg_sym) %>%
    filter(trace_type %in% c("No tracing", "Manual + digital")) %>%
    group_by(r0_base, backtrace_distance, contact_limit_manual,
             trace_neg_symptomatic, scenario_type) %>%
    spread(trace_type, effective_r0_mean) %>%
    summarise(eff_r0_no = sum(`No tracing`, na.rm = TRUE),
              eff_r0_combo  = sum(`Manual + digital`, na.rm = TRUE)) %>%
    mutate(rel_eff_r0 = eff_r0_combo/eff_r0_no) %>% ungroup %>%
    mutate(comparison = "Hybrid vs no tracing")
}

process_data_mvno <- function(data, p_traced = 0.9, trace_neg_sym = TRUE){
  # Prepare data for hybrid vs no tracing
  data %>% prepare_data(p_traced, trace_neg_sym) %>%
    filter(trace_type %in% c("No tracing", "Manual only")) %>%
    group_by(r0_base, backtrace_distance, contact_limit_manual,
             trace_neg_symptomatic, scenario_type) %>%
    spread(trace_type, effective_r0_mean) %>%
    summarise(eff_r0_no = sum(`No tracing`, na.rm = TRUE),
              eff_r0_manual = sum(`Manual only`, na.rm = TRUE)) %>%
    mutate(rel_eff_r0 = eff_r0_manual/eff_r0_no) %>% ungroup %>%
    mutate(comparison = "Manual vs no tracing")
}

data_vmanual <- process_data_vmanual(data)
data_vno <- process_data_vno(data)
data_mvno <- process_data_mvno(data)
data_all <- bind_rows(data_vmanual, data_vno, data_mvno)

#==============================================================================
# Make plots
#==============================================================================

make_rel_reff_plot <- function(data, p_traced = 0.9, 
                               trace_neg_sym = TRUE,
                               coord_ratio = 4){
  data %>%
    ggplot(aes(x=r0_base, y=rel_eff_r0,
               colour = factor(backtrace_distance),
               linetype = factor(contact_limit_manual),
               shape = factor(contact_limit_manual))) +
    geom_line() + geom_point(size=2) +
    geom_hline(yintercept = 1, linetype = "dotted", size = 1) +
    scale_x_continuous(name = expression(R[0])) +
    scale_y_continuous(name = expression(paste("Relative ",R[eff])),
                       breaks = seq(0, 10, 0.2),
                       minor_breaks = seq(0, 10, 0.1)) +
    facet_grid(comparison ~ scenario_type) +
    scale_colour_brewer(palette = "Dark2", labels = label_backtrace, name = NULL) +
    scale_linetype_discrete(name = NULL, labels = label_limits_manual) +
    scale_shape_discrete(name = NULL, labels = label_limits_manual) +
    coord_fixed(ratio = coord_ratio) +
    guides(colour=guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2),
           shape = guide_legend(ncol = 1, order = 2)) +
    theme_base
}

rel_reff_plot <- make_rel_reff_plot(data_all)

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio, device = "png"){
  ggsave(filename=paste0(path_prefix, path_suffix, ".", device), plot = plot,
         device = device, width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/rel_reff"
row_height <- 8
plot_scale <- row_height * 3.5
plot_ratio <- 1

save_fig(path_prefix, "", rel_reff_plot, plot_scale, plot_ratio)