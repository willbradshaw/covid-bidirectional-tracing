library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

data_path <- "figure_data/si/fig3_baseline_explore_1k_scenario.tsv.gz"
data <- suppressMessages(read_tsv(data_path))

digital_path <- "figure_data/si/fig3_baseline_explore_digital_1k_scenario.tsv.gz"
data_digital <- suppressMessages(read_tsv(digital_path))

data <- data %>% filter(!(p_traced_manual == 0 & p_traced_auto == -1)) %>%
  bind_rows(data_digital)

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

add_trace_type <- function(data, ttype_exclude = NULL){
  data_ttype <- mutate(data, trace_type = ifelse(p_traced_auto == 0,
    ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
    ifelse(p_traced_manual == 0, "Digital only", "Manual + digital")))
  if (!is.null(ttype_exclude)){
    data_ttype <- filter(data_ttype, !(trace_type %in% ttype_exclude))
  }
  return(data_ttype)
}

data_ttype <- data %>%
  add_trace_type(ttype_exclude = c("Digital only", "Manual only")) %>%
  filter(p_traced_auto == 0 | p_traced_auto == p_traced_manual | p_traced_auto == -1)

data_processed <- data_ttype %>%
  group_by(r0_base, backtrace_distance, contact_limit_manual, trace_neg_symptomatic) %>%
  mutate(eff_r0_notrace = effective_r0_mean[trace_type == "No tracing"][1]) %>%
  rename(eff_r0_combo = effective_r0_mean) %>%
  filter(trace_type != "No tracing") %>%
  mutate(rel_eff_r0 = eff_r0_combo/eff_r0_notrace)

#==============================================================================
# Make plots
#==============================================================================

rel_reff_plot <- data_processed %>% ungroup %>%
  mutate(p_traced_manual = factor(p_traced_manual, levels = c(0.9, 0.8)),
         trace_neg_symptomatic = factor(trace_neg_symptomatic,
                                        levels=c(TRUE, FALSE))) %>%
  ggplot(aes(x=r0_base, y=rel_eff_r0,
             colour = factor(backtrace_distance),
             linetype = factor(contact_limit_manual),
             shape = factor(contact_limit_manual))) +
  geom_line() + geom_point(size=2) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 1) +
  scale_x_continuous(name = expression(R[0])) +
  scale_y_continuous(Relative~R[eff]~"(combined vs no tracing)") +
  scale_colour_brewer(palette = "Dark2", labels = label_backtrace, name = NULL) +
  scale_linetype_discrete(name = NULL, labels = label_limits_manual) +
  scale_shape_discrete(name = NULL, labels = label_limits_manual) +
  facet_grid(p_traced_manual~trace_neg_symptomatic,
             labeller = labeller(p_traced_manual = label_trace,
                                 trace_neg_symptomatic = label_policy)) +
  coord_fixed(ratio = 4) +
  guides(colour=guide_legend(ncol = 1, order = 1),
         linetype = guide_legend(ncol = 1, order = 2),
         shape = guide_legend(ncol = 1, order = 2)) +
  theme_base

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = "svg", width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures/si/rel_reff_baseline_vs_no"
row_height <- 10

save_fig(path_prefix, ".svg", rel_reff_plot, row_height*2, 3/2 * 0.8)