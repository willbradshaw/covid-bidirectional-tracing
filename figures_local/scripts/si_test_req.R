library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Universal coverage
data_path <- "figures_local/data/fig2_test_req_equiv_500_scenario.tsv.gz"
data <- suppressMessages(read_tsv(data_path)) %>% 
  filter(test_sensitivity == 0.8)

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
  ifelse(as.numeric(x) == Inf, "No data limit", paste0(x, "-day data limit"))
}

label_vars <- function(x){
  ifelse(x=="p_data_sharing_auto", "% of cases sharing data",
         "% of cases with\ntrace-enabled smartphones")
}

label_pc <- function(x) round(as.numeric(x)*100)

label_sens <- function(x){
  paste0(as.numeric(x)*100, "% sensitivity")
}

#==============================================================================
# Make plots
#==============================================================================

legend_vspace = -1

# Forward + reverse

plot_ctrl_combo <- data %>% filter(p_traced_auto < 0) %>%
  ggplot(aes(x=p_traced_manual, y=p_controlled, 
             ymin = p_controlled_lower, ymax = p_controlled_upper,
             colour = factor(backtrace_distance, levels=c(0,Inf)),
             linetype = factor(contact_limit_manual, levels=c(7,2)),
             shape = factor(contact_limit_manual, levels=c(7,2)))) +
  geom_errorbar(width=0.05) + geom_line() + geom_point(size=2) +
  # facet_grid(test_sensitivity~., 
  #            labeller=labeller(test_sensitivity=label_sens)) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = label_pc) +
  scale_x_continuous(name="% of non-environmental\ncontacts traced",
                     limits=c(0,1.05), breaks = seq(0,1,0.2), labels=label_pc) +
  scale_colour_brewer(labels=label_backtrace, name=NULL, palette = "Dark2") +
  scale_linetype_discrete(name=NULL, labels = label_limits_manual) +
  scale_shape_discrete(name=NULL, labels = label_limits_manual) +
  coord_fixed(ratio = 1) +
  guides(colour=guide_legend(ncol = 1, order = 1),
         linetype = guide_legend(ncol = 1, order = 2),
         shape = guide_legend(ncol = 1, order = 2)) +
  theme_base + theme(
    legend.justification = c("left","top"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.position = c(0.05,1),
    legend.spacing.y = unit(legend_vspace, "mm"),
  )

plot_reff_combo <- data %>% filter(p_traced_auto < 0) %>%
  ggplot(aes(x=p_traced_manual, y=effective_r0_mean, 
             colour = factor(backtrace_distance, levels=c(0,Inf)),
             linetype = factor(contact_limit_manual, levels=c(7,2)),
             shape = factor(contact_limit_manual, levels=c(7,2)))) +
  geom_line() + geom_point(size=2) +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  geom_hline(yintercept = 2.5, linetype = "dotted", size=1,
                      colour = "red") +
  annotate("text", x=0.98, y=2.37, hjust=0.5, vjust = 0.5,
           colour = "red", size=4.6, label=expression(R[0])) +
  # facet_grid(test_sensitivity~., 
  #            labeller=labeller(test_sensitivity=label_sens)) +
  scale_y_continuous(name = "Mean effective reproduction number", limits = c(0,2.5),
                     breaks = seq(0,10,0.5)) +
  scale_x_continuous(name="% of non-environmental\ncontacts traced",
                     limits=c(0,1), breaks = seq(0,1,0.2), labels=label_pc) +
  scale_colour_brewer(labels=label_backtrace, name=NULL, palette = "Dark2") +
  scale_linetype_discrete(name=NULL, labels = label_limits_manual) +
  scale_shape_discrete(name=NULL, labels = label_limits_manual) +
  coord_fixed(ratio = 1/2.5) +
  guides(colour=guide_legend(ncol = 1, order = 1),
         linetype = guide_legend(ncol = 1, order = 2),
         shape = guide_legend(ncol = 1, order = 2)) +
  theme_base + theme(
    legend.justification = c("left","bottom"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.position = c(0.02,0.01),
    legend.spacing.y = unit(legend_vspace, "mm"),
  )


# Combine together
plot_combo <- plot_grid(plot_ctrl_combo, plot_reff_combo, 
                         labels = "auto", nrow = 1, ncol = 2, align = "hv",
                         axis = "l", label_size = fontsize_base * fontscale_label,
                         label_fontfamily = titlefont, label_colour = "black")


# Forward only
plot_ctrl_manual <- data %>% filter(p_traced_auto == 0) %>%
  ggplot(aes(x=p_traced_manual, y=p_controlled, 
             ymin = p_controlled_lower, ymax = p_controlled_upper,
             colour = factor(backtrace_distance, levels=c(0,Inf)),
             linetype = factor(contact_limit_manual, levels=c(7,2)),
             shape = factor(contact_limit_manual, levels=c(7,2)))) +
  geom_errorbar(width=0.05) + geom_line() + geom_point(size=2) +
  facet_grid(test_sensitivity~., 
             labeller=labeller(test_sensitivity=label_sens)) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = label_pc) +
  scale_x_continuous(name="% of non-environmental\ncontacts traced",
                     limits=c(0,1.05), breaks = seq(0,1,0.2), labels=label_pc) +
  scale_colour_brewer(labels=label_backtrace, name=NULL, palette = "Dark2") +
  scale_linetype_discrete(name=NULL, labels = label_limits_manual) +
  scale_shape_discrete(name=NULL, labels = label_limits_manual) +
  coord_fixed(ratio = 1) +
  guides(colour=guide_legend(ncol = 1, order = 1),
         linetype = guide_legend(ncol = 1, order = 2),
         shape = guide_legend(ncol = 1, order = 2)) +
  theme_base + theme(
    legend.justification = c("left","top"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.position = c(0.05,1),
    legend.spacing.y = unit(legend_vspace, "mm"),
  )

plot_reff_manual <- data %>% filter(p_traced_auto == 0) %>%
  ggplot(aes(x=p_traced_manual, y=effective_r0_mean, 
             colour = factor(backtrace_distance, levels=c(0,Inf)),
             linetype = factor(contact_limit_manual, levels=c(7,2)),
             shape = factor(contact_limit_manual, levels=c(7,2)))) +
  geom_line() + geom_point(size=2) +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  geom_hline(yintercept = 2.5, linetype = "dotted", size=1,
             colour = "red") +
  annotate("text", x=0.98, y=2.37, hjust=0.5, vjust = 0.5,
           colour = "red", size=4.6, label=expression(R[0])) +
  facet_grid(test_sensitivity~., 
             labeller=labeller(test_sensitivity=label_sens)) +
  scale_y_continuous(name = "Mean effective reproduction number", limits = c(0,2.5),
                     breaks = seq(0,10,0.5)) +
  scale_x_continuous(name="% of non-environmental\ncontacts traced",
                     limits=c(0,1), breaks = seq(0,1,0.2), labels=label_pc) +
  scale_colour_brewer(labels=label_backtrace, name=NULL, palette = "Dark2") +
  scale_linetype_discrete(name=NULL, labels = label_limits_manual) +
  scale_shape_discrete(name=NULL, labels = label_limits_manual) +
  coord_fixed(ratio = 1/2.5) +
  guides(colour=guide_legend(ncol = 1, order = 1),
         linetype = guide_legend(ncol = 1, order = 2),
         shape = guide_legend(ncol = 1, order = 2)) +
  theme_base + theme(
    legend.justification = c("left","bottom"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.position = c(0.02,0.01),
    legend.spacing.y = unit(legend_vspace, "mm"),
  )


# Combine together
plot_manual <- plot_grid(plot_ctrl_manual, plot_reff_manual, 
                        labels = "auto", nrow = 1, ncol = 2, align = "hv",
                        axis = "l", label_size = fontsize_base * fontscale_label,
                        label_fontfamily = titlefont, label_colour = "black")

#==============================================================================
# Save output
#==============================================================================

save_fig <- function(path_prefix, path_suffix, plot, plot_scale,
                     plot_ratio){
  ggsave(filename=paste0(path_prefix, path_suffix), plot = plot,
         device = "png", width = plot_scale * plot_ratio,
         height = plot_scale, units = "cm", dpi = 320, limitsize=FALSE)
}

path_prefix <- "figures_local/img/si_test_req_"
row_height <- 14

save_fig(path_prefix, "combo.png", plot_combo, row_height, 0.8*2)
save_fig(path_prefix, "manual.png", plot_manual, row_height, 0.8*2)

