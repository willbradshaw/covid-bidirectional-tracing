library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

dispersion_path <- "figure_data/fig3/dispersion_200_sc.tsv.gz"
dispersion_data <- suppressMessages(read_tsv(dispersion_path))

#==============================================================================
# Process data
#==============================================================================

dispersion_data_processed <- dispersion_data %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Automated only", 
                                    "Manual + automated")))

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

# Control plots

control_plot_90 <- dispersion_data_processed %>% filter(p_traced_auto != 0.8,
                                           p_traced_manual != 0.8) %>%
  ggplot(aes(x=dispersion, y=p_controlled, 
             ymin=p_controlled_lower, ymax=p_controlled_upper,
             colour = trace_type,
             linetype = factor(backtrace_distance),
             shape = factor(backtrace_distance))) +
  geom_errorbar(width=0.05) +
  geom_line() + geom_point(size=2) +
  facet_grid(.~contact_limit_manual,
             labeller=labeller(backtrace_distance=label_backtrace,
                               contact_limit_manual=label_limits_manual)) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_log10(name = "Overdispersion (k)") +
  scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
  scale_linetype_discrete(name=NULL, labels=label_backtrace) +
  scale_shape_discrete(name=NULL, labels=label_backtrace) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(nrow=2, order=2),
         shape=guide_legend(nrow=2, order=2)
  ) +
  coord_fixed(ratio=1) +
  theme_base +
  ggtitle("Control rates with varying overdispersion\n(90% contacts traced, 200 runs)")
  
control_plot_80 <- dispersion_data_processed %>% filter(p_traced_auto != 0.9,
                                                        p_traced_manual != 0.9) %>%
  ggplot(aes(x=dispersion, y=p_controlled, 
             ymin=p_controlled_lower, ymax=p_controlled_upper,
             colour = trace_type,
             linetype = factor(backtrace_distance),
             shape = factor(backtrace_distance))) +
  geom_errorbar(width=0.05) +
  geom_line() + geom_point(size=2) +
  facet_grid(.~contact_limit_manual,
             labeller=labeller(backtrace_distance=label_backtrace,
                               contact_limit_manual=label_limits_manual)) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_log10(name = "Overdispersion (k)") +#, limits = c(0,0.5),
  #breaks = seq(0,1,0.1), labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
  scale_linetype_discrete(name=NULL, labels=label_backtrace) +
  scale_shape_discrete(name=NULL, labels=label_backtrace) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(nrow=2, order=2),
         shape=guide_legend(nrow=2, order=2)
  ) +
  coord_fixed(ratio=1) +
  theme_base +
  ggtitle("Control rates with varying overdispersion\n(80% contacts traced, 200 runs)")

# R_eff plots

r_eff_plot_90 <- dispersion_data_processed %>% filter(p_traced_auto != 0.8,
                                                p_traced_manual != 0.8) %>%
  ggplot(aes(x=dispersion, y=effective_r0_mean, 
             colour = trace_type,
             linetype = factor(backtrace_distance),
             shape = factor(backtrace_distance),
             ymin = effective_r0_mean-effective_r0_sd,
             ymax = effective_r0_mean+effective_r0_sd)) +
  geom_line() + geom_point(size=2) +
  geom_errorbar(width=0.025) +
  facet_grid(.~contact_limit_manual,
             labeller=labeller(backtrace_distance=label_backtrace,
                               contact_limit_manual=label_limits_manual)) +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  geom_hline(yintercept=2.5, linetype="dotted", size=1, colour="red") +
  scale_y_continuous(name = "Avg effective reprod. number",
                     limits=c(0,NA), breaks = seq(0,10,0.5)) +
  scale_x_log10(name = "Overdispersion (k)") +
  scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
  scale_linetype_discrete(name=NULL, labels=label_backtrace) +
  scale_shape_discrete(name=NULL, labels=label_backtrace) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(nrow=2, order=2),
         shape=guide_legend(nrow=2, order=2)
  ) +
  coord_fixed(ratio=0.5) +
  theme_base +
  ggtitle("Effective R with varying overdispersion\n(90% contacts traced, 200 runs)")


r_eff_plot_80 <- dispersion_data_processed %>% filter(p_traced_auto != 0.9,
                                                      p_traced_manual != 0.9) %>%
  ggplot(aes(x=dispersion, y=effective_r0_mean, 
             colour = trace_type,
             linetype = factor(backtrace_distance),
             shape = factor(backtrace_distance),
             ymin = effective_r0_mean-effective_r0_sd,
             ymax = effective_r0_mean+effective_r0_sd)) +
  geom_line() + geom_point(size=2) +
  geom_errorbar(width=0.025) +
  facet_grid(.~contact_limit_manual,
             labeller=labeller(backtrace_distance=label_backtrace,
                               contact_limit_manual=label_limits_manual)) +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  geom_hline(yintercept=2.5, linetype="dotted", size=1, colour="red") +
  scale_y_continuous(name = "Avg effective reprod. number",
                     limits=c(0,NA), breaks = seq(0,10,0.5)) +
  scale_x_log10(name = "Overdispersion (k)") +
  scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
  scale_linetype_discrete(name=NULL, labels=label_backtrace) +
  scale_shape_discrete(name=NULL, labels=label_backtrace) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(nrow=2, order=2),
         shape=guide_legend(nrow=2, order=2)
  ) +
  coord_fixed(ratio=0.5) +
  theme_base +
  ggtitle("Effective R with varying overdispersion\n(80% contacts traced, 200 runs)")


#==============================================================================
# Save output
#==============================================================================

plot_scale <- 10
plot_prefix <- "figures/fig3/dispersion_"
ggsave(filename=paste0(plot_prefix, "control_90.png"), plot = control_plot_90,
       device = "png", width = plot_scale*3,
       height = plot_scale*2, units = "cm",
       dpi = 320, limitsize = FALSE)
ggsave(filename=paste0(plot_prefix, "control_80.png"), plot = control_plot_80,
       device = "png", width = plot_scale*3,
       height = plot_scale*2, units = "cm",
       dpi = 320, limitsize = FALSE)
ggsave(filename=paste0(plot_prefix, "r_eff_90.png"), plot = r_eff_plot_90,
       device = "png", width = plot_scale*3,
       height = plot_scale*2, units = "cm",
       dpi = 320, limitsize = FALSE)
ggsave(filename=paste0(plot_prefix, "r_eff_80.png"), plot = r_eff_plot_80,
       device = "png", width = plot_scale*3,
       height = plot_scale*2, units = "cm",
       dpi = 320, limitsize = FALSE)