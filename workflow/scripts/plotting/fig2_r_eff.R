library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Fig. 2a - universal coverage
univ_path <- "saved_data/fig2b_scenario.tsv.gz"
univ_data <- suppressMessages(read_tsv(fig2b_scenario_path))

# Fig. 2b - data sharing
sharing_path <- "saved_data/fig2c_sharing_scenario.tsv.gz"
sharing_data <- suppressMessages(read_tsv(fig2c_sharing_scenario_path))

# Fig. 2c - smartphone usage
phone_path <- "saved_data/fig2c_smartphones_scenario.tsv.gz"
phone_data <- suppressMessages(read_tsv(fig2c_smartphones_scenario_path))

# Fig. 2d - manual only
manual_path <- "saved_data/fig2e_scenario.tsv.gz"
manual_data <- suppressMessages(read_tsv(fig2e_scenario_path))

# Fig. 2e - manual + automated
combo_path <- "saved_data/fig2d_scenario.tsv.gz"
combo_data <- suppressMessages(read_tsv(fig2d_scenario_path))

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

#==============================================================================
# Make plots
#==============================================================================

# Data sharing

sharing_r_eff_plot <- ggplot(sharing_data,
  aes(x=p_data_sharing_auto, y=effective_r0_mean,
      ymin = effective_r0_mean-effective_r0_sd,
      ymax = effective_r0_mean+effective_r0_sd,
      colour = factor(backtrace_distance))) +
  geom_errorbar(width=0.05) +
  geom_line() + geom_point(size=2) +
  ggtitle("Automated tracing\n(partial data sharing)") +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  geom_hline(yintercept = combo_data$r0_base[1], linetype = "dotted", size=1,
             colour = "red") +
  annotate("text", x=0.98, y=2.37, hjust=0.5, vjust = 0.5,
           colour = "red", size=4.6, label=expression(R[0])) +
  scale_y_continuous(name = "Mean effective\nreproduction number",
                     limits=c(0,NA), breaks = seq(0,2.5,0.5)) +
  scale_x_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  guides(colour=guide_legend(nrow=2)) +
  coord_fixed(ratio=0.3) +
  theme_base + theme(
    legend.position=c(0.05,0.2),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm"))
  )

phone_r_eff_plot <- ggplot(phone_data,
                             aes(x=p_smartphone_overall, y=effective_r0_mean,
                                 ymin = effective_r0_mean-effective_r0_sd,
                                 ymax = effective_r0_mean+effective_r0_sd,
                                 colour = factor(backtrace_distance))) +
  geom_errorbar(width=0.05) +
  geom_line() + geom_point(size=2) +
  ggtitle("Automated tracing\n(partial data sharing)") +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  geom_hline(yintercept = combo_data$r0_base[1], linetype = "dotted", size=1,
             colour = "red") +
  annotate("text", x=0.98, y=2.37, hjust=0.5, vjust = 0.5,
           colour = "red", size=4.6, label=expression(R[0])) +
  scale_y_continuous(name = "Mean effective\nreproduction number",
                     limits=c(0,NA), breaks = seq(0,2.5,0.5)) +
  scale_x_continuous(name = "% of cases with chirping\nsmartphones", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  guides(colour=guide_legend(nrow=2)) +
  coord_fixed(ratio=0.3) +
  theme_base + theme(
    legend.position=c(0.05,0.2),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm"))
  )

#==============================================================================
# Save output
#==============================================================================

plot_scale <- 10
plot_prefix <- "figures/fig2/"
ggsave(filename=paste0(plot_prefix, "sharing_r_eff.png"), plot = sharing_r_eff_plot,
       device = "png", width = plot_scale*2,
       height = plot_scale*2, units = "cm",
       dpi = 320, limitsize = FALSE)
ggsave(filename=paste0(plot_prefix, "phones_r_eff.png"), plot = phone_r_eff_plot,
       device = "png", width = plot_scale*2,
       height = plot_scale*2, units = "cm",
       dpi = 320, limitsize = FALSE)