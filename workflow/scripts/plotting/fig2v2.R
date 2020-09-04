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

# Universal coverage
univ_plot <- ggplot(univ_data,
                aes(x=p_traced_auto, y=p_controlled,
                    colour = factor(backtrace_distance),
                    linetype = factor(contact_limit_auto, levels=c(Inf, 2)),
                    shape = factor(contact_limit_auto, levels=c(Inf, 2)))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  ggtitle("Automated tracing\n(universal coverage)") +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  scale_linetype_discrete(name=NULL, labels=label_limits) +
  scale_shape_discrete(name=NULL, labels=label_limits) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(order=2),
         shape=guide_legend(order=2)
  ) +
  coord_fixed() +
  theme_base + theme(
    legend.position=c(0.05,0.75),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.spacing.y = unit(0, "mm"),
  )

# Data sharing
sharing_plot <- ggplot(sharing_data,
                        aes(x=p_data_sharing_auto, y=p_controlled, 
                            colour = factor(backtrace_distance))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  ggtitle("Automated tracing\n(partial data sharing)") +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  guides(colour=guide_legend(nrow=2)) +
  coord_fixed() +
  theme_base + theme(
    legend.position=c(0.05,0.85),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm"))
  )

# Smartphone coverage
phone_plot <- ggplot(phone_data,
                            aes(x=p_smartphone_overall, y=p_controlled, 
                                colour = factor(backtrace_distance))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  ggtitle("Automated tracing\n(partial smartphone usage)") +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of cases with chirping\nsmartphones", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  guides(colour=guide_legend(nrow=2)) +
  coord_fixed() +
  theme_base + theme(
    legend.position=c(0.05,0.85),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm"))
  )

# Manual tracing only
manual_plot <- ggplot(manual_data,
                aes(x=p_traced_manual, y=p_controlled,
                    colour = factor(backtrace_distance),
                    linetype = factor(contact_limit_manual),
                    shape = factor(contact_limit_manual))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  ggtitle("Manual tracing\nonly") +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced manually", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  scale_linetype_discrete(name=NULL, labels=label_limits_manual) +
  scale_shape_discrete(name=NULL, labels=label_limits_manual) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(order=2),
         shape=guide_legend(order=2)
  ) +
  coord_fixed() +
  theme_base + theme(
    legend.position=c(0.05,0.75),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.spacing.y = unit(-1, "mm"),
  )

# Combined tracing
combo_plot <- ggplot(combo_data,
                aes(x=p_traced_manual, y=p_controlled,
                    colour = factor(backtrace_distance),
                    linetype = factor(contact_limit_manual),
                    shape = factor(contact_limit_manual))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  ggtitle("Manual + automated\ntracing") +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced manually", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  scale_linetype_discrete(name=NULL, labels=label_limits_manual) +
  scale_shape_discrete(name=NULL, labels=label_limits_manual) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(order=2),
         shape=guide_legend(order=2)
  ) +
  coord_fixed() +
  theme_base + theme(
    legend.position=c(0.05,0.75),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.spacing.y = unit(-1, "mm"),
  )

# R0 plot
r0_plot <- bind_rows(combo_data %>% mutate(trace_type = "Manual + automated"),
                     manual_data %>% mutate(trace_type = "Manual tracing only")) %>%
  filter(backtrace_distance == Inf) %>%
  ggplot(aes(x=p_traced_manual, y=effective_r0_mean,
                    colour = trace_type,
                    linetype = factor(contact_limit_manual),
                    shape = factor(contact_limit_manual))) +
  geom_line() + geom_point(size=2) +
  ggtitle(expression(bold(atop("Effect on"~R[eff],"(forward + reverse tracing)")))) +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  geom_hline(yintercept = combo_data$r0_base[1], linetype = "dotted", size=1,
             colour = "red") +
  annotate("text", x=0.98, y=2.37, hjust=0.5, vjust = 0.5,
           colour = "red", size=4.6, label=expression(R[0])) +
  scale_y_continuous(name = "Avg effective reprod. number",
                     limits=c(0,2.5), breaks = seq(0,2.5,0.5)) +
  scale_x_continuous(name = "% of contacts traced manually", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Set1", name=NULL) +
  scale_linetype_discrete(name=NULL, labels=label_limits_manual) +
  scale_shape_discrete(name=NULL, labels=label_limits_manual) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(order=2),
         shape=guide_legend(order=2)
  ) +
  coord_fixed(ratio=0.4) +
  theme_base + theme(
    legend.position=c(0.01,0.21),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.spacing.y = unit(-1, "mm"),
  )

#==============================================================================
# Combine and label
#==============================================================================

# Configure size
row_height <- 14
width_ratio <- 0.8

# Six-panel version
col_widths <- c(1,1,1)
fig2 <- plot_grid(univ_plot, sharing_plot, phone_plot,
                  manual_plot, combo_plot, r0_plot,
                  labels = "auto", nrow=2, ncol=3,
                  rel_widths = col_widths, align="hv", axis="l",
                  label_size = fontsize_base * fontscale_label,
                  label_fontfamily = titlefont, label_colour = "black")

ggsave(filename="figures/fig2.png", plot = fig2,
       device = "png", width = row_height*sum(col_widths)*width_ratio,
       height = row_height*2, units = "cm",
       dpi = 320, limitsize = FALSE)