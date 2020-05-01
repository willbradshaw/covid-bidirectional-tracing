library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

# Fig. 2a
fig2a_scenario_path <- "saved_data/fig2a_scenario.tsv.gz"
fig2a_scenario_data <- suppressMessages(read_tsv(fig2a_scenario_path))

# Fig. 2b
fig2b_scenario_path <- "saved_data/fig2b_scenario.tsv.gz"
fig2b_scenario_data <- suppressMessages(read_tsv(fig2b_scenario_path))

# Fig. 2c
fig2c_sharing_scenario_path <- "saved_data/fig2c_sharing_scenario.tsv.gz"
fig2c_smartphones_scenario_path <- "saved_data/fig2c_smartphones_scenario.tsv.gz"
fig2c_sharing_scenario_data <- suppressMessages(read_tsv(fig2c_sharing_scenario_path))
fig2c_smartphones_scenario_data <- suppressMessages(read_tsv(fig2c_smartphones_scenario_path))
fig2c_heatmap_scenario_path <- "saved_data/fig2c_heatmap_scenario.tsv.gz"
fig2c_heatmap_scenario_data <- suppressMessages(read_tsv(fig2c_heatmap_scenario_path))

# Fig. 2d
fig2d_scenario_path <- "saved_data/fig2d_scenario.tsv.gz"
fig2d_scenario_data <- suppressMessages(read_tsv(fig2d_scenario_path))

# Fig. 2e
fig2e_scenario_path <- "saved_data/fig2e_scenario.tsv.gz"
fig2e_scenario_data <- suppressMessages(read_tsv(fig2e_scenario_path))

#==============================================================================
# Transform data where needed
#==============================================================================

fig2c_data_sharing_processed <- fig2c_sharing_scenario_data %>%
  select(scenario, p_environmental, p_data_sharing_manual,
         data_limit_auto, backtrace_distance,
         p_controlled, p_controlled_upper, p_controlled_lower,
         p_data_sharing_auto) %>%
  gather(variable, p, -(scenario:p_controlled_lower))
fig2c_data_smartphones_processed <- fig2c_smartphones_scenario_data %>%
  select(scenario, p_environmental, p_data_sharing_manual,
         data_limit_auto, backtrace_distance,
         p_controlled, p_controlled_upper, p_controlled_lower,
         p_smartphone_overall) %>%
  gather(variable, p, -(scenario:p_controlled_lower))
fig2c_data_processed <- bind_rows(fig2c_data_sharing_processed,
                                  fig2c_data_smartphones_processed)

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
  ifelse(as.numeric(x) == Inf, "No limit", paste0(x, "-day limit"))
}

label_limits_manual <- function(x){
  ifelse(as.numeric(x) == Inf, "No limit", paste0(x, "-day manual limit"))
}


label_vars <- function(x){
  ifelse(x=="p_data_sharing_auto", "% of cases sharing data",
         "% of cases with\ntrace-enabled smartphones")
}

#==============================================================================
# Make plots
#==============================================================================

# Fig. 2a
fig2a <- ggplot(fig2a_scenario_data,
                aes(x=p_traced_auto, y=p_controlled,
                    colour = factor(backtrace_distance))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  #ggtitle("Automated tracing\n(universal coverage)") +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2),
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

# Fig. 2b
fig2b <- ggplot(fig2b_scenario_data,
                aes(x=p_traced_auto, y=p_controlled,
                    colour = factor(backtrace_distance),
                    linetype = factor(contact_limit_auto),
                    shape = factor(contact_limit_auto))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  #ggtitle("Automated tracing\n(universal coverage)") +
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

# Fig. 2c (various options)
fig2c_sharing <- ggplot(fig2c_sharing_scenario_data,
                        aes(x=p_data_sharing_auto, y=p_controlled, 
                            colour = factor(backtrace_distance))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of cases with smartphones\nsharing data", breaks = seq(0,1,0.2),
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

fig2c_smartphones <- ggplot(fig2c_smartphones_scenario_data,
                        aes(x=p_smartphone_overall, y=p_controlled, 
                            colour = factor(backtrace_distance))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of cases with trace-enabled\nsmartphones", breaks = seq(0,1,0.2),
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

fig2c_combo <- ggplot(fig2c_data_processed,
                aes(x=p, y=p_controlled,
                    colour = factor(backtrace_distance),
                    linetype = variable, shape = variable)) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "%", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=label_backtrace,
                      name=NULL) +
  scale_linetype_discrete(name = NULL, labels = label_vars) +
  scale_shape_discrete(name = NULL, labels = label_vars) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(order=2),
         shape=FALSE
  ) +
  coord_fixed() +
  theme_base + theme(
    legend.position=c(0.05,0.75),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.spacing.y = unit(-1, "mm"),
  )

fig2c_heatmap <- ggplot(fig2c_heatmap_scenario_data,
                        aes(y=p_smartphone_overall, x=p_data_sharing_auto,
                            fill = p_controlled)) +
  geom_tile() +
  scale_x_continuous(name = "% of cases sharing data",
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_y_continuous(name = "% of cases with\ntrace-enabled smartphones", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_fill_continuous(type = "viridis", name = "% of outbreaks\ncontrolled",
                        labels = function(x) round(x*100)) +
  coord_fixed() +
  theme_base + theme(
    legend.position = "right",
    plot.margin = margin(l=0.5, unit="cm"),
    # axis.title.x = element_text(size = fontsize_base * 1.2),
    # axis.title.y = element_text(size = fontsize_base * 1.2),
  )
  

# Fig. 2d
fig2d <- ggplot(fig2d_scenario_data,
                aes(x=p_traced_manual, y=p_controlled,
                    colour = factor(backtrace_distance),
                    linetype = factor(contact_limit_manual),
                    shape = factor(contact_limit_manual))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
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

# Fig. 2e
fig2e <- ggplot(fig2e_scenario_data %>% filter(p_smartphone_overall == 0),
                aes(x=p_traced_manual, y=p_controlled,
                    colour = factor(backtrace_distance),
                    linetype = factor(contact_limit_manual),
                    shape = factor(contact_limit_manual))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
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

fig2e_7panel <- fig2e +
  scale_x_continuous(name = "% of contacts traced\nmanually", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100))
  

# Fig2f
fig2f <- ggplot(bind_rows(fig2d_scenario_data %>% mutate(trace_type = "Manual + automated"),
                          fig2e_scenario_data %>% mutate(trace_type = "Manual tracing only")) %>%
                  filter(backtrace_distance == Inf),
                aes(x=p_traced_manual, y=effective_r0_mean,
                    colour = trace_type,
                    linetype = factor(contact_limit_manual),
                    shape = factor(contact_limit_manual))) +
  geom_line() + geom_point(size=2) +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  scale_y_continuous(name = "Avg effective reprod. number",#\n(mean ± S.E.M.)",
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

fig2f_6panel <- fig2f +
  scale_y_continuous(name = "Average effective\nreproduction number",
                     limits=c(0,2.5), breaks = seq(0,2.5,0.5)) +
  theme_base + theme(
    legend.position=c(0.01,0.21),
    legend.justification = c("left","center"),
    legend.background = element_rect(fill=alpha("white", 0.5)),
    legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
    legend.spacing.y = unit(-1, "mm"),
    plot.margin = margin(r=4, unit="cm")
  )

#==============================================================================
# Combine and label
#==============================================================================

# Configure size
row_height <- 11.7

# Six-panel version
col_widths <- c(1,1,1.4)
fig2_6panel <- plot_grid(fig2a, fig2b, fig2c_heatmap,
                         fig2d, fig2e, fig2f_6panel,
                         labels = "auto", nrow=2, ncol=3,
                         rel_widths = col_widths, align="hv", axis="l",
                         label_size = fontsize_base * fontscale_label,
                         label_fontfamily = titlefont, label_colour = "black")

ggsave(filename="figures/fig2_6panel.png", plot = fig2_6panel,
       device = "png", width = row_height*sum(col_widths),
       height = row_height*2, units = "cm",
       dpi = 320, limitsize = FALSE)

# Seven-panel version
fig2_row1 <- plot_grid(fig2a + coord_fixed(ratio=0.5),
                       fig2b + coord_fixed(ratio=0.5),
                       labels = c("a", "b"), nrow=1, ncol=2,
                       label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")
fig2_row2 <- plot_grid(fig2c_smartphones, fig2c_sharing, fig2e_7panel,
                       labels = c("c", "d", "e"), nrow=1, ncol=3,
                       label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")
fig2_row3 <- plot_grid(fig2d + coord_fixed(ratio=0.5),
                       fig2f + coord_fixed(ratio=0.2),
                       labels = c("f", "g"), nrow=1, ncol=2,
                       label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")
fig2_7panel <- plot_grid(fig2_row1, fig2_row2, fig2_row3, nrow=3, ncol=1)
ggsave(filename="figures/fig2_7panel.png", plot = fig2_7panel,
       device = "png", width = row_height*3,
       height = row_height*3, units = "cm",
       dpi = 320, limitsize = FALSE)





# 
# ggsave(filename="figures/fig2a.png",
#        plot = fig2a, device="png",
#        width=20, height=15, units="cm", dpi=320, limitsize=FALSE)
# ggsave(filename="figures/fig2b.png",
#        plot = fig2b, device="png",
#        width=20, height=15, units="cm", dpi=320, limitsize=FALSE)
# ggsave(filename="figures/fig2c.png",
#        plot = fig2c, device="png",
#        width=20, height=15, units="cm", dpi=320, limitsize=FALSE)
