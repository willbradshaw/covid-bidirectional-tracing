library(tidyverse)
library(cowplot)

#==============================================================================
# Read in data
#==============================================================================

data_path <- "figures_local/data/si_presym_500_scenario.tsv.gz"
data <- suppressMessages(read_tsv(data_path))

#==============================================================================
# Process data
#==============================================================================

alpha_to_presym <- function(alpha){
  presym <- numeric(length(alpha))
  presym[alpha == 1.38]  <- 0.2
  presym[alpha == 0.73]  <- 0.3
  presym[alpha == 0.325]  <- 0.4
  presym[alpha == 0]    <- 0.5
  presym[alpha == -0.325] <- 0.6
  presym[alpha == -0.73] <- 0.7
  presym[alpha == -1.38] <- 0.8
  return(presym)
}

data_processed <- data %>%
  mutate(trace_type = ifelse(p_traced_auto == 0,
                             ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                             ifelse(p_traced_manual == 0, "Digital only", 
                                    "Manual + digital")),
         p_presymptomatic = alpha_to_presym(generation_alpha),
         backtrace_distance = factor(backtrace_distance, levels = c(Inf,0)))

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

label_rel_r0 <- function(x){
  paste("Rel. R0:", x)
}

label_p_asym <- function(x){
  paste0(round(as.numeric(x)*100), "% asym. carriers")
}

#==============================================================================
# Make plots
#==============================================================================

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

# Control plots

control_plot_90 <- data_processed %>%
  ggplot(aes(x=p_presymptomatic, y=p_controlled, 
             ymin=p_controlled_lower, ymax=p_controlled_upper,
             colour = factor(trace_type, levels=ttype_levels),
             linetype = factor(backtrace_distance),
             shape = factor(backtrace_distance))) +
  geom_errorbar(width=0.05) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% presymptomatic transmission") +
  scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
  scale_linetype_discrete(name=NULL, labels=label_backtrace) +
  scale_shape_discrete(name=NULL, labels=label_backtrace) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(nrow=2, order=2),
         shape=guide_legend(nrow=2, order=2)
  ) +
  coord_fixed(ratio=0.9*0.7) +
  theme_base


# R_eff plots
r_eff_plot_90 <- data_processed %>%
  ggplot(aes(x=p_presymptomatic, y=effective_r0_mean, 
             colour = factor(trace_type, levels=ttype_levels),
             linetype = factor(backtrace_distance),
             shape = factor(backtrace_distance),
             ymin = effective_r0_mean-effective_r0_sd,
             ymax = effective_r0_mean+effective_r0_sd)) +
  geom_line() + geom_point(size=2) +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  scale_y_continuous(name = "Avg effective reprod. number",
                     limits=c(0,NA), breaks = seq(0,10,0.5)) +
  scale_x_continuous(name = "% presymptomatic transmission") +
  scale_colour_brewer(type = "div", palette = "Set1", name="Trace type") +
  scale_linetype_discrete(name=NULL, labels=label_backtrace) +
  scale_shape_discrete(name=NULL, labels=label_backtrace) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(nrow=2, order=2),
         shape=guide_legend(nrow=2, order=2)
  ) +
  coord_fixed(ratio=0.3*0.7) +
  theme_base

out_plot <- plot_grid(control_plot_90, r_eff_plot_90,
                      labels = "auto", nrow = 1, ncol = 2, align = "hv",
                      axis = "l", label_size = fontsize_base * fontscale_label,
                      label_fontfamily = titlefont, label_colour = "black")


#==============================================================================
# Save output
#==============================================================================

plot_scale <- 14
plot_prefix <- "figures_local/img/si_presymptomatic"
ggsave(filename=paste0(plot_prefix, ".svg"), plot = out_plot,
       device = "svg", width = plot_scale*2.6,
       height = plot_scale, units = "cm",
       dpi = 320, limitsize = FALSE)