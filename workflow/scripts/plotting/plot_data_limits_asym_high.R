library(tidyverse)

# Read in data
scenario_path <- "saved_data/data_limits_asym_high_scenario.tsv.gz"
scenario_data <- suppressMessages(read_tsv(scenario_path))

# Make plot
label_limit <- function(x) {
  x1 <- ifelse(as.numeric(x) == Inf, "No",
               paste0(x, "-day"))
  return(paste(x1, "limit"))
}
label_gen <- function(x) paste0(x, "-generation delay")
g <- ggplot(scenario_data, aes(x=p_traced_auto, y=p_controlled,
                               colour = factor(backtrace_distance))) +
  geom_errorbar(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
                width=0.05) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  facet_grid(rollout_delay_gen~data_limit_auto,
             labeller = labeller(data_limit_auto = label_limit,
                                 rollout_delay_gen = label_gen)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  theme_bw() + theme(
    legend.position = "right",
    panel.spacing.x = unit(0.3, "cm")
  )

ggsave(filename="saved_data/plot_data_limits_asym_high.png",
        plot = g, device="png",
        width=22, height=12, units="cm", dpi=320, limitsize=FALSE)
