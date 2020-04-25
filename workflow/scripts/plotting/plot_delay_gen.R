library(tidyverse)

# Read in data
scenario_path <- "saved_data/delay_gen_scenario.tsv.gz"
scenario_data <- suppressMessages(read_tsv(scenario_path))

# Make plot
label_asym <- function(x) {
  paste0(round(as.numeric(x)*100), "% asymptomatic")
}
label_delay <- function(y){
  paste0(round(as.numeric(y)), " gen.")
}
g <- ggplot(scenario_data, aes(x=p_traced_auto, y=p_controlled,
                               colour = factor(backtrace_distance),
                               fill = factor(backtrace_distance))) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  facet_grid(p_asymptomatic~rollout_delay_gen,
             labeller = labeller(
               p_asymptomatic = label_asym,
               rollout_delay_gen = label_delay
               )) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  scale_fill_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  theme_bw() + theme(
    legend.position = "right",
    panel.spacing.x = unit(0.4, "cm")
  )

ggsave(filename="saved_data/test_delay_gen.png",
        plot = g, device="png",
        width=22, height=12, units="cm", dpi=320, limitsize=FALSE)
