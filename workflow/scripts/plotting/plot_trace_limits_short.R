library(tidyverse)

# Read in data
scenario_0gen_path <- "saved_data/trace_limit_0gen_short_scenario.tsv.gz"
scenario_4gen_path <- "saved_data/trace_limit_4gen_short_scenario.tsv.gz"
scenario_0gen_data <- suppressMessages(read_tsv(scenario_0gen_path))
scenario_4gen_data <- suppressMessages(read_tsv(scenario_4gen_path))
scenario_data <- bind_rows(scenario_0gen_data, scenario_4gen_data)

# Make plot
label_limit <- function(x) {
  x1 <- ifelse(as.numeric(x) == Inf, "No",
               paste0(x, "-day"))
  return(paste(x1, "limit"))
}
label_gen <- function(x) paste0(x, "-generation delay")
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
  facet_grid(rollout_delay_gen~contact_limit_auto,
             labeller = labeller(contact_limit_auto = label_limit,
                                 rollout_delay_gen = label_gen)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  scale_fill_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  theme_bw() + theme(
    legend.position = "right",
    panel.spacing.x = unit(0.3, "cm")
  )

ggsave(filename="saved_data/plot_trace_limits_short.png",
       plot = g, device="png",
       width=22, height=12, units="cm", dpi=320, limitsize=FALSE)
