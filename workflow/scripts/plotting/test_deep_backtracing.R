library(tidyverse)

# Read in data
scenario_path <- "output_files/data/reprod_scenario.tsv.gz"
#scenario_file <- gzfile(scenario_path, "rb")
scenario_data <- suppressMessages(read_tsv(scenario_path))
#close(scenario_file)

# Make plot
g_deep_backtracing <- ggplot(scenario_data,
                                  aes(x=p_traced_auto, y=p_controlled,
                                      colour = factor(backtrace_distance))) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper,
                  fill = factor(backtrace_distance)), alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  facet_grid(paste0(p_asymptomatic*100,"% asymptomatic")~paste(rollout_delay_gen, "gen.")) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  scale_fill_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  theme_bw() + theme(
    legend.position = "right",
    panel.spacing.x = unit(0.3, "cm")
  )

ggsave(filename="output_files/figures/test_deep_backtracing.png",
        plot = g_deep_backtracing, device="png", 
        width=22, height=12, units="cm", dpi=320, limitsize=FALSE)
