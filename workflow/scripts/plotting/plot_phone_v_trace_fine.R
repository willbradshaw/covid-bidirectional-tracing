library(tidyverse)
library(cowplot)

# Read in data
scenario_path_bt <- "saved_data/phone_v_trace_bt_scenario.tsv.gz"
scenario_path_nobt <- "saved_data/phone_v_trace_nobt_scenario.tsv.gz"
scenario_path_pfine <- "saved_data/phone_v_trace_pfine_scenario.tsv.gz"
scenario_data_bt <- suppressMessages(read_tsv(scenario_path_bt))
scenario_data_nobt <- suppressMessages(read_tsv(scenario_path_nobt))
scenario_data_pfine <- suppressMessages(read_tsv(scenario_path_pfine))

scenario_data <- bind_rows(scenario_data_bt, scenario_data_nobt,
                           scenario_data_pfine) %>%
  filter(p_traced_auto >= 0.4, p_smartphone_overall >= 0.6)
scenario_data_cols <- scenario_data %>% 
  group_by(p_smartphone_overall, p_traced_auto, backtrace_distance, p_controlled,
           p_controlled_lower, p_controlled_upper) %>% summarise

# Make plots
g <- ggplot(scenario_data_cols, aes(x=p_traced_auto, y=p_controlled,
                                     colour = factor(backtrace_distance),
                                     fill = factor(backtrace_distance))) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  facet_grid(~p_smartphone_overall, labeller = labeller(
    p_smartphone_overall = function(x) paste0(as.numeric(x)*100, "% smartphone cov.")
  )) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of smartphone-linked contacts traced", 
                     breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                      name="Max. backtrace\ndistance") +
  scale_fill_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                    name="Max. backtrace\ndistance") +
  theme_bw() + theme(
    legend.position = "bottom",
    panel.spacing.x = unit(0.3, "cm")
  )

ggsave(filename="saved_data/plot_phone_v_trace_fine.png",
       plot = g, device="png",
       width=20, height=12, units="cm", dpi=320, limitsize=FALSE)
