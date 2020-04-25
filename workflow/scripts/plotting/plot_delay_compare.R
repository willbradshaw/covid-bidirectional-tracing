library(tidyverse)

# Read in and process

# Generation delay
scenario_path_gen <- "saved_data/delay_gen_scenario.tsv.gz"
scenario_data_gen <- suppressMessages(read_tsv(scenario_path_gen)) %>%
  rename(generations = rollout_delay_gen) %>%
  select(p_traced_auto, backtrace_distance, p_asymptomatic,
         generations, p_controlled,
         p_controlled_upper, p_controlled_lower) %>%
  gather(delay_type, delay, -p_traced_auto, -backtrace_distance,
         -p_asymptomatic, -matches("p_controlled"))

# Time delay
scenario_path_days <- "saved_data/day_delay_scenario.tsv.gz"
scenario_data_days <- suppressMessages(read_tsv(scenario_path_days))
scenario_data_weeks <- scenario_data_days %>% select(-rollout_delay_gen) %>%
  mutate(weeks = round(rollout_delay_days/7)) %>%
  select(p_traced_auto, backtrace_distance, p_asymptomatic,
         weeks, p_controlled,
         p_controlled_upper, p_controlled_lower) %>%
  gather(delay_type, delay, -p_traced_auto, -backtrace_distance,
         -p_asymptomatic, -matches("p_controlled"))

# Combine
scenario_data <- bind_rows(scenario_data_gen, scenario_data_weeks)

# Make plot
label_asym <- function(x) {
  paste0(round(as.numeric(x)*100), "% asymptomatic")
}
label_delay <- function(y){
  paste0(round(as.numeric(y)), " gen.")
}
g_asym_25 <- scenario_data %>% filter(p_asymptomatic == 0.25) %>%
  ggplot(aes(x=p_traced_auto, y=p_controlled,
             colour = delay_type, fill = delay_type)) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  facet_grid(backtrace_distance ~ delay, labeller = labeller(
    backtrace_distance = function(x) paste("BT distance:", x),
    delay = function(x) paste("delay:", x)
  )) + ggtitle("25% asymptomatic carriers") +
  scale_colour_brewer(type = "div", palette = "Set1",
                      name="Delay type") +
  scale_fill_brewer(type = "div", palette = "Set1",
                      name="Delay type") +
  theme_bw() + theme(
    legend.position = "right",
    panel.spacing.x = unit(0.4, "cm")
  )


g_asym_50 <- scenario_data %>% filter(p_asymptomatic == 0.5) %>%
  ggplot(aes(x=p_traced_auto, y=p_controlled,
             colour = delay_type, fill = delay_type)) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  facet_grid(backtrace_distance ~ delay, labeller = labeller(
    backtrace_distance = function(x) paste("BT distance:", x),
    delay = function(x) paste("delay:", x)
  )) + ggtitle("50% asymptomatic carriers") +
  scale_colour_brewer(type = "div", palette = "Set1",
                      name="Delay type") +
  scale_fill_brewer(type = "div", palette = "Set1",
                    name="Delay type") +
  theme_bw() + theme(
    legend.position = "right",
    panel.spacing.x = unit(0.4, "cm")
  )

g_both <- plot_grid(g_asym_25, g_asym_50, ncol=1)

ggsave(filename="saved_data/plot_delay_compare.png",
        plot = g_both, device="png",
        width=18, height=20, units="cm", dpi=320, limitsize=FALSE)
