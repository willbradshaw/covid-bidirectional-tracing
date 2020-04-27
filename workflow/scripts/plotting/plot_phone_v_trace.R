library(tidyverse)
library(cowplot)

# Read in data
scenario_path_bt <- "saved_data/phone_v_trace_bt_scenario.tsv.gz"
scenario_path_nobt <- "saved_data/phone_v_trace_nobt_scenario.tsv.gz"
scenario_data_bt <- suppressMessages(read_tsv(scenario_path_bt))
scenario_data_nobt <- suppressMessages(read_tsv(scenario_path_nobt))
scenario_data <- bind_rows(scenario_data_bt, scenario_data_nobt)
scenario_data_cols <- scenario_data %>% 
  group_by(p_smartphone_overall, p_traced_auto, backtrace_distance, p_controlled,
           p_controlled_lower, p_controlled_upper) %>%
  summarise %>% filter(p_traced_auto >= 0.4, p_smartphone_overall >= 0.4)

# Calculate control ratios
ratio_tab <- expand_grid(x = seq(0,1,0.2), y=seq(0,1,0.2),
                         bt = c(0,Inf)) %>%
  mutate(x = round(x,1), y = round(y,1)) %>%
  mutate(p_st = sapply(1:nrow(.), function(m) scenario_data %>% filter(p_smartphone_overall == x[m],
                                                                     p_traced_auto == y[m],
                                                                     backtrace_distance == bt[m]) %>% pull(p_controlled)),
         p_ts = sapply(1:nrow(.), function(m) scenario_data %>% filter(p_smartphone_overall == y[m],
                                                                     p_traced_auto == x[m],
                                                                     backtrace_distance == bt[m]) %>% pull(p_controlled))
  ) %>% filter(p_st != 0 , p_ts != 0) %>%
  mutate(p_ratio = p_ts/p_st) %>% group_by(bt)

# Make plots
g0 <- ggplot(scenario_data_cols, aes(x=p_traced_auto, y=p_controlled,
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
h0 <- ggplot(scenario_data_cols, aes(x=p_smartphone_overall, y=p_controlled,
                                     colour = factor(backtrace_distance),
                                     fill = factor(backtrace_distance))) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  facet_grid(~p_traced_auto, labeller = labeller(
    p_traced_auto = function(x) paste0(as.numeric(x)*100, "% contacts_traced")
  )) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% smartphone coverage", 
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

# label_limit <- function(x) {
#   x1 <- ifelse(as.numeric(x) == Inf, "No",
#                paste0(x, "-day"))
#   return(paste(x1, "limit"))
# }
# label_gen <- function(x) paste0(x, "-generation delay")
g <- ggplot(scenario_data, aes(x=p_traced_auto, y=p_controlled,
                               colour = factor(p_smartphone_overall),
                               fill = factor(p_smartphone_overall))) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  facet_grid(~backtrace_distance) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of smartphone-linked contacts traced", 
                     breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                      name="% of individuals with\ntrace-enabled smartphones") +
  scale_fill_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                      name="% of individuals with\ntrace-enabled smartphones") +
  theme_bw() + theme(
    legend.position = "bottom",
    panel.spacing.x = unit(0.3, "cm")
  )
h <- ggplot(scenario_data, aes(x=p_smartphone_overall, y=p_controlled,
                               colour = factor(p_traced_auto),
                               fill = factor(p_traced_auto))) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  facet_grid(~backtrace_distance) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of individuals with trace-enabled smartphones", 
                     breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                      name="% of smartphone-linked\ncontacts traced") +
  scale_fill_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                    name="% of smartphone-linked\ncontacts traced") +
  theme_bw() + theme(
    legend.position = "bottom",
    panel.spacing.x = unit(0.3, "cm")
  )

i <- plot_grid(g,h)

ggsave(filename="saved_data/plot_phone_v_trace.png",
       plot = i, device="png",
       width=22, height=12, units="cm", dpi=320, limitsize=FALSE)
