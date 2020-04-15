#library(ringbp)
library(tidyverse)
library(data.table)
source("R/parameter_sweep.R")

scenarios <- tidyr::expand_grid(
  # Varying parameters
  p_traced = seq(0, 1, 0.2),
  quarantine = c(TRUE, FALSE),
  sero_test = c(TRUE, FALSE),
  backtrace = c(TRUE, FALSE),
  # Set parameters
  run_new = TRUE,
  p_asymptomatic = 0.3, #seq(0,0.5,0.1),
  cap_max_days = 365,
  cap_cases = 5000,
  r0isolated = 0,
  disp.iso = 1,
  disp.com = 0.16,
  delay_group = "SARS",
  delay_shape = 1.651524,
  delay_scale = 4.287786,
  theta = "30%",
  k = 0.7,
  r0community = 2.5,
  num.initial.cases = 20) %>%
  dplyr::mutate(scenario = 1:dplyr::n())

## Set up multicore if using see ?future::plan for details
## Use the workers argument to control the number of cores used.
#future::plan("multiprocess")


## Run parameter sweep
sweep_results_backtracing <- parameter_sweep(scenarios, n_rep = 100,
                                             use_future = FALSE)

## Save sweep results
sweep_results_backtracing_out <- sweep_results_backtracing %>%
  mutate(cases_per_gen = sapply(cases_per_gen, function(x) paste(x, collapse="|")))
write_tsv(sweep_results_backtracing_out, "results/test_backtracing.tsv")

## Calculate control rates
sweep_results_backtracing_out <- suppressMessages(read_tsv("results/test_backtracing.tsv"))
outcomes_backtracing <- sweep_results_backtracing_out %>% 
  group_by(p_traced, quarantine, sero_test, backtrace, sim, final_total_cases) %>%
  filter(between(week, 12, 16)) %>% summarise(trial_cases = sum(weekly_cases)) %>%
  mutate(controlled = (trial_cases == 0) & final_total_cases < 5000)
outcomes_backtracing_final <- outcomes_backtracing %>% 
  group_by(p_traced, quarantine, sero_test, backtrace) %>% 
  summarise(p_controlled = mean(controlled))

## Plot results

facet_names <- c(`FALSE`="No testing", `TRUE`="Perfect testing")
g_backtracing <- ggplot(outcomes_backtracing_final, aes(x=p_traced, y=p_controlled, colour=quarantine)) +
  geom_line(aes(linetype = backtrace)) + geom_point(aes(shape = backtrace), size=2) +
  facet_wrap(~ sero_test, labeller = labeller(sero_test = facet_names)) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2)) +
  scale_colour_discrete(labels = c("No quarantine", "Perfect quarantine")) +
  scale_linetype_discrete(labels = c("No backtracing", "Backtracing")) +
  scale_shape_discrete(labels = c("No backtracing", "Backtracing")) +
  theme_light() + theme(
    legend.position = "right",
    legend.title = element_blank()
  )

ggsave(filename="test_backtracing.png", plot=g_backtracing, device="png", 
       width=15, height=9, units="cm", dpi=320, limitsize=FALSE)
       