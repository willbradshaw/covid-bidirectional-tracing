#library(ringbp)
library(tidyverse)
library(data.table)
source("R/parameter_sweep.R")

scenarios_overlap2 <- tidyr::expand_grid(
  # Varying parameters
  p_asymptomatic = seq(0,0.3,0.1),
  p_traced = seq(0, 1, 0.2),
  run_new = c(TRUE, FALSE),
  # Set parameters
  quarantine = FALSE,
  sero_test = FALSE, #c(TRUE, FALSE),
  backtrace = FALSE, #c(TRUE, FALSE),
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
sweep_results_overlap2 <- parameter_sweep(scenarios_overlap2, n_rep = 100,
                                 use_future = FALSE)

## Save sweep results
sweep_results_overlap2_out <- sweep_results_overlap2 %>%
  mutate(cases_per_gen = sapply(cases_per_gen, function(x) paste(x, collapse="|")))
write_tsv(sweep_results_overlap2_out, "results/test_overlap2.tsv")

## Calculate control rates
sweep_results_overlap2_out <- suppressMessages(read_tsv("results/test_overlap2.tsv"))
outcomes_overlap2 <- sweep_results_overlap2_out %>% 
  group_by(p_traced, p_asymptomatic, run_new, sim, final_total_cases) %>%
  filter(between(week, 12, 16)) %>% summarise(trial_cases = sum(weekly_cases)) %>%
  mutate(controlled = (trial_cases == 0) & final_total_cases < 5000)
outcomes_overlap2_final <- outcomes_overlap2 %>% 
  group_by(p_traced, p_asymptomatic, run_new) %>% 
  summarise(p_controlled = mean(controlled))

## Plot results
g_overlap2 <- ggplot(outcomes_overlap2_final, 
                     aes(x=p_traced, y=p_controlled, colour=factor(p_asymptomatic))) +
  geom_line(aes(linetype = run_new)) + geom_point(aes(shape = run_new), size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2)) +
  scale_colour_brewer(name = "% Asymptomatic carriers",
                      type="qual", palette="Dark2") +
  scale_linetype_discrete(name = "Codebase", labels = c("Original code", "Modified code")) +
  scale_shape_discrete(name = "Codebase", labels = c("Original code", "Modified code")) +
  theme_light() + theme(
    legend.position = "right"
  )
ggsave(filename="test_overlap2.png", plot=g_overlap2, device="png", 
       width=15, height=9, units="cm", dpi=320, limitsize=FALSE)
