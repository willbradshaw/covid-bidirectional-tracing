#library(ringbp)
library(tidyverse)
library(data.table)
source("R/wrappers.R")

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------

n_iterations <- 100
test_name <- "backtracing_2"

var_params <- c("p_traced", "sero_test", "backtrace")
max_weeks_run <- 52
max_weeks_test <- 12
cap_cases <- 5000
cap_max_generations <- 100

scenarios <- tidyr::expand_grid(
  # Varying parameters
  p_traced = seq(0, 1, 0.2),
  sero_test = c(TRUE, FALSE),
  backtrace = c(TRUE, FALSE),
  # Set parameters
  quarantine = FALSE, # Currently redundant with sero_test
  r0_base = 2.5,
  r0_asymptomatic = 2.5,
  p_asymptomatic = 0.3,
  p_isolation = 1, # Perfect isolation
  dispersion = 0.16,
  cap_max_weeks = max_weeks_run,
  cap_max_generations = cap_max_generations,
  cap_cases = cap_cases,
  n_initial_cases = 20,
  delay_shape = 1.651524, # Hellewell short delay case
  delay_scale = 4.287786,
  incubation_shape = 2.322737, # Hellewell default case
  incubation_scale = 6.492272,
  generation_omega = 2, # Hellewell default case
  generation_k = 0.7, # Hellewell 30% pre-symptomatic transmission case
  ) %>% dplyr::mutate(scenario = 1:dplyr::n())

#------------------------------------------------------------------------------
# Run sweep and save results
#------------------------------------------------------------------------------

assign(paste0("sweep_results_", test_name),
       parameter_sweep(scenarios = scenarios, n_iterations = n_iterations,
                       show_progress = FALSE, use_future = FALSE,
                       report = TRUE))

write_tsv(get(paste0("sweep_results_", test_name)),
          paste0("results/test_", test_name, ".tsv"))

#------------------------------------------------------------------------------
# Read results and compute control rates
#------------------------------------------------------------------------------

# Read results
assign(paste0("sweep_results_", test_name, "_out"),
       suppressMessages(read_tsv(paste0("results/test_", test_name, ".tsv"))))

# Determine outcome for each iteration
assign(paste0("outcomes_", test_name),
       get(paste0("sweep_results_", test_name, "_out")) %>%
         group_by_at(vars(one_of(var_params))) %>%
         group_by(sim, final_total_cases, add = TRUE) %>%
         filter(week >= max_weeks_test) %>%
         summarise(trial_cases = sum(weekly_cases)) %>%
         mutate(controlled = (trial_cases == 0) &
                  (final_total_cases < cap_cases)))

# Calculate control rate for each scenario
assign(paste0("outcomes_", test_name, "_final"),
       get(paste0("outcomes_", test_name)) %>%
         group_by_at(vars(one_of(var_params))) %>%
         summarise(p_controlled = mean(controlled)))
#          
# # TODO: Write processed results
# #------------------------------------------------------------------------------
# # Plot results
# #------------------------------------------------------------------------------
# # TODO: Follow up on generation counts
# 
g_backtracing_2 <- ggplot(outcomes_backtracing_2_final,
                          aes(x=p_traced, y=p_controlled, colour = sero_test)) +
  geom_line(aes(linetype = backtrace)) + geom_point(aes(shape = backtrace), size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2)) +
  scale_colour_discrete(labels = c("No testing", "Perfect testing")) +
  scale_linetype_discrete(labels = c("No backtracing", "Backtracing")) +
  scale_shape_discrete(labels = c("No backtracing", "Backtracing")) +
  theme_light() + theme(
    legend.position = "right",
    legend.title = element_blank()
  )
# 
# # ggsave(filename="test_backtracing.png", plot=g_backtracing, device="png", 
# #        width=15, height=9, units="cm", dpi=320, limitsize=FALSE)
#        
