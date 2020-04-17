#library(ringbp)
library(tidyverse)
library(data.table)
source("R/wrappers.R")

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------

n_iterations <- 100
test_name <- "backtracing_n_big_gen"

var_params <- c("backtrace_distance", "p_traced", "rollout_delay_generations",
                "p_asymptomatic")
cap_max_weeks <- 52
cap_cases <- 5000
cap_max_generations <- 100

scenarios <- tidyr::expand_grid(
  # Varying parameters
  p_traced = seq(0.4, 1, 0.2),
  backtrace_distance = c(0,1,Inf),
  p_asymptomatic = c(0.25, 0.5),
  rollout_delay_generations = 0:4,
  # Set parameters
  rollout_delay_days = 0,
  n_initial_cases = 20,
  sero_test = TRUE,
  quarantine = FALSE, # Currently redundant with sero_test
  r0_base = 2.5,
  r0_asymptomatic = 2.5,
  p_isolation = 1, # Perfect isolation
  dispersion = 0.16,
  cap_max_weeks = cap_max_weeks,
  cap_max_generations = cap_max_generations,
  cap_cases = cap_cases,
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
# Count additions at each generation
#------------------------------------------------------------------------------

assign(paste0("gen_counts_", test_name),
       get(paste0("sweep_results_", test_name)) %>%
         group_by_at(vars(one_of(var_params))) %>%
         group_by(sim, cases_per_gen, add = TRUE) %>%
         summarise %>% mutate(cases_per_gen = str_split(cases_per_gen, "\\|")) %>%
         unnest(cases_per_gen) %>% mutate(generation = row_number()-1) %>%
         mutate(cases_per_gen = as.integer(cases_per_gen)))

#------------------------------------------------------------------------------
# Read results and compute control rates
#------------------------------------------------------------------------------

# Read results
assign(paste0("sweep_results_", test_name),
       suppressMessages(read_tsv(paste0("results/test_", test_name, ".tsv"))))

# Determine outcome for each iteration
assign(paste0("outcomes_", test_name),
       get(paste0("sweep_results_", test_name)) %>%
         group_by_at(vars(one_of(var_params))) %>%
         group_by(sim, final_total_cases, outbreak_extinct, add = TRUE) %>%
         summarise() %>% 
         mutate(controlled = outbreak_extinct &
                  (final_total_cases < cap_cases)))

# Calculate control rate for each scenario
assign(paste0("outcomes_", test_name, "_final"),
       get(paste0("outcomes_", test_name)) %>%
         group_by_at(vars(one_of(var_params))) %>%
         summarise(p_controlled = mean(controlled)))

# TODO: Write processed results
#------------------------------------------------------------------------------
# Plot results
#------------------------------------------------------------------------------
# TODO: Follow up on generation counts

g_backtracing_n_big_gen <- ggplot(outcomes_backtracing_n_big_gen_final,
                              aes(x=p_traced, y=p_controlled,
                                  colour = factor(backtrace_distance))) +
  geom_line() +
  geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2)) +
  scale_x_continuous(name = "% of contacts traced", breaks = seq(0,1,0.2)) +
  facet_grid(paste0(p_asymptomatic*100,"% asymptomatic")~paste(rollout_delay_generations, "gen.")) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=c("0","1","∞"),
                      name="Max. backtrace\ndistance") +
  theme_bw() + theme(
    legend.position = "right",
    #legend.title = element_blank()
  )

ggsave(filename="figures/test_backtracing_n_big_gen.png",
       plot = g_backtracing_n_big_gen, device="png", 
       width=22, height=12, units="cm", dpi=320, limitsize=FALSE)
