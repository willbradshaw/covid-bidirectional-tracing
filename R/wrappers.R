#' Run a specified number of simulations with identical parameters
#' @author Joel Hellewell (original version), Will Bradshaw (modified version)

source("R/outbreak.R")

library(tidyverse)
library(data.table)

scenario_sim <- function(n_iterations = NULL, r0_base = NULL, dispersion = NULL,
                         cap_max_weeks = NULL, cap_max_generations = NULL,
                         r0_asymptomatic = NULL, p_traced = NULL, p_asymptomatic = NULL,
                         backtrace_distance = NULL, quarantine = NULL, sero_test = NULL,
                         p_isolation = NULL, delay_shape = NULL,
                         delay_scale = NULL, incubation_shape = NULL,
                         incubation_scale = NULL, generation_omega = NULL,
                         generation_k = NULL, report = NULL,
                         n_initial_cases = NULL, cap_cases = NULL){
    #' Run a specified number of simulations with identical parameters
    if (!is.null(report) & !is.na(report)){cat(report, ": ", sep="")}
    # Run each iteration
    iter_out <- purrr::map(.x = 1:n_iterations,
                           ~ outbreak_model(n_initial_cases = n_initial_cases,
                                            r0_base = r0_base, dispersion = dispersion,
                                            cap_max_weeks = cap_max_weeks,
                                            cap_max_generations = cap_max_generations,
                                            r0_asymptomatic = r0_asymptomatic,
                                            p_traced = p_traced, p_asymptomatic = p_asymptomatic,
                                            backtrace_distance = backtrace_distance, 
                                            quarantine = quarantine,
                                            sero_test = sero_test, p_isolation = p_isolation,
                                            delay_shape = delay_shape, delay_scale = delay_scale,
                                            incubation_shape = incubation_shape,
                                            incubation_scale = incubation_scale,
                                            generation_omega = generation_omega,
                                            generation_k = generation_k, cap_cases = cap_cases))
    # Label and concatenate
    results <- lapply(1:n_iterations, function(n) iter_out[[n]][, sim:= n]) %>%
        data.table::rbindlist(fill = TRUE)
  if (!is.null(report) & !is.na(report)){cat(timestamp(quiet=TRUE), "\n")}
  return(results)
}

parameter_sweep <- function(scenarios = NULL, n_iterations = NULL,
                            show_progress = NULL, use_future = NULL,
                            report = NULL){
    #' Run one set of simulations for each scenario in a table of scenarios
    # Define sweep function
    safe_sim_fn <- scenario_sim #purrr::safely(scenario_sim) # TODO: Figure out how to capture errors
    map_fn <- ifelse(use_future,
                     function(...) furrr::future_map(..., .progress = show_progress),
                     purrr::map) # TODO: Debug future version
    # Nest scenarios into sub-tables
    scenario_data <- scenarios %>% mutate(report = scenario) %>%
        dplyr::group_by(scenario) %>%
        tidyr::nest() %>%
        dplyr::ungroup() %>%
        ##Randomise the order of scenarios - helps share the load across cores
        {if (use_future) dplyr::sample_frac(., size = 1, replace = FALSE) else .}
    # Perform sweep
    scenario_sims <- dplyr::mutate(scenario_data, sims = map_fn(data, ~ safe_sim_fn(
        n_iterations = n_iterations, r0_base = .$r0_base, dispersion = .$dispersion,
        cap_max_weeks = .$cap_max_weeks, cap_max_generations = .$cap_max_generations,
        r0_asymptomatic = .$r0_asymptomatic, p_traced = .$p_traced,
        p_asymptomatic = .$p_asymptomatic, n_initial_cases = .$n_initial_cases,
        backtrace_distance = .$backtrace_distance, quarantine = .$quarantine, 
        sero_test = .$sero_test,
        p_isolation = .$p_isolation, delay_shape = .$delay_shape, cap_cases = .$cap_cases,
        delay_scale = .$delay_scale, incubation_shape = .$incubation_shape,
        incubation_scale = .$incubation_scale, generation_omega = .$generation_omega,
        generation_k = .$generation_k, report = ifelse(report, .$report, NA))))
    # Unnest for analysis
    scenario_sims_out <- scenario_sims %>% tidyr::unnest(c("data", "sims"))
    # TODO: Debug future case, set up multiprocessing
    return(scenario_sims_out)
}
