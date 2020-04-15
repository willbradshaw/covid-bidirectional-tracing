#' Sweep across parameters
#' @author Sam Abbott

source("R/scenario_sim.R")

parameter_sweep <- function(scenarios = NULL, n_rep = 1,
                            show_progress = TRUE, use_future = TRUE) {

  safe_sim_fn <- ifelse(use_future, purrr::safely(scenario_sim), scenario_sim)
  
  scenario_data <- scenarios %>% mutate(scenario_report = scenario) %>%
    dplyr::group_by(scenario) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    ##Randomise the order of scenarios - helps share the load across cores
    {if (use_future) dplyr::sample_frac(., size = 1, replace = FALSE) else .}
  
  map_fn <- ifelse(use_future, function(...) furrr::future_map(..., .progress = show_progress)[[1]],
                   purrr::map)
  
  scenario_sims <- dplyr::mutate(scenario_data, sims = map_fn(data, ~ safe_sim_fn(
    n.sim = n_rep, p_asymptomatic = .$p_asymptomatic, p_traced = .$p_traced,
    cap_max_days = .$cap_max_days, cap_cases = .$cap_cases, r0isolated = .$r0isolated,
    r0community = .$r0community, disp.iso = .$disp.iso, disp.com = .$disp.com,
    delay_shape = .$delay_shape, delay_scale = .$delay_scale, k = .$k,
    num.initial.cases = .$num.initial.cases, quarantine = .$quarantine,
    backtrace = .$backtrace, sero_test = .$sero_test, run_new = .$run_new,
    report = .$scenario_report)))
  
  if (use_future){
    scenario_sims_out <- scenario_sims %>% tidyr::unnest(cols = "data")
  } else {
    scenario_sims_out <- scenario_sims %>% tidyr::unnest(cols = c("data", "sims"))
  }
  
  return(scenario_sims_out)
}
