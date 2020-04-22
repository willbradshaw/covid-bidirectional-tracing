#' Run a specified number of simulations with identical parameters
#' @author Joel Hellewell (original version), Will Bradshaw (modified version)

source("R/outbreak.R")

library(tidyverse)
library(data.table)

#------------------------------------------------------------------------------
# Run parameter sweep
#------------------------------------------------------------------------------

parameter_sweep <- function(scenarios = NULL, n_iterations = NULL,
                            parallel = NULL, show_progress = NULL,
                            report = FALSE){
  #' Run one set of simulations for each scenario in a table of scenarios
  # Nest scenarios into sub-tables
  scenario_data <- scenarios %>% mutate(report = scenario) %>%
    dplyr::group_by(scenario) %>%
    tidyr::nest() %>%
    dplyr::ungroup()
  if (parallel) future::plan(future::multiprocess)
  map_fn <- ifelse(parallel,
                   function(...) furrr::future_map(..., .progress = show_progress),
                   purrr::map) # TODO: Debug future version
  scenario_sims <- dplyr::mutate(scenario_data, sims = map_fn(data, ~ scenario_sim(
    n_iterations = n_iterations, dispersion = .$dispersion, r0_base = .$r0_base, 
    r0_asymptomatic = .$r0_asymptomatic, p_asymptomatic = .$p_asymptomatic,
    generation_omega = .$generation_omega, generation_alpha = .$generation_alpha,
    recovery_quantile = .$recovery_quantile, incubation_time = .$incubation_time,
    test_time = .$test_time, trace_time_auto = .$trace_time_auto,
    trace_time_manual = .$trace_time_manual, delay_time = .$delay_time,
    n_initial_cases = .$n_initial_cases, test_sensitivity = .$test_sensitivity,
    test_serological = .$test_serological, p_ident_sym = .$p_ident_sym,
    p_smartphone_overall = .$p_smartphone_overall, p_smartphone_link = .$p_smartphone_link,
    trace_neg_symptomatic = .$trace_neg_symptomatic, p_traced_auto = .$p_traced_auto,
    p_traced_manual = .$p_traced_manual, data_limit_auto = .$data_limit_auto,
    data_limit_manual = .$data_limit_manual, 
    contact_limit_auto_asym = .$contact_limit_auto_asym,
    contact_limit_auto_sym = .$contact_limit_auto_sym,
    contact_limit_manual_asym = .$contact_limit_manual_asym,
    contact_limit_manual_sym = .$contact_limit_manual_sym,
    rollout_delay_gen = .$rollout_delay_gen, rollout_delay_days = .$rollout_delay_days,
    p_blocked_isolation = .$p_blocked_isolation,
    p_blocked_quarantine = .$p_blocked_quarantine,
    cap_max_generations = .$cap_max_generations,
    cap_max_weeks = .$cap_max_weeks, cap_cases = .$cap_cases,
    backtrace_distance = .$backtrace_distance,
    p_environmental = .$p_environmental, report = ifelse(report, .$report, NA)))) %>%
    tidyr::unnest(c("data", "sims"))
  return(scenario_sims)
}

#------------------------------------------------------------------------------
# Summarise data
#------------------------------------------------------------------------------

summarise_by_week <- function(case_data, group_vars){
  #' Convert a database of case reports into one of weekly case counts
  group_vars <- c("scenario", "run", group_vars[! group_vars %in% c("scenario", "run")])
  data_out <- case_data %>% as.data.table %>% copy %>%
    .[, week := floor(exposure/7)] %>%
    .[, .(cases = .N), by = c("week", group_vars)] %>%
    .[order(week), cumulative_cases := cumsum(cases), by = group_vars]
  return(data_out)
}

summarise_by_generation <- function(case_data, group_vars){
  #' Convert a database of case reports into one of per-generation case counts
  group_vars <- c("scenario", "run", group_vars[! group_vars %in% c("scenario", "run")])
  data_out <- case_data %>% as.data.table %>% copy %>%
    .[, .(cases = .N), by = c("generation", group_vars)] %>%
    .[order(generation), cumulative_cases := cumsum(cases), by=group_vars] %>%
    .[order(generation), effective_r0 := lead(cases)/cases, by=group_vars]
  return(data_out)
}

summarise_by_run <- function(case_data, group_vars){
  #' Convert a database of case reports into one of per-run outcomes
  group_vars <- c("scenario", "run", group_vars[! group_vars %in% c("scenario", "run")])
  # Get per-generation summaries
  db_gens  <- summarise_by_generation(case_data, group_vars)
  # Get key summary statistics
  case_stats <- case_data %>% as.data.table %>% copy %>%
    .[, .(last_exposure_days = max(exposure),
          extinct = all(processed)), by = group_vars] %>%
    .[, `:=`(last_exposure_weeks = last_exposure_days/7,
             last_exposure_days = NULL)]
  gen_stats <- db_gens %>% copy %>%
    .[, .(avg_effective_r0 = sum(effective_r0*cases, na.rm = TRUE)/
            sum(cases + (effective_r0*0), na.rm = TRUE), # Weight each generation's R0 by the number of cases in that generation
          total_cases = sum(cases),
          max_generation = max(generation)), by = group_vars]
  run_stats <- case_stats[gen_stats, on=group_vars] %>%
    .[, controlled := extinct & 
        max_generation < cap_max_generations &
        last_exposure_weeks < cap_max_weeks &
        total_cases < cap_cases]
  return(run_stats)
}

summarise_by_scenario <- function(case_data, group_vars, ci_width = 0.95,
                                  alpha_prior = 1, beta_prior = 1){
  #' Convert a database of case reports
  group_vars <- c("scenario", group_vars[! group_vars %in% c("scenario", "run")])
  beta_bound <- function(q,x,y) qbeta(q, alpha_prior+x, beta_prior+y)
  beta_lower <- function(x,y) beta_bound((1-ci_width)/2, x, y)
  beta_upper <- function(x,y) beta_bound(1-(1-ci_width)/2, x, y)
  db_runs <- summarise_by_run(case_data, group_vars)
  scenario_stats <- db_runs %>% copy %>%
    .[, .(n_runs = length(run), p_extinct = mean(extinct),
          p_exceed_case_cap = mean(total_cases > cap_cases),
          p_exceed_week_cap = mean(last_exposure_weeks > cap_max_weeks),
          p_exceed_gen_cap = mean(max_generation > cap_max_generations),
          n_controlled = sum(controlled),
          n_uncontrolled = sum(!controlled),
          p_controlled = mean(controlled)), by = group_vars] %>%
    .[, `:=`(p_controlled_lower = beta_lower(n_controlled, n_uncontrolled),
             p_controlled_upper = beta_upper(n_controlled, n_uncontrolled))]
  return(scenario_stats)
}

#------------------------------------------------------------------------------
# Write output
#------------------------------------------------------------------------------

write_dfile <- function(data, path_prefix, compress = TRUE){
  if (!compress) f <- file else f <- gzfile
  if (!compress) ext <- ".tsv" else ext <- ".tsv.gz"
  dpath <- paste0(path_prefix, ext)
  dfile <- f(dpath, "w")
  write_tsv(data, dfile)
  close(dfile)
}

write_data <- function(case_data, path_prefix, write_raw = FALSE,
                       write_weekly = FALSE, write_generational = FALSE,
                       write_run = FALSE, write_scenario = FALSE,
                       compress = TRUE, group_vars = NULL, ci_width = NULL,
                       alpha_prior = NULL, beta_prior = NULL){
  #' Write simulated case data and summaries to file
  if (write_raw){
    write_dfile(case_data, paste0(path_prefix, "_raw"), compress)
  }
  if (write_weekly){
    week_data <- summarise_by_week(case_data, group_vars)
    write_dfile(week_data, paste0(path_prefix, "_weekly"), compress)
  }
  if (write_generational){
    gen_data <- summarise_by_generation(case_data, group_vars)
    write_dfile(gen_data, paste0(path_prefix, "_gen"), compress)
  }
  if (write_run){
    run_data <- summarise_by_run(case_data, group_vars)
    write_dfile(run_data, paste0(path_prefix, "_run"), compress)
  }
  if (write_scenario){
    scenario_data <- summarise_by_scenario(case_data, group_vars, ci_width,
                                           alpha_prior, beta_prior)
    write_dfile(scenario_data, paste0(path_prefix, "_scenario"), compress)
  }
}

#------------------------------------------------------------------------------
# Overall wrapper function
#------------------------------------------------------------------------------

simulate_process <- function(scenario_parameters = NULL, n_iterations = NULL,
                             parallel = NULL, report = NULL,
                             show_progress = NULL, path_prefix = NULL,
                             write_raw = NULL, write_weekly = NULL,
                             write_generational = NULL, write_run = NULL,
                             write_scenario = NULL, compress_output = NULL,
                             log_path = NULL, ci_width = NULL,
                             alpha_prior = NULL, beta_prior = NULL){
  #' Generate a table of scenarios, run parameter sweep, and write output
  if (!is.null(log_path) && !is.na(log_path)){
    log_file <- file(log_path, "wt")
    sink(log_file)
    sink(log_file, type = "message")
  }
  # 1. Generate scenario table from input list of name/value combinations
  scenarios <- scenario_parameters %>% 
    expand.grid(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
    as_tibble %>% dplyr::mutate(scenario = 1:dplyr::n())
  # 2. Run parameter sweep
  sweep <- parameter_sweep(scenarios = scenarios, n_iterations = n_iterations,
                           parallel = parallel, report = report,
                           show_progress = show_progress)
  # 3. Write data and summaries
  write_data(case_data = sweep, path_prefix = path_prefix,
             write_raw = write_raw, write_weekly = write_weekly,
             write_generational = write_generational,
             write_run = write_run, write_scenario = write_scenario,
             compress = compress_output, group_vars = colnames(scenarios),
             ci_width = ci_width, alpha_prior = alpha_prior,
             beta_prior = beta_prior)
  if (!is.null(log_path) && !is.na(log_path)){
    sink()
    sink(type = "message")
    close(log_file)
  }
}