#' Wrapper code for branching process model

#------------------------------------------------------------------------------
# Run parameter sweep
#------------------------------------------------------------------------------

parameter_sweep <- function(scenarios = NULL, n_iterations = NULL,
                            threads = NULL, show_progress = NULL, report = NULL,
                            write_raw = NULL, write_weekly = NULL,
                            write_generational = NULL, write_run = NULL,
                            path_prefix = NULL, compress = NULL, ci_width = NULL,
                            alpha_prior = NULL, beta_prior = NULL){
    #' Run one set of simulations for each scenario in a table of scenarios
    # Nest scenarios into sub-tables
    sim_fn <- function(n) scenario_sim(scenario = scenarios[n,],
                                       n_iterations = n_iterations,
                                       report = report,
                                       write_raw = write_raw,
                                       write_weekly = write_weekly,
                                       write_generational = write_generational,
                                       write_run = write_run,
                                       path_prefix = path_prefix,
                                       compress = compress,
                                       ci_width = ci_width,
                                       alpha_prior = alpha_prior,
                                       beta_prior = beta_prior)
    sim_data <- mclapply(1:nrow(scenarios), sim_fn, mc.cores = threads)
    #print(sim_data)
    return(sim_data %>% rbindlist)
}

#------------------------------------------------------------------------------
# Overall wrapper function
#------------------------------------------------------------------------------

simulate_process <- function(scenario_parameters = NULL, n_iterations = NULL,
                             report = NULL, threads = NULL,
                             show_progress = NULL, path_prefix = NULL,
                             write_raw = NULL, write_weekly = NULL,
                             write_generational = NULL, write_run = NULL,
                             compress_output = NULL,
                             log_path = NULL, ci_width = NULL,
                             alpha_prior = NULL, beta_prior = NULL){
  #' Generate a table of scenarios, run parameter sweep, and write output
  if (!is.null(log_path) && !is.na(log_path)){
    log_file <- file(log_path, "wt")
    sink(log_file)
    sink(log_file, type = "message")
    start <- proc.time()
    cat("Beginning simulation: ", date(), "\n", sep="")
  }
  # 1. Generate scenario table from input list of name/value combinations
  scenarios <- scenario_parameters %>% 
    expand.grid(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
    as_tibble %>% dplyr::mutate(scenario = 1:dplyr::n())
  cat("Number of scenarios:", nrow(scenarios), "\n")
  cat("Requested cores:", threads, "\n")
  cat("Available cores:", detectCores(), "\n")
  threads <- max(threads, detectCores())
  # 2. Run parameter sweep
  sweep <- parameter_sweep(scenarios = scenarios, n_iterations = n_iterations,
                           report = report, threads = threads,
                           show_progress = show_progress,
                           path_prefix = path_prefix,
                           write_raw = write_raw, write_weekly = write_weekly,
                           write_generational = write_generational,
                           write_run = write_run,
                           compress = compress_output,
                           ci_width = ci_width,
                           alpha_prior = alpha_prior,
                           beta_prior = beta_prior)
  # 3. Write summarised scenario data
  write_dfile(sweep, paste0(path_prefix, "_scenario"), compress_output)
  if (!is.null(log_path) && !is.na(log_path)){
    cat("End of simulation: ", date(), " (", timetaken(start), ")\n", sep="")
    sink()
    sink(type = "message")
    close(log_file)
  }
}
