#' Run a specified number of simulations with identical parameters
#' @author Joel Hellewell

source("R/outbreak_model.R")

library(tidyverse)
library(data.table)


scenario_sim <- function(n.sim = NULL, p_traced = NULL, cap_max_days = NULL, cap_cases = NULL,
                         r0isolated = NULL, r0community = NULL, disp.iso = NULL, disp.com = NULL, k = NULL,
                         delay_shape = NULL, delay_scale = NULL, num.initial.cases = NULL, p_asymptomatic = NULL,
                         quarantine = NULL, backtrace = NULL, sero_test = NULL, report = NULL,
                         run_new = NULL) {
  
  if (!is.null(report)){cat(report, ": ", sep="")}

  # Run n.sim number of model runs and put them all together in a big data.frame
  res_list <- purrr::map(.x = 1:n.sim, ~ outbreak_model(num.initial.cases = num.initial.cases,
                                             p_traced = p_traced,
                                             cap_max_days = cap_max_days,
                                             cap_cases = cap_cases,
                                             r0isolated = r0isolated,
                                             r0community = r0community,
                                             disp.iso = disp.iso,
                                             disp.com = disp.com,
                                             delay_shape = delay_shape,
                                             delay_scale = delay_scale,
                                             k = k,
                                             p_asymptomatic = p_asymptomatic,
                                             quarantine = quarantine,
                                             backtrace = backtrace,
                                             sero_test = sero_test,
                                             run_new = run_new
                                             ))

  # bind output together and add simulation index
  res <- data.table::rbindlist(lapply(1:n.sim, function(n) res_list[[n]] <- res_list[[n]][, sim := n]))
  # res <- data.table::rbindlist(res_list)
  # res[, sim := rep(1:n.sim, rep(floor(cap_max_days / 7) + 1, n.sim)), ] # WJB: This might win my award for kludge of the year
  
  if (!is.null(report)){cat(timestamp(quiet=TRUE), "\n")}
  
  return(res)
}
