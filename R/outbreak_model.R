#' Run a single instance of the branching process model

source("R/outbreak_setup_old.R")
source("R/outbreak_setup_new.R")
source("R/outbreak_step_old.R")
source("R/outbreak_step_new.R")
source("R/aux_functions.R")

library(tidyverse)
library(data.table)

outbreak_model <- function(num.initial.cases = NULL, p_traced = NULL,
                           cap_max_days = NULL, cap_cases = NULL,
                           r0isolated = NULL, r0community = NULL,
                           disp.iso = NULL, disp.com = NULL,
                           k = NULL, delay_shape = NULL,
                           delay_scale = NULL, p_asymptomatic = NULL,
                           quarantine = NULL, backtrace = NULL,
                           sero_test = NULL, run_new = NULL) {

  # Set up functions to sample from distributions
  # incubation period sampling function
  incfn <- dist_setup(dist_shape = 2.322737,
                      dist_scale = 6.492272)
  # incfn <- dist_setup(dist_shape = 3.303525,dist_scale = 6.68849) # incubation function for ECDC run
  # onset to isolation delay sampling function
  delayfn <- dist_setup(delay_shape,
                        delay_scale)

  # Set initial values for loop indices
  total.cases <- num.initial.cases
  latest.exposure <- 0
  extinct <- FALSE

  # Initial setup
  if (!run_new){
    case_data <- outbreak_setup_old(num.initial.cases = num.initial.cases,
                                    incfn = incfn,
                                    p_asymptomatic = p_asymptomatic,
                                    delayfn = delayfn,
                                    k = k)
  } else {
    case_data <- outbreak_setup_new(num.initial.cases = num.initial.cases,
                                    incfn = incfn,
                                    p_asymptomatic = p_asymptomatic,
                                    delayfn = delayfn,
                                    k = k)
  }
  
  initial_asym <- sum(case_data$asym)

  # Preallocate
  effective_r0_vect <- c()
  cases_in_gen_vect <- c()


  # Model loop
  while (latest.exposure < cap_max_days & total.cases < cap_cases & !extinct) {

    if (!run_new){
      out <- outbreak_step(case_data = case_data,
                           disp.iso = disp.iso,
                           disp.com = disp.com,
                           r0isolated = r0isolated,
                           r0community = r0community,
                           incfn = incfn,
                           delayfn = delayfn,
                           p_traced = p_traced,
                           k = k,
                           quarantine = quarantine,
                           p_asymptomatic = p_asymptomatic)
    } else {
      out <- outbreak_step_backtrace(case_data = case_data,
                                     disp.iso = disp.iso,
                                     disp.com = disp.com,
                                     r0isolated = r0isolated,
                                     r0community = r0community,
                                     incfn = incfn,
                                     delayfn = delayfn,
                                     p_traced = p_traced,
                                     k = k,
                                     quarantine = quarantine,
                                     p_asymptomatic = p_asymptomatic,
                                     backtrace = backtrace,
                                     sero_test = sero_test)
    }

    case_data <- out[[1]]
    effective_r0_vect <- c(effective_r0_vect, out[[2]])
    cases_in_gen_vect <- c(cases_in_gen_vect, out[[3]])
    total.cases <- nrow(case_data)
    latest.exposure <- max(case_data$exposure)
    extinct <- all(case_data$isolated)
  }

  # Prepare output, group into weeks
  weekly_cases <- case_data[, week := floor(exposure / 7) # WJB: Changed onset to exposure here so as not to neglect asymptomatic cases in case count
                            ][, .(weekly_cases = .N), by = week
                              ]
  # maximum outbreak week
  max_week <- floor(cap_max_days / 7)
  # weeks with 0 cases in 0:max_week
  missing_weeks <- (0:max_week)[!(0:max_week %in% weekly_cases$week)]

  # add in missing weeks if any are missing
  if (length(missing_weeks > 0)) {
    weekly_cases <- data.table::rbindlist(list(weekly_cases,
                                               data.table(week = missing_weeks,
                                                          weekly_cases = 0)))
  }
  # order and sum up
  weekly_cases <- weekly_cases[order(week)
                               ][, cumulative := cumsum(weekly_cases)]
  # cut at max_week
  weekly_cases <- weekly_cases[week <= max_week]

  # Add effective R0
  weekly_cases <- weekly_cases[, `:=`(effective_r0 = mean(effective_r0_vect,
                                                          na.rm = TRUE),
                                      cases_per_gen = list(cases_in_gen_vect),
                                      final_total_cases = total.cases,
                                      final_latest_exposure = latest.exposure,
                                      final_extinct = extinct,
                                      initial_asym = initial_asym)]
  # return
  return(weekly_cases)
}
