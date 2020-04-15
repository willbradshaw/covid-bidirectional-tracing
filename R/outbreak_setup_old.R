#' Set up initial cases for branching process (old code)
#' @author Joel Hellewell

outbreak_setup_old <- function(num.initial.cases, incfn, delayfn, k, p_asymptomatic) {
  # Set up table of initial cases
  inc_samples <- incfn(num.initial.cases)

  case_data <- data.table(exposure = rep(0, num.initial.cases), # Exposure time of 0 for all initial cases
                          asym = purrr::rbernoulli(num.initial.cases, p_asymptomatic),
                          caseid = 1:(num.initial.cases), # set case id
                          infector = 0,
                          missed = TRUE,
                          onset = inc_samples,
                          new_cases = NA)

  # set isolation time for cluster to minimum time of onset of symptoms + draw from delay distribution
  case_data <- case_data[, isolated_time := onset + delayfn(1)
                         ][, isolated := FALSE]

  # Asymptomatic cases have infinite isolation time (TODO: Revert this after fixing incubation function)
  case_data$isolated_time[case_data$asym] <- Inf

  # return
  return(case_data)
}
