#' Move forward one generation in the branching process (new code)
#'
#' @author Joel Hellewell (original version), Will Bradshaw (modified version)

outbreak_step_backtrace <- function(case_data = NULL, disp.iso = NULL, disp.com = NULL,
                                    r0isolated = NULL, r0community = NULL, p_asymptomatic = NULL,
                                    incfn = NULL, delayfn = NULL, p_traced = NULL,
                                    k = NULL, quarantine = NULL, backtrace = NULL,
                                    sero_test = NULL) {
  
  # A vectorised version of isTRUE
  vect_isTRUE <- function(x) {
    purrr::map_lgl(x, isTRUE)
  }
  
  vect_max <- function(x, y) {
    purrr::map2_dbl(x, y, max)
  }
  
  vect_min <- function(x, y) {
    purrr::map2_dbl(x, y, min)
  }
  
  # For each case in case_data, draw new_cases from a negative binomial distribution
  # with an R0 and dispersion dependent on if isolated=TRUE
  case_data[, new_cases := purrr::map2_dbl(
    ifelse(vect_isTRUE(isolated), disp.iso, disp.com),
    ifelse(vect_isTRUE(isolated),
           r0isolated,
           r0community),
    ~ rnbinom(1, size = .x, mu = .y))
    ]
  
  # Select cases that have generated any new cases
  new_case_data <- case_data[new_cases > 0]
  # The total new cases generated
  total_new_cases <- case_data[, sum(new_cases), ]
  
  # If no new cases drawn, outbreak is over so return case_data
  if (total_new_cases == 0) {
    # If everyone is isolated it means that either control has worked or everyone has had a chance to infect but didn't
    case_data$isolated <- TRUE
    
    effective_r0 <- 0
    cases_in_gen <- 0
    out <- list(case_data, effective_r0, cases_in_gen)
    names(out) <- c("cases", "effective_r0", "cases_in_gen")
    
    return(out)
  }
  
  #----------------------------------------------------------------------------
  # Generate (putative) secondary cases
  #----------------------------------------------------------------------------
  
  # Compile a data.table for all new cases, new_cases is the amount of people that each infector has infected
  inc_samples <- incfn(total_new_cases)
  
  prob_samples <- data.table(
    # time when new cases were exposed, a draw from serial interval based on infector's onset
    exposure = unlist(purrr::map2(new_case_data$new_cases, new_case_data$onset,
                                  function(x, y) {
                                    inf_fn(rep(y, x), k)
                                  })),
    # records the infector of each new person
    infector = unlist(purrr::map2(new_case_data$caseid, new_case_data$new_cases,
                                  function(x, y) {
                                    rep(as.integer(x), as.integer(y))
                                  })),
    # records when infector developed symptoms
    infector_onset = unlist(purrr::map2(new_case_data$onset, new_case_data$new_cases,
                                        function(x, y) {
                                          rep(x, as.integer(y))
                                        })),
    # records when infector was isolated
    infector_iso_time = unlist(purrr::map2(new_case_data$isolated_time, new_case_data$new_cases,
                                           function(x, y) {
                                             rep(x, as.integer(y))
                                           })),
    # records if infector asymptomatic
    infector_asym = unlist(purrr::map2(new_case_data$asym, new_case_data$new_cases,
                                       function(x, y) {
                                         rep(x, y)
                                       })),
    # draws a sample to see if this person is asymptomatic
    asym = purrr::rbernoulli(n = total_new_cases, p = p_asymptomatic),
    # draws a sample to see if this person is traced
    missed = purrr::rbernoulli(n = total_new_cases, p = 1 - p_traced),
    # sample from the incubation period for each new person
    incubfn_sample = inc_samples,
    isolated = FALSE,
    new_cases = NA
  )
  
  prob_samples <- prob_samples[exposure < infector_iso_time][, # filter out new cases prevented by isolation
                                                             `:=`(# onset of new case is exposure + incubation period sample
                                                               onset = exposure + incubfn_sample)]
  
  # Cases whose parent was never isolated are automatically missed (until symptoms or backtracing)
  prob_samples$missed[prob_samples$infector_iso_time == Inf] <- TRUE
  
  # asymptomatic cases have infinite onset and are missed if testing is unavailable
  # prob_samples$onset[vect_isTRUE(prob_samples$asym)] <- Inf (TODO: Restore this after fixing incubation function)
  # if (!sero_test){
  #   prob_samples$missed[vect_isTRUE(prob_samples$asym)] <- TRUE
  # }
  
  #----------------------------------------------------------------------------
  # Determine isolation time
  #----------------------------------------------------------------------------
  
  n_asym <- sum(prob_samples$asym)
  n_sym  <- nrow(prob_samples) - n_asym
  
  # Asymptomatic case:
  # - If no quarantine and no testing, or if missed, never isolated (no symptoms)
  # - If quarantine or testing, and not missed, isolated upon tracing (i.e. infector symptom onset)
  # NB: Quarantine and testing cases will diverge once we implement testing delays
  #     or imperfect testing, but are equivalent for now
  
  if (n_asym > 0){ # TODO: Just functionalise this already
    prob_samples_asym <- prob_samples[asym==TRUE] %>%
      .[, iso_time_untraced := Inf] %>%
      .[, isolated_time := ifelse(!quarantine && !sero_test | missed,
                                  iso_time_untraced, infector_iso_time)]
  } else {prob_samples_asym <- head(prob_samples, 0)}
  
  # Symptomatic case
  # - If missed, you are isolated at symptom onset (plus delay)
  # - If not missed, you are isolated at symptom onset (plus delay) or 
  #   post-tracing (see next point) whichever comes first
  # - Post-tracing response depends on quarantine and testing status:
  #   * If quarantine active, you are isolated immediately upon tracing
  #   * If quarantine inactive but testing is in place, you are isolated immediately upon tracing
  #   * Otherwise, you are isolated immediately upon symptom onset
  # NB: Currently quarantine and tracing options are redundant;
  #     this will change once we implement testing delays or imperfect testing
  
  if (n_sym > 0){
    prob_samples_sym <- prob_samples[asym==FALSE] %>%
      .[, iso_time_untraced := onset + delayfn(1)] %>%
      .[, isolated_time := ifelse(missed, iso_time_untraced,
                                  pmin(iso_time_untraced,
                                       ifelse(!quarantine && !sero_test,
                                              pmax(onset, infector_iso_time), infector_iso_time)))]
  } else {prob_samples_sym <- head(prob_samples, 0)[, isolated_time := NA]}

    
  prob_samples_iso <- rbind(prob_samples_asym, prob_samples_sym, fill=TRUE)
  
  #----------------------------------------------------------------------------
  # Backtrace to parents, then retrace
  #----------------------------------------------------------------------------
  
  if (backtrace){
    # Asymptomatic infectors are considered (for backtracing) to have infinite onset
    # WJB: This is only necessary as long as the onset itself can't be infinite; fix after fixing generation time distr
    prob_samples_iso$infector_onset[vect_isTRUE(prob_samples_iso$infector_asym)] <- Inf
    
    # Identify cases whose isolation time precedes that of parent and test
    # successful contact tracing of parent
    prob_samples_pre <- prob_samples_iso[isolated_time < infector_iso_time] %>%
      .[, infector_traced := purrr::rbernoulli(n = nrow(.), p = p_traced)] %>%
      .[infector_traced==TRUE]
    
    # Determine new isolation time of traced parents (NB: can ignore cases where infector missed, since these filtered out)
    n_asym <- sum(prob_samples_pre$infector_asym)
    n_sym  <- nrow(prob_samples_pre) - n_asym
    
    if (n_asym > 0){
      prob_samples_pre_asym <- prob_samples_pre[infector_asym==TRUE] %>%
        .[, infector_iso_new := ifelse(!quarantine && !sero_test, infector_iso_time, isolated_time)]
    } else {prob_samples_pre_asym <- head(prob_samples_pre,0)}
    
    if (n_sym > 0){
      prob_samples_pre_sym <- prob_samples_pre[infector_asym==FALSE] %>%
        .[, infector_iso_new := pmin(infector_iso_time,
                                     ifelse(!quarantine && !sero_test,
                                            pmax(infector_onset, isolated_time), isolated_time))]
    } else {prob_samples_pre_sym <- head(prob_samples_pre,0)[, infector_iso_new := NA]}
    
    # Combine, then filter out cases where "new" infector isolation time exceeds old
    prob_samples_pre_iso <- rbind(prob_samples_pre_asym, prob_samples_pre_sym,
                                  fill = TRUE) %>%
      .[infector_iso_new < infector_iso_time]

    if (nrow(prob_samples_pre_iso) > 0){
      
      # New parent isolation time becomes min across all successful backtraces
      parents_backtraced <- prob_samples_pre_iso[, .(infector_iso_new = min(infector_iso_new)), by = "infector"]
      
      # Separate children of backtraced parents from others (no change to others) and filter out
      # children from latter that were infected after new isolation time
      prob_samples_not_backtraced <- prob_samples_iso[! infector %in% parents_backtraced$infector]
      prob_samples_backtraced <- prob_samples_iso[parents_backtraced, on="infector"] %>%
        .[, infector_iso_time := pmin(infector_iso_time, infector_iso_new, na.rm=TRUE)] %>%
        .[exposure < infector_iso_time] %>%
        .[, missed := purrr::rbernoulli(n = nrow(.), p = 1 - p_traced)]
      
      # Repeat contact tracing for remaining children of backtraced parents
      # (Don't think I need to worry about which child originated the backtrace, since
      # their isolation time should be unchanged) # TODO: Check this more thoroughly
      n_asym <- sum(prob_samples_backtraced$asym)
      n_sym  <- nrow(prob_samples_backtraced) - n_asym
      
      if (n_asym > 0){
        prob_samples_backtraced_retraced_asym <- prob_samples_backtraced[asym==TRUE] %>%
          .[, isolated_time := ifelse(!quarantine && !sero_test | missed,
                                      iso_time_untraced, infector_iso_time)]
      } else {prob_samples_backtraced_retraced_asym <- head(prob_samples_backtraced,0)}
      
      if (n_sym > 0){
        prob_samples_backtraced_retraced_sym <- prob_samples_backtraced[asym==FALSE] %>%
          .[, isolated_time := ifelse(missed, iso_time_untraced,
                                      pmin(iso_time_untraced,
                                           ifelse(!quarantine && !sero_test,
                                                  pmax(onset, infector_iso_time), infector_iso_time)))]
      } else {prob_samples_backtraced_retraced_sym <- head(prob_samples_backtraced,0)[, isolated_time := NA]}

      
      prob_samples_out <- rbind(prob_samples_not_backtraced,
                                prob_samples_backtraced_retraced_asym,
                                prob_samples_backtraced_retraced_sym, fill = TRUE)
    } else {
      prob_samples_out <- prob_samples_iso
    }
  } else {
    prob_samples_out <- prob_samples_iso
  }
  
  #----------------------------------------------------------------------------
  # Finalise and return output
  #----------------------------------------------------------------------------
  
  # Chop out unneeded sample columns
  #prob_samples[, c("incubfn_sample", "infector_iso_time", "infector_asym") := NULL]
  
  # Set new case ids for new people
  prob_samples_out$caseid <- (nrow(case_data) + 1):(nrow(case_data) + nrow(prob_samples_out))
  
  ## Number of new cases
  cases_in_gen <- nrow(prob_samples_out)
  
  ## Estimate the effective r0
  effective_r0 <- nrow(prob_samples_out) / nrow(case_data[!vect_isTRUE(case_data$isolated)])
  
  # Everyone in case_data so far has had their chance to infect and are therefore considered isolated
  case_data$isolated <- TRUE
  
  # bind original cases + new secondary cases
  case_data_out <- data.table::rbindlist(list(case_data, prob_samples_out),
                                         use.names = TRUE, fill=TRUE)
  
  # Return
  out <- list(cases=case_data_out, effective_r0=effective_r0, cases_in_gen=cases_in_gen)
  
  return(out)
}
