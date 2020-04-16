#' Core functionality for branching-chain outbreak model
#'
#' @author Joel Hellewell (original version), Will Bradshaw (modified version)

library(tidyverse)
library(data.table)
library(sn)

#----------------------------------------------------------------------------
# Auxiliary functions
#----------------------------------------------------------------------------

draw_new_cases <- function(parents, r0_sym, r0_asym, dispersion){
    #' Calculate the number of new cases produced by each individual in a
    #' parent case dataset by sampling from a negative binomial distribution,
    #' with mean=r0 depending on whether the individual is asymptomatic.
    return(parents %>% .[, n_children := rnbinom(n=nrow(.), size = dispersion,
                                          mu = ifelse(asym, r0_asym, r0_sym))])
} # TODO: Remember to filter out older generations before running this

map2_tab <- function(tab, col1, col2, fn){
    unlist(purrr::map2(tab[[col1]], tab[[col2]], function(x, y) fn(x,y)))
}
rep_new_cases <- function(tab, col){
    map2_tab(tab, col, "n_children", function(x,y) rep(x, as.integer(y)))
}

generate_new_cases <- function(parents, p_asymptomatic, p_traced, k,
                               generation_time, incubation_time, delay_time,
                               p_isolation){
    #' Generate a table of new child cases from a table of parent cases
    data.table(infector = rep_new_cases(parents, "case_id"), # Parent case
               infector_onset_gen = rep_new_cases(parents, "onset_gen"), # Symptom onset of parent
               infector_onset_true = rep_new_cases(parents, "onset_true"), # Symptom onset of parent
               infector_iso_time = rep_new_cases(parents, "isolation_time"), # Parent's iso time
               infector_asym = rep_new_cases(parents, "asym"), # Is parent asymptomatic?
               generation = rep_new_cases(parents, "generation") + 1, # New generation number
               processed = FALSE, n_children = NA) %>%
    .[,`:=`(exposure = generation_time(infector_onset_gen),
            # TODO: Refactor generation time to depend on parent's exposure, rather than onset?
            asym = purrr::rbernoulli(n = nrow(.), p = p_asymptomatic),
            # Successfully traced from parent?
            traced_fwd = purrr::rbernoulli(n = nrow(.), p = p_traced),
            # Successful backtracing to parent (if used)?
            traced_rev = purrr::rbernoulli(n = nrow(.), p = p_traced),
            escapes_isolation = purrr::rbernoulli(n = nrow(.), p = 1-p_isolation))] %>%
    # Symptom onset for determining generation time # TODO: Refactor gentime so can drop
    .[, onset_gen := exposure + incubation_time(nrow(.))] %>%
    # True symptom onset (infinite for asymptomatic individuals)
    .[, onset_true := ifelse(asym, Inf, onset_gen)] %>%
    # Maximum isolation time (if untraced)
    .[, iso_time_untraced := onset_true + delay_time(nrow(.))]
}  # TODO: Pass a distribution of compliance rates to be sampled by each parent case, rather than a fixed point value?

terminate_extinct <- function(cases){
    #' End-of-step cleanup for extinction case
    return(list(cases = cases %>% .[, processed := TRUE],
                effective_r0 = 0, cases_in_generation = 0))
}

filter_new_cases <- function(cases){
    #' Filter new cases based on isolation time of parent
    # TODO: Add quarantine time here as well
    # TODO: Pass parent and child tables separately, merge in-function
    # TODO: Generalise filter column (currently only infector_iso_time)
    # Separate out previously-escaped cases and cases that precede parent's isolation
    cases[exposure < infector_iso_time | escapes_isolation == TRUE]
}

trace_forward <- function(cases, quarantine, sero_test){
    #' Determine isolation times for secondary cases based on forward tracing
    # TODO: Pass parent and child tables separately, merge in-function
    # TODO: Generalise infector_iso_time to general fwd_tracing_time (to account for delays etc)
    # TODO: Distinguish quarantine time from isolation time; tracing determines the former
    cases %>% .[, isolation_time := ifelse(!traced_fwd, iso_time_untraced,
                                           pmin(iso_time_untraced,
                                                ifelse(rep(!quarantine && !sero_test, nrow(.)),
                                                       pmax(onset_true, infector_iso_time),
                                                       infector_iso_time)))]
    # - Default untraced isolation time depends on syptomatic status
    # - If missed in forward tracing, isolation time defaults to untraced time
    # - Otherwise, isolation time is the minimum of the untraced and traced times
    # - If quarantine in place, traced isolation time equals infector isolation time
    # - Otherwise, traced isolation time is the maximum of infector isolation time and symptom onset
    # - For asymptomatic individuals, this means isolation time is infinite unless quarantined/tested
    # - For individuals with untraced asymptomatic parents (i.e. parents with infinite onset),
    #   this produces the untraced time in all cases (TODO: test this)
}

run_backtrace_nstep <- function(cases, quarantine, sero_test,
                                max_distance){
    #' Perform backtracing up to some specified number of generations before
    #' the most recent generation.
    # Determine which generations to act on and separate out others
    gen_latest <- max(cases$generation)
    gen_earliest <- max(0, gen_latest - max_distance)
    cases_old <- cases[generation < gen_earliest] # Keep unchanged for later
    cases_act <- cases[generation >= gen_earliest]
    # Split remaining cases by generation, then copy for initial upwards trace
    cases_split <- split(cases_act, by="generation") # NB: generation labels are *strings*
    cases_up <- cases_split
    # Upwards trace: starting with latest generation, identify cases whose 
    # isolation time precedes that of parent and filter by backtrace success
    for (k in seq(gen_latest, gen_earliest+1)){
        cases_pre_filtered <- cases_up[[as.character(k)]] %>%
            .[traced_rev == TRUE] %>% .[isolation_time < infector_iso_time]
        if (nrow(cases_pre_filtered) == 0){k <- k + 1; break}
            # Determine new isolation times of traced parents for each backtrace option
        cases_iso_new <- cases_pre_filtered %>%
            .[, infector_iso_new := pmin(infector_iso_time,
                                         ifelse(rep(!quarantine && !sero_test, nrow(.)),
                                            pmax(infector_onset_true, isolation_time),
                                            isolation_time))] %>%
            .[infector_iso_new < infector_iso_time] # Fairly sure this filter is pointless, but JIC
        if (nrow(cases_iso_new) == 0){k <- k + 1; break}
        # New parent isolation time becomes min across all successful backtraces
        parents_backtraced <- cases_iso_new[, .(infector_iso_new = min(infector_iso_new)),
                                            by = "infector"]
        # Copy new parental isolation time to children (other copy), then to parents (don't filter yet)
        cases_split[[as.character(k)]] <- merge(cases_split[[as.character(k)]], parents_backtraced, 
                                                by="infector", all=TRUE) %>%
            .[, infector_iso_time := pmin(infector_iso_time, infector_iso_new, na.rm=TRUE)] %>%
            .[,infector_iso_new:=NULL]
        cases_split[[as.character(k-1)]] <- parents_backtraced %>% 
            setnames("infector", "case_id") %>%
            setnames("infector_iso_new", "isolation_time_new") %>% 
            merge(cases_up[[as.character(k-1)]], by="case_id", all=TRUE) %>%
            .[, isolation_time := pmin(isolation_time, isolation_time_new, na.rm=TRUE)]
        cases_up[[as.character(k-1)]] <- cases_split[[as.character(k-1)]][!is.na(isolation_time_new)]
        cases_split[[as.character(k-1)]][,isolation_time_new:=NULL]
        cases_up[[as.character(k-1)]][,isolation_time_new:=NULL]
    }
    if (k > gen_latest) return(cases) # Aborted at first generation; nothing to be done
    # Downwards trace: starting with generation k, prune cases that occur after parental 
    # isolation or whose parents were pruned, then repeat forward tracing
    for (l in seq(k, gen_latest)){
        parent_ids <- cases_split[[as.character(l-1)]]$case_id
        cases_split[[as.character(l)]] <- cases_split[[as.character(l)]] %>%
            .[infector %in% parent_ids] %>% filter_new_cases() %>%
            trace_forward(quarantine, sero_test)
    }
    return(rbind(cases_old, data.table::rbindlist(cases_split, fill=TRUE), fill=TRUE))
}

run_backtrace_1step <- function(cases, quarantine, sero_test){
    #' Perform backtracing for a single generation (i.e. trace parents, but not grandparents etc)
    # Identify cases whose isolation time precedes that of parent and filter by backtrace success
    cases_pre_filtered <- cases[traced_rev == TRUE] %>% .[isolation_time < infector_iso_time]
    if (nrow(cases_pre_filtered) == 0) return(cases)
    # Determine new isolation times of traced parents for each backtrace option
    cases_iso_new <- cases_pre_filtered %>%
        .[, infector_iso_new := pmin(infector_iso_time,
                                     ifelse(rep(!quarantine && !sero_test, nrow(.)),
                                            pmax(infector_onset_true, isolation_time),
                                            isolation_time))] %>%
        .[infector_iso_new < infector_iso_time] # Fairly sure this filter is pointless, but JIC
    if (nrow(cases_iso_new) == 0) return(cases)
    # New parent isolation time becomes min across all successful backtraces
    parents_backtraced <- cases_iso_new[, .(infector_iso_new = min(infector_iso_new)),
                                        by = "infector"]
    # Separate out children from backtraced parents (no change to others), then filter out
    # those that were infected after new isolation time
    cases_not_backtraced <- cases[! infector %in% parents_backtraced$infector]
    cases_backtraced <- cases[parents_backtraced, on="infector"] %>%
        .[, infector_iso_time := pmin(infector_iso_time, infector_iso_new, na.rm=TRUE)] %>%
        filter_new_cases()
    if (nrow(cases_backtraced) == 0) return(cases_not_backtraced)
    # Repeat forward tracing on remaining child cases
    cases_backtraced_retraced <- trace_forward(cases_backtraced, quarantine, sero_test)
    return(rbind(cases_not_backtraced, cases_backtraced_retraced, fill = TRUE))
}

terminate_open <- function(cases){
    #' End-of-step cleanup for non-extinction case
    latest_generation <- max(cases$generation)
    n_latest <- nrow(cases[generation == latest_generation])
    n_prev <- nrow(cases[generation == latest_generation-1])
    return(list(cases = cases,
                cases_in_generation = n_latest,
                effective_r0 = n_latest/n_prev))
}

compute_weekly_cases <- function(case_data, max_weeks){
    #' Convert a database of case reports into one of weekly case counts
    # Compute weekly case counts
    weekly_cases <- case_data[, week := floor(exposure / 7)] %>%
        .[, .(weekly_cases = .N), by = week]
    # Add missing weeks
    missing_weeks <- (0:max_weeks)[!(0:max_weeks %in% weekly_cases$week)]
    if (length(missing_weeks > 0)) {
        missing_db <- data.table(week = missing_weeks, weekly_cases = 0)
        weekly_cases <- rbind(weekly_cases, missing_db, fill=TRUE)
    }
    # Cut at max week (for some reason)
    weekly_cases <- weekly_cases[week <= max_weeks]
    return(weekly_cases)
}

compute_symptomatic_r0 <- function(r0_base, r0_asymptomatic, p_asymptomatic){
    (r0_base - r0_asymptomatic*p_asymptomatic)/(1-p_asymptomatic)
}

#----------------------------------------------------------------------------
# Core functions
#----------------------------------------------------------------------------

outbreak_setup <- function(n_initial_cases, p_asymptomatic,
                           incubation_time, delay_time, p_isolation){
    #' Set up a table of initial cases
    data.table(infector = 0, infector_onset_gen = NA, infector_onset_true = NA,
               infector_iso_time = NA, infector_asym = NA, generation = 0,
               processed = FALSE, n_children = NA, exposure = 0,
               case_id = 1:n_initial_cases,
               traced_fwd = FALSE, traced_rev = FALSE) %>%
    .[, `:=`(asym = purrr::rbernoulli(nrow(.), p_asymptomatic), # Asymptomatic?
             escapes_isolation = purrr::rbernoulli(n = nrow(.), p = 1-p_isolation),
             # Symptom onset for determining generation time (TODO: refactor gentime so can drop)
             onset_gen = incubation_time(nrow(.)))] %>%
    # True symptom onset (infinite for asymptomatic individuals)
    .[, onset_true := ifelse(asym, Inf, onset_gen)] %>%
    # Maximum isolation time (if untraced)
    .[, iso_time_untraced := onset_true + delay_time(nrow(.))] %>%
    # Initial isolation time = maximum
    .[, isolation_time := iso_time_untraced]
}

outbreak_step_backtrace <- function(case_data = NULL, dispersion = NULL,
                                    r0_symptomatic = NULL, r0_asymptomatic = NULL,
                                    p_asymptomatic = NULL, p_traced = NULL,
                                    generation_time = NULL, incubation_time = NULL,
                                    delay_time = NULL, backtrace_distance = NULL,
                                    quarantine = NULL, sero_test = NULL,
                                    p_isolation = NULL){
    #' Move forward one generation in the branching process
    # Separate current parents from previous generations and draw new cases
    case_data_old <- case_data %>% .[processed == TRUE]
    case_data_new <- case_data %>% .[processed == FALSE] %>%
        draw_new_cases(r0_symptomatic, r0_asymptomatic, dispersion)
    # If no new cases, terminate step
    if (sum(case_data_new$n_children) == 0){
        return(terminate_extinct(rbind(case_data_old, case_data_new)))
    }
    # Generate putative secondary cases
    children_putative <- generate_new_cases(case_data_new, p_asymptomatic,
                                            p_traced, k, generation_time,
                                            incubation_time, delay_time, p_isolation)
    # Filter secondary cases based on isolation time of parents # TODO: Add quarantine time
    children_filtered <- filter_new_cases(children_putative)
    # Determine isolation time of remaining secondary cases
    children_iso <- trace_forward(children_filtered, quarantine, sero_test)
    # Assign case IDs to child cases and mark parents as processed
    if (nrow(children_iso) > 0){
        children_iso[["case_id"]] <- 1:nrow(children_iso) + max(case_data[["case_id"]])
    }
    case_data_new$processed <- TRUE
    # Perform backtracing, if applicable
    cases_iso <- data.table::rbindlist(list(case_data_old, case_data_new, children_iso),
                                       fill=TRUE)
    if (backtrace_distance > 0) {
        cases_out <- run_backtrace_nstep(cases_iso, quarantine, sero_test, 
                                         backtrace_distance)
    } else {
        cases_out <- cases_iso
    }
    return(terminate_open(cases_out))
}

outbreak_model <- function(n_initial_cases = NULL, r0_base = NULL,
                           dispersion = NULL, cap_max_weeks = NULL,
                           cap_max_generations = NULL, r0_asymptomatic = NULL,
                           p_asymptomatic = NULL, p_traced = NULL,
                           backtrace_distance = NULL, quarantine = NULL,
                           sero_test = NULL, p_isolation = NULL,
                           delay_shape = NULL, delay_scale = NULL,
                           incubation_shape = NULL, incubation_scale = NULL,
                           generation_omega = NULL, generation_k = NULL,
                           cap_cases = NULL){
    #' Run a single complete instance of the branching-process model
    # Set up incubation, generation-time, and delay distributions
    # TODO: Generalise these / replace with Ferretti functions
    delay_time <- function(n) rweibull(n, shape=delay_shape, scale=delay_scale)
    incubation_time <- function(n) rweibull(n, shape=incubation_shape,
                                            scale=incubation_scale)
    generation_time <- function(onsets){
        sn::rsn(n = length(onsets), xi = onsets, omega = generation_omega,
                alpha = generation_k) %>%  ifelse(. < 1, 1, .)
    } # TODO: Not sure why constraining to 1 or greater here
    # Compute symptomatic R0
    r0_symptomatic <- compute_symptomatic_r0(r0_base, r0_asymptomatic, p_asymptomatic)
    # Set initial values for loop indices and preallocate space for metrics
    # TODO: Change from cumulative case threshold to simultaneous
    total_cases <- n_initial_cases
    latest_exposure_days <- 0
    latest_exposure_weeks <- 0
    outbreak_extinct <- FALSE
    effective_r0_vect <- numeric(cap_max_generations)
    cases_in_gen_vect <- numeric(cap_max_generations)
    latest_generation <- 0
    # Set up initial cases
    case_data <- outbreak_setup(n_initial_cases, p_asymptomatic,
                                incubation_time, delay_time, p_isolation)
    # Run outbreak loop
    while (latest_exposure_weeks < cap_max_weeks & total_cases < cap_cases &
           !outbreak_extinct & latest_generation <= cap_max_generations){
        out <- outbreak_step_backtrace(case_data = case_data, dispersion = dispersion,
                                       r0_symptomatic = r0_symptomatic,
                                       r0_asymptomatic = r0_asymptomatic,
                                       p_asymptomatic = p_asymptomatic, p_traced = p_traced,
                                       generation_time = generation_time,
                                       incubation_time = incubation_time,
                                       delay_time = delay_time, p_isolation = p_isolation,
                                       quarantine = quarantine, sero_test = sero_test,
                                       backtrace_distance = backtrace_distance)
        case_data <- out[["cases"]]
        latest_generation <- latest_generation + 1
        effective_r0_vect[latest_generation] <- out[["effective_r0"]]
        cases_in_gen_vect[latest_generation] <- out[["cases_in_generation"]]
        total_cases <- nrow(case_data)
        latest_exposure_days <- max(case_data$exposure)
        latest_exposure_weeks <- floor(latest_exposure_days / 7)
        outbreak_extinct <- all(case_data$processed)
    }
    # Compute weekly cases and add final data
    avg_effective_r0 <- sum(effective_r0_vect*cases_in_gen_vect)/
        sum(cases_in_gen_vect)
    weekly_cases <- compute_weekly_cases(case_data, cap_max_weeks) %>%
        .[, `:=`(avg_effective_r0 = avg_effective_r0,
                 cases_per_gen = paste(cases_in_gen_vect[1:latest_generation], collapse = "|"),
                 final_total_cases = total_cases,
                 outbreak_extinct = outbreak_extinct,
                 final_generation = latest_generation)]
    return(weekly_cases)
}

