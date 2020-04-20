#' Core functionality for branching-chain outbreak model
#'
#' @author Joel Hellewell (original version), Will Bradshaw (modified version)

library(tidyverse)
library(data.table)
library(sn)

#----------------------------------------------------------------------------
# Create new cases
#----------------------------------------------------------------------------

set_primary_immutables_index <- function(n_initial_cases, p_asymptomatic,
                                         test_time, p_ident_sym){
    #' Set primary immutable keys for index cases
    data.table(case_id = 1:n_initial_cases,
               asym = purrr::rbernoulli(n_initial_cases, p_asymptomatic),
               blocked_isolation = NA, # These two are for child cases only
               blocked_quarantine = NA,
               test_delay = test_time(n_initial_cases),
               ident_sym = purrr::rbernoulli(n_initial_cases, p_ident_sym),
               generation = 0, infector_id = 0, infector_onset_gen = NA,
               infector_onset_true = NA, infector_asym = NA, 
               infector_has_smartphone = NA)
}

map2_tab <- function(tab, col1, col2, fn){
    unlist(purrr::map2(tab[[col1]], tab[[col2]], function(x, y) fn(x,y)))
}
rep_new_cases <- function(tab, col){
    map2_tab(tab, col, "n_children", function(x,y) rep(x, as.integer(y)))
}

set_primary_immutables_child <- function(parents, p_asymptomatic, p_blocked_isolation,
                                         p_blocked_quarantine, test_time, p_ident_sym){
    #' Set primary immutable keys for index cases
    n_children_total <- sum(parents$n_children)
    case_ids <- max(parents$case_id)+1:n_children_total
    child_cases <- data.table(case_id = case_ids,
                              asym = purrr::rbernoulli(n_children_total, p_asymptomatic),
                              blocked_isolation = purrr::rbernoulli(n_children_total, p_blocked_isolation),
                              blocked_quarantine = purrr::rbernoulli(n_children_total, p_blocked_quarantine),
                              test_delay = test_time(n_children_total),
                              ident_sym = purrr::rbernoulli(n_children_total, p_ident_sym),
                              generation = rep_new_cases(parents, "generation") + 1, # New generation number
                              infector_id = rep_new_cases(parents, "case_id"), # Parent case
                              infector_onset_gen = rep_new_cases(parents, "onset_gen"), # Symptom onset of parent
                              infector_onset_true = rep_new_cases(parents, "onset_true"), # Symptom onset of parent
                              infector_asym = rep_new_cases(parents, "asym"), # Is parent asymptomatic?
                              infector_has_smartphone = rep_new_cases(parents, "has_smartphone")) # Is parent asymptomatic?
    return(child_cases)
}

create_child_cases <- function(parents, p_asymptomatic, p_blocked_isolation,
                               p_blocked_quarantine, test_time, p_ident_sym,
                               ...){
    #' Generate a table of new child cases from a table of parent cases
    cases <- set_primary_immutables_child(parents, p_asymptomatic,
                                          p_blocked_isolation,
                                          p_blocked_quarantine,
                                          test_time, p_ident_sym)
    # TODO: Set secondary immutables, parental mutables, nonparental mutables
    
}


generate_new_cases <- function(parents, p_asymptomatic, p_traced, k,
                               generation_time, incubation_time, delay_time,
                               p_isolation, rollout_delay_days,
                               rollout_delay_generations){
    data.table(                              infector_iso_time = rep_new_cases(parents, "isolation_time"), # Parent's iso time

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
        .[, isolation_time := pmax(onset_true, rollout_delay_days) +
              delay_time(nrow(.))] %>%
        .[, isolation_time := ifelse(generation < rep(rollout_delay_generations, nrow(.)),
                                     Inf, isolation_time)]
}  # TODO: Pass a distribution of compliance rates to be sampled by each parent case, rather than a fixed point value?


set_secondary_immutables <- function(cases, index, p_smartphone_overall,
                                     p_smartphone_parent_yes, p_smartphone_parent_no,
                                     generation_time,
                                     dispersion, r0_asymptomatic, r0_symptomatic,
                                     incubation_time, p_traced_auto, p_traced_manual,
                                     trace_time_auto, trace_time_manual,
                                     recovery_time, test_function){
    #' Set secondary immutable keys
    cases %>% .[, `:=`(has_smartphone = ifelse(rep(index, nrow(.)), p_smartphone_overall,
                                               ifelse(infector_has_smartphone,
                                                      p_smartphone_parent_yes,
                                                      p_smartphone_parent_no)),
                       exposure = ifelse(rep(index, nrow(.)), 0,
                                         generation_time(infector_onset_gen)),
                       n_children = rnbinom(n=nrow(.), size = dispersion,
                                            mu = ifelse(asym, r0_asymptomatic, 
                                                        r0_symptomatic)),
                       processed = FALSE)] %>% # (Not technically immutable, but independent of tracing)
        .[, `:=`(auto_traced = infector_has_smartphone & has_smartphone,
                 onset_gen := exposure + incubation_time(nrow(.)))] %>%
        .[, `:=`(traceable_fwd = purrr::rbernoulli(nrow(.), ifelse(auto_traced, p_traced_auto, p_traced_manual)),
                 traceable_rev = purrr::rbernoulli(nrow(.), ifelse(auto_traced, p_traced_auto, p_traced_manual)),
                 trace_delay = ifelse(auto_traced, trace_time_auto(nrow(.)), trace_time_manual(nrow(.))),
                 onset_true = ifelse(asym, Inf, onset_gen),
                 recovery = recovery_time(onset_gen))] %>%
        .[, test_positive = test_function(recovery)]
}

set_nonparental_mutables <- function(cases, trace_neg_symptomatic){
    #' Set values of nonparental mutable keys (except identification time)
    cases %>% .[, `:=`(quarantine_time = identification_time,
                       isolation_time = ifelse(!asym, quarantine_time,
                                               ifelse(test_positive, quarantine_time + test_delay,
                                                      Inf)),
                       trace_init_time = ifelse(!asym & rep(trace_neg_symptomatic, nrow(.)),
                                                quarantine_time,
                                                ifelse(test_positive, quarantine_time + test_delay,
                                                       Inf)))]
}

create_index_cases <- function(n_initial_cases, p_asymptomatic, test_time,
                               p_ident_sym, p_smartphone_overall,
                               r0_asymptomatic, r0_symptomatic,
                               incubation_time, p_traced_auto,
                               p_traced_manual, trace_time_auto,
                               trace_time_manual, recovery_time,
                               test_function, rollout_delay_gen,
                               rollout_delay_days, delay_time,
                               trace_neg_symptomatic){
    #' Set up a table of index cases
    # Set immutable keys
    cases <- set_primary_immutable_index(n_initial_cases, p_asymptomatic,
                                         test_time, p_ident_sym)
    cases <- set_secondary_immutable(cases, TRUE, p_smartphone_overall,
                                     NA, NA, NA,
                                     dispersion, r0_asymptomatic, r0_symptomatic,
                                     incubation_time, p_traced_auto, p_traced_manual,
                                     trace_time_auto, trace_time_manual,
                                     recovery_time, test_function)
    # Set identification time (inc. rollout delay) and parent mutables
    cases <- cases %>% .[, `:=`(parent_ident_time = NA,
                                parent_quar_time = NA,
                                parent_iso_time = NA,
                                parent_trace_init_time = NA,
                                identification_time = ifelse(rep(rollout_delay_gen > 0, nrow(.)), Inf,
                                                             ifelse(recovery < rollout_delay_days, Inf,
                                                                    onset_true + delay_time(nrow(.)))))]
    # Set other mutable keys
    cases <- set_nonparental_mutables(cases, trace_neg_symptomatic)
    return(cases)
}











#----------------------------------------------------------------------------
# Auxiliary functions
#----------------------------------------------------------------------------

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
    return(cases[exposure < infector_iso_time | escapes_isolation == TRUE])
}

trace_forward <- function(cases, quarantine, sero_test){
    #' Determine isolation times for secondary cases based on forward tracing
    # TODO: Pass parent and child tables separately, merge in-function
    # TODO: Generalise infector_iso_time to general fwd_tracing_time (to account for delays etc)
    # TODO: Distinguish quarantine time from isolation time; tracing determines the former
    cases %>% copy %>% .[, isolation_time := ifelse(!traced_fwd, isolation_time,
                                                    pmin(isolation_time,
                                                         ifelse(rep(!quarantine && !sero_test, nrow(.)),
                                                         pmax(onset_true, infector_iso_time),
                                                         infector_iso_time)))] %>% return
    # - Default untraced isolation time depends on syptomatic status
    # - If missed in forward tracing, isolation time defaults to untraced time
    # - Otherwise, isolation time is the minimum of the untraced and traced times
    # - If quarantine in place, traced isolation time equals infector isolation time
    # - Otherwise, traced isolation time is the maximum of infector isolation time and symptom onset
    # - For asymptomatic individuals, this means isolation time is infinite unless quarantined/tested
    # - For individuals with untraced asymptomatic parents (i.e. parents with infinite onset),
    #   this produces the untraced time in all cases (TODO: test this)
    # - If an individual has already been traced in the past (e.g. through backtracing),
    #   re-running trace_forward does not reset its isolation time to the untraced time (TODO: Consider this more thoroughly)
}

run_backtrace_update_parent_isolation <- function(cases){
    #' Update parent isolation times in case table
    # Separate out individuals whose parents aren't in table (e.g. gen-0-ers)
    cases_no_parent <- cases[!(infector %in% cases$case_id)]
    cases_parent <- cases[infector %in% cases$case_id]
    # Get updated parental isolation times for latter
    parent_iso_indices <- match(cases_parent$infector, cases$case_id)
    parent_isos <- cases$isolation_time[parent_iso_indices]
    cases_parent[, infector_iso_time := parent_isos]
    return(rbind(cases_no_parent, cases_parent, fill=TRUE))
}

run_backtrace_update_filter <- function(cases_merged, gen_earliest,
                                        quarantine, sero_test){
    #' Perform update, filter and forward tracing steps during backtracing
    cases_updated <- run_backtrace_update_parent_isolation(cases_merged)
    # Remove entries whose exposure now succeeds parent isolation time
    cases_filtered <- filter_new_cases(cases_updated)
    # Remove entries whose parents are no longer in dataset
    nodes_orphaned <- cases_filtered %>% 
        .[generation > gen_earliest & !(infector %in% .$case_id)] %>% .$case_id
    while (length(nodes_orphaned) > 0){
        cases_filtered <- cases_filtered[!(case_id %in% nodes_orphaned)]
        nodes_orphaned <- cases_filtered %>% 
            .[generation > gen_earliest & !(infector %in% .$case_id)] %>% .$case_id
    }
    # Forward trace remaining entries
    cases_retraced <- trace_forward(cases_filtered, quarantine, sero_test)
    setindex(cases_retraced, "infector")
    # Test for inequality with input and return
    changed <- !isTRUE(all.equal(cases_merged, cases_retraced))
    return(list(cases=cases_retraced, changed=changed))
}

run_backtrace_iter <- function(cases_iter, gen_earliest, quarantine, sero_test){
    #' Perform backtracing on all cases in a (pre-filtered) dataframe
    # 1. Look for everyone preceding their parent in isolation time
    #    and filter by backtrace success (and whether parent is in cases at all)
    cases_pre_filtered <- cases_iter %>% 
        .[traced_rev == TRUE & isolation_time < infector_iso_time &
              infector %in% .$case_id]
    if (nrow(cases_pre_filtered) == 0) return(list(cases=cases_iter, changed=FALSE))
    # 2. Determine new putative parental isolation times
    cases_iso_new <- cases_pre_filtered %>%
        .[, infector_iso_new := pmin(infector_iso_time,
                                     ifelse(rep(!quarantine && !sero_test, nrow(.)),
                                            pmax(infector_onset_true, isolation_time),
                                            isolation_time))] %>%
        .[infector_iso_new < infector_iso_time] # Fairly sure this filter is pointless, but JIC
    if (nrow(cases_iso_new) == 0) return(list(cases=cases_iter, changed=FALSE))
    # 3. New parent isolation time becomes min across all successful backtraces
    parents_backtraced <- cases_iso_new[, .(infector_iso_new = min(infector_iso_new)),
                                        by = "infector"]
    parents_backtraced_renamed <- parents_backtraced %>% copy %>%
        setnames(c("infector", "infector_iso_new"), c("case_id", "isolation_time_new"))
    # 4. Update isolation times in parent entries
    cases_merged <- merge(cases_iter, parents_backtraced_renamed, by="case_id", all=TRUE) %>%
        .[, isolation_time := pmin(isolation_time, isolation_time_new, na.rm=TRUE)] %>%
        .[, isolation_time_new := NULL]
    # 5. Iterate update, filter and forward-trace steps until convergence
    cases_updated <- run_backtrace_update_filter(cases_merged, gen_earliest, 
                                                 quarantine, sero_test)
    while(cases_updated$changed){
        cases_updated <- run_backtrace_update_filter(cases_updated$cases, gen_earliest,
                                                     quarantine, sero_test)
    }
    # Test for inequality with input and return
    setindex(cases_updated$cases, "infector")
    changed <- !isTRUE(all.equal(cases_iter, cases_updated$cases))
    return(list(cases=cases_updated$cases, changed=changed))
}

run_backtrace_nstep <- function(cases_iso, quarantine, sero_test,
                                backtrace_distance){
    #' Perform backtracing up to some specified number of generations before
    #' the most recent generation.
    if (backtrace_distance == 0) return(cases_iso)
    # Determine which generations to act on and separate out others
    gen_latest <- max(cases_iso$generation)
    gen_earliest <- max(0, gen_latest - backtrace_distance)
    cases_old <- cases_iso[generation < gen_earliest] # Keep unchanged for later
    cases_act <- cases_iso[generation >= gen_earliest]
    # Run first iteration of backtracing
    setindex(cases_act, "infector")
    cases_iter <- run_backtrace_iter(cases_act, gen_earliest, quarantine, sero_test)
    while(cases_iter$changed){
        cases_iter <- run_backtrace_iter(cases_iter$cases, gen_earliest, quarantine, sero_test)
    }
    return(rbind(cases_old, cases_iter$cases, fill=TRUE))
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
} # TODO: Collapse terminate functions into one function with a switch

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
    ifelse(p_asymptomatic==1, 0,
           (r0_base - r0_asymptomatic*p_asymptomatic)/(1-p_asymptomatic))
}

#----------------------------------------------------------------------------
# Core functions
#----------------------------------------------------------------------------

outbreak_setup <- function(n_initial_cases, p_asymptomatic,
                           incubation_time, delay_time, p_isolation,
                           rollout_delay_generations, rollout_delay_days){
    #' Set up a table of initial cases
    cases <- data.table(infector = 0, infector_onset_gen = Inf, 
                        infector_onset_true = Inf,
                        infector_iso_time = Inf, infector_asym = NA, 
                        generation = 0, processed = FALSE, n_children = NA,
                        exposure = 0, case_id = 1:n_initial_cases,
                        traced_fwd = FALSE, traced_rev = FALSE) %>%
        .[, `:=`(asym = purrr::rbernoulli(nrow(.), p_asymptomatic), # Asymptomatic?
                 escapes_isolation = purrr::rbernoulli(n = nrow(.), p = 1-p_isolation),
                 # Symptom onset for determining generation time (TODO: refactor gentime so can drop)
                 onset_gen = incubation_time(nrow(.)))] %>%
        # True symptom onset (infinite for asymptomatic individuals)
        .[, onset_true := ifelse(asym, Inf, onset_gen)] %>%
        # Maximum isolation time (if untraced)
        .[, isolation_time := pmax(onset_true, rollout_delay_days) +
              delay_time(nrow(.))]
    if (rollout_delay_generations > 0) cases$isolation_time <- Inf
    return(cases)
}

outbreak_step <- function(case_data = NULL, dispersion = NULL,
                          r0_symptomatic = NULL, r0_asymptomatic = NULL,
                          p_asymptomatic = NULL, p_traced = NULL,
                          generation_time = NULL, incubation_time = NULL,
                          delay_time = NULL, backtrace_distance = NULL,
                          quarantine = NULL, sero_test = NULL,
                          p_isolation = NULL, rollout_delay_generations = NULL,
                          rollout_delay_days = NULL){
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
                                            incubation_time, delay_time, p_isolation,
                                            rollout_delay_days,
                                            rollout_delay_generations)
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
    cases_out <- run_backtrace_nstep(cases_iso, quarantine, sero_test, 
                                     backtrace_distance)
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
                           cap_cases = NULL, rollout_delay_generations = NULL,
                           rollout_delay_days = NULL){
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
                                incubation_time, delay_time, p_isolation,
                                rollout_delay_generations,
                                rollout_delay_days)
    # Run outbreak loop
    while (latest_exposure_weeks < cap_max_weeks & total_cases < cap_cases &
           !outbreak_extinct & latest_generation <= cap_max_generations){
        out <- outbreak_step(case_data = case_data, dispersion = dispersion,
                                       r0_symptomatic = r0_symptomatic,
                                       r0_asymptomatic = r0_asymptomatic,
                                       p_asymptomatic = p_asymptomatic, p_traced = p_traced,
                                       generation_time = generation_time,
                                       incubation_time = incubation_time,
                                       delay_time = delay_time, p_isolation = p_isolation,
                                       quarantine = quarantine, sero_test = sero_test,
                                       backtrace_distance = backtrace_distance,
                             rollout_delay_generations = rollout_delay_generations,
                             rollout_delay_days = rollout_delay_days)
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

