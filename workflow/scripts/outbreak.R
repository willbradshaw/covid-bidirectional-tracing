#' Core functionality for branching-chain outbreak model

#----------------------------------------------------------------------------
# Create new cases
#----------------------------------------------------------------------------

map2_tab <- function(tab, col1, col2, fn){
    unlist(purrr::map2(tab[[col1]], tab[[col2]], function(x, y) fn(x,y)))
}

rep_new_cases <- function(tab, col){
    map2_tab(tab, col, "n_children", function(x,y) rep(x, as.integer(y)))
}

set_primary_immutables_index <- function(n_initial_cases, p_asymptomatic,
                                         test_time, test_sensitivity,
                                         test_serological, p_ident_sym,
                                         p_compliance_isolation,
                                         p_data_sharing_auto, p_data_sharing_manual
                                         ){
    #' Set primary immutable keys for index cases
    data.table(case_id = 1:n_initial_cases,
               asym = purrr::rbernoulli(n_initial_cases, p_asymptomatic),
               environmental = TRUE, # Index cases can never be traced to/from their parents (who don't exist)
               sero_test = test_serological,
               test_delay = test_time(n_initial_cases), # TODO: Attach testing to child? (would allow multiple testing during reverse tracing)
               testable = purrr::rbernoulli(n_initial_cases, test_sensitivity), # TODO: See above
               ident_sym = purrr::rbernoulli(n_initial_cases, p_ident_sym),
               isolates = purrr::rbernoulli(n_initial_cases, p_compliance_isolation),
               shares_data_auto = purrr::rbernoulli(n_initial_cases, p_data_sharing_auto),
               shares_data_manual = purrr::rbernoulli(n_initial_cases, p_data_sharing_manual),
               generation = 0, infector_id = 0, infector_onset_gen = Inf,
               infector_exposure = Inf,
               infector_onset_true = Inf, infector_asym = NA, infector_recovery = Inf,
               infector_has_smartphone = NA, infector_shares_data_auto = FALSE,
               infector_shares_data_manual = FALSE, infector_isolates = FALSE)
}

set_primary_immutables_child <- function(parents, p_asymptomatic, p_compliance_isolation,
                                         test_time, test_serological,
                                         test_sensitivity, p_ident_sym, p_environmental,
                                         p_data_sharing_auto, p_data_sharing_manual
                                         ){
    #' Set primary immutable keys for index cases
    n_children_total <- sum(parents$n_children)
    case_ids <- max(parents$case_id) + 1:n_children_total
    child_cases <- data.table(case_id = case_ids,
                              asym = purrr::rbernoulli(n_children_total, p_asymptomatic),
                              environmental = purrr::rbernoulli(n_children_total, p_environmental),
                              test_delay = test_time(n_children_total),
                              sero_test = test_serological,
                              testable = purrr::rbernoulli(n_children_total, test_sensitivity),
                              ident_sym = purrr::rbernoulli(n_children_total, p_ident_sym),
                              isolates = purrr::rbernoulli(n_children_total, p_compliance_isolation),
                              shares_data_auto = purrr::rbernoulli(n_children_total, p_data_sharing_auto),
                              shares_data_manual = purrr::rbernoulli(n_children_total, p_data_sharing_manual),
                              generation = rep_new_cases(parents, "generation") + 1, # New generation number
                              infector_id = rep_new_cases(parents, "case_id"), # Parent case
                              infector_exposure = rep_new_cases(parents, "exposure"), # Exposure of parent (not used, but useful to know)
                              infector_onset_gen = rep_new_cases(parents, "onset_gen"), # Symptom onset of parent
                              infector_onset_true = rep_new_cases(parents, "onset_true"), # Symptom onset of parent
                              infector_asym = rep_new_cases(parents, "asym"), # Is parent asymptomatic?
                              infector_recovery = rep_new_cases(parents, "recovery"), # Recovery time of parent
                              infector_has_smartphone = rep_new_cases(parents, "has_smartphone"), # Does parent have smartphone?
                              infector_shares_data_auto = rep_new_cases(parents, "shares_data_auto"), # Does parent case share their contact data for tracing?
                              infector_shares_data_manual = rep_new_cases(parents, "shares_data_manual"), # Does parent case share their contact data for tracing?
                              infector_isolates = rep_new_cases(parents, "isolates")) # Does parent comply with isolation?
    return(child_cases)
}

set_secondary_immutables <- function(cases, index, p_smartphone_overall,
                                     p_smartphone_infector_yes, p_smartphone_infector_no,
                                     generation_time, n_children_fn,
                                     trace_neg_symptomatic,
                                     incubation_time, p_traced_auto, p_traced_manual,
                                     trace_time_auto, trace_time_manual,
                                     recovery_time, data_limit_auto, data_limit_manual,
                                     contact_limit_auto, contact_limit_manual){
    #' Set secondary immutable keys
    cases %>% .[, `:=`(has_smartphone = purrr::rbernoulli(nrow(.),
                                                          ifelse(rep(index, nrow(.)), p_smartphone_overall,
                                                                 ifelse(infector_has_smartphone,
                                                                        p_smartphone_infector_yes,
                                                                        p_smartphone_infector_no))),
                       exposure = ifelse(rep(index, nrow(.)), 0,
                                         generation_time(infector_onset_gen)),
                       n_children = n_children_fn(asym),
                       trace_if_neg = ifelse(asym, FALSE, trace_neg_symptomatic),
                       processed = FALSE)] %>% # (Not actually immutable, but independent of tracing)
        .[, `:=`(auto_traced_fwd = rep(!index, nrow(.)) & infector_has_smartphone & has_smartphone & 
                     infector_shares_data_auto & !environmental &
                     purrr::rbernoulli(nrow(.), p_traced_auto),
                 auto_traced_rev = rep(!index, nrow(.)) & infector_has_smartphone & has_smartphone & 
                     shares_data_auto & !environmental &
                     purrr::rbernoulli(nrow(.), p_traced_auto),
                 onset_gen = exposure + incubation_time(nrow(.)))] %>%
        .[, `:=`(manual_traced_fwd = (rep(!index, nrow(.)) & !auto_traced_fwd & !environmental &
                     infector_shares_data_manual & purrr::rbernoulli(nrow(.), p_traced_manual)),
                 manual_traced_rev = (rep(!index, nrow(.)) & !auto_traced_rev & !environmental &
                     shares_data_manual & purrr::rbernoulli(nrow(.), p_traced_manual)))] %>%
        .[, `:=`(traceable_fwd = (auto_traced_fwd | manual_traced_fwd),
                 traceable_rev = (auto_traced_rev | manual_traced_rev),
                 trace_delay_fwd = ifelse(auto_traced_fwd, trace_time_auto(nrow(.)),
                                          ifelse(manual_traced_fwd, trace_time_manual(nrow(.)), Inf)),
                 trace_delay_rev = ifelse(auto_traced_rev, trace_time_auto(nrow(.)),
                                          ifelse(manual_traced_rev, trace_time_manual(nrow(.)), Inf)),
                 onset_true = ifelse(asym, Inf, onset_gen),
                 recovery = recovery_time(onset_gen),
                 data_limit_fwd = ifelse(auto_traced_fwd, data_limit_auto, data_limit_manual),
                 data_limit_rev = ifelse(auto_traced_rev, data_limit_auto, data_limit_manual),
                 contact_limit_fwd = ifelse(auto_traced_fwd, contact_limit_auto, contact_limit_manual),
                 contact_limit_rev = ifelse(auto_traced_rev, contact_limit_auto, contact_limit_manual)
        )]
}

initialise_parental_mutables <- function(cases, parents){
    #' Set initial values of parental mutables for child cases
    cases %>% .[, `:=`(infector_ident_time = rep_new_cases(parents, "identification_time"),
                       infector_quar_time = rep_new_cases(parents, "quarantine_time"),
                       infector_iso_time = rep_new_cases(parents, "isolation_time"),
                       infector_trace_init_time = rep_new_cases(parents, "trace_init_time")
                       )] %>%
        .[, `:=`(in_data_threshold_fwd = (exposure + data_limit_fwd >= infector_trace_init_time),
                 in_contact_threshold_fwd = (exposure + contact_limit_fwd >= pmin(infector_onset_true,
                                                                                  infector_trace_init_time))
                 )]
}

initialise_identification_time <- function(cases, rollout_delay_gen,
                                           rollout_delay_days, delay_time){
    #' Compute initial identification time of cases, before any tracing
    #' (Depends on ident_sym, rollout delay, symptomatic status)
    n_cases <- nrow(cases)
    cases_identified <- cases$ident_sym &
        (cases$generation >= rollout_delay_gen) &
        (cases$recovery >= rollout_delay_days)
    cases <- cases[, identification_time := ifelse(!cases_identified, Inf,
                                                   onset_true + delay_time(n_cases))]
    return(cases)
}

set_nonparental_mutables <- function(cases){
    #' Set (or update) values of nonparental mutable keys (except identification time)
    cases %>% .[, `:=`(quarantine_time = identification_time,
                       test_positive = testable & (identification_time > exposure) &
                           (sero_test | (identification_time < recovery)))] %>%
        .[, `:=`(isolation_time = ifelse(onset_true <= quarantine_time, quarantine_time, # If already showing symptoms, isolate immediately
                                         ifelse(test_positive, pmin(onset_true, quarantine_time + test_delay), # Else if test positive, isolate on test result (or symptom onset, if that comes first)
                                                onset_true)), # Otherwise, isolate on symptom onset (which might be never)
                 trace_init_time = ifelse(trace_if_neg, # If tracing symptomatics regardless of test result, then...
                                          ifelse(onset_true <= quarantine_time, quarantine_time, # If already showing symptoms, trace immediately
                                                 ifelse(test_positive, pmin(onset_true, quarantine_time + test_delay), # Otherwise, trace as soon as either symptoms appear or test comes back positive
                                                        onset_true)), # If test comes back negative, trace on symptom onset
                                          ifelse(test_positive, quarantine_time + test_delay, Inf)))] %>% # If only tracing symptomatics that test positive, trace on positive test result, otherwise never
        .[, `:=`(in_data_threshold_rev = (exposure + data_limit_rev >= trace_init_time),
                 in_contact_threshold_rev = exposure + contact_limit_rev >= pmin(onset_true,
                                                                                 trace_init_time)
                 )]
}

create_index_cases <- function(n_initial_cases, p_asymptomatic, test_time,
                               test_sensitivity, test_serological, p_ident_sym, 
                               p_smartphone_overall, n_children_fn,
                               trace_neg_symptomatic,
                               incubation_time, p_traced_auto,
                               p_traced_manual, trace_time_auto,
                               trace_time_manual, recovery_time,
                               data_limit_auto, data_limit_manual,
                               contact_limit_auto, contact_limit_manual,
                               rollout_delay_gen, p_data_sharing_auto,
                               p_data_sharing_manual, p_compliance_isolation,
                               rollout_delay_days, delay_time){
    #' Set up a table of index cases
    # Set immutable keys
    cases <- set_primary_immutables_index(n_initial_cases = n_initial_cases,
                                          p_asymptomatic = p_asymptomatic,
                                          test_time = test_time,
                                          test_sensitivity = test_sensitivity,
                                          test_serological = test_serological,
                                          p_ident_sym = p_ident_sym,
                                          p_compliance_isolation = p_compliance_isolation,
                                          p_data_sharing_auto = p_data_sharing_auto,
                                          p_data_sharing_manual = p_data_sharing_manual)
    cases <- set_secondary_immutables(cases = cases, index = TRUE, p_smartphone_overall = p_smartphone_overall,
                                      p_smartphone_infector_yes = NA, p_smartphone_infector_no = NA,
                                      generation_time = generation_time, n_children_fn = n_children_fn,
                                      trace_neg_symptomatic = trace_neg_symptomatic,
                                      incubation_time = incubation_time, p_traced_auto = p_traced_auto, 
                                      p_traced_manual = p_traced_manual, trace_time_auto = trace_time_auto,
                                      trace_time_manual = trace_time_manual, recovery_time = recovery_time,
                                      data_limit_auto = data_limit_auto, data_limit_manual = data_limit_manual,
                                      contact_limit_auto = contact_limit_auto,
                                      contact_limit_manual = contact_limit_manual)
    # Set identification time (dependent on rollout delay, ident_sym) and parent mutables
    cases <- cases %>% .[, `:=`(infector_ident_time = Inf,
                                infector_quar_time = Inf,
                                infector_iso_time = Inf,
                                infector_trace_init_time = Inf)]
    cases <- initialise_identification_time(cases = cases, delay_time = delay_time,
                                            rollout_delay_gen = rollout_delay_gen, 
                                            rollout_delay_days = rollout_delay_days)
    # Set other mutable keys
    cases <- set_nonparental_mutables(cases)
    return(cases)
}

create_child_cases <- function(parents, p_asymptomatic, p_compliance_isolation,
                               test_time, test_sensitivity, test_serological,
                               p_ident_sym, p_smartphone_infector_yes, p_smartphone_infector_no,
                               generation_time, n_children_fn, trace_neg_symptomatic,
                               incubation_time, p_traced_auto,
                               p_traced_manual, trace_time_auto, trace_time_manual,
                               recovery_time, data_limit_auto, data_limit_manual,
                               contact_limit_auto, contact_limit_manual,
                               rollout_delay_gen, p_data_sharing_auto,
                               p_data_sharing_manual,
                               rollout_delay_days, delay_time, p_environmental){
    #' Generate a table of new child cases from a table of parent cases
    # Set immutable keys
    cases <- set_primary_immutables_child(parents = parents, p_asymptomatic = p_asymptomatic,
                                          p_compliance_isolation = p_compliance_isolation,
                                          test_time = test_time, test_serological = test_serological,
                                          test_sensitivity = test_sensitivity, p_ident_sym = p_ident_sym,
                                          p_environmental = p_environmental,
                                          p_data_sharing_auto = p_data_sharing_auto,
                                          p_data_sharing_manual = p_data_sharing_manual)
    cases <- set_secondary_immutables(cases = cases, index = FALSE, p_smartphone_overall = NA,
                                      p_smartphone_infector_yes = p_smartphone_infector_yes,
                                      p_smartphone_infector_no = p_smartphone_infector_no,
                                      generation_time = generation_time, n_children_fn = n_children_fn,
                                      trace_neg_symptomatic = trace_neg_symptomatic,
                                      incubation_time = incubation_time, p_traced_auto = p_traced_auto, 
                                      p_traced_manual = p_traced_manual, trace_time_auto = trace_time_auto,
                                      trace_time_manual = trace_time_manual, recovery_time = recovery_time,
                                      data_limit_auto = data_limit_auto, data_limit_manual = data_limit_manual,
                                      contact_limit_auto = contact_limit_auto,
                                      contact_limit_manual = contact_limit_manual)
    # Initialise parental mutables
    cases <- initialise_parental_mutables(cases, parents)
    # Initialise identification time
    cases <- initialise_identification_time(cases = cases, delay_time = delay_time,
                                            rollout_delay_gen = rollout_delay_gen, 
                                            rollout_delay_days = rollout_delay_days)
    # Initialise other nonparental mutables
    cases <- set_nonparental_mutables(cases)
    # Remove (hopefully very few) cases that were exposed after parental recovery time
    cases <- cases[exposure < infector_recovery]
    return(cases)
}

filter_new_cases <- function(cases){
    #' Filter new cases based on quarantine/isolation times of parent
    # Determine which cases were exposed during parental quarantine/isolation
    exposure_in_isolation <- cases$exposure > cases$infector_iso_time
    exposure_in_quarantine <- (cases$exposure > cases$infector_quar_time) &
        !exposure_in_isolation
    # Determine which such cases are blocked
    blocked_by_isolation <- exposure_in_isolation & cases$infector_isolates
    blocked_by_quarantine <- exposure_in_quarantine & cases$infector_isolates
    # TODO: Distinguish compliance with isolation from compliance with quarantine?
    # Return non-blocked cases
    return(cases[(!blocked_by_isolation) & (!blocked_by_quarantine)])
}

#----------------------------------------------------------------------------
# Contact tracing
#----------------------------------------------------------------------------

trace_forward <- function(cases){
    #' Update mutable keys of cases based on forward tracing
    # Determine which contacts are forward traced (based on baseline
    # traceability and data limits)
    traced_fwd <- cases$traceable_fwd & 
        cases$in_data_threshold_fwd  &
        cases$in_contact_threshold_fwd
    # Update identification time of traced contacts
    cases_traced <- cases %>% copy %>%
        .[, identification_time := ifelse(!traced_fwd, identification_time,
                                          pmin(identification_time, 
                                               infector_trace_init_time + trace_delay_fwd))]
    # Update other mutables based on identification time
    cases_updated <- set_nonparental_mutables(cases_traced)
    return(cases_updated)
}

trace_reverse <- function(cases, backtrace_distance){
    #' Perform reverse tracing up to some specified number of generations before
    #' the most recent generation.
    if (backtrace_distance == 0) return(cases)
    #if (backtrace_distance == 1) return(trace_reverse_1step(cases))
    # Determine which cases to act on and separate out others
    gen_latest <- max(cases$generation)
    gen_earliest <- max(0, gen_latest - backtrace_distance)
    cases_old <- cases[generation < gen_earliest] # Keep unchanged for later
    cases_act <- cases[generation >= gen_earliest]
    # Run first iteration of reverse tracing
    setindex(cases_act, "infector_id")
    cases_iter <- trace_reverse_iter(cases_act, gen_earliest)
    # Repeat until convergence
    while (cases_iter$changed){
        cases_iter <- trace_reverse_iter(cases_iter$cases, gen_earliest)
    }
    # Combine with untraced generations and return
    return(rbind(cases_old, cases_iter$cases, fill=TRUE))
}

trace_reverse_iter <- function(cases, gen_earliest){
    #' Perform one iteration of reverse tracing on all cases in a 
    #' pre-filtered dataset
    # 1. Determine which contacts are meaningfully reverse-traced (based on baseline
    #    traceability, whether the trace occurs before the parent is already
    #    identified, data limits, and whether parent is in the dataset at all)
    # (NB: If we later attach testing in reverse-tracing to the child case, 
    #      will need to remove fourth condition here)
    traced_rev <- cases$traceable_rev & 
        cases$in_data_threshold_rev  &
        cases$in_contact_threshold_rev &
        (cases$trace_init_time + cases$trace_delay_rev < cases$infector_ident_time) &
        (cases$infector_id %in% cases$case_id)
    if (sum(traced_rev) == 0) return(list(cases=cases, changed=FALSE))
    # 2. Filter to successful reverse traces and determine new putative
    #    parental identification times
    # (NB: If we later attach testing in reverse-tracing to the child case,
    #      changes will be needed here)
    cases_ident_new <- cases %>% copy %>% .[traced_rev == TRUE] %>%
        .[, infector_ident_new := pmin(infector_ident_time,
                                       trace_init_time + trace_delay_rev)]
    # 3. New parent identification time becomes min across all successful backtraces
    parents_backtraced <- cases_ident_new[, .(infector_ident_new = min(infector_ident_new)),
                                        by = "infector_id"]
    parents_backtraced_renamed <- parents_backtraced %>% copy %>%
        setnames(c("infector_id", "infector_ident_new"),
                 c("case_id", "identification_time_new"))
    # 4. Update identification times in parent entries, then update other mutables
    cases_merged <- merge(cases, parents_backtraced_renamed, by="case_id", all=TRUE) %>%
        .[, identification_time := pmin(identification_time, identification_time_new, 
                                        na.rm=TRUE)] %>%
        .[, identification_time_new := NULL]
    cases_parents_updated <- set_nonparental_mutables(cases_merged)
    # 5. Iterate child-update, filter and forward-trace steps until convergence
    cases_updated <- trace_reverse_update_filter(cases_parents_updated, gen_earliest)
    while (cases_updated$changed){
        cases_updated <- trace_reverse_update_filter(cases_updated$cases, gen_earliest)
    }
    # 6. Test for inequality with input and return
    setindex(cases_updated$cases, "infector_id")
    changed <- !isTRUE(all.equal(cases, cases_updated$cases))
    return(list(cases=cases_updated$cases, changed=changed))
}

trace_reverse_update_filter <- function(cases, gen_earliest){
    #' Perform update, filter and forward tracing steps during backtracing
    # TODO: Some sort of early break?
    # 1. Update parental mutables
    cases_updated <- update_parental_mutables(cases)
    # 2. Filter child cases (based on parental quarantine and isolation)
    cases_filtered <- filter_new_cases(cases_updated)
    # 3. Remove entries whose parents are no longer in the dataset
    cases_non_orphaned <- filter_orphaned_cases(cases_filtered, gen_earliest)
    # 4. Forward-trace remaining entries
    cases_retraced <- trace_forward(cases_non_orphaned)
    # Test for inequality with input and return
    setindex(cases_retraced, "infector_id")
    changed <- !isTRUE(all.equal(cases, cases_retraced))
    return(list(cases=cases_retraced, changed=changed))
}

filter_orphaned_cases <- function(cases, gen_earliest){
    #' Remove child cases whose ancestors are no longer in the dataset
    orphaned <- (cases$generation > gen_earliest) &
        !(cases$infector_id %in% cases$case_id)
    while (sum(orphaned) > 0){
        cases <- cases[!orphaned] # Remove orphaned nodes
        orphaned <- (cases$generation > gen_earliest) &
            !(cases$infector_id %in% cases$case_id) # Re-test orphan status
    }
    return(cases)
}

update_parental_mutables <- function(cases){
    #' Update parental mutable keys for child cases during reverse-tracing
    # 1. Separate out individuals whose parents aren't in data
    parent_in_data <- cases$infector_id %in% cases$case_id
    if (sum(parent_in_data) == 0) return(cases)
    cases_no_parent <- cases[parent_in_data==FALSE]
    cases_parent <- cases[parent_in_data==TRUE]
    # 2. Find indices of parent cases in dataset
    parent_indices <- match(cases_parent$infector_id, cases$case_id)
    # 3. Update primary parental mutables
    cases_parent[, `:=`(infector_ident_time = cases$identification_time[parent_indices],
                        infector_quar_time = cases$quarantine_time[parent_indices],
                        infector_iso_time = cases$isolation_time[parent_indices],
                        infector_trace_init_time = cases$isolation_time[parent_indices])]
    # 4. Update secondary parental mutables
    cases_parent[, `:=`(in_data_threshold_fwd = (exposure + data_limit_fwd >= infector_trace_init_time),
                        in_contact_threshold_fwd = exposure + contact_limit_fwd >= pmin(infector_onset_true,
                                                                                        infector_trace_init_time))]
    
    # 5. Combine with parentless cases and return
    return(rbind(cases_no_parent, cases_parent, fill=TRUE))
}

#----------------------------------------------------------------------------
# Other auxiliary functions
#----------------------------------------------------------------------------

compute_symptomatic_r0 <- function(r0_base, r0_asymptomatic, p_asymptomatic){
    ifelse(p_asymptomatic==1, 0,
           (r0_base - r0_asymptomatic*p_asymptomatic)/(1-p_asymptomatic))
}

compute_p_smartphone_infector_no <- function(p_smartphone_overall,
                                           p_smartphone_infector_yes){
    #' Compute the probability that a contact of a non-smartphone-haver
    #' has a smartphone (using Markov chain assumption)
    if(p_smartphone_infector_yes==1) return(0)
    return((p_smartphone_overall/(1-p_smartphone_overall)) * 
               (1-p_smartphone_infector_yes))
}

#----------------------------------------------------------------------------
# Outer simulation functions
#----------------------------------------------------------------------------

outbreak_step <- function(case_data = NULL,
                          child_case_fn = NULL,
                          backtrace_distance = NULL){
    #' Move forward one generation in the branching process
    # 1. Separate current parents from previous generations and check for new cases
    case_data_old <- case_data[processed == TRUE]
    case_data_new <- case_data[processed == FALSE]
    case_data_new[, processed := TRUE] # After split, set processed to TRUE
    if (sum(case_data_new$n_children) == 0) return(rbind(case_data_old, case_data_new))
    # 2. Generate putative secondary cases
    children_putative <- child_case_fn(case_data_new)
    # 3. Filter secondary cases based on isolation time of parents
    children_filtered <- filter_new_cases(children_putative)
    # 4. Perform forward tracing on secondary cases
    children_traced <- trace_forward(children_filtered)
    # 5. Perform reverse tracing on entire case dataset, then return
    cases_all <- data.table::rbindlist(list(case_data_old, case_data_new, children_traced),
                                       fill=TRUE)
    cases_out <- trace_reverse(cases_all, backtrace_distance)
    #gc(verbose=FALSE, full=TRUE)
    return(cases_out)
}

run_outbreak <- function(index_case_fn = NULL, child_case_fn = NULL,
                         cap_max_generations = NULL, cap_max_weeks = NULL,
                         cap_cases = NULL, backtrace_distance = NULL){
    #' Run a single complete instance of the branching-process model
    # Set up initial cases
    case_data <- index_case_fn()
    # Run outbreak loop
    while (max(case_data$exposure)/7 < cap_max_weeks & # Time limit
           nrow(case_data) < cap_cases & # Cumulative case limit (TODO: Change to simultaneous/generational limit?)
           max(case_data$generation) < cap_max_generations & # Generation limit
           any(!case_data$processed)){ # Extinction condition
        case_data <- outbreak_step(case_data = case_data,
                                   child_case_fn = child_case_fn,
                                   backtrace_distance = backtrace_distance)
        #gc(verbose=FALSE, full=TRUE)
    }
    gc(verbose=FALSE, full=TRUE)
    return(case_data)
}

scenario_sim <- function(n_iterations = NULL, dispersion = NULL, r0_base = NULL,
                         rel_r0_asymptomatic = NULL, p_asymptomatic = NULL,
                         generation_omega = NULL, generation_alpha = NULL,
                         recovery_quantile = NULL, incubation_time = NULL,
                         test_time = NULL, trace_time_auto = NULL,
                         trace_time_manual = NULL,
                         delay_time = NULL, n_initial_cases = NULL,
                         test_sensitivity = NULL, test_serological = NULL, 
                         p_ident_sym = NULL, p_smartphone_overall = NULL,
                         p_smartphone_link = NULL, trace_neg_symptomatic = NULL,
                         p_traced_auto = NULL, p_traced_manual = NULL,
                         data_limit_auto = NULL, data_limit_manual = NULL,
                         contact_limit_auto = NULL,
                         contact_limit_manual = NULL,
                         rollout_delay_gen = NULL, rollout_delay_days = NULL,
                         p_compliance_isolation = NULL,
                         cap_max_generations = NULL, cap_max_weeks = NULL,
                         cap_cases = NULL, backtrace_distance = NULL,
                         p_environmental = NULL, report = NULL,
                         p_data_sharing_auto = NULL,
                         p_data_sharing_manual = NULL,
                         scenario = NULL
                         ){
    #' Run a specified number of outbreaks with identical parameters
    if (report){
        start <- proc.time()
        cat("Scenario ", scenario, " beginning at: ", date(), "\n", sep="")
    }
    # Compute auxiliary parameters
    r0_asymptomatic <- rel_r0_asymptomatic * r0_base
    r0_symptomatic <- compute_symptomatic_r0(r0_base, r0_asymptomatic, p_asymptomatic)
    p_smartphone_infector_yes <- ifelse(p_smartphone_link < 0, p_smartphone_overall,
                                        p_smartphone_link)
    p_smartphone_infector_no <- compute_p_smartphone_infector_no(p_smartphone_overall,
                                                             p_smartphone_infector_yes)
    # Prepare fixed-shape distributions
    n_children_fn <- function(asym) rnbinom(n=length(asym), size=dispersion,
                                            mu=ifelse(asym, r0_asymptomatic,
                                                      r0_symptomatic))
    generation_time <- function(onsets){
        sn::rsn(n = length(onsets), xi = onsets, omega = generation_omega,
                alpha = generation_alpha) %>%  ifelse(. < 0, 0, .)
    } # "generation_alpha" previously known as "k" or "generation_k"
    recovery_time <- function(onsets){
        sn::qsn(recovery_quantile, xi = onsets, omega = generation_omega,
                alpha = generation_alpha)
    }
    # Parse flexible-shape distributions
    incubation_time <- eval(parse(text=incubation_time))
    test_time <- eval(parse(text=test_time))
    trace_time_auto <- eval(parse(text=trace_time_auto))
    trace_time_manual <- eval(parse(text=trace_time_manual))
    delay_time <- eval(parse(text=delay_time))
    # Prepare case creation functions
    index_case_fn <- purrr::partial(create_index_cases, n_initial_cases = n_initial_cases,
                                    p_asymptomatic = p_asymptomatic,
                                    test_time = test_time, test_sensitivity = test_sensitivity,
                                    test_serological = test_serological, p_ident_sym = p_ident_sym,
                                    p_smartphone_overall = p_smartphone_overall,
                                    n_children_fn = n_children_fn,
                                    trace_neg_symptomatic = trace_neg_symptomatic,
                                    incubation_time = incubation_time,
                                    p_traced_auto = p_traced_auto,
                                    p_traced_manual = p_traced_manual,
                                    trace_time_auto = trace_time_auto,
                                    trace_time_manual = trace_time_manual,
                                    recovery_time = recovery_time,
                                    data_limit_auto = data_limit_auto,
                                    data_limit_manual = data_limit_manual,
                                    contact_limit_auto = contact_limit_auto,
                                    contact_limit_manual = contact_limit_manual,
                                    rollout_delay_gen = rollout_delay_gen,
                                    rollout_delay_days = rollout_delay_days,
                                    delay_time = delay_time,
                                    p_data_sharing_auto = p_data_sharing_auto,
                                    p_data_sharing_manual = p_data_sharing_manual,
                                    p_compliance_isolation = p_compliance_isolation
                                    )
    child_case_fn <- purrr::partial(create_child_cases,
                                    p_asymptomatic = p_asymptomatic,
                                    p_compliance_isolation = p_compliance_isolation,
                                    test_time = test_time, test_sensitivity = test_sensitivity,
                                    test_serological = test_serological, p_ident_sym = p_ident_sym,
                                    p_smartphone_infector_yes = p_smartphone_infector_yes,
                                    p_smartphone_infector_no = p_smartphone_infector_no,
                                    generation_time = generation_time,
                                    n_children_fn = n_children_fn,
                                    trace_neg_symptomatic = trace_neg_symptomatic,
                                    incubation_time = incubation_time,
                                    p_traced_auto = p_traced_auto,
                                    p_traced_manual = p_traced_manual,
                                    trace_time_auto = trace_time_auto,
                                    trace_time_manual = trace_time_manual,
                                    recovery_time = recovery_time,
                                    data_limit_auto = data_limit_auto,
                                    data_limit_manual = data_limit_manual,
                                    contact_limit_auto = contact_limit_auto,
                                    contact_limit_manual = contact_limit_manual,
                                    rollout_delay_gen = rollout_delay_gen,
                                    rollout_delay_days = rollout_delay_days,
                                    delay_time = delay_time,
                                    p_environmental = p_environmental,
                                    p_data_sharing_auto = p_data_sharing_auto,
                                    p_data_sharing_manual = p_data_sharing_manual
                                    )
    # Execute runs
    iter_out <- purrr::map(.x = 1:n_iterations,
                           ~ run_outbreak(index_case_fn = index_case_fn,
                                          child_case_fn = child_case_fn,
                                          cap_max_generations = cap_max_generations,
                                          cap_max_weeks = cap_max_weeks,
                                          cap_cases = cap_cases,
                                          backtrace_distance = backtrace_distance))
    gc(verbose=FALSE, full=TRUE)
    # Label and concatenate
    results_raw <- lapply(1:n_iterations, function(n) iter_out[[n]][, run := n]) %>%
        data.table::rbindlist(fill = TRUE) %>% .[, scenario := scenario]
    gc(verbose=FALSE, full=TRUE)
    if (report){
        cat("Scenario ", scenario, " concluding at: ", date(), " (", timetaken(start), ")\n", sep="")
        cat("Scenario ", scenario, " final case count across all runs: ", nrow(results_raw), "\n", sep="")
    }
    return(results_raw)
}
