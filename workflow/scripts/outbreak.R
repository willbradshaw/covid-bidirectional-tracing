# Core functionality for branching-chain outbreak model
# Author: William J. Bradshaw (heavily modified from code by Joel Hellewell
#                              and Sam Abbott)
# Date: 5 May 2020

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
                                         p_data_sharing_auto, p_data_sharing_manual,
                                         p_smartphone_listens, p_smartphone_chirps
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
               smartphone_listens = purrr::rbernoulli(n_initial_cases, p_smartphone_listens),
               smartphone_chirps  = purrr::rbernoulli(n_initial_cases, p_smartphone_chirps),
               generation = 0, infector_id = 0, infector_onset_gen = Inf,
               infector_exposure = Inf,
               infector_onset_true = Inf, infector_asym = NA, infector_recovery = Inf,
               infector_has_smartphone = NA, infector_shares_data_auto = FALSE,
               infector_shares_data_manual = FALSE, infector_isolates = FALSE,
               infector_smartphone_listens = NA, infector_smartphone_chirps = NA)
}

set_primary_immutables_child <- function(parents, p_asymptomatic, p_compliance_isolation,
                                         test_time, test_serological,
                                         test_sensitivity, p_ident_sym, p_environmental,
                                         p_data_sharing_auto, p_data_sharing_manual,
                                         p_smartphone_listens, p_smartphone_chirps
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
                              smartphone_listens = purrr::rbernoulli(n_children_total, p_smartphone_listens),
                              smartphone_chirps  = purrr::rbernoulli(n_children_total, p_smartphone_chirps),
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
                              infector_isolates = rep_new_cases(parents, "isolates"),
                              infector_smartphone_listens = rep_new_cases(parents, "smartphone_listens"), # Does parent's smartphone (if any) listen for chirps?
                              infector_smartphone_chirps = rep_new_cases(parents, "smartphone_chirps"), # Does parent's smartphone (if any) broadcast chirps?
                              )
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
        .[, `:=`(smartphone_match_fwd_1way = infector_has_smartphone & has_smartphone & 
                    infector_smartphone_chirps & smartphone_listens, # Compatible phones? (1-way)
                smartphone_match_fwd_2way = ifelse(rep(!reciprocal_chirping, nrow(.)), FALSE,
                    infector_has_smartphone & has_smartphone &
                    smartphone_chirps & infector_smartphone_listens), # Compatible phones? (2-way)
                 smartphone_match_rev_1way = infector_has_smartphone & has_smartphone & 
                    smartphone_chirps & infector_smartphone_listens, # Compatible phones? (1-way)
                smartphone_match_rev_2way = ifelse(rep(!reciprocal_chirping, nrow(.)), FALSE,
                    infector_has_smartphone & has_smartphone &
                    infector_smartphone_chirps & smartphone_listens), # Compatible phones? (2-way)
                smartphone_match_fwd = smartphone_match_fwd_1way | smartphone_match_fwd_2way,
                smartphone_match_rev = smartphone_match_rev_1way | smartphone_match_rev_2way)
        ] %>%
        .[, `:=`(auto_traced_fwd = rep(!index, nrow(.)) & # Not an index case
                     smartphone_match_fwd & # Compatible smartphones
                     infector_shares_data_auto & # Infector shares data
                     !environmental & # Contact not environmental
                     purrr::rbernoulli(nrow(.), p_traced_auto), # Contact is recorded
                 auto_traced_rev = rep(!index, nrow(.)) & # Not an index case
                     smartphone_match_rev & # Compatible smartphones
                     shares_data_auto & # This case shares data
                     !environmental & # Contact not environmental
                     purrr::rbernoulli(nrow(.), p_traced_auto), # Contact is recorded
                 onset_gen = exposure + incubation_time(nrow(.)))] %>%
        .[, `:=`(manual_traced_fwd = (rep(!index, nrow(.)) & # Not an index case
                     !auto_traced_fwd & # Not digitally traced
                     !environmental & # Contact non environmental
                     infector_shares_data_manual & # Infector shares data
                     purrr::rbernoulli(nrow(.), p_traced_manual)), # Contact is recorded
                 manual_traced_rev = (rep(!index, nrow(.)) & # Not an index case
                     !auto_traced_rev & # Not digitally traced
                     !environmental & # Contact not environmental
                     shares_data_manual & # This case shares data
                     purrr::rbernoulli(nrow(.), p_traced_manual)))] %>% # Contact is recorded
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
                                                                                  infector_quar_time))
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
                                                                                 quarantine_time)
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
                               p_smartphone_listens, p_smartphone_chirps,
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
                                          p_data_sharing_manual = p_data_sharing_manual,
                                          p_smartphone_listens = p_smartphone_listens,
                                          p_smartphone_chirps = p_smartphone_chirps)
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
                               p_smartphone_listens, p_smartphone_chirps,
                               rollout_delay_days, delay_time, p_environmental){
    #' Generate a table of new child cases from a table of parent cases
    # Set immutable keys
    cases <- set_primary_immutables_child(parents = parents, p_asymptomatic = p_asymptomatic,
                                          p_compliance_isolation = p_compliance_isolation,
                                          test_time = test_time, test_serological = test_serological,
                                          test_sensitivity = test_sensitivity,
                                          p_ident_sym = p_ident_sym,
                                          p_environmental = p_environmental,
                                          p_data_sharing_auto = p_data_sharing_auto,
                                          p_data_sharing_manual = p_data_sharing_manual,
                                          p_smartphone_listens = p_smartphone_listens,
                                          p_smartphone_chirps = p_smartphone_chirps)
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
                                                                                        infector_quar_time))]
    
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
        .[, .(cases = .N, extinct = all(processed)), by = c("generation", group_vars)] %>%
        .[order(generation), cumulative_cases := cumsum(cases), by=group_vars] %>%
        .[order(generation), effective_r0 := ifelse(extinct & generation == max(generation),
                                                    0, lead(cases)/cases), by=group_vars]
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

summarise_by_scenario <- function(run_data, group_vars, ci_width = 0.95,
                                  alpha_prior = 1, beta_prior = 1){
    #' Convert a database of case reports
    group_vars <- c("scenario", group_vars[! group_vars %in% c("scenario", "run")])
    beta_bound <- function(q,x,y) qbeta(q, alpha_prior+x, beta_prior+y)
    beta_lower <- function(x,y) beta_bound((1-ci_width)/2, x, y)
    beta_upper <- function(x,y) beta_bound(1-(1-ci_width)/2, x, y)
    scenario_stats <- run_data %>% copy %>%
        .[, .(n_runs = length(run), p_extinct = mean(extinct),
              effective_r0_mean = mean(avg_effective_r0),
              effective_r0_sd = sd(avg_effective_r0),
              avg_total_cases_total = mean(total_cases),
              avg_total_cases_controlled = sum(total_cases * controlled)/sum(controlled),
              avg_total_cases_uncontrolled = sum(total_cases * !controlled)/sum(!controlled),
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

write_dfile <- function(data, path_prefix, compress = TRUE){
    if (!compress) f <- file else f <- gzfile
    if (!compress) ext <- ".tsv" else ext <- ".tsv.gz"
    dpath <- paste0(path_prefix, ext)
    dfile <- f(dpath, "w")
    write_tsv(data, dfile)
    close(dfile)
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
                         write_raw = NULL, write_weekly = NULL,
                         write_generational = NULL, path_prefix = NULL,
                         compress = NULL, scenario_data = NULL){
    #' Run a single complete instance of the branching-process model
    # Set up initial cases
    case_data <- index_case_fn()
    # Run outbreak loop
    while (max(case_data$exposure)/7 < scenario_data$cap_max_weeks & # Time limit
           nrow(case_data) < scenario_data$cap_cases & # Cumulative case limit
           max(case_data$generation) < scenario_data$cap_max_generations & # Generation limit
           any(!case_data$processed)){ # Extinction condition
        case_data <- outbreak_step(case_data = case_data,
                                   child_case_fn = child_case_fn,
                                   backtrace_distance = scenario_data$backtrace_distance)
        #gc(verbose=FALSE, full=TRUE)
    }
    # Summarise and write
    case_data[, `:=`(scenario=scenario_data$scenario, run=scenario_data$run)]
    cases_out <- as.data.table(scenario_data)[case_data, on=c("scenario", "run")]
    path_prefix <- paste(path_prefix, scenario_data$scenario, scenario_data$run,
                         sep="_")
    if (write_raw){
        write_dfile(cases_out, paste0(path_prefix, "_raw"), compress)
    }
    if (write_weekly){
        week_data <- summarise_by_week(cases_out, colnames(scenario_data))
        write_dfile(week_data, paste0(path_prefix, "_weekly"), compress)
    }
    if (write_generational){
        gen_data <- summarise_by_generation(cases_out, colnames(scenario_data))
        write_dfile(gen_data, paste0(path_prefix, "_gen"), compress)
    }
    run_data <- summarise_by_run(cases_out, colnames(scenario_data))
    gc(verbose=FALSE, full=TRUE)
    return(run_data)
}

scenario_sim <- function(scenario = NULL, n_iterations = NULL, report = NULL,
                         write_raw = NULL, write_weekly = NULL,
                         write_generational = NULL, write_run = NULL,
                         path_prefix = NULL, compress = NULL,
                         ci_width = NULL,
                         alpha_prior = NULL, beta_prior = NULL){
    #' Run a specified number of outbreaks with identical parameters
    if (report){
        start <- proc.time()
        cat("Scenario ", scenario$scenario, " beginning at: ", date(), "\n", sep="")
    }
    # Compute auxiliary parameters
    r0_asymptomatic <- scenario$rel_r0_asymptomatic * scenario$r0_base
    r0_symptomatic <- compute_symptomatic_r0(scenario$r0_base, r0_asymptomatic,
                                             scenario$p_asymptomatic)
    p_smartphone_infector_yes <- ifelse(scenario$p_smartphone_link < 0,
                                        scenario$p_smartphone_overall,
                                        scenario$p_smartphone_link)
    p_smartphone_infector_no <- compute_p_smartphone_infector_no(scenario$p_smartphone_overall,
                                                             p_smartphone_infector_yes)
    p_traced_manual <- ifelse(scenario$p_traced_manual < 0,
                              scenario$p_traced_auto,
                              scenario$p_traced_manual)
    p_traced_auto <- ifelse(scenario$p_traced_auto < 0,
                              scenario$p_traced_manual,
                              scenario$p_traced_auto)
    # Prepare fixed-shape distributions
    n_children_fn <- function(asym) rnbinom(n=length(asym), size=scenario$dispersion,
                                            mu=ifelse(asym, r0_asymptomatic,
                                                      r0_symptomatic))
    generation_time <- function(onsets){
        sn::rsn(n = length(onsets), xi = onsets, omega = scenario$generation_omega,
                alpha = scenario$generation_alpha) %>%  ifelse(. < 0, 0, .)
    } # "generation_alpha" previously known as "k" or "generation_k"
    recovery_time <- function(onsets){
        sn::qsn(scenario$recovery_quantile, xi = onsets, omega = scenario$generation_omega,
                alpha = scenario$generation_alpha)
    }
    # Parse flexible-shape distributions
    incubation_time <- eval(parse(text=scenario$incubation_time))
    test_time <- eval(parse(text=scenario$test_time))
    trace_time_auto <- eval(parse(text=scenario$trace_time_auto))
    trace_time_manual <- eval(parse(text=scenario$trace_time_manual))
    delay_time <- eval(parse(text=scenario$delay_time))
    # Prepare case creation functions
    index_case_fn <- purrr::partial(create_index_cases, n_initial_cases = scenario$n_initial_cases,
                                    p_asymptomatic = scenario$p_asymptomatic,
                                    test_time = test_time, test_sensitivity = scenario$test_sensitivity,
                                    test_serological = scenario$test_serological,
                                    p_ident_sym = scenario$p_ident_sym,
                                    p_smartphone_overall = scenario$p_smartphone_overall,
                                    n_children_fn = n_children_fn,
                                    trace_neg_symptomatic = scenario$trace_neg_symptomatic,
                                    incubation_time = incubation_time,
                                    p_traced_auto = p_traced_auto,
                                    p_traced_manual = p_traced_manual,
                                    trace_time_auto = trace_time_auto,
                                    trace_time_manual = trace_time_manual,
                                    recovery_time = recovery_time,
                                    data_limit_auto = scenario$data_limit_auto,
                                    data_limit_manual = scenario$data_limit_manual,
                                    contact_limit_auto = scenario$contact_limit_auto,
                                    contact_limit_manual = scenario$contact_limit_manual,
                                    rollout_delay_gen = scenario$rollout_delay_gen,
                                    rollout_delay_days = scenario$rollout_delay_days,
                                    delay_time = delay_time,
                                    p_data_sharing_auto = scenario$p_data_sharing_auto,
                                    p_data_sharing_manual = scenario$p_data_sharing_manual,
                                    p_compliance_isolation = scenario$p_compliance_isolation,
                                    p_smartphone_listens = scenario$p_smartphone_listens,
                                    p_smartphone_chirps = scenario$p_smartphone_chirps
                                    )
    child_case_fn <- purrr::partial(create_child_cases,
                                    p_asymptomatic = scenario$p_asymptomatic,
                                    p_compliance_isolation = scenario$p_compliance_isolation,
                                    test_time = test_time, test_sensitivity = scenario$test_sensitivity,
                                    test_serological = scenario$test_serological,
                                    p_ident_sym = scenario$p_ident_sym,
                                    p_smartphone_infector_yes = p_smartphone_infector_yes,
                                    p_smartphone_infector_no = p_smartphone_infector_no,
                                    generation_time = generation_time,
                                    n_children_fn = n_children_fn,
                                    trace_neg_symptomatic = scenario$trace_neg_symptomatic,
                                    incubation_time = incubation_time,
                                    p_traced_auto = p_traced_auto,
                                    p_traced_manual = p_traced_manual,
                                    trace_time_auto = trace_time_auto,
                                    trace_time_manual = trace_time_manual,
                                    recovery_time = recovery_time,
                                    data_limit_auto = scenario$data_limit_auto,
                                    data_limit_manual = scenario$data_limit_manual,
                                    contact_limit_auto = scenario$contact_limit_auto,
                                    contact_limit_manual = scenario$contact_limit_manual,
                                    rollout_delay_gen = scenario$rollout_delay_gen,
                                    rollout_delay_days = scenario$rollout_delay_days,
                                    delay_time = delay_time,
                                    p_environmental = scenario$p_environmental,
                                    p_data_sharing_auto = scenario$p_data_sharing_auto,
                                    p_data_sharing_manual = scenario$p_data_sharing_manual,
                                    p_smartphone_listens = scenario$p_smartphone_listens,
                                    p_smartphone_chirps = scenario$p_smartphone_chirps
                                    )
    # Execute runs
    iter_out <- purrr::map(.x = 1:n_iterations,
                           ~ run_outbreak(index_case_fn = index_case_fn,
                                          child_case_fn = child_case_fn,
                                          write_raw = write_raw, write_weekly = write_weekly,
                                          write_generational = write_generational,
                                          path_prefix = path_prefix, 
                                          compress = compress,
                                          scenario_data = scenario %>% mutate(run=.x)))
    # Concatenate and write
    run_data <- iter_out %>% rbindlist
    if (write_run){
        path_prefix <- paste(path_prefix, scenario$scenario, sep="_")
        write_dfile(run_data, paste0(path_prefix, "_run"), compress)
    }
    scenario_data <- summarise_by_scenario(run_data, colnames(scenario), ci_width,
                                           alpha_prior, beta_prior)
    gc(verbose=FALSE, full=TRUE)
    if (report){
        cat("Scenario ", scenario$scenario, " concluding at: ", date(), " (", timetaken(start), ")\n", sep="")
    }
    return(scenario_data)
}
