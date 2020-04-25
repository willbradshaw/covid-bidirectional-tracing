library(tidyverse)
library(data.table)
library(sn)

#==============================================================================
# Setup
#==============================================================================

# Set random seed
set.seed(09122)

# Set background parameters
n_iterations = 100;
dispersion = 0.16;
r0_base = 2.5;
r0_asymptomatic = r0_base;
generation_omega = 2; 
generation_alpha = 0.7;
recovery_quantile = 0.999;
incubation_time = "function(n) rweibull(n=n, scale=6.492272, shape=2.322737)"
test_time = "function(n) rep(0, n)";
trace_time_auto = "function(n) rep(0, n)";
trace_time_manual = "function(n) rep(0, n)";
delay_time = "function(n) rweibull(n=n, scale = 4.287786, shape = 1.651524)";
n_initial_cases = 20;
test_sensitivity = 1;
test_serological = TRUE;
p_smartphone_overall = 1;
p_smartphone_link = 1;
trace_neg_symptomatic = TRUE;
p_traced_manual = 0;
data_limit_auto = Inf;
data_limit_manual = Inf;
contact_limit_manual = Inf;
p_blocked_isolation = 1;
p_blocked_quarantine = 1;
cap_max_generations = Inf;
cap_max_weeks = 52;
cap_cases = 5000;
p_environmental = 0;

# Set parameters of interest
p_traced_auto = 0.8;
rollout_delay_gen = Inf;
rollout_delay_days = Inf;
backtrace_distance = Inf;
p_asymptomatic = 0.3;
p_ident_sym = 0;
n_gen_run <- 8;
contact_limit_auto = Inf;


# Define auxiliary functions (from scenario_sim)
r0_symptomatic <- compute_symptomatic_r0(r0_base, r0_asymptomatic, p_asymptomatic)
p_smartphone_infector_yes <- ifelse(p_smartphone_link == Inf, p_smartphone_overall,
                                    p_smartphone_link)
p_smartphone_infector_no <- compute_p_smartphone_infector_no(p_smartphone_overall,
                                                             p_smartphone_infector_yes)
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
incubation_time <- eval(parse(text=incubation_time))
test_time <- eval(parse(text=test_time))
trace_time_auto <- eval(parse(text=trace_time_auto))
trace_time_manual <- eval(parse(text=trace_time_manual))
delay_time <- eval(parse(text=delay_time))
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
                                delay_time = delay_time)
child_case_fn <- purrr::partial(create_child_cases,
                                p_asymptomatic = p_asymptomatic,
                                p_blocked_isolation = p_blocked_isolation,
                                p_blocked_quarantine = p_blocked_quarantine,
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
                                p_environmental = p_environmental)

#==============================================================================
# Null run (no isolation or tracing)
#==============================================================================

# Initialise new cases
new_cases_null <- index_case_fn()
case_data_null <- new_cases_null %>% copy

# Run loops
for (n in 1:n_gen_run){
  case_data_null <- outbreak_step(case_data = case_data_null,
                             child_case_fn = child_case_fn,
                             backtrace_distance = backtrace_distance)
}

#==============================================================================
# Analysis
#==============================================================================

#------------------------------------------------------------------------------
# 1. Distribution of recovery times
#------------------------------------------------------------------------------

# Mean ± SD
recovery_tab <- case_data_null %>% group_by(generation) %>%
  mutate(recovery_weeks = recovery/7) %>%
  summarise(recovery_weeks_mean = mean(recovery_weeks),
            recovery_weeks_sd = sd(recovery_weeks))
g_recovery <- ggplot(recovery_tab) +
  geom_ribbon(aes(x=generation, ymin=recovery_weeks_mean-recovery_weeks_sd,
                  ymax = recovery_weeks_mean + recovery_weeks_sd),
              alpha = 0.3) +
  geom_line(aes(x=generation, y=recovery_weeks_mean)) +
  geom_abline(intercept=1, colour = "red", linetype = "dashed") +
  scale_y_continuous(name="Recovery time (mean ± s.d.)", limits=c(0,NA),
                     breaks = seq(0,n_gen_run*2, 2)) +
  scale_x_continuous(name="Generation") +
  theme_bw()

# Histogram
g_recovery_hist <- ggplot(case_data_null) +
  geom_histogram(aes(x=recovery), binwidth=1) +
  geom_vline(xintercept = 28, colour="red") +
  scale_x_continuous(name = "Recovery time (weeks)", 
                     breaks = seq(0, max(case_data_null$recovery), 7),
                     labels = function(x) round(x/7)) +
  facet_grid(generation~., labeller = labeller(generation = function(x) paste("gen.", x)),
             scales = "free_y") +
  theme_bw()

#------------------------------------------------------------------------------
# 2. Number/proportion of hidden cases
#------------------------------------------------------------------------------

# Number of individuals blocked by response delay
delay_tab <- tibble(delay = 1:n_gen_run-1) %>%
  mutate(generations = sapply(delay, function(d) sum(case_data_null$generation < d)),
         weeks = sapply(delay, function(d) sum(case_data_null$recovery < d * 7))) %>%
  gather(delay_type, n_blocked, -delay) %>% 
  mutate(p_blocked = n_blocked / nrow(case_data_null))

g_n_hidden <- ggplot(delay_tab, aes(x=delay, y=n_blocked, colour=delay_type)) +
  geom_line() + geom_point(size=2) +
  scale_x_continuous(name="Delay", breaks = seq(0,n_gen_run,2)) +
  scale_y_log10(name="# cases hidden by delay") +
  theme_bw() + ggtitle(paste("R0 =", r0_base))

# Proportion of individuals in each generation hidden by response delay
p_hidden_tab_gen <- expand_grid(delay = 1:n_gen_run-1,
                                generation=0:max(case_data_null$generation)) %>%
  mutate(p_hidden = (generation < delay) * 1, delay_type = "generations")
p_hidden_tab_weeks <- expand_grid(delay = 1:n_gen_run-1,
                                  generation=0:max(case_data_null$generation)) %>%
  inner_join(case_data_null, by="generation") %>% group_by(delay, generation) %>%
  summarise(p_hidden = mean(recovery < delay*7), delay_type = "weeks")
p_hidden_tab <- bind_rows(p_hidden_tab_gen, p_hidden_tab_weeks)

g_p_hidden <- ggplot(p_hidden_tab_weeks, aes(x=delay, y=p_hidden, 
                                             colour=factor(generation))) +
  geom_line() + geom_point(size=2) +
  scale_x_continuous(name="Delay (weeks)", breaks = seq(0,n_gen_run,2)) +
  scale_y_continuous(name="Fraction of cases hidden by delay", breaks = seq(0,1,0.2)) +
  scale_color_discrete(name="Generation") +
  theme_bw()

#------------------------------------------------------------------------------
# 3. Time from hypothetical trace time to exposure (for data limits)
#------------------------------------------------------------------------------

# Get database of trace windows for different trace types
exposure_to_trace_fwd <- case_data_null[traceable_fwd==TRUE] %>% copy %>%
  .[, trace_window := infector_onset_true + delay_time(nrow(.)) - exposure] %>%
  .[, trace_type := "forward tracing"]
exposure_to_trace_rev_1gen <- case_data_null[traceable_rev==TRUE] %>% copy %>%
  .[, trace_window := onset_true + delay_time(nrow(.)) - exposure] %>%
  .[, trace_type := "reverse tracing (1 gen.)"]
exposure_to_trace_rev_2gen <- case_data_null[traceable_rev==TRUE] %>% copy %>%
  .[, trace_window := onset_true + delay_time(nrow(.)) - infector_exposure] %>%
  .[, trace_type := "reverse tracing (2 gen.)"]
trace_windows <- rbindlist(list(exposure_to_trace_fwd, exposure_to_trace_rev_1gen,
                                exposure_to_trace_rev_2gen))
trace_windows_finite <- trace_windows[trace_window < Inf]

g_trace_window <- ggplot(trace_windows_finite) +
  geom_histogram(aes(x=trace_window), binwidth=1) +
  geom_vline(xintercept = c(2,7,14), linetype="dashed", colour="red") +
  facet_grid(trace_type~.) +
  scale_x_continuous(name="Time from exposure to trace initiation",
                     breaks = c(2,seq(0,max(trace_windows_finite$trace_window), 7))) +
  theme_bw()

# As before, but taking minimum successful reverse trace for each parent
exposure_to_trace_rev_1gen_min <- case_data_null[traceable_rev==TRUE] %>% copy %>%
  .[, trace_init_time := onset_true + delay_time(nrow(.))] %>%
  .[, .SD[which.min(trace_init_time)], by = "infector_id"] %>%
  .[, `:=`(trace_window = trace_init_time - exposure,
           trace_type = "reverse tracing (1 gen.)")]
exposure_to_trace_rev_2gen_min <- case_data_null[traceable_rev==TRUE] %>% copy %>%
  .[, trace_init_time := onset_true + delay_time(nrow(.))] %>%
  .[, .SD[which.min(trace_init_time)], by = "infector_id"] %>%
  .[, `:=`(trace_window = trace_init_time - infector_exposure,
           trace_type = "reverse tracing (2 gen.)")]
trace_windows_min <- rbindlist(list(exposure_to_trace_fwd,
                                    exposure_to_trace_rev_1gen_min,
                                    exposure_to_trace_rev_2gen_min), fill=TRUE)
trace_windows_min_finite <- trace_windows_min[trace_window < Inf]

g_trace_window_min <- ggplot(trace_windows_min_finite) +
  geom_histogram(aes(x=trace_window), binwidth=1, boundary=7) +
  geom_vline(xintercept = c(2,7,14), linetype="dashed", colour="red") +
  facet_grid(trace_type~., scales="free_y") +
  scale_x_continuous(name="Time from exposure to trace initiation",
                     breaks = c(2,seq(0,max(trace_windows_finite$trace_window), 7))) +
  theme_bw()


  
#------------------------------------------------------------------------------
# 3. Chain length
#------------------------------------------------------------------------------
# 
# # For an N-generation delay, all identified individuals have a chain of
# # at least N hidden ancestors that could be uncovered by multigenerational
# # (but not 1-generational) backtracing
# 
# list_descendents <- function(cases, ids, include_orig = TRUE){
#   # Get a list of all descendents of a set of cases
#   ids_out <- numeric(0)
#   for (m in 1:n_gen_run){
#     new_ids <- cases[infector_id %in% c(ids, ids_out)]$case_id
#     ids_out <- c(ids_out, new_ids)
#   }
#   if (include_orig) ids_out <- c(ids, ids_out)
#   return(unique(ids_out))
# }
# 
# purge_descendents <- function(cases, ids){
#   # Remove all cases descended from a list of case IDs
#   for (m in 1:n_gen_run){
#     new_ids <- cases[infector_id %in% ids]$case_id
#     ids <- unique(c(ids, new_ids))
#   }
#   return(cases[!case_id %in% ids])
# }
# 
# gen_filter <- function(n){
#   # Goal: find all cases that would be the *first* identified member in their
#   # lineage given an n-generation response delay
#   # 1. Filter out any cases hidden by response delay
#   cases <- case_data_null[generation >= n]
#   # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
#   #    and discard their descendents
#   harbingers_sym <- list(); length(harbingers_sym) <- n_gen_run
#   for (m in 1:(n_gen_run - n + 1)){
#     # a. Find and separate out all symptomatic individuals with no parents
#     harbingers_sym[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
#     # b. Discard their descendents
#     cases <- purge_descendents(cases, harbingers_sym[[m]]$case_id)
#     # c. Discard any asymptomatic individuals with no parents, then repeat cycle
#     cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
#   }
#   # 3. Combine into one dataframe
#   cases <- rbindlist(harbingers_sym)
#   # 4. Count generations and return
#   cases_out <- cases[, .(n_cases = .N), by="generation"][
#     , `:=`(p_cases = n_cases/sum(n_cases), delay = n, delay_type = "generations")
#     ]
#   return(cases_out)
# }
# week_filter <- function(n){
#   # Goal: find all cases that would be the *first* identified member in their
#   # lineage given an n-week response delay
#   # 1. Filter out any cases hidden by response delay
#   cases <- case_data_null[recovery >= n*7]
#   # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
#   #    and discard their descendents
#   harbingers_sym <- list(); length(harbingers_sym) <- n_gen_run
#   for (m in 1:n_gen_run){
#     # a. Find and separate out all symptomatic individuals with no parents
#     harbingers_sym[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
#     # b. Discard their descendents
#     cases <- purge_descendents(cases, harbingers_sym[[m]]$case_id)
#     # c. Discard any asymptomatic individuals with no parents, then repeat cycle
#     cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
#   }
#   # 3. Combine into one dataframe
#   cases <- rbindlist(harbingers_sym)
#   # 4. Count generations and return
#   cases_out <- cases[, .(n_cases = .N), by="generation"][
#     , `:=`(p_cases = n_cases/sum(n_cases), delay = n, delay_type = "weeks")
#     ]
#   return(cases_out)
# }
# 
# gen_delays <- lapply(1:n_gen_run-1, function(n) gen_filter(n)) %>% rbindlist
# week_delays <- lapply(1:n_gen_run-1, function(n) week_filter(n)) %>% rbindlist
# filters <- rbind(gen_delays, week_delays)
# 
# g_filters <- ggplot(filters, aes(x=generation, y=p_cases, colour=delay_type)) +
#   geom_line() + geom_point(size=2) +
#   facet_grid(delay ~ ., scales="free_y") +
#   scale_y_continuous(limits=c(0,1)) +
#   theme_bw()
# 
# # ...
# # What % of cases are accessible to 1-generational backtracing after an n-week delay?
# # What % of cases are descended from cases that
# # - (i) are non-hidden
# # - (ii) are hidden, but whose direct descendents are non-hidden
# 
# p_cases_accessible_weeks <- function(delay){
#   # 1. Filter out any cases hidden by response delay
#   cases <- case_data_null[recovery >= delay*7]
#   # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
#   #    and discard their descendents
#   harbingers <- list(); length(harbingers) <- n_gen_run
#   for (m in 1:n_gen_run){
#     # a. Find and separate out all symptomatic individuals with no parents
#     harbingers[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
#     # b. Discard their descendents
#     cases <- purge_descendents(cases, harbingers[[m]]$case_id)
#     # c. Discard any asymptomatic individuals with no parents, then repeat cycle
#     cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
#   }
#   # 3. Combine into one dataframe
#   harbingers_all <- rbindlist(harbingers)
#   # 4. Count proportion of cases that are descended from all harbingers (and their parents)
#   desc_harb <- list_descendents(case_data_null, harbingers_all$case_id, TRUE)
#   desc_harb_inf <- c(desc_harb, 
#                      list_descendents(case_data_null, harbingers_all$infector_id, TRUE)) %>%
#     unique %>% .[. != 0]
#   t_desc <- tibble(delay = delay, delay_type = "weeks",
#                    n_desc_harb = length(desc_harb),
#                    n_desc_harb_inf = length(desc_harb_inf))
#   return(t_desc)
# }
# 
# p_cases_accessible_gen <- function(delay){
#   # 1. Filter out any cases hidden by response delay
#   cases <- case_data_null[generation >= delay]
#   # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
#   #    and discard their descendents
#   harbingers <- list(); length(harbingers) <- n_gen_run
#   for (m in 1:n_gen_run){
#     # a. Find and separate out all symptomatic individuals with no parents
#     harbingers[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
#     # b. Discard their descendents
#     cases <- purge_descendents(cases, harbingers[[m]]$case_id)
#     # c. Discard any asymptomatic individuals with no parents, then repeat cycle
#     cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
#   }
#   # 3. Combine into one dataframe
#   harbingers_all <- rbindlist(harbingers)
#   # 4. Count proportion of cases that are descended from all harbingers (and their parents)
#   desc_harb <- list_descendents(case_data_null, harbingers_all$case_id, TRUE)
#   desc_harb_inf <- c(desc_harb, 
#                      list_descendents(case_data_null, harbingers_all$infector_id, TRUE)) %>%
#     unique %>% .[. != 0]
#   t_desc <- tibble(delay = delay, delay_type = "generations",
#                    n_desc_harb = length(desc_harb),
#                    n_desc_harb_inf = length(desc_harb_inf))
#   return(t_desc)
# }
# 
# p_acc_gen <- lapply(1:n_gen_run-1, function(n) p_cases_accessible_gen(n)) %>% rbindlist
# p_acc_weeks <- lapply(1:n_gen_run-1, function(n) p_cases_accessible_weeks(n)) %>% rbindlist
# p_acc <- rbind(p_acc_gen, p_acc_weeks)
# 
# g_acc <- p_acc %>% gather(backtrace_distance, n_hit, -delay, -delay_type) %>%
#   mutate(backtrace_distance = ifelse(backtrace_distance == "n_desc_harb_inf", 1, 0),
#          n_missed = nrow(case_data_null) - n_hit, 
#          p_missed = n_missed / nrow(case_data_null)) %>% 
#   ggplot(aes(x=delay, y=n_missed, colour=factor(backtrace_distance))) + 
#   geom_line(aes(linetype = delay_type)) + geom_point(size=2) +
#   theme_bw()
# # There might be something salvageable here, but we're not there yet