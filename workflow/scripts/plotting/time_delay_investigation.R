library(tidyverse)
library(data.table)
library(sn)

#==============================================================================
# Setup
#==============================================================================

# Set random seed
#set.seed(09122)
set.seed(132421)
#set.seed(555555)

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

print(nrow(case_data_null))

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

get_trace_windows <- function(cases, backtrace_distance, correct_min){
  # Filter by traceability (depends on p_trace)
  if (backtrace_distance > 0){
    cases_filtered <- cases[traceable_rev == TRUE] %>% copy
  } else {
    cases_filtered <- cases[traceable_fwd == TRUE] %>% copy
  }
  # Determine trace initiation times
  if (backtrace_distance == 0){
    cases_init <- cases_filtered %>% copy %>%
      .[, infector_trace_init_time := infector_onset_true + delay_time(nrow(.))]
  } else {
    cases_init <- cases_filtered %>% copy %>%
      .[, trace_init_time := onset_true + delay_time(nrow(.))]
  }
  # If taking minimum across traces, do that now
  if (correct_min & backtrace_distance > 0){
    cases_min <- cases_init %>% copy %>%
      .[, .SD[which.min(trace_init_time)], by = "infector_id"]
  } else {
    cases_min <- cases_init %>% copy
  }

  # Compare to relevant exposure and determine trace windows
  if (backtrace_distance == 0){
    cases_window <- cases_min %>% copy %>%
      .[, `:=`(trace_window = infector_trace_init_time - exposure,
               trace_type = "forward tracing",
               corrected = correct_min)]
  } else if (backtrace_distance == 1){
    cases_window <- cases_min %>% copy %>%
      .[, `:=`(trace_window = trace_init_time - exposure,
               trace_type = paste0("reverse tracing (1 gen.)"),
               corrected = correct_min)]
  } else {
    cases_window <- cases_min %>% copy %>%
      .[, `:=`(trace_window = trace_init_time - infector_exposure,
               trace_type = paste0("reverse tracing (", backtrace_distance, " gen.)"),
               corrected = correct_min)]
  }
  return(cases_window %>% select(trace_window, trace_type, corrected))
}

get_trace_windows_ptrace <- function(p_trace, correct_min){
  cases_ptrace <- case_data_null[generation > 0] %>% copy %>%
    .[, `:=`(traceable_fwd = purrr::rbernoulli(nrow(.), p_trace),
             traceable_rev = purrr::rbernoulli(nrow(.), p_trace))]
  windows_fwd <- get_trace_windows(cases_ptrace, 0, correct_min)
  windows_rev1 <- get_trace_windows(cases_ptrace, 1, correct_min)
  windows_rev2 <- get_trace_windows(cases_ptrace, 2, correct_min)
  windows_out <- rbind(windows_fwd, windows_rev1, windows_rev2) %>%
    mutate(p_traced = p_trace)
  return(windows_out)
}

trace_windows_uncorrected <- lapply(seq(0.2, 1, 0.2), function(p)
  get_trace_windows_ptrace(p, FALSE)) %>% rbindlist
trace_windows_corrected <- lapply(seq(0.2, 1, 0.2), function(p)
  get_trace_windows_ptrace(p, TRUE)) %>% rbindlist
trace_windows_all <- rbind(trace_windows_uncorrected,
                           trace_windows_corrected)
x_max <- trace_windows_all %>% filter(trace_window < Inf) %>%
  pull(trace_window) %>% max
x_min <- trace_windows_all %>% filter(trace_window < Inf) %>%
  pull(trace_window) %>% min


g_trace_window_uncorrected <- ggplot(trace_windows_uncorrected %>% filter(p_traced == 0.8)) +
  geom_histogram(aes(x=trace_window), binwidth=1, 
                 boundary = 7) +
  geom_vline(xintercept = c(2,7,14,28), linetype="dashed", colour="red") +
  facet_grid(trace_type~., scales="free_y") +
  scale_x_continuous(name="Time from tracee exposure to trace initiation",
                     breaks = c(2,seq(0,x_max,7)), limits=c(x_min, x_max)) +
  scale_fill_brewer(name="% contacts\ntraced",
                    labels = function(x) round(as.numeric(x)*100),
                    type = "div", palette = "Set1") +
  theme_bw()
g_trace_window_corrected_080 <- ggplot(trace_windows_corrected %>% filter(p_traced == 0.8)) +
  geom_histogram(aes(x=trace_window), binwidth=1, 
                 boundary = 7) +
  geom_vline(xintercept = c(2,7,14,28), linetype="dashed", colour="red") +
  facet_grid(trace_type~., scales="free_y") +
  scale_x_continuous(name="Time from tracee exposure to trace initiation",
                     breaks = c(2,seq(0,x_max,7)), limits=c(x_min, x_max)) +
  scale_fill_brewer(name="% contacts\ntraced",
                    labels = function(x) round(as.numeric(x)*100),
                    type = "div", palette = "Set1") +
  theme_bw() # TODO: Pick a p_traced to use for this

# Mean and SD
trace_windows_summary <- trace_windows_all %>% 
  group_by(corrected, p_traced, trace_type) %>%
  summarise(window_mean = mean(trace_window[trace_window < Inf]),
            window_sd = sd(trace_window[trace_window < Inf]),
            p_inf = mean(trace_window == Inf))

g_windows_summary <- ggplot(trace_windows_summary, 
                            aes(x=p_traced, y=window_mean, 
                                colour=trace_type, fill = trace_type)) +
  geom_ribbon(aes(ymin=window_mean-window_sd,
                  ymax=window_mean+window_sd), alpha=0.3, colour=NA) +
  geom_line() + geom_point(size=2) +
  facet_grid(.~corrected, 
             labeller=labeller(corrected = setNames(c("Uncorrected", "Corrected"),
                                                    c(FALSE, TRUE)))) +
  scale_y_continuous(name="Mean time from tracee exposure to trace initiation",
                     breaks = seq(0,20,4)) +
  scale_x_continuous(name="% of all contacts traced",
                     labels = function(x) x*100) +
  theme_bw()

# % falling under different thresholds
trace_windows_thresholds <- trace_windows_all %>% 
  group_by(corrected, p_traced, trace_type) %>%
  summarise(p_under_2 = mean(trace_window <= 2),
            p_under_7 = mean(trace_window <= 7),
            p_under_5= mean(trace_window <= 5),
            p_under_12 = mean(trace_window <= 12),
            p_under_14 = mean(trace_window <= 14)) %>%
  gather(threshold, p_under_threshold, -(corrected:trace_type)) %>%
  mutate(threshold = sub("p_under_", "", threshold))

g_windows_thresholds <- ggplot(trace_windows_thresholds, 
                            aes(x=p_traced, y=p_under_threshold, 
                                colour=threshold)) +
  geom_line() + geom_point(size=2) +
  facet_grid(corrected~trace_type,
             labeller=labeller(corrected = setNames(c("Uncorrected", "Corrected"),
                                                    c(FALSE, TRUE)))) +
  scale_y_continuous(name="% of trace windows within threshold",
                     labels = function(x) x*100,
                     breaks = seq(0,1,0.2), limits=c(0,1)) +
  scale_colour_discrete(name="Threshold\n(days)") +
  scale_x_continuous(name="% of contacts traced", labels = function(x) x*100) +
  theme_bw()

#------------------------------------------------------------------------------
# 4. Raw distribution of delay times
#------------------------------------------------------------------------------

tab_delays <- tibble(x=delay_time(1e6))

g_delays <- ggplot(tab_delays) + 
  geom_histogram(aes(x=x), boundary=0, binwidth=1) + 
  scale_y_continuous(name="Frequency", labels = NULL) + 
  scale_x_continuous(name="Delay between symptom onset and isolation/tracing (days)") +
  theme_bw()

#------------------------------------------------------------------------------
# 4. Proportion of hiddens with non-hidden siblings
#------------------------------------------------------------------------------

count_hidden_with_siblings <- function(p_traced, p_asymptomatic, delay_weeks){
  # Create dataset
  cases <- case_data_null %>% copy %>%
    .[, `:=`(p_traced_auto = purrr::rbernoulli(nrow(.), p_traced),
             asym = purrr::rbernoulli(nrow(.), p_asymptomatic),
             hidden = (recovery < delay_weeks*7))]
  # Count siblings in each class
  cases_sibs <- cases %>% group_by(infector_id) %>%
    summarise(n_total = n(), n_hidden_valid = sum(hidden & traceable_fwd),
              n_hidden_invalid = sum(hidden & !traceable_fwd),
              n_unhidden_valid = sum(!hidden & traceable_rev & !asym),
              n_unhidden_invalid = sum(!hidden & !(traceable_rev & !asym))) %>%
    mutate(n_hidden = n_hidden_invalid + n_hidden_valid,
           n_unhidden = n_unhidden_invalid + n_unhidden_valid) %>%
    mutate(n_hidden_count = n_hidden_valid * (n_unhidden_valid > 0))
  # Summarise and return
  cases_count <- cases_sibs %>% ungroup %>%
    summarise(n_hidden_count = sum(n_hidden_count),
              n_hidden_valid = sum(n_hidden_valid),
              n_hidden_total = sum(n_hidden)) %>%
    mutate(p_hidden_count_valid = n_hidden_count/n_hidden_valid,
           p_hidden_count_total = n_hidden_count/n_hidden_total,
           delay = delay_weeks,
           p_asymptomatic = p_asymptomatic,
           p_traced = p_traced)
  return(cases_count)
}

p_hidden_count <- lapply(c(0.25, 0.5), function(p_asym)
  lapply(c(0.8, 1), function(p_traced)
    lapply(1:n_gen_run-1, function(d) count_hidden_with_siblings(p_traced, p_asym, d)) %>%
      rbindlist) %>% rbindlist) %>% rbindlist

g_hidden_count <- ggplot(p_hidden_count %>% filter(!is.nan(p_hidden_count_total)),
                         aes(x=delay, y=p_hidden_count_total,
                             colour = factor(p_asymptomatic))) +
  geom_line(aes(linetype=factor(p_traced))) + geom_point(size=2) +
  scale_x_continuous(name="Delay (weeks)") +
  scale_y_continuous(name="% of hidden cases accessible from an unhidden sibling",
                     labels = function(y) y*100, limits = c(0,1),
                     breaks = seq(0,1,0.2)) +
  scale_colour_brewer(type = "div", palette = "Set1", 
                      name="% asymptomatic\ncarriers",
                      labels = function(x) as.numeric(x)*100) +
  scale_linetype_discrete(name="% contacts traced",
                          labels = function(x) as.numeric(x)*100) +
  theme_bw()

#------------------------------------------------------------------------------
# 3. Chain length
#------------------------------------------------------------------------------

# For an N-generation delay, all identified individuals have a chain of
# at least N hidden ancestors that could be uncovered by multigenerational
# (but not 1-generational) backtracing

list_descendents <- function(cases, ids, include_orig = TRUE){
  # Get a list of all descendents of a set of cases
  ids_out <- numeric(0)
  for (m in 1:n_gen_run){
    new_ids <- cases[infector_id %in% c(ids, ids_out)]$case_id
    ids_out <- c(ids_out, new_ids)
  }
  if (include_orig) ids_out <- c(ids, ids_out)
  return(unique(ids_out))
}

purge_descendents <- function(cases, ids){
  # Remove all cases descended from a list of case IDs
  for (m in 1:n_gen_run){
    new_ids <- cases[infector_id %in% ids]$case_id
    ids <- unique(c(ids, new_ids))
  }
  return(cases[!case_id %in% ids])
}

gen_filter <- function(n){
  # Goal: find all cases that would be the *first* identified member in their
  # lineage given an n-generation response delay
  # 1. Filter out any cases hidden by response delay
  cases <- case_data_null[generation >= n]
  # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
  #    and discard their descendents
  harbingers_sym <- list(); length(harbingers_sym) <- n_gen_run
  for (m in 1:(n_gen_run - n + 1)){
    # a. Find and separate out all symptomatic individuals with no parents
    harbingers_sym[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
    # b. Discard their descendents
    cases <- purge_descendents(cases, harbingers_sym[[m]]$case_id)
    # c. Discard any asymptomatic individuals with no parents, then repeat cycle
    cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
  }
  # 3. Combine into one dataframe
  cases <- rbindlist(harbingers_sym)
  # 4. Count generations and return
  cases_out <- cases[, .(n_cases = .N), by="generation"][
    , `:=`(p_cases = n_cases/sum(n_cases), delay = n, delay_type = "generations")
    ]
  return(cases_out)
}
week_filter <- function(n){
  # Goal: find all cases that would be the *first* identified member in their
  # lineage given an n-week response delay
  # 1. Filter out any cases hidden by response delay
  cases <- case_data_null[recovery >= n*7]
  # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
  #    and discard their descendents
  harbingers_sym <- list(); length(harbingers_sym) <- n_gen_run
  for (m in 1:n_gen_run){
    # a. Find and separate out all symptomatic individuals with no parents
    harbingers_sym[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
    # b. Discard their descendents
    cases <- purge_descendents(cases, harbingers_sym[[m]]$case_id)
    # c. Discard any asymptomatic individuals with no parents, then repeat cycle
    cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
  }
  # 3. Combine into one dataframe
  cases <- rbindlist(harbingers_sym)
  # 4. Count generations and return
  cases_out <- cases[, .(n_cases = .N), by="generation"][
    , `:=`(p_cases = n_cases/sum(n_cases), delay = n, delay_type = "weeks")
    ]
  return(cases_out)
}

gen_delays <- lapply(1:n_gen_run-1, function(n) gen_filter(n)) %>% rbindlist
week_delays <- lapply(1:n_gen_run-1, function(n) week_filter(n)) %>% rbindlist
filters <- rbind(gen_delays, week_delays)

g_filters <- ggplot(filters, aes(x=generation, y=p_cases, colour=delay_type)) +
  geom_line() + geom_point(size=2) +
  facet_grid(delay ~ ., scales="free_y") +
  scale_y_continuous(limits=c(0,1), name="Fraction of chains") +
  scale_x_continuous(name="Hidden chain length") +
  scale_color_discrete(name="Delay type") +
  theme_bw()

# Average chain length
cl_avg <- filters %>% group_by(delay, delay_type) %>% 
  summarise(chain_length_avg = sum(generation*p_cases))
g_cl_avg <- ggplot(cl_avg, aes(x=delay, y=chain_length_avg, colour=delay_type)) +
  geom_line() + geom_point(size=2) +
  scale_x_continuous(name="Delay") +
  scale_y_continuous(name="Average hidden chain length") +
  scale_color_discrete(name="Delay type") +
  theme_bw()

# # ...
# # What % of cases are accessible to 1-generational backtracing after an n-week delay?
# # What % of cases are descended from cases that
# # - (i) are non-hidden
# # - (ii) are hidden, but whose direct descendents are non-hidden
# 
p_cases_accessible_weeks <- function(delay){
  # 1. Filter out any cases hidden by response delay
  cases <- case_data_null[recovery >= delay*7]
  # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
  #    and discard their descendents
  harbingers <- list(); length(harbingers) <- n_gen_run
  for (m in 1:n_gen_run){
    # a. Find and separate out all symptomatic individuals with no parents
    harbingers[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
    # b. Discard their descendents
    cases <- purge_descendents(cases, harbingers[[m]]$case_id)
    # c. Discard any asymptomatic individuals with no parents, then repeat cycle
    cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
  }
  # 3. Combine into one dataframe
  harbingers_all <- rbindlist(harbingers)
  # 4. Count proportion of cases that are descended from all harbingers (and their parents)
  desc_harb <- list_descendents(case_data_null, harbingers_all$case_id, TRUE)
  desc_harb_inf <- c(desc_harb,
                     list_descendents(case_data_null, harbingers_all$infector_id, TRUE)) %>%
    unique %>% .[. != 0]
  t_desc <- tibble(delay = delay, delay_type = "weeks",
                   n_desc_harb = length(desc_harb),
                   n_desc_harb_inf = length(desc_harb_inf))
  return(t_desc)
}

p_cases_accessible_gen <- function(delay){
  # 1. Filter out any cases hidden by response delay
  cases <- case_data_null[generation >= delay]
  # 2. Find all symptomatic individuals with no symptomatic parents, separate them out
  #    and discard their descendents
  harbingers <- list(); length(harbingers) <- n_gen_run
  for (m in 1:n_gen_run){
    # a. Find and separate out all symptomatic individuals with no parents
    harbingers[[m]] <- cases[!(infector_id %in% cases$case_id) & asym==FALSE]
    # b. Discard their descendents
    cases <- purge_descendents(cases, harbingers[[m]]$case_id)
    # c. Discard any asymptomatic individuals with no parents, then repeat cycle
    cases <- cases[(infector_id %in% cases$case_id) | asym==FALSE]
  }
  # 3. Combine into one dataframe
  harbingers_all <- rbindlist(harbingers)
  # 4. Count proportion of cases that are descended from all harbingers (and their parents)
  desc_harb <- list_descendents(case_data_null, harbingers_all$case_id, TRUE)
  desc_harb_inf <- c(desc_harb,
                     list_descendents(case_data_null, harbingers_all$infector_id, TRUE)) %>%
    unique %>% .[. != 0]
  t_desc <- tibble(delay = delay, delay_type = "generations",
                   n_desc_harb = length(desc_harb),
                   n_desc_harb_inf = length(desc_harb_inf))
  return(t_desc)
}

p_acc_gen <- lapply(1:n_gen_run-1, function(n) p_cases_accessible_gen(n)) %>% rbindlist
p_acc_weeks <- lapply(1:n_gen_run-1, function(n) p_cases_accessible_weeks(n)) %>% rbindlist
p_acc <- rbind(p_acc_gen, p_acc_weeks)

g_acc <- p_acc %>% gather(backtrace_distance, n_hit, -delay, -delay_type) %>%
  mutate(backtrace_distance = ifelse(backtrace_distance == "n_desc_harb_inf", 1, 0),
         n_missed = nrow(case_data_null) - n_hit,
         p_missed = n_missed / nrow(case_data_null)) %>%
  ggplot(aes(x=delay, y=n_missed, colour=factor(backtrace_distance))) +
  geom_line(aes(linetype = delay_type)) + geom_point(size=2) +
  theme_bw()
# There might be something salvageable here, but we're not there yet