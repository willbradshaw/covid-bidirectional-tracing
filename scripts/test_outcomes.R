
n.sim = 100
num.initial.cases = 20
prop.ascertain = 0
cap_max_days = 365
cap_cases = 5000
r0isolated = 0
r0community = 2.5
disp.iso = 1
disp.com = 0.16
k = 0.7
delay_shape = 1.651524
delay_scale = 2.305172
prop.asym = 0.3
quarantine = FALSE
sero_test = FALSE
backtrace = FALSE
  

# Outcome for old code
s_old <- scenario_sim(n.sim = n.sim, num.initial.cases = num.initial.cases, prop.ascertain = prop.ascertain, 
                      cap_max_days = cap_max_days, cap_cases = cap_cases, r0isolated = r0isolated,
                      r0community = r0community, disp.iso = disp.iso, disp.com = disp.com, k = k, 
                      delay_shape = delay_shape, delay_scale = delay_scale, prop.asym = prop.asym, 
                      quarantine = quarantine, sero_test = sero_test, backtrace = backtrace,
                      run_old_code = TRUE)

# Outcome for new code
s_new <- scenario_sim(n.sim = n.sim, num.initial.cases = num.initial.cases, prop.ascertain = prop.ascertain, 
                      cap_max_days = cap_max_days, cap_cases = cap_cases, r0isolated = r0isolated,
                      r0community = r0community, disp.iso = disp.iso, disp.com = disp.com, k = k, 
                      delay_shape = delay_shape, delay_scale = delay_scale, prop.asym = prop.asym, 
                      quarantine = quarantine, sero_test = sero_test, backtrace = backtrace,
                      run_old_code = FALSE)

# Process and compare
s_both <- bind_rows(s_old %>% mutate(code="old"), s_new %>% mutate(code="new"))
outcomes <- s_both %>% group_by(code, sim, final_total_cases, final_latest_exposure, 
                                final_extinct, initial_asym) %>%
  filter(between(week, 12, 16)) %>% summarise(trial_cases = sum(weekly_cases)) %>%
  mutate(controlled = (trial_cases == 0) & final_total_cases < 5000)
outcomes_final <- outcomes %>% group_by(code) %>% summarise(p_controlled = mean(controlled))

print(outcomes_final)