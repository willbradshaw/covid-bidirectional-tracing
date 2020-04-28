#' Run simulation within a Snakemake workflow

# Libraries and functions
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(sn))
suppressMessages(library(parallel))
source("workflow/scripts/wrappers.R")
source("workflow/scripts/outbreak.R")

# Process scenario parameters for infinities
sp <- snakemake@params[["scenario_parameters"]]
numeric_keys <- c("r0_base", "rel_r0_asymptomatic", "dispersion", "generation_omega",
                  "generation_alpha", "recovery_quantile", "n_initial_cases",
                  "test_sensitivity", "data_limit_auto", "data_limit_manual",
                  "contact_limit_auto", "contact_limit_manual",
                  "rollout_delay_gen", "rollout_delay_days",
                  "cap_max_generations", "cap_max_weeks", "cap_cases",
                  "backtrace_distance")
for (k in numeric_keys){
    sp[[k]] <- as.numeric(unlist(sp[[k]]))
}

# Run simulation
simulate_process(scenario_parameters = sp,
                 threads = snakemake@threads,
                 n_iterations = snakemake@params[["n_iterations"]],
                 report = snakemake@params[["report"]],
                 show_progress = snakemake@params[["show_progress"]],
                 path_prefix = snakemake@params[["path_prefix"]],
                 write_raw = snakemake@params[["output_parameters"]][["write_raw"]],
                 write_weekly = snakemake@params[["output_parameters"]][["write_weekly"]],
                 write_generational = snakemake@params[["output_parameters"]][["write_gen"]],
                 write_run = snakemake@params[["output_parameters"]][["write_run"]],
                 write_scenario = TRUE, compress_output = TRUE,
                 log_path = snakemake@log[[1]],
                 ci_width = snakemake@params[["ci_parameters"]][["ci_width"]],
                 alpha_prior = snakemake@params[["ci_parameters"]][["alpha_prior"]],
                 beta_prior = snakemake@params[["ci_parameters"]][["beta_prior"]])

