# Snakefile for running branching-process model with bidirectional tracing

import math
import numpy as np

configfile: "config.yaml"

# Pre-emptively compute number of scenarios (for thread allocation)

spv = config["scenario_parameters"].values()
n_scenarios = 1

for p in spv:
    if not (isinstance(p, int) or isinstance(p, str) or isinstance(p, float)):
        n_scenarios *= len(p)

n_threads = int(np.min([n_scenarios, config["n_threads_max"]]))

# Run rule

rule run_simulation:
    """Run simulation with specified parameter combinations."""
    output: os.path.join("output_files/data", config["simulation_name"]) + \
        "_scenario.tsv.gz"
    log: "log_files/run_simulation_" + config["simulation_name"] + ".log"
    threads: n_threads
    params:
        scenario_parameters = config["scenario_parameters"],
        n_iterations = config["n_iterations"],
        path_prefix = os.path.join("output_files/data", \
                                   config["simulation_name"]),
        output_parameters = config["output_parameters"],
        ci_parameters = config["ci_parameters"],
        report = config["print_reports"],
    conda: config["env_path"]
    script: "workflow/scripts/run_simulation_snakemake.R"
