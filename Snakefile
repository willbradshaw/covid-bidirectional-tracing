# Snakefile for running branching-process model with bidirectional tracing

import math

configfile: "config.yaml"

rule run_simulation:
    """Run simulation with specified parameter combinations."""
    output: os.path.join("output_files/data", config["simulation_name"]) + \
        "_scenario.tsv.gz"
    log: "log_files/run_simulation_" + config["simulation_name"] + ".log"
    threads: config["n_threads"]
    params:
        scenario_parameters = config["scenario_parameters"],
        n_iterations = config["n_iterations"],
        parallel = config["parallel"],
        path_prefix = os.path.join("output_files/data", \
                                   config["simulation_name"]),
        output_parameters = config["output_parameters"],
        ci_parameters = config["ci_parameters"],
    conda: config["env_path"]
    script: "workflow/scripts/run_simulation_snakemake.R"
