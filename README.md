# Bidirectional contact tracing is required for reliable COVID-19 control

William J. Bradshaw, Ethan C. Alley, Jonathan H. Huggins, Alun L. Lloyd, and Kevin M. Esvelt

DOI: https://doi.org/10.1101/2020.05.06.20093369 (medRxiv preprint)

Adapted from code originally published [here](https://github.com/cmmid/ringbp). 
 
## Abstract 

Contact tracing is critical to limiting the spread of pandemics such as COVID-19, but most protocols only “forward-trace” to notify people who were recently exposed. Using a stochastic branching process model, we find that “bidirectional” tracing to identify infector individuals robustly outperforms forward-only approaches across a wide range of scenarios. Rapid smartphone-based exposure notification underperformed conventional manual tracing unless uptake of the digital system was near-universal; however, the combination of manual and digital approaches outperforms either approach used in isolation. Taken together, the combination of manual, digital, and bidirectional tracing more than doubles the probability of controlling outbreaks across three epidemiological scenarios, but only when exposure events can be detected by nearly all smartphones. Implementing combined bidirectional tracing may be critical to controlling COVID-19 without more costly interventions.

## Usage

To run the branching-process model for a given set of assumptions:

- [Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on your system.
- Navigate to the home directory of this project.
- Edit `config.yaml` to include the values of interest; the model will run all combinations of values across all `scenario_parameters` keys. Example config files can be found in `figures/configs`.
- To the model with all available cores, run:
```snakemake --use-conda --cores```
- To specify a number of cores, run:
```snakemake --use-conda --cores <n_cores>```

To generate figures from existing data:

- [Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on your system.
- Navigate to the `figures` directory.
- Our data files for every figure are included in `figures/data`. If you have run your own simulations and want to use those instead, update the paths to data files in `figures/config.yaml`.
- From within the `figures` directory, run Snakemake using the same commands described above. Output PNG images will be written to `figures/output_files`.
