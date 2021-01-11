# Bidirectional contact tracing could dramatically improve COVID-19 control

William J. Bradshaw, Ethan C. Alley, Jonathan H. Huggins, Alun L. Lloyd, and Kevin M. Esvelt

DOI: https://doi.org/10.1038/s41467-020-20325-7 (*Nature Communications* paper)

Adapted from code originally published [here](https://github.com/cmmid/ringbp). 
 
## Abstract 

Contact tracing is critical to controlling COVID-19, but most protocols only “forward-trace” to notify people who were recently exposed. Using a stochastic branching-process model, we find that “bidirectional” tracing to identify infector individuals and their other infectees robustly improves outbreak control. In our model, bidirectional tracing more than doubles the reduction in effective reproduction number (R_eff) achieved by forward-tracing alone, while dramatically increasing resilience to low case ascertainment and test sensitivity. The greatest gains are realised by expanding the manual tracing window from 2 to 6 days pre-symptom-onset or, alternatively, by implementing high-uptake smartphone-based exposure notification; however, to achieve the performance of the former approach, the latter requires nearly all smartphones to detect exposure events. With or without exposure notification, our results suggest that implementing bidirectional tracing could dramatically improve COVID-19 control.

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
