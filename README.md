# MRGI-simulation

Stochastic simulations comparing the **desired gain index (DGI)** and the **minimum required gain index (MRGI)** in a doubled-haploid wheat breeding program, run with [AlphaSimR](https://cran.r-project.org/package=AlphaSimR).

This repository contains the code used for the simulation study in:

> Martini, Vitale, Cerón-Rojas, Montesinos-López, Crossa, Ortiz. *Minimum Required Gains: A Practical Reframing of Desired-Gain Selection Indices.* (in preparation)

## Overview

- **DGI** (Pesek and Baker, 1969): weights `b = G⁻¹d`, so the expected response is proportional to the desired gains `d`.
- **MRGI** (Itoh and Yamada, 1988): weights that minimise `b'Pb` subject to `Gb ≥ M`, so each trait gains *at least* its minimum required gain `M`. Traits to be decreased are sign-reversed before solving; traits to be maintained enter as equality constraints (`Gb = 0`).

Both indices are applied at every selection step of the same breeding pipeline, with the same number of lines selected at each step. The index only ranks the candidates, so the two methods differ only in which lines are selected.

### Breeding pipeline (every year)

| Step | Lines | Selected |
|---|---|---|
| Crossing block | 100 biparental crosses among 60 parents | |
| DH lines | 100 per F1 (10,000) | |
| Headrows | 3,000 | random (30%) |
| Stage 1 yield trial | 1,000 | 33.3% |
| Stage 2 yield trial | 200 | 20% |
| Stage 3 yield trial | 50 | 25% |
| New parents | 60 from stage 1 and 2 lines | 5% |

Pipeline filling and burn-in: 6 cycles and 5 years of random selection. Index selection: 25 years. 5 replicates per scenario (`set.seed(REP)`, so DGI and MRGI share the same burn-in population and random seed in each replicate).

### Scenarios

| Scenario | Traits | Target gains per selection step (genetic SD) |
|---|---|---|
| 1. Equal | four positively correlated generic traits (T1 to T4) | 1, 1, 1, 1 |
| 2. Uneven | same traits and burn-in as scenario 1 | 1, 0.1, 0.8, 0.1 |
| 3. Realistic | grain yield (GY), disease severity (DS), lodging (LOD), days to anthesis (DTA) | GY +0.60, DS −0.20, LOD −0.10, DTA 0 (maintained) |

Full settings (founder means and variances, heritabilities per stage, genetic correlations) are in `scripts/sim_param_scenarios_1_2.R`, `scripts/sim_param_scenario_3.R` and Supplementary Table S1 of the paper.

## Repository structure

```
scripts/
  MRGI_functions.R               functions only (sourced): DGI and MRGI weights, selection intensity <-> proportion
  sim_param_scenarios_1_2.R      parameters only (sourced): scenarios 1 and 2
  sim_param_scenario_3.R         parameters only (sourced): scenario 3 (GY, DS, LOD, DTA)
  01_burnin_scenarios_1_2.R      pipeline filling + burn-in (shared by scenarios 1 and 2)
  02_burnin_scenario_3.R         pipeline filling + burn-in, scenario 3
  03_DGI_scenario_1_equal.R      25 years of DGI selection, scenario 1
  04_MRGI_scenario_1_equal.R     25 years of MRGI selection, scenario 1
  05_DGI_scenario_2_uneven.R     DGI, scenario 2
  06_MRGI_scenario_2_uneven.R    MRGI, scenario 2
  07_DGI_scenario_3_realistic.R  DGI, scenario 3
  08_MRGI_scenario_3_realistic.R MRGI, scenario 3
  09_figures_scenarios_1_2.R     F6 mean and variance over time + intensity table, scenarios 1 and 2
  10_figures_scenario_3.R        same, scenario 3
  11_table_S2.R                  genetic gain and required selection intensity (Table S2)
run_all.R                        runs the numbered scripts in order
```

The numbered scripts are run; the three unnumbered files are only sourced by them. Simulation outputs (.RData, .rds) are written to `output/` (not tracked by git) and figures and tables to `results/`; both folders are created when needed.

## How to run

Set the working directory to the repository root and run:

```r
install.packages(c("AlphaSimR", "quadprog", "dplyr", "tidyr", "ggplot2"))
source("run_all.R")
```

or run the scripts one at a time in numerical order. All paths are relative to the repository root. The burn-in scripts must be run before the index-selection scripts, and the index-selection scripts before scripts 09 to 11.

### Outputs

For each scenario and index, `output/` gets, per replicate:

- `<scenario>_<REP>.rds`: mean and variance of the genetic values per stage, year and trait;
- `<scenario>_selection_log_<REP>.rds`: for every selection step and year, the number of candidates and selected lines, the applied selection intensity, and the intensity each index would require to reach its target;
- `<scenario>_sessionInfo.txt`: R and package versions.

## Software

Results in the paper were obtained with R 4.6.1, AlphaSimR 2.1.0 and quadprog 1.5-8. Simulation results depend on the AlphaSimR version, so use the same version to reproduce the exact numbers.

## Notes

- The indices use the true genetic (co)variance matrix `G` of the candidates at each step, available in simulation. In practice `G` would be estimated, so the results show the behaviour of the indices without estimation error.
- The burn-in scripts do not set a random seed; rerunning them gives new burn-in populations, and results will differ slightly from those in the paper.
- The breeding program structure is adapted from the AlphaSimR line breeding example by C. Werner (CIMMYT).

## Citation

If you use this code, please cite the paper above and AlphaSimR:

> Gaynor, R.C., Gorjanc, G., Hickey, J.M. (2021). AlphaSimR: an R package for breeding program simulations. *G3 Genes|Genomes|Genetics*, 11(2), jkaa017. https://doi.org/10.1093/g3journal/jkaa017

## License

MIT (see `LICENSE`).
