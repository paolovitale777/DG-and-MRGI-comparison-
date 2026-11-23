## Runs the full study from the repository root: burn-ins, the six index-selection
## simulations (3 scenarios x 2 indices, 5 replicates each), figures and Table S2.
## Each script runs in its own R session (Rscript), as the burn-in scripts use save.image().

scripts <- c(
  "scripts/01_burnin_scenarios_1_2.R",
  "scripts/02_burnin_scenario_3.R",
  "scripts/03_DGI_scenario_1_equal.R",
  "scripts/04_MRGI_scenario_1_equal.R",
  "scripts/05_DGI_scenario_2_uneven.R",
  "scripts/06_MRGI_scenario_2_uneven.R",
  "scripts/07_DGI_scenario_3_realistic.R",
  "scripts/08_MRGI_scenario_3_realistic.R",
  "scripts/09_figures_scenarios_1_2.R",
  "scripts/10_figures_scenario_3.R",
  "scripts/11_table_S2.R")

dir.create("output",  showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
for (s in scripts) {
  message("Running ", s)
  status <- system2("Rscript", s)
  if (status != 0) stop("Script failed: ", s)
}
